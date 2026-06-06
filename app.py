from flask import Flask, render_template, request, jsonify
import numpy as np

app = Flask(__name__)

# ─── Molecular descriptor computation (pure Python, no RDKit needed) ───────────
# We use Lipinski/Veber-style rules based on known SMILES patterns.
# For a full deployment with RDKit, see README_DEPLOY.md

ATOM_WEIGHTS = {
    'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007,
    'S': 32.06,  'F': 18.998, 'Cl': 35.45, 'Br': 79.904,
    'I': 126.904, 'P': 30.974
}

def parse_smiles_descriptors(smiles: str) -> dict:
    """
    Extract molecular descriptors from SMILES string.
    Uses pattern-matching for a lightweight, dependency-free approximation
    matching the feature set from the original notebook.
    """
    s = smiles.strip()
    if not s:
        raise ValueError("Empty SMILES string")

    # Molecular Weight approximation
    mw = 0.0
    i = 0
    while i < len(s):
        if i+1 < len(s) and s[i:i+2] in ATOM_WEIGHTS:
            mw += ATOM_WEIGHTS[s[i:i+2]]
            i += 2
        elif s[i] in ATOM_WEIGHTS:
            mw += ATOM_WEIGHTS[s[i]]
            i += 1
        else:
            i += 1
    # Add implicit hydrogens (rough estimate)
    c_count = s.count('C') - s.lower().count('cl')
    mw += c_count * 1.5  # average implicit H

    # LogP (Wildman-Crippen approximation)
    logp = 0.0
    logp += s.count('C') * 0.53
    logp += s.count('F') * 0.14
    logp += s.count('Cl') * 0.60
    logp += s.count('Br') * 0.88
    logp += s.count('I') * 1.12
    logp -= s.count('O') * 0.67
    logp -= s.count('N') * 0.35
    logp -= s.count('S') * 0.09
    logp -= s.count('[OH]') * 0.4
    logp -= s.count('(=O)') * 0.5

    # H-bond donors (OH, NH patterns)
    hbd = s.count('[OH]') + s.count('[NH]') + s.count('[NH2]') + s.count('O') // 4

    # H-bond acceptors (O, N atoms)
    hba = s.count('O') + s.count('N') + s.count('n')

    # Rotatable bonds (single bonds not in ring, rough count)
    rotatable = max(0, s.count('-') + s.count('C') // 3 - s.count('1') - s.count('2') - 2)

    # TPSA approximation
    tpsa = 0.0
    tpsa += s.count('O') * 9.23
    tpsa += s.count('N') * 12.36
    tpsa += s.count('[OH]') * 4.5
    tpsa += s.count('[NH]') * 3.2
    tpsa += s.count('(=O)') * 8.6

    # Ring count
    ring_count = max(s.count('1'), s.count('2'), s.count('3'), s.count('4'))

    # Aromatic ring indicator
    aromatic = 1 if ('c' in s or 'n' in s or 'o' in s) else 0

    return {
        'MW': round(mw, 2),
        'LogP': round(logp, 3),
        'HBD': hbd,
        'HBA': hba,
        'RotatableBonds': rotatable,
        'TPSA': round(tpsa, 2),
        'RingCount': ring_count,
        'Aromatic': aromatic,
    }


def predict_toxicity(descriptors: dict) -> dict:
    """
    Rule-based toxicity model based on:
    - Lipinski's Rule of Five
    - Veber's rules
    - Tox21 / QSAR literature thresholds
    
    Returns a toxicity probability and verdict.
    """
    score = 0.0
    flags = []

    mw     = descriptors['MW']
    logp   = descriptors['LogP']
    hbd    = descriptors['HBD']
    hba    = descriptors['HBA']
    tpsa   = descriptors['TPSA']
    rot    = descriptors['RotatableBonds']
    rings  = descriptors['RingCount']

    # --- Lipinski violations (each violation ~+15% toxicity risk) ---
    if mw > 500:
        score += 0.15
        flags.append(f"High MW ({mw:.0f} > 500) – poor absorption")
    if logp > 5:
        score += 0.15
        flags.append(f"High LogP ({logp:.2f} > 5) – lipophilicity concern")
    if hbd > 5:
        score += 0.10
        flags.append(f"High H-bond donors ({hbd} > 5)")
    if hba > 10:
        score += 0.10
        flags.append(f"High H-bond acceptors ({hba} > 10)")

    # --- Veber rules ---
    if rot > 10:
        score += 0.10
        flags.append(f"Many rotatable bonds ({rot} > 10) – poor oral bioavailability")
    if tpsa > 140:
        score += 0.15
        flags.append(f"High TPSA ({tpsa:.0f} > 140 Å²) – poor membrane permeability")

    # --- Toxicophore heuristics ---
    if logp < -3:
        score += 0.10
        flags.append(f"Very low LogP ({logp:.2f}) – renal toxicity risk")
    if mw < 100:
        score += 0.08
        flags.append(f"Very low MW ({mw:.0f}) – reactive small molecule")
    if rings > 5:
        score += 0.10
        flags.append(f"Many ring systems ({rings} > 5) – complex/flat molecule")

    # Positive factors (reduce score)
    if 150 <= mw <= 500 and 0 <= logp <= 4 and tpsa <= 90:
        score -= 0.05
        flags.append("✓ Drug-like physicochemical profile (low risk)")

    # Clamp between 0 and 1
    score = max(0.0, min(1.0, score))

    verdict = "TOXIC" if score >= 0.50 else "NON-TOXIC"
    risk_level = (
        "High Risk"   if score >= 0.65 else
        "Moderate Risk" if score >= 0.40 else
        "Low Risk"
    )

    return {
        'toxicity_probability': round(score * 100, 1),
        'verdict': verdict,
        'risk_level': risk_level,
        'flags': flags,
    }


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    smiles = data.get('smiles', '').strip()

    if not smiles:
        return jsonify({'error': 'Please enter a SMILES string.'}), 400

    try:
        descriptors = parse_smiles_descriptors(smiles)
        result = predict_toxicity(descriptors)
        return jsonify({
            'descriptors': descriptors,
            'result': result,
            'smiles': smiles,
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
