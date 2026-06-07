from flask import Flask, render_template, request, jsonify
import numpy as np

app = Flask(__name__)

ATOM_WEIGHTS = {
    'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007,
    'S': 32.06,  'F': 18.998, 'Cl': 35.45, 'Br': 79.904,
    'I': 126.904, 'P': 30.974
}

def parse_smiles_descriptors(smiles: str) -> dict:
    s = smiles.strip()
    if not s:
        raise ValueError("Empty SMILES string")

    # Molecular Weight
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
    c_count = s.count('C') - s.lower().count('cl')
    mw += c_count * 1.5

    # LogP
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

    # H-bond donors
    hbd = s.count('[OH]') + s.count('[NH]') + s.count('[NH2]') + s.count('O') // 4

    # H-bond acceptors
    hba = s.count('O') + s.count('N') + s.count('n')

    # Rotatable bonds
    rotatable = max(0, s.count('-') + s.count('C') // 3 - s.count('1') - s.count('2') - 2)

    # TPSA
    tpsa = 0.0
    tpsa += s.count('O') * 9.23
    tpsa += s.count('N') * 12.36
    tpsa += s.count('[OH]') * 4.5
    tpsa += s.count('[NH]') * 3.2
    tpsa += s.count('(=O)') * 8.6

    # Ring count
    ring_count = max(s.count('1'), s.count('2'), s.count('3'), s.count('4'))

    # Aromatic
    aromatic = 1 if ('c' in s or 'n' in s or 'o' in s) else 0

    # Heavy atom count
    heavy_atoms = sum(1 for c in s if c.isupper())

    return {
        'MW': round(mw, 2),
        'LogP': round(logp, 3),
        'HBD': hbd,
        'HBA': hba,
        'RotatableBonds': rotatable,
        'TPSA': round(tpsa, 2),
        'RingCount': ring_count,
        'Aromatic': aromatic,
        'HeavyAtoms': heavy_atoms,
    }


def predict_toxicity(descriptors: dict, smiles: str) -> dict:
    score = 0.0
    flags = []

    mw          = descriptors['MW']
    logp        = descriptors['LogP']
    hbd         = descriptors['HBD']
    hba         = descriptors['HBA']
    tpsa        = descriptors['TPSA']
    rot         = descriptors['RotatableBonds']
    rings       = descriptors['RingCount']
    heavy_atoms = descriptors['HeavyAtoms']
    s           = smiles

    # ── Lipinski Rule of Five ─────────────────────────────────────────
    if mw > 500:
        score += 0.15
        flags.append(f"High MW ({mw:.0f} > 500) – poor absorption")
    if logp > 5:
        score += 0.15
        flags.append(f"High LogP ({logp:.2f} > 5) – lipophilicity concern")
    if logp > 8:
        score += 0.10
        flags.append(f"Very high LogP ({logp:.2f} > 8) – severe accumulation risk")
    if hbd > 5:
        score += 0.10
        flags.append(f"High H-bond donors ({hbd} > 5)")
    if hba > 10:
        score += 0.10
        flags.append(f"High H-bond acceptors ({hba} > 10)")

    # ── Veber Rules ───────────────────────────────────────────────────
    if rot > 10:
        score += 0.10
        flags.append(f"Many rotatable bonds ({rot} > 10) – poor oral bioavailability")
    if tpsa > 140:
        score += 0.15
        flags.append(f"High TPSA ({tpsa:.0f} > 140 Å²) – poor membrane permeability")

    # ── Toxicophore Patterns ──────────────────────────────────────────
    toxic_patterns = [
        # Classic dangerous groups
        ('[N+](=O)[O-]', 0.25, 'Nitro group – mutagenicity risk'),
        ('C(=O)Cl',      0.25, 'Acid chloride – highly reactive'),
        ('[N-]=[N+]',    0.20, 'Diazo group – carcinogenicity risk'),
        ('SS',           0.15, 'Disulfide – hepatotoxicity risk'),
        ('C(F)(F)F',     0.15, 'Trifluoromethyl – metabolic concern'),
        # Heavy metals
        ('[Hg]',         0.40, 'Mercury – heavy metal toxicity'),
        ('[As]',         0.40, 'Arsenic – heavy metal toxicity'),
        ('[Pb]',         0.40, 'Lead – heavy metal toxicity'),
        ('[Cd]',         0.40, 'Cadmium – heavy metal toxicity'),
        # Reactive groups
        ('N=N',          0.18, 'Azo group – potential carcinogen'),
        ('C=C=C',        0.15, 'Allene group – reactive intermediate'),
        ('OO',           0.18, 'Peroxide bond – oxidative stress risk'),
        ('[N+]',         0.10, 'Quaternary nitrogen – membrane disruption'),
        ('C#N',          0.12, 'Nitrile group – cyanide release risk'),
        ('S(=O)(=O)Cl',  0.22, 'Sulfonyl chloride – highly reactive'),
        # Aromatic toxicity
        ('c1ccccc1',     0.12, 'Benzene ring – aromatic toxicity risk'),
        ('c1ccc2ccccc2c1', 0.20, 'Polycyclic aromatic – carcinogenicity risk'),
        # Receptor activity / bioactivity flags
        ('OCCN',         0.12, 'Aminoethanol group – membrane disruption'),
        ('C=CC=C',       0.15, 'Conjugated diene – reactive metabolite risk'),
        ('[nH]',         0.10, 'Aromatic NH – bioactivity flag'),
        # Halogenated compounds
        ('CBr',          0.12, 'Alkyl bromide – alkylating agent risk'),
        ('CI',           0.14, 'Alkyl iodide – highly reactive'),
        ('CCl',          0.10, 'Alkyl chloride – potential hepatotoxin'),
    ]

    for pattern, penalty, reason in toxic_patterns:
        if pattern.lower() in s.lower():
            score += penalty
            flags.append(f"⚠ {reason}")

    # ── Heavy Atom Checks ─────────────────────────────────────────────
    if heavy_atoms > 40:
        score += 0.10
        flags.append(f"Large molecule ({heavy_atoms} heavy atoms) – clearance concern")
    if heavy_atoms < 5:
        score += 0.12
        flags.append(f"Very small molecule ({heavy_atoms} atoms) – high reactivity risk")

    # ── Additional Heuristics ─────────────────────────────────────────
    if logp < -3:
        score += 0.10
        flags.append(f"Very low LogP ({logp:.2f}) – renal toxicity risk")
    if mw < 100:
        score += 0.08
        flags.append(f"Very low MW ({mw:.0f}) – reactive small molecule")
    if rings > 5:
        score += 0.10
        flags.append(f"Many ring systems ({rings} > 5) – complex/flat molecule")

    # CYP450 inhibition risk
    if hba > 8 and logp > 3:
        score += 0.08
        flags.append("High HBA + LogP – possible CYP450 inhibition risk")

    # ── Positive Factors (reduce score) ───────────────────────────────
    if 150 <= mw <= 500 and 0 <= logp <= 4 and tpsa <= 90:
        score -= 0.02
        flags.append("✓ Drug-like physicochemical profile (low risk)")

    if 5 <= heavy_atoms <= 30 and logp <= 3:
        score -= 0.02
        flags.append("✓ Appropriate size and polarity")

    # Clamp between 0 and 1
    score = max(0.0, min(1.0, score))

    verdict = "TOXIC" if score >= 0.25 else "NON-TOXIC"

    risk_level = (
        "High Risk"     if score >= 0.55 else
        "Moderate Risk" if score >= 0.35 else
        "Low Risk"      if score >= 0.20 else
        "Very Low Risk"
    )

    confidence = (
        "High"   if len(flags) >= 4 else
        "Medium" if len(flags) >= 2 else
        "Low"
    )

    return {
        'toxicity_probability': round(score * 100, 1),
        'verdict': verdict,
        'risk_level': risk_level,
        'confidence': confidence,
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
        result = predict_toxicity(descriptors, smiles)
        return jsonify({
            'descriptors': descriptors,
            'result': result,
            'smiles': smiles,
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
