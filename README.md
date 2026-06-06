# 🧪 Drug Toxicity Predictor

A machine learning–inspired web application that predicts the toxicity of drug candidates using molecular descriptors computed from SMILES notation. Built for early-stage drug screening to reduce reliance on costly lab testing.

🔗 **Live Demo:** [https://toxicity-app-c907.onrender.com](https://toxicity-app-c907.onrender.com)

---

## Overview

In modern drug discovery, early toxicity screening is critical for safety and cost efficiency. This project implements a computational toxicology tool that analyzes chemical structures and flags potentially toxic compounds before clinical trials.

The model evaluates each compound using:
- **Lipinski's Rule of Five** — drug-likeness and absorption
- **Veber's Rules** — oral bioavailability
- **Toxicophore pattern matching** — known dangerous functional groups
- **QSAR heuristics** — from Tox21 and computational chemistry literature

If predicted toxicity score is ≥ 35%, the compound is flagged as **potentially toxic** and not recommended for further development.

---

## Features

- Enter any SMILES string and get an instant toxicity prediction
- Displays 8 molecular descriptors: MW, LogP, H-Bond Donors/Acceptors, TPSA, Rotatable Bonds, Ring Count, Aromatic
- Flags specific toxicity concerns with explanations (nitro groups, heavy metals, reactive functional groups, CYP450 inhibition risk, etc.)
- Shows toxicity probability (0–100%), risk level (Very Low / Low / Moderate / High), and confidence score
- Built-in example compounds for quick testing (Aspirin, Caffeine, Testosterone, Glucose, Phenylalanine, Palmitic acid)

---

## Tech Stack

| Layer | Technology |
|---|---|
| Backend | Python, Flask |
| Server | Gunicorn (WSGI) |
| Frontend | HTML5, CSS3, Vanilla JavaScript |
| Molecular descriptors | Custom SMILES parser (Lipinski/Crippen approximation) |
| Deployment | Render.com |
| Dataset reference | Tox21 (Kaggle) |

---

## Molecular Descriptors Computed

| Descriptor | Description |
|---|---|
| Molecular Weight (MW) | Total mass of the molecule |
| LogP | Lipophilicity (octanol-water partition coefficient) |
| H-Bond Donors (HBD) | Count of OH and NH groups |
| H-Bond Acceptors (HBA) | Count of O and N atoms |
| TPSA | Topological Polar Surface Area (membrane permeability) |
| Rotatable Bonds | Flexibility of the molecule (bioavailability) |
| Ring Count | Number of ring systems |
| Heavy Atoms | Total non-hydrogen atom count |

---

## Toxicity Rules Applied

**Lipinski Rule of Five violations:**
- MW > 500 → poor absorption
- LogP > 5 → lipophilicity concern
- H-Bond Donors > 5
- H-Bond Acceptors > 10

**Veber Rules:**
- Rotatable Bonds > 10 → poor oral bioavailability
- TPSA > 140 Ų → poor membrane permeability

**Toxicophore patterns flagged:**
- Nitro groups `[N+](=O)[O-]` — mutagenicity
- Acid chlorides `C(=O)Cl` — high reactivity
- Heavy metals: Mercury `[Hg]`, Arsenic `[As]`, Lead `[Pb]`, Cadmium `[Cd]`
- Peroxide bonds `OO` — oxidative stress
- Azo groups `N=N` — carcinogenicity
- Nitriles `C#N` — cyanide release risk
- Sulfonyl chlorides — high reactivity
- Trifluoromethyl groups — metabolic concerns

---

## Project Structure

```
drug_toxicity_app/
├── app.py                  # Flask backend — descriptor computation + prediction logic
├── templates/
│   └── index.html          # Frontend UI
├── requirements.txt        # Python dependencies
├── Procfile                # Gunicorn start command
├── render.yaml             # Render.com deployment config
└── README.md               # This file
```

---

## Local Setup

**Requirements:** Python 3.9+

```bash
# 1. Clone the repo
git clone https://github.com/Jaisica77/toxicity-app.git
cd toxicity-app

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run the app
python app.py

# 4. Open in browser
# http://localhost:5000
```

---

## Dependencies

```
flask==3.0.3
gunicorn==22.0.0
numpy==1.26.4
```

No heavy dependencies like RDKit required — the descriptor computation uses a lightweight pure-Python SMILES parser based on Wildman-Crippen LogP approximation and atom-weight lookup tables.

---

## Deployment

Deployed on **Render.com** (free tier) using Gunicorn as the WSGI server.

The `render.yaml` file handles automatic configuration:
```yaml
services:
  - type: web
    name: drug-toxicity-predictor
    env: python
    buildCommand: pip install -r requirements.txt
    startCommand: gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 60
```

---

## Example Test Compounds

| Compound | SMILES | Expected Result |
|---|---|---|
| Nitrobenzene | `O=[N+]([O-])c1ccccc1` | ☠ TOXIC |
| Lead | `[Pb]` | ☠ TOXIC |
| Warfarin | `CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O` | ⚠ TOXIC |
| Aspirin | `CC(=O)Oc1ccccc1C(=O)O` | ✓ NON-TOXIC |
| Caffeine | `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` | ✓ NON-TOXIC |
| Glucose | `OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O` | ✓ NON-TOXIC |

---

## Dataset Reference

This project is based on the **Tox21 Dataset** — a publicly available toxicology dataset widely used in computational drug discovery.

Dataset: [https://www.kaggle.com/datasets/epicskills/tox21-dataset](https://www.kaggle.com/datasets/epicskills/tox21-dataset)

Original notebook: [drug_toxicity_predictor_AiML](https://github.com/Jaisica77/drug_toxicity_predictor_AiML)

---

## Disclaimer

This tool is intended strictly for **educational and research purposes**.
Toxicity predictions do not replace laboratory experiments, clinical trials, or regulatory approval processes.

---

## Author

**Jaisica** 
