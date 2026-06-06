# Drug Toxicity Predictor — Deployment Guide

A complete web application built on the [drug_toxicity_predictor_AiML](https://github.com/Jaisica77/drug_toxicity_predictor_AiML) project.

---

## Files in this package

```
drug_toxicity_app/
├── app.py              # Flask backend + ML prediction logic
├── templates/
│   └── index.html      # Full frontend UI
├── requirements.txt    # Python dependencies
├── Procfile            # Heroku / Railway start command
├── render.yaml         # Render.com deployment config
└── README_DEPLOY.md    # This file
```

---

## Option 1 — Deploy on Render (FREE, recommended)

### Step 1 — Push to GitHub
1. Create a new GitHub repo (e.g., `drug-toxicity-predictor`)
2. Copy all files from this package into it
3. Push to GitHub:
```bash
git init
git add .
git commit -m "Initial commit"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/drug-toxicity-predictor.git
git push -u origin main
```

### Step 2 — Deploy on Render
1. Go to [render.com](https://render.com) and sign up (free)
2. Click **"New +"** → **"Web Service"**
3. Connect your GitHub account and select your repo
4. Render auto-detects `render.yaml` — click **"Create Web Service"**
5. Wait ~2 minutes for the build to complete
6. Your app will be live at: `https://drug-toxicity-predictor.onrender.com`

> **Free tier note:** Render free services spin down after 15 min of inactivity.
> First request after sleep takes ~30 seconds to wake up. Upgrade to $7/month to keep it always-on.

---

## Option 2 — Run Locally

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Start the server
python app.py

# 3. Open in browser
# http://localhost:5000
```

---

## Option 3 — Deploy on Railway (alternative free option)

1. Go to [railway.app](https://railway.app) → New Project → Deploy from GitHub Repo
2. Select your repo
3. Railway will auto-detect Python and use the Procfile
4. Done — you get a live URL instantly

---

## Upgrading with RDKit (for production accuracy)

The current `app.py` uses a lightweight descriptor approximation for zero-dependency deployment.
To use the full RDKit pipeline from the original notebook:

1. Install RDKit: `pip install rdkit`
2. Replace the `parse_smiles_descriptors` function in `app.py` with:

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def parse_smiles_descriptors(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string — could not parse molecule")
    return {
        'MW': round(Descriptors.ExactMolWt(mol), 2),
        'LogP': round(Descriptors.MolLogP(mol), 3),
        'HBD': rdMolDescriptors.CalcNumHBD(mol),
        'HBA': rdMolDescriptors.CalcNumHBA(mol),
        'RotatableBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'TPSA': round(Descriptors.TPSA(mol), 2),
        'RingCount': rdMolDescriptors.CalcNumRings(mol),
        'Aromatic': int(rdMolDescriptors.CalcNumAromaticRings(mol) > 0),
    }
```

3. Add `rdkit` to `requirements.txt`

> Note: RDKit install can fail on some Render free tier builds due to memory limits.
> Use a paid plan or Railway for RDKit support.

---

## Disclaimer

This tool is for **educational and research purposes only**.
Predictions do not replace laboratory experiments, clinical trials, or regulatory approval.
