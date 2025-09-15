---
title: ARR-MEDIC CYP3A4 Demo
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: gradio
app_file: app.py
pinned: false
license: mit
short_description: Interactive CYP3A4 drug interaction prediction demo
tags:
- drug-discovery
- molecular-prediction
- cyp3a4
- pharmacology
- cheminformatics
- educational
---

# 🧬 ARR-MEDIC CYP3A4 Opensource Demo

**Interactive CYP3A4 Drug Interaction Prediction System**

This demo allows you to predict CYP3A4 inhibition for small molecules using SMILES notation. It's designed for research and educational purposes.

## 🚀 How to Use

1. **Enter SMILES**: Input a molecular structure in SMILES format
2. **Add ID** (optional): Give your compound a name or identifier
3. **Predict**: Click the prediction button to analyze the compound
4. **View Results**: See the prediction, confidence, and molecular visualization

## 🧪 Example Compounds

Try these example compounds to get started:

- **Ethanol**: `CCO` (simple alcohol)
- **Caffeine**: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` (stimulant)
- **Ibuprofen**: `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O` (NSAID)

## 🔬 About CYP3A4

CYP3A4 is the most important drug-metabolizing enzyme in humans:

- **Location**: Primarily in liver and intestines
- **Function**: Metabolizes ~50% of all medications
- **Importance**: CYP3A4 inhibition can cause drug-drug interactions

## 📊 Model Information

- **Type**: Rule-based molecular descriptor analysis
- **Accuracy**: ~70% (educational baseline)
- **Features**: Molecular weight, LogP, TPSA, hydrogen bonds
- **Visualization**: 2D molecular structure rendering

## ⚠️ Important Disclaimers

🔴 **NOT FOR CLINICAL USE**: This tool is intended **only for research and educational purposes**. Do not use for:
- Clinical decision-making
- Patient care
- Drug development
- Regulatory submissions

The predictions are based on a simplified model and should not be considered clinically validated.

## 🛠️ Technical Details

- **Framework**: Gradio + FastAPI backend
- **Molecular Descriptors**: RDKit (with fallback mode)
- **Prediction Method**: Rule-based scoring system
- **Visualization**: RDKit molecule rendering

## 🔗 Related Links

- **GitHub Repository**: [ARR-MEDIC CYP3A4 Opensource](https://github.com/your-org/arr-medic-cyp3a4-opensource)
- **Documentation**: Full installation and API guide
- **License**: MIT License (Open Source)

## 📚 Learn More

### What is CYP3A4 Inhibition?

When a compound inhibits CYP3A4, it can:
- Slow down metabolism of other drugs
- Increase blood levels of co-administered medications
- Potentially cause adverse effects or toxicity

### Common CYP3A4 Inhibitors

- **Strong**: Ketoconazole, Ritonavir, Clarithromycin
- **Moderate**: Erythromycin, Diltiazem, Verapamil
- **Weak**: Cimetidine, Ranitidine

### Molecular Properties Associated with Inhibition

- Higher molecular weight (300-600 Da)
- Moderate lipophilicity (LogP 2-5)
- Aromatic rings and heteroatoms
- Specific functional groups (azoles, macrolides)

---

**Built with ❤️ for the Open Source Community**

*This demo showcases the power of open-source drug discovery tools and educational resources.*