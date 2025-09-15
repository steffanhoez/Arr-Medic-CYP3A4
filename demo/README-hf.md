---
title: ARR-MEDIC CYP3A4 Demo
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: Interactive CYP3A4 drug interaction prediction with molecular visualization
tags:
- drug-discovery
- molecular-prediction
- cyp3a4
- pharmacology
- cheminformatics
- educational
---

# 🧬 ARR-MEDIC CYP3A4 Interactive Demo

**Real-time CYP3A4 Drug Interaction Prediction with Molecular Visualization**

This demo allows you to predict CYP3A4 inhibition for small molecules using SMILES notation, complete with molecular structure visualization.

## 🚀 How to Use

1. **Enter SMILES**: Input a molecular structure in SMILES format
2. **Add Name** (optional): Give your compound an identifier
3. **Predict**: Click the prediction button to analyze CYP3A4 inhibition potential
4. **View Results**: See prediction probability, confidence score, and molecular structure

## 🧪 Try These Examples

Click the example buttons or copy these SMILES:

- **Ethanol**: `CCO` (simple alcohol, likely non-inhibitor)
- **Caffeine**: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` (mild stimulant)
- **Ibuprofen**: `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O` (common NSAID)
- **Celecoxib**: `CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F`

## 🔬 About CYP3A4

CYP3A4 is the most important drug-metabolizing enzyme, responsible for metabolizing over 50% of clinically used drugs. Understanding CYP3A4 inhibition helps predict:

- **Drug-drug interactions** - preventing dangerous combinations
- **Patient safety** - especially important for elderly patients on multiple medications
- **Dosing adjustments** - when co-administered drugs affect metabolism
- **Drug development** - early screening to avoid late-stage failures

## 📊 Model Information

- **Type**: Rule-based molecular descriptor analysis
- **Accuracy**: ~70% (educational baseline for learning purposes)
- **Features**: Molecular weight, LogP, TPSA, hydrogen bonds, aromatic rings
- **Visualization**: 2D molecular structure rendering with RDKit

## ⚠️ Important Disclaimers

🔴 **NOT FOR CLINICAL USE**: This tool is intended **only for research and educational purposes**.

- ❌ Not validated for clinical decision-making
- ❌ Not for patient care or treatment planning
- ❌ Not for drug development or regulatory submissions
- ✅ Perfect for learning drug metabolism concepts
- ✅ Great for educational demonstrations and research

## 🔗 Links

- **GitHub Repository**: [ARR-MEDIC CYP3A4 Opensource](https://github.com/Flamehaven/Arr-Medic-CYP3A4)
- **Documentation**: Full installation guides and API documentation
- **License**: MIT License (Open Source)

---

*Built with Gradio + RDKit + Docker for reliable molecular visualization*