# ARR-MEDIC CYP3A4 Opensource

## ⚠️ Disclaimer

🔴 **Not for Clinical or Diagnostic Use**
This project is intended **only for research and educational purposes**.
Do **not** use it in clinical decision-making, patient care, or regulatory submissions.

---

🧬 **Open Source CYP3A4 Drug Interaction Prediction System**

> MIT License | For Academic Research and Learning

## 🎯 Overview

ARR-MEDIC CYP3A4 Opensource is a simplified, community-driven version of the CYP3A4 drug metabolic prediction system. This version provides basic drug interaction prediction capabilities with ~70% accuracy, perfect for:

- 🎓 **Academic Research**
- 🔬 **Educational Purposes** 
- 🚀 **Proof of Concept Development**
- 🌐 **Community Contributions**

## ⚡ Quick Start

### Option 1: Docker (Recommended)
```bash
# Clone repository
git clone https://github.com/your-org/arr-medic-cyp3a4-opensource
cd arr-medic-cyp3a4-opensource

# Start with Docker Compose
docker-compose up -d

# Check API health
curl http://localhost:8000/health
```

### Option 2: Local Development

#### Install with pip (simplified)
```bash
# Clone repository
git clone https://github.com/your-org/arr-medic-cyp3a4-opensource
cd arr-medic-cyp3a4-opensource

# Set up backend
cd backend
pip install -r requirements.txt
cp .env.example .env

# Start API server
uvicorn main:app --reload --host 0.0.0.0 --port 8000

# API will be available at http://localhost:8000
# Documentation at http://localhost:8000/docs
```

#### Install with conda (recommended for RDKit)
```bash
# Clone repository
git clone https://github.com/your-org/arr-medic-cyp3a4-opensource
cd arr-medic-cyp3a4-opensource

# Create conda environment
conda create -n arr-medic python=3.10
conda activate arr-medic

# Install RDKit via conda (recommended)
conda install -c conda-forge rdkit

# Install remaining dependencies
cd backend
pip install -r requirements.txt
cp .env.example .env

# Start API server
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

## 🏗️ Architecture

```
arr-medic-cyp3a4-opensource/
├── backend/              # FastAPI prediction service
├── frontend/             # React web interface  
├── models/              # Basic ML models
├── data/                # Sample datasets
├── docker/              # Container configs
└── docs/                # Documentation
```

⚠️ **Disclaimer**: This project is for **research and educational purposes only**.
It is **not intended for clinical or diagnostic use**.

## Roadmap
- v1.0: Rule-based baseline (~70% accuracy)
- v2.0: Add RDKit descriptors + ML (RandomForest, XGBoost)
- v3.0: GNN / Transformer support

For advanced **clinical-grade Pro version (90%+ accuracy, Flamehaven)**,
please see commercial offerings. Migration toolkit (logs/DB transfer) will be provided.

---

## Example: Run in Google Colab

You can try the predictor interactively in a Jupyter/Colab notebook.

```python
!pip install arr-medic-cyp3a4-opensource fastapi uvicorn

from predictor import predict

print(predict("CCN"))  # simple demo
```

## 🔬 Features

✅ **Core Features**:
- Basic CYP3A4 inhibition prediction
- Single & batch prediction endpoints
- REST API with OpenAPI documentation
- Async SQLite database storage
- Docker containerization support
- Simple molecular descriptor analysis
- Comprehensive error handling

❌ **Not Included**:
- Advanced clinical ethics integration
- Flame-based emotional therapies
- Hospital EMR/FHIR connectivity
- Real-time patient data processing
- Premium ML models

## 📊 Performance

- **Accuracy**: ~70% on standard datasets
- **Throughput**: 100-500 predictions/minute
- **Latency**: < 2 seconds per prediction
- **Memory**: < 2GB RAM required
- **Batch Processing**: Up to 100 compounds per request

## 🔧 API Endpoints

- `GET /` - API information
- `GET /health` - Health check
- `POST /predict` - Single compound prediction
- `POST /predict/batch` - Batch prediction (up to 100 compounds)
- `GET /docs` - Interactive API documentation

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🔗 Related Projects

- **ARR-MEDIC Professional**: Commercial version with 90%+ accuracy
- **Flamehaven Platform**: Enterprise medical AI ecosystem

---

**⚠️ Disclaimer**: This opensource version is for research purposes only. Not intended for clinical use without proper validation and regulatory approval.