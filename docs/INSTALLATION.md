# Installation Guide

This guide provides detailed instructions for installing ARR-MEDIC CYP3A4 Opensource on different platforms and environments.

## Quick Start (Automated)

### Linux/macOS
```bash
# Make the script executable and run
chmod +x scripts/install.sh
./scripts/install.sh
```

### Windows
```cmd
# Run the installation script
scripts\install.bat
```

## Manual Installation

### Option 1: Conda Installation (Recommended)

This is the recommended approach as RDKit works best when installed via conda.

#### Step 1: Install Miniconda/Anaconda
If you don't have conda installed:
- **Miniconda** (lightweight): https://docs.conda.io/en/latest/miniconda.html
- **Anaconda** (full package): https://www.anaconda.com/download

#### Step 2: Create Environment
```bash
# Create a new environment with Python 3.10
conda create -n arr-medic python=3.10
conda activate arr-medic
```

#### Step 3: Install RDKit
```bash
# Install RDKit via conda-forge (this is crucial for molecular descriptors)
conda install -c conda-forge rdkit
```

#### Step 4: Install Other Dependencies
```bash
cd backend
pip install -r requirements.txt
```

#### Step 5: Verify Installation
```bash
python -c "
from rdkit import Chem
from predictor import CYP3A4BasicPredictor
predictor = CYP3A4BasicPredictor()
print('✅ Installation successful!')
print('RDKit available:', predictor.is_ready())
"
```

### Option 2: Pip-Only Installation (Simplified)

⚠️ **Warning**: This method doesn't include RDKit, so the predictor will use simplified molecular descriptors with reduced accuracy.

```bash
cd backend
pip install -r requirements.txt
```

**What you'll miss without RDKit:**
- Accurate molecular descriptors (MW, LogP, TPSA, etc.)
- SMILES validation
- Full molecular property analysis

The predictor will fall back to string-based heuristics.

### Option 3: Docker Installation

Perfect for production deployment or if you want everything containerized.

```bash
# Build and start services
docker-compose up -d

# Check if running
curl http://localhost:8000/health
```

## Troubleshooting

### Common Issues

#### 1. RDKit Installation Fails
```bash
# Error: "Could not find a version that satisfies the requirement rdkit"
# Solution: Use conda instead of pip for RDKit
conda install -c conda-forge rdkit
```

#### 2. Import Error: "No module named 'rdkit'"
```bash
# Check if you're in the right environment
conda info --envs
conda activate arr-medic

# Verify RDKit installation
python -c "from rdkit import Chem; print('RDKit OK')"
```

#### 3. Permission Errors (Linux/macOS)
```bash
# Make scripts executable
chmod +x scripts/install.sh

# If pip install fails
pip install --user -r requirements.txt
```

#### 4. Windows Path Issues
```cmd
REM Use forward slashes or escaped backslashes
python "backend\main.py"

REM Or use the Windows-specific paths
cd backend && python main.py
```

#### 5. Port Already in Use
```bash
# Check what's using port 8000
lsof -i :8000  # macOS/Linux
netstat -ano | findstr :8000  # Windows

# Use a different port
uvicorn main:app --host 0.0.0.0 --port 8080
```

### Environment-Specific Notes

#### macOS (Apple Silicon M1/M2)
```bash
# RDKit works well with conda on Apple Silicon
conda install -c conda-forge rdkit

# If you get architecture issues:
arch -arm64 conda install -c conda-forge rdkit
```

#### Windows
```cmd
REM Use Anaconda Prompt for best conda experience
REM Install Microsoft Visual C++ Redistributable if needed

REM For PowerShell users:
conda init powershell
```

#### Linux
```bash
# Install system dependencies if needed (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install python3-dev build-essential

# For CentOS/RHEL
sudo yum groupinstall "Development Tools"
sudo yum install python3-devel
```

## Performance Optimization

### Development Mode
```bash
# Fast reload for development
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

### Production Mode
```bash
# Multi-worker production setup
gunicorn main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000
```

### Memory Optimization
```bash
# Set environment variables for lower memory usage
export PRED_PROB_FLOOR=0.1
export PRED_PROB_CEIL=0.9
```

## Testing Installation

### Quick Test
```python
import requests

# Test health endpoint
response = requests.get("http://localhost:8000/health")
print(response.json())

# Test prediction
response = requests.post(
    "http://localhost:8000/predict",
    json={"smiles": "CCO", "compound_id": "ethanol"}
)
print(response.json())
```

### Run Test Suite
```bash
cd backend
pytest tests/ -v
```

## Next Steps

After successful installation:

1. **Start the API server**: `uvicorn main:app --reload`
2. **Access documentation**: http://localhost:8000/docs
3. **Try the Colab notebook**: Open `notebooks/demo_colab.ipynb`
4. **Read the API guide**: See main README.md for usage examples

## Getting Help

If you encounter issues not covered here:

1. Check the [GitHub Issues](../../issues) for known problems
2. Create a new issue using our bug report template
3. Join the discussion in our community channels

## Uninstallation

### Conda Environment
```bash
conda deactivate
conda env remove -n arr-medic
```

### Docker
```bash
docker-compose down --volumes
docker rmi arr-medic-cyp3a4-opensource_backend
```

### Pip Installation
```bash
pip uninstall -r requirements.txt
```