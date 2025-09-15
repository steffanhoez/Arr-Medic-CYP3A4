#!/bin/bash
# ARR-MEDIC CYP3A4 Installation Script
# Supports both conda and pip installation paths

set -e

echo "🧬 ARR-MEDIC CYP3A4 Opensource Installation"
echo "=========================================="

# Check if conda is available
if command -v conda &> /dev/null; then
    CONDA_AVAILABLE=true
    echo "✅ Conda detected"
else
    CONDA_AVAILABLE=false
    echo "❌ Conda not found"
fi

# Check if python is available
if ! command -v python &> /dev/null && ! command -v python3 &> /dev/null; then
    echo "❌ Python not found. Please install Python 3.8+ first."
    exit 1
fi

PYTHON_CMD=$(command -v python3 2>/dev/null || command -v python)
echo "✅ Python found: $PYTHON_CMD"

# Function to install with conda
install_conda() {
    echo ""
    echo "🐍 Installing with Conda (Recommended for RDKit)"
    echo "================================================"

    # Create conda environment
    read -p "Enter environment name [arr-medic]: " ENV_NAME
    ENV_NAME=${ENV_NAME:-arr-medic}

    echo "Creating conda environment: $ENV_NAME"
    conda create -n "$ENV_NAME" python=3.10 -y

    echo "Activating environment..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"

    echo "Installing RDKit via conda..."
    conda install -c conda-forge rdkit -y

    echo "Installing remaining dependencies..."
    cd backend
    pip install -r requirements.txt

    echo ""
    echo "✅ Installation completed!"
    echo "To activate environment: conda activate $ENV_NAME"
    echo "To start server: uvicorn main:app --reload"
}

# Function to install with pip only
install_pip() {
    echo ""
    echo "📦 Installing with pip (Simplified mode)"
    echo "========================================"
    echo "⚠️  Warning: RDKit will not be available"
    echo "    The predictor will use simplified descriptors"

    read -p "Continue with pip-only installation? [y/N]: " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Installation cancelled."
        exit 1
    fi

    echo "Installing dependencies..."
    cd backend
    pip install -r requirements.txt

    echo ""
    echo "✅ Installation completed (simplified mode)!"
    echo "To start server: uvicorn main:app --reload"
    echo ""
    echo "📝 Note: For full functionality with RDKit, use conda installation"
}

# Function to install with Docker
install_docker() {
    echo ""
    echo "🐳 Installing with Docker"
    echo "========================"

    if ! command -v docker &> /dev/null; then
        echo "❌ Docker not found. Please install Docker first."
        return 1
    fi

    echo "Building Docker image..."
    docker-compose build

    echo "Starting services..."
    docker-compose up -d

    echo ""
    echo "✅ Docker installation completed!"
    echo "API available at: http://localhost:8000"
    echo "To stop: docker-compose down"
}

# Main installation menu
echo ""
echo "Choose installation method:"
echo "1) Conda + RDKit (Recommended)"
echo "2) pip only (Simplified)"
echo "3) Docker"
echo "4) Exit"

read -p "Enter choice [1-4]: " choice

case $choice in
    1)
        if [ "$CONDA_AVAILABLE" = true ]; then
            install_conda
        else
            echo "❌ Conda not available. Please install Miniconda/Anaconda first."
            echo "   Download from: https://docs.conda.io/en/latest/miniconda.html"
            exit 1
        fi
        ;;
    2)
        install_pip
        ;;
    3)
        install_docker
        ;;
    4)
        echo "Installation cancelled."
        exit 0
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo ""
echo "🎉 Setup complete! Check the documentation for next steps."
echo "📖 Docs: http://localhost:8000/docs (after starting the server)"