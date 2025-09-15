#!/bin/bash
# Deploy ARR-MEDIC CYP3A4 Demo to Hugging Face Spaces

set -e

USERNAME="$1"
PROJECT_NAME="arr-medic-cyp3a4"
DEMO_NAME="$PROJECT_NAME-demo"

if [ -z "$USERNAME" ]; then
    echo "Usage: $0 <HF_USERNAME>"
    echo "Example: $0 Flamehaven"
    exit 1
fi

REPO_URL="https://huggingface.co/spaces/$USERNAME/$DEMO_NAME"

echo "🚀 Deploying ARR-MEDIC CYP3A4 Demo to Hugging Face Spaces"
echo "========================================================="
echo "Username: $USERNAME"
echo "Demo URL: $REPO_URL"
echo ""

# Check if demo directory exists
if [ ! -d "demo" ]; then
    echo "❌ Error: demo directory not found"
    echo "Please run this script from the project root directory"
    exit 1
fi

# Create temporary deployment directory
DEPLOY_DIR="temp_deploy"
rm -rf "$DEPLOY_DIR"
mkdir -p "$DEPLOY_DIR"
cd "$DEPLOY_DIR"

echo "📥 Cloning HF Space repository..."
git clone "$REPO_URL" "$DEMO_NAME" 2>/dev/null || {
    echo "📝 Repository doesn't exist, will create new one"
    mkdir "$DEMO_NAME"
    cd "$DEMO_NAME"
    git init
    git remote add origin "$REPO_URL"
}

cd "$DEMO_NAME"

echo "📋 Copying demo files..."
cp ../../demo/app.py ./
cp ../../demo/requirements.txt ./
cp ../../demo/README.md ./
cp ../../demo/predictor.py ./
cp ../../demo/models.py ./
cp ../../demo/.env ./

# Create .gitignore
cat > .gitignore << 'EOF'
__pycache__/
*.pyc
*.pyo
*.pyd
.Python
env/
venv/
.venv/
pip-log.txt
pip-delete-this-directory.txt
.tox/
.coverage
.coverage.*
.cache
nosetests.xml
coverage.xml
*.log
.env
.DS_Store
EOF

echo "📝 Setting up git configuration..."
git add .
git status

# Create commit message
COMMIT_MSG="Deploy ARR-MEDIC CYP3A4 interactive demo

Features:
- SMILES-based CYP3A4 inhibition prediction
- Interactive Gradio interface
- Molecular visualization (when RDKit available)
- Educational examples and documentation

🔴 Research/Educational use only - Not for clinical use

Built with Gradio for Hugging Face Spaces"

echo "💾 Committing changes..."
git commit -m "$COMMIT_MSG" || echo "No changes to commit"

echo "🚀 Pushing to Hugging Face Spaces..."
git push -u origin main

echo ""
echo "✅ Deployment completed!"
echo "🌐 Demo URL: $REPO_URL"
echo "⏳ It may take a few minutes to build and deploy"
echo ""
echo "📊 Next steps:"
echo "1. Visit the demo URL to see your Space"
echo "2. Check the build logs on HF Spaces"
echo "3. Test the demo functionality"
echo "4. Share with the community!"

# Cleanup
cd ../../..
rm -rf "$DEPLOY_DIR"

echo ""
echo "🎉 ARR-MEDIC CYP3A4 demo is now live on Hugging Face Spaces!"