@echo off
REM ARR-MEDIC CYP3A4 Installation Script for Windows
setlocal enabledelayedexpansion

echo 🧬 ARR-MEDIC CYP3A4 Opensource Installation
echo ==========================================

REM Check if conda is available
conda --version >nul 2>&1
if %errorlevel%==0 (
    set CONDA_AVAILABLE=true
    echo ✅ Conda detected
) else (
    set CONDA_AVAILABLE=false
    echo ❌ Conda not found
)

REM Check if python is available
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ❌ Python not found. Please install Python 3.8+ first.
    pause
    exit /b 1
)

echo ✅ Python found
python --version

echo.
echo Choose installation method:
echo 1) Conda + RDKit (Recommended)
echo 2) pip only (Simplified)
echo 3) Docker
echo 4) Exit

set /p choice=Enter choice [1-4]:

if "%choice%"=="1" goto conda_install
if "%choice%"=="2" goto pip_install
if "%choice%"=="3" goto docker_install
if "%choice%"=="4" goto exit_script
goto invalid_choice

:conda_install
if "%CONDA_AVAILABLE%"=="false" (
    echo ❌ Conda not available. Please install Miniconda/Anaconda first.
    echo    Download from: https://docs.conda.io/en/latest/miniconda.html
    pause
    exit /b 1
)

echo.
echo 🐍 Installing with Conda (Recommended for RDKit)
echo ================================================

set /p env_name=Enter environment name [arr-medic]:
if "%env_name%"=="" set env_name=arr-medic

echo Creating conda environment: %env_name%
call conda create -n %env_name% python=3.10 -y

echo Activating environment...
call conda activate %env_name%

echo Installing RDKit via conda...
call conda install -c conda-forge rdkit -y

echo Installing remaining dependencies...
cd backend
pip install -r requirements.txt

echo.
echo ✅ Installation completed!
echo To activate environment: conda activate %env_name%
echo To start server: uvicorn main:app --reload
goto end

:pip_install
echo.
echo 📦 Installing with pip (Simplified mode)
echo ========================================
echo ⚠️  Warning: RDKit will not be available
echo     The predictor will use simplified descriptors
echo.

set /p confirm=Continue with pip-only installation? [y/N]:
if /i not "%confirm%"=="y" (
    echo Installation cancelled.
    goto end
)

echo Installing dependencies...
cd backend
pip install -r requirements.txt

echo.
echo ✅ Installation completed (simplified mode)!
echo To start server: uvicorn main:app --reload
echo.
echo 📝 Note: For full functionality with RDKit, use conda installation
goto end

:docker_install
echo.
echo 🐳 Installing with Docker
echo ========================

docker --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ❌ Docker not found. Please install Docker first.
    goto end
)

echo Building Docker image...
docker-compose build

echo Starting services...
docker-compose up -d

echo.
echo ✅ Docker installation completed!
echo API available at: http://localhost:8000
echo To stop: docker-compose down
goto end

:invalid_choice
echo Invalid choice. Exiting.
goto end

:exit_script
echo Installation cancelled.
goto end

:end
echo.
echo 🎉 Setup complete! Check the documentation for next steps.
echo 📖 Docs: http://localhost:8000/docs (after starting the server)
pause