\# Changelog

All notable changes to this project will be documented in this file.  

The format is based on \[Keep a Changelog](https://keepachangelog.com/en/1.1.0/),  

and this project adheres to \[Semantic Versioning](https://semver.org/).



\## \[1.1.0] - 2025-09-15

\### Added 🌐
\*\*Multilingual Interface Support\*\*
\- \*\*Korean-English Language Toggle\*\* - Real-time interface switching in Gradio demo
\- \*\*Complete UI Translation\*\* - All labels, buttons, messages, and results support both languages
\- \*\*Localized Prediction Results\*\* - Molecular descriptors and warnings display in selected language
\- \*\*Educational Enhancement\*\* - Better accessibility for Korean-speaking researchers and students

🚀 \*\*Demo Improvements\*\*
\- \*\*Enhanced Gradio Interface\*\* - Improved user experience with language selector
\- \*\*HF Spaces Deployment\*\* - Successfully deployed to https://huggingface.co/spaces/Flamehaven/arr-medic-cyp3a4-demo
\- \*\*Mobile-Friendly Design\*\* - Responsive interface works on all devices
\- \*\*Pre-loaded Examples\*\* - Quick testing with Ethanol, Caffeine, and Ibuprofen compounds

🔧 \*\*Technical Updates\*\*
\- \*\*Fixed RDKit Dependencies\*\* - Stable conda-based installation approach
\- \*\*Improved Docker Configuration\*\* - Updated Dockerfile with optimized conda/pip setup
\- \*\*Enhanced Requirements\*\* - Compatible gradio~=3.50 with rdkit-pypi integration
\- \*\*Docker Compose Enhancement\*\* - Added multilingual demo service

\### Changed
\- \*\*Updated Requirements.txt\*\* - Switched to gradio~=3.50 for better Pydantic v1 compatibility
\- \*\*Enhanced README.md\*\* - Added comprehensive multilingual feature documentation
\- \*\*Improved Docker Setup\*\* - Optimized container build process for scientific packages

\### Fixed
\- \*\*IndentationError Resolution\*\* - Fixed Python syntax issues in deployed app.py
\- \*\*HF Spaces Compatibility\*\* - Resolved dependency conflicts for stable deployment
\- \*\*RDKit Integration\*\* - Stable molecular visualization across platforms

\## \[1.0.1] - 2025-09-14

\### Added
\- Automated installation scripts (`scripts/install.sh`, `scripts/install.bat`)
\- Comprehensive installation guide (`docs/INSTALLATION.md`) with troubleshooting
\- Enhanced requirements.txt with detailed RDKit installation guidance
\- Interactive Jupyter notebook demo (`notebooks/demo_colab.ipynb`)
\- GitHub issue templates for bug reports and feature requests
\- Security middleware: GZip compression and TrustedHost validation
\- Multiple installation pathways: Conda+RDKit (recommended), pip-only, Docker
\- Cross-platform installation support (Windows/Linux/macOS)

\### Changed
\- Enhanced README with automated installer instructions and roadmap
\- Updated disclaimer section with clinical usage warning
\- Improved installation documentation with conda/pip options
\- Enhanced requirements.txt with clear dependency categories and comments

\### Fixed
\- Removed `.pytest_cache/` from version control
\- Improved dependency management for RDKit installation challenges



\## \[1.0.0] - 2025-09-15

\### Added

\- Initial open-source release of \*\*ARR-MEDIC CYP3A4 Opensource\*\*

\- FastAPI backend with REST API endpoints:

&nbsp; - `/predict` for single compound prediction

&nbsp; - `/predict/batch` for batch predictions

&nbsp; - `/health` system health check

&nbsp; - `/stats` basic API usage statistics

\- Rule-based CYP3A4 predictor (~70% accuracy, educational baseline)

\- SQLite database integration for prediction history

\- Basic CI workflow with pytest and Bandit security scan

\- Dockerfile and docker-compose configuration

\- Example `.env.example` for configuration

\- Pytest-based API tests



\### Documentation

\- Initial `README.md` with installation, usage, and API examples

\- Added \*\*⚠️ Research/Education disclaimer\*\* (not for clinical or diagnostic use)

\- Contribution guidelines (`CONTRIBUTING.md`)



\### Security

\- Integrated Bandit security scan into GitHub Actions CI workflow





