\# Changelog

All notable changes to this project will be documented in this file.  

The format is based on \[Keep a Changelog](https://keepachangelog.com/en/1.1.0/),  

and this project adheres to \[Semantic Versioning](https://semver.org/).



\## \[Unreleased]

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





