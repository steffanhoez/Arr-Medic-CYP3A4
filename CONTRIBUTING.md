# Contributing to ARR-MEDIC CYP3A4 Opensource

We welcome contributions to the ARR-MEDIC CYP3A4 Opensource project! This guide will help you get started.

## 🚀 Quick Start

1. **Fork the repository**
2. **Clone your fork**
   ```bash
   git clone https://github.com/your-username/arr-medic-cyp3a4-opensource.git
   cd arr-medic-cyp3a4-opensource
   ```
3. **Set up development environment**
   ```bash
   cd backend
   pip install -r requirements.txt
   cp .env.example .env
   ```

## 🛠️ Development Setup

### Backend Development
```bash
cd backend

# Install dependencies
pip install -r requirements.txt

# Run tests
pytest tests/ -v

# Start development server
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

### Using Docker
```bash
# Build and run with Docker Compose
docker-compose up -d

# View logs
docker-compose logs -f backend
```

## 🧪 Testing

### Running Tests
```bash
cd backend
pytest tests/ -v

# With coverage
pip install pytest-cov
pytest tests/ --cov=. --cov-report=html
```

### Manual API Testing
```bash
# Health check
curl http://localhost:8000/health

# Single prediction
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "compound_id": "ethanol"}'

# Batch prediction
curl -X POST http://localhost:8000/predict/batch \
  -H "Content-Type: application/json" \
  -d '{"compounds": [{"smiles": "CCO", "compound_id": "ethanol"}]}'
```

## 📝 Code Style

### Python Code Style
- Follow PEP 8
- Use type hints
- Maximum line length: 100 characters
- Use docstrings for all functions and classes

### Example
```python
async def predict_compound(
    smiles: str,
    compound_id: Optional[str] = None
) -> Dict[str, Any]:
    """
    Predict CYP3A4 interaction for a compound.

    Args:
        smiles: SMILES string representation
        compound_id: Optional compound identifier

    Returns:
        Prediction result dictionary
    """
    # Implementation here
    pass
```

## 🐛 Bug Reports

When reporting bugs, please include:

1. **Clear description** of the issue
2. **Steps to reproduce** the problem
3. **Expected vs actual behavior**
4. **Environment details** (Python version, OS, etc.)
5. **Error messages** or logs if applicable

### Bug Report Template
```markdown
**Bug Description**
A clear description of what the bug is.

**To Reproduce**
1. Go to '...'
2. Click on '....'
3. See error

**Expected Behavior**
What you expected to happen.

**Environment**
- OS: [e.g. Ubuntu 20.04]
- Python: [e.g. 3.10.0]
- Version: [e.g. 1.0.0]
```

## ✨ Feature Requests

For feature requests, please:

1. **Check existing issues** to avoid duplicates
2. **Describe the use case** clearly
3. **Explain the benefit** to users
4. **Consider backwards compatibility**

## 🔄 Pull Request Process

### Before Submitting
1. **Create an issue** first to discuss major changes
2. **Fork and create a branch** from `main`
3. **Write tests** for new functionality
4. **Update documentation** if needed
5. **Ensure tests pass** locally

### PR Guidelines
- **Clear title** and description
- **Reference related issues** (e.g., "Fixes #123")
- **Keep PRs focused** - one feature per PR
- **Include tests** for new code
- **Update documentation** if needed

### PR Template
```markdown
**Description**
Brief description of changes.

**Type of Change**
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement

**Testing**
- [ ] Tests pass locally
- [ ] Added tests for new functionality
- [ ] Manual testing completed

**Checklist**
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
```

## 📚 Areas for Contribution

### 🔬 Science & Algorithms
- Improve prediction accuracy
- Add new molecular descriptors
- Implement ensemble methods
- Optimize model performance

### 🛠️ Engineering
- API improvements
- Database optimizations
- Containerization enhancements
- CI/CD pipeline improvements

### 📖 Documentation
- API documentation
- Tutorial content
- Code examples
- Performance benchmarks

### 🧪 Testing
- Unit test coverage
- Integration tests
- Performance tests
- Security testing

## 📄 Licensing

By contributing, you agree that your contributions will be licensed under the MIT License.

## 💬 Community

- **Issues**: Use GitHub Issues for bug reports and feature requests
- **Discussions**: Use GitHub Discussions for questions and ideas
- **Email**: opensource@arr-medic.com for sensitive issues

## 🙏 Recognition

Contributors will be recognized in:
- `CONTRIBUTORS.md` file
- Release notes for significant contributions
- Project documentation

---

Thank you for contributing to ARR-MEDIC CYP3A4 Opensource! 🚀