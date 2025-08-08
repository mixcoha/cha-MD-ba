# Contributing to CHA-MD-BA

Thank you for your interest in contributing to CHA-MD-BA! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

By participating in this project, you agree to abide by our Code of Conduct. Please be respectful and considerate of others.

## How to Contribute

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## Development Setup

1. Clone your fork:
```bash
git clone https://github.com/yourusername/cha-md-ba.git
cd cha-md-ba
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install development dependencies:
```bash
pip install -e ".[dev]"
```

## Coding Standards

- Follow PEP 8 style guide
- Use type hints
- Write docstrings for all public functions and classes
- Keep functions focused and small
- Write tests for new features

## Testing

Run tests with:
```bash
pytest
```

## Documentation

- Update documentation when adding new features
- Keep docstrings up to date
- Add examples for new features
- Update README if necessary

## Pull Request Process

1. Ensure your code passes all tests
2. Update documentation
3. Add appropriate tests
4. Submit pull request with description of changes

## Issue Reporting

When reporting issues, please include:
1. Description of the problem
2. Steps to reproduce
3. Expected behavior
4. Actual behavior
5. Environment details

## Feature Requests

For feature requests, please:
1. Describe the feature
2. Explain the use case
3. Provide examples if possible

## Release Process

1. Update version number
2. Update changelog
3. Create release tag
4. Build and upload to PyPI

## Questions?

Feel free to open an issue or contact the maintainers. 