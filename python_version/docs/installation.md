# Installation Guide

## Prerequisites

Before installing CHA-MD-BA, ensure you have the following prerequisites:

- Python 3.8 or higher
- GROMACS 2022 or higher
- pip (Python package installer)

## Installation Steps

1. Clone the repository:
```bash
git clone https://github.com/yourusername/cha-md-ba.git
cd cha-md-ba
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the package:
```bash
pip install -e .
```

## Verifying the Installation

To verify that CHA-MD-BA is installed correctly, run:
```bash
cha-md-ba --version
```

## Troubleshooting

### Common Issues

1. **GROMACS not found**
   - Ensure GROMACS is installed and in your PATH
   - Set the GROMACS environment variables:
     ```bash
     source /path/to/gromacs/bin/GMXRC
     ```

2. **Python dependencies**
   - If you encounter dependency issues, try:
     ```bash
     pip install -r requirements.txt
     ```

3. **Virtual environment issues**
   - If you have problems with the virtual environment, try:
     ```bash
     python -m pip install --upgrade pip
     python -m pip install --upgrade virtualenv
     ```

## Getting Help

If you encounter any issues during installation, please:
1. Check the [FAQ](faq.md)
2. Open an issue on GitHub
3. Contact the development team 