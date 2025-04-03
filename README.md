# CHA-MD-BA: Molecular Dynamics Simulation Pipeline

A comprehensive Python package for molecular dynamics simulations, focusing on protein structure analysis and dynamics.

## Overview

CHA-MD-BA provides a complete pipeline for molecular dynamics simulations, including:
- Protein structure preparation
- Energy minimization
- NVT equilibration
- NPT equilibration
- Production runs
- Analysis tools

## Installation

```bash
pip install cha-md-ba
```

## Requirements

- Python 3.8+
- GROMACS 2022+
- NumPy
- MDAnalysis
- Rich

## Usage

The package provides both a command-line interface and a Python API:

### Command Line Interface

```bash
cha-md-ba prepare --input protein.pdb
cha-md-ba minimize
cha-md-ba nvt --force-constant 1000
cha-md-ba npt
cha-md-ba production
```

### Python API

```python
from cha_md_ba import prepare, minimize, nvt, npt

# Prepare the system
preparator = prepare.SystemPreparator("protein.pdb")
preparator.prepare()

# Minimize the system
minimizer = minimize.Minimizer()
minimizer.minimize()

# NVT equilibration
nvt_equilibrator = nvt.NVTEquilibrator()
nvt_equilibrator.equilibrate(force_constant=1000)

# NPT equilibration
npt_equilibrator = npt.NPTEquilibrator()
npt_equilibrator.equilibrate()
```

## Documentation

Detailed documentation is available in the `docs` directory:
- [Installation Guide](docs/installation.md)
- [User Guide](docs/user_guide.md)
- [API Reference](docs/api.md)
- [Examples](docs/examples.md)

## Contributing

Contributions are welcome! Please read our [Contributing Guidelines](CONTRIBUTING.md) before submitting pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- GROMACS development team
- MDAnalysis developers
- All contributors to this project

Edgar Mixcoha
    
    
    
    
    
    
    
   
    
    
