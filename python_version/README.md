# CHA-MD-BA (Python Version)

This is the new version of CHA-MD-BA implemented in Python, offering a more robust and modular interface for automating molecular dynamics simulations.

## Features

- Complete MD simulation pipeline
- Intuitive command-line interface
- Automated trajectory analysis
- GROMACS integration
- Comprehensive documentation

## Installation

```bash
pip install cha-md-ba
```

## Usage

```bash
# System preparation
cha-md-ba prepare 1tim.pdb output_dir --forcefield amber99sb-ildn --water-model tip3p

# Minimization
cha-md-ba minimize system.gro system.top min_dir

# NVT equilibration
cha-md-ba nvt system.gro system.top nvt_dir

# NPT equilibration
cha-md-ba npt system.gro system.top npt_dir

# Production
cha-md-ba production system.gro system.top prod_dir --num-runs 3

# Analysis
cha-md-ba analyze trajectory.xtc topology.tpr --selection "protein"
```

## Requirements

- Python 3.8+
- GROMACS
- Dependencies listed in `requirements.txt`

## Documentation

Complete documentation is available in the `docs/` folder and online at [Documentation URL].

## Development

To contribute to development:

1. Clone the repository
2. Install development dependencies: `pip install -e ".[dev]"`
3. Run tests: `pytest`
4. Follow contribution guidelines in `CONTRIBUTING.md`

## License

This project is licensed under the MIT License - see the `LICENSE` file for details. 