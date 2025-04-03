# User Guide

## Introduction

CHA-MD-BA is a comprehensive pipeline for molecular dynamics simulations. This guide will walk you through the main features and workflows.

## Basic Workflow

The typical workflow consists of the following steps:

1. System Preparation
2. Energy Minimization
3. NVT Equilibration
4. NPT Equilibration
5. Production Run
6. Analysis

## System Preparation

### Input Requirements

- PDB file of your protein structure
- Force field parameters
- Solvent model selection

### Command Line Usage

```bash
cha-md-ba prepare --input protein.pdb --force-field amber99sb-ildn --water-model tip3p
```

### Python API Usage

```python
from cha_md_ba import prepare

preparator = prepare.SystemPreparator(
    input_file="protein.pdb",
    force_field="amber99sb-ildn",
    water_model="tip3p"
)
preparator.prepare()
```

## Energy Minimization

### Purpose

Energy minimization ensures the system is in a stable configuration before starting dynamics.

### Command Line Usage

```bash
cha-md-ba minimize --max-steps 50000 --emtol 1000
```

### Python API Usage

```python
from cha_md_ba import minimize

minimizer = minimize.Minimizer(
    max_steps=50000,
    emtol=1000
)
minimizer.minimize()
```

## NVT Equilibration

### Purpose

NVT equilibration brings the system to the desired temperature.

### Command Line Usage

```bash
cha-md-ba nvt --force-constant 1000 --temperature 300
```

### Python API Usage

```python
from cha_md_ba import nvt

nvt_equilibrator = nvt.NVTEquilibrator(
    force_constant=1000,
    temperature=300
)
nvt_equilibrator.equilibrate()
```

## NPT Equilibration

### Purpose

NPT equilibration brings the system to the desired pressure.

### Command Line Usage

```bash
cha-md-ba npt --pressure 1.0 --temperature 300
```

### Python API Usage

```python
from cha_md_ba import npt

npt_equilibrator = npt.NPTEquilibrator(
    pressure=1.0,
    temperature=300
)
npt_equilibrator.equilibrate()
```

## Production Run

### Purpose

The production run generates the actual trajectory for analysis.

### Command Line Usage

```bash
cha-md-ba production --time 100 --dt 0.002
```

### Python API Usage

```python
from cha_md_ba import production

prod = production.ProductionRunner(
    time=100,
    dt=0.002
)
prod.run()
```

## Analysis

### Available Analysis Tools

- RMSD calculation
- RMSF analysis
- Secondary structure analysis
- Hydrogen bond analysis
- Cluster analysis

### Command Line Usage

```bash
cha-md-ba analyze --type rmsd --reference reference.pdb --trajectory traj.xtc
```

### Python API Usage

```python
from cha_md_ba import analysis

analyzer = analysis.Analyzer(
    reference="reference.pdb",
    trajectory="traj.xtc"
)
rmsd = analyzer.calculate_rmsd()
```

## Best Practices

1. Always check the system preparation step carefully
2. Monitor the energy minimization convergence
3. Verify temperature and pressure equilibration
4. Save checkpoints regularly during production runs
5. Validate analysis results with multiple methods

## Common Issues and Solutions

1. **System instability**
   - Check force field parameters
   - Verify solvent model compatibility
   - Adjust minimization parameters

2. **Temperature/pressure not equilibrating**
   - Check thermostat/barostat settings
   - Verify system size and composition
   - Consider longer equilibration times

3. **Analysis errors**
   - Ensure trajectory files are complete
   - Check reference structure alignment
   - Verify analysis parameters

## Getting Help

For additional support:
1. Check the [FAQ](faq.md)
2. Review the [API Reference](api.md)
3. Open an issue on GitHub 