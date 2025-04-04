# Examples

## Basic Usage

### System Preparation

```python
from cha_md_ba import prepare

# Initialize the preparator
preparator = prepare.SystemPreparator(
    input_file="protein.pdb",
    force_field="amber99sb-ildn",
    water_model="tip3p"
)

# Prepare the system
preparator.prepare()
```

### Energy Minimization

```python
from cha_md_ba import minimize

# Initialize the minimizer
minimizer = minimize.Minimizer(
    max_steps=50000,
    emtol=1000.0
)

# Perform minimization
minimizer.minimize()
```

### NVT Equilibration

```python
from cha_md_ba import nvt

# Initialize the NVT equilibrator
nvt_equilibrator = nvt.NVTEquilibrator(
    force_constant=1000.0,
    temperature=300.0
)

# Perform NVT equilibration
nvt_equilibrator.equilibrate()
```

### NPT Equilibration

```python
from cha_md_ba import npt

# Initialize the NPT equilibrator
npt_equilibrator = npt.NPTEquilibrator(
    pressure=1.0,
    temperature=300.0
)

# Perform NPT equilibration
npt_equilibrator.equilibrate()
```

### Production Run

```python
from cha_md_ba import production

# Initialize the production runner
prod = production.ProductionRunner(
    time=100.0,  # 100 ns
    dt=0.002
)

# Run the production simulation
prod.run()
```

### Analysis

```python
from cha_md_ba import analysis

# Initialize the analyzer
analyzer = analysis.Analyzer(
    reference="reference.pdb",
    trajectory="traj.xtc"
)

# Calculate RMSD
rmsd = analyzer.calculate_rmsd()

# Calculate RMSF
rmsf = analyzer.calculate_rmsf()

# Analyze secondary structure
ss = analyzer.analyze_secondary_structure()

# Analyze hydrogen bonds
hbonds = analyzer.analyze_hbonds()

# Perform cluster analysis
clusters = analyzer.cluster_analysis()
```

## Advanced Usage

### Custom Force Field

```python
from cha_md_ba import prepare

# Initialize with custom force field
preparator = prepare.SystemPreparator(
    input_file="protein.pdb",
    force_field="charmm36",
    water_model="tip3p"
)

# Add custom parameters
preparator.add_custom_parameters("custom.itp")
preparator.prepare()
```

### Restrained Simulation

```python
from cha_md_ba import nvt

# Initialize with position restraints
nvt_equilibrator = nvt.NVTEquilibrator(
    force_constant=1000.0,
    temperature=300.0,
    restraints=["backbone", "CA"]
)

# Perform restrained equilibration
nvt_equilibrator.equilibrate()
```

### Temperature Ramp

```python
from cha_md_ba import nvt

# Initialize with temperature ramp
nvt_equilibrator = nvt.NVTEquilibrator(
    force_constant=1000.0,
    temperature=300.0,
    temperature_ramp={
        "start": 100.0,
        "end": 300.0,
        "steps": 50000
    }
)

# Perform temperature ramp
nvt_equilibrator.equilibrate()
```

### Custom Analysis

```python
from cha_md_ba import analysis

# Initialize with custom selections
analyzer = analysis.Analyzer(
    reference="reference.pdb",
    trajectory="traj.xtc",
    selections={
        "protein": "protein",
        "backbone": "backbone",
        "sidechains": "not backbone"
    }
)

# Calculate custom properties
properties = analyzer.calculate_properties(
    properties=["radius_of_gyration", "sasa"],
    selections=["protein", "backbone"]
)
```

## Command Line Examples

### System Preparation

```bash
cha-md-ba prepare \
    --input protein.pdb \
    --force-field amber99sb-ildn \
    --water-model tip3p \
    --box-type dodecahedron \
    --box-distance 1.0
```

### Energy Minimization

```bash
cha-md-ba minimize \
    --max-steps 50000 \
    --emtol 1000.0 \
    --emstep 0.01
```

### NVT Equilibration

```bash
cha-md-ba nvt \
    --force-constant 1000.0 \
    --temperature 300.0 \
    --dt 0.002 \
    --nsteps 50000
```

### NPT Equilibration

```bash
cha-md-ba npt \
    --pressure 1.0 \
    --temperature 300.0 \
    --dt 0.002 \
    --nsteps 50000
```

### Production Run

```bash
cha-md-ba production \
    --time 100.0 \
    --dt 0.002 \
    --temperature 300.0 \
    --pressure 1.0
```

### Analysis

```bash
cha-md-ba analyze \
    --type rmsd \
    --reference reference.pdb \
    --trajectory traj.xtc \
    --selection protein
``` 