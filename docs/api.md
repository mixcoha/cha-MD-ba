# API Reference

## Package Structure

```
cha_md_ba/
├── __init__.py
├── prepare.py
├── minimize.py
├── nvt.py
├── npt.py
├── production.py
├── analysis.py
└── cli.py
```

## Core Classes

### SystemPreparator

```python
class SystemPreparator:
    def __init__(
        self,
        input_file: str,
        force_field: str = "amber99sb-ildn",
        water_model: str = "tip3p",
        box_type: str = "dodecahedron",
        box_distance: float = 1.0,
        pname: str = "NA",
        nname: str = "CL",
        conc: float = 0.15
    )
```

Prepares the system for molecular dynamics simulations.

**Parameters:**
- `input_file`: Path to input PDB file
- `force_field`: Force field to use
- `water_model`: Water model to use
- `box_type`: Type of simulation box
- `box_distance`: Minimum distance between protein and box edge
- `pname`: Positive ion name
- `nname`: Negative ion name
- `conc`: Ion concentration

**Methods:**
- `prepare()`: Prepares the system
- `generate_topology()`: Generates topology file
- `solvate()`: Solvates the system
- `add_ions()`: Adds ions to the system

### Minimizer

```python
class Minimizer:
    def __init__(
        self,
        max_steps: int = 50000,
        emtol: float = 1000.0,
        emstep: float = 0.01
    )
```

Performs energy minimization of the system.

**Parameters:**
- `max_steps`: Maximum number of steps
- `emtol`: Energy tolerance
- `emstep`: Energy step size

**Methods:**
- `minimize()`: Performs energy minimization
- `check_convergence()`: Checks if minimization converged

### NVTEquilibrator

```python
class NVTEquilibrator:
    def __init__(
        self,
        force_constant: float = 1000.0,
        temperature: float = 300.0,
        dt: float = 0.002,
        nsteps: int = 50000
    )
```

Performs NVT equilibration.

**Parameters:**
- `force_constant`: Position restraint force constant
- `temperature`: Target temperature
- `dt`: Time step
- `nsteps`: Number of steps

**Methods:**
- `equilibrate()`: Performs NVT equilibration
- `check_equilibration()`: Checks if system is equilibrated

### NPTEquilibrator

```python
class NPTEquilibrator:
    def __init__(
        self,
        pressure: float = 1.0,
        temperature: float = 300.0,
        dt: float = 0.002,
        nsteps: int = 50000
    )
```

Performs NPT equilibration.

**Parameters:**
- `pressure`: Target pressure
- `temperature`: Target temperature
- `dt`: Time step
- `nsteps`: Number of steps

**Methods:**
- `equilibrate()`: Performs NPT equilibration
- `check_equilibration()`: Checks if system is equilibrated

### ProductionRunner

```python
class ProductionRunner:
    def __init__(
        self,
        time: float = 100.0,
        dt: float = 0.002,
        temperature: float = 300.0,
        pressure: float = 1.0
    )
```

Runs the production simulation.

**Parameters:**
- `time`: Total simulation time (ns)
- `dt`: Time step
- `temperature`: Target temperature
- `pressure`: Target pressure

**Methods:**
- `run()`: Runs the production simulation
- `save_checkpoint()`: Saves simulation checkpoint

### Analyzer

```python
class Analyzer:
    def __init__(
        self,
        reference: str,
        trajectory: str,
        topology: str = None
    )
```

Analyzes simulation trajectories.

**Parameters:**
- `reference`: Path to reference structure
- `trajectory`: Path to trajectory file
- `topology`: Path to topology file

**Methods:**
- `calculate_rmsd()`: Calculates RMSD
- `calculate_rmsf()`: Calculates RMSF
- `analyze_secondary_structure()`: Analyzes secondary structure
- `analyze_hbonds()`: Analyzes hydrogen bonds
- `cluster_analysis()`: Performs cluster analysis

## Utility Functions

### File Management

```python
def create_directory(path: str) -> None
def remove_directory(path: str) -> None
def copy_file(src: str, dst: str) -> None
```

### System Checks

```python
def check_gromacs() -> bool
def check_dependencies() -> bool
def check_system_requirements() -> bool
```

### Logging

```python
def setup_logging(level: str = "INFO") -> None
def log_message(message: str, level: str = "INFO") -> None
```

## Error Handling

The package uses custom exceptions for error handling:

```python
class CHAMDBAError(Exception)
class PreparationError(CHAMDBAError)
class MinimizationError(CHAMDBAError)
class EquilibrationError(CHAMDBAError)
class ProductionError(CHAMDBAError)
class AnalysisError(CHAMDBAError)
```

## Configuration

The package can be configured using environment variables:

- `CHA_MD_BA_GMXBIN`: Path to GROMACS binaries
- `CHA_MD_BA_FORCEFIELD`: Default force field
- `CHA_MD_BA_WATERMODEL`: Default water model
- `CHA_MD_BA_TEMP`: Default temperature
- `CHA_MD_BA_PRESSURE`: Default pressure 