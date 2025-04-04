# CHA-MD-BA (Bash Version)

This is the original version of CHA-MD-BA implemented in bash, designed to automate molecular dynamics simulations using GROMACS.

## Script Structure

- `analisis.sh`: Script for trajectory analysis
- `centrar.sh`: Molecule centering
- `centrar_prot.sh`: Protein-specific centering
- `ext_carpeta_cola_proc_cg.sh`: Coarse-grain simulation extension
- `extend_carpeta_cola_proc_ns.sh`: Simulation extension
- `kk.sh`: Utility script
- `minimiza.sh`: Energy minimization
- `minimiza_atomico.sh`: Atomic system minimization
- `npt.sh`: NPT simulation
- `nvt2.sh`: NVT simulation
- `plot_several.sh`: Plot generation

## Usage

To use this version, ensure you have GROMACS installed and configured on your system. The scripts are designed to be executed in sequential order:

1. System preparation
2. Minimization
3. NVT equilibration
4. NPT equilibration
5. Production
6. Analysis

## Requirements

- GROMACS
- Bash shell
- Basic Unix utilities

## Documentation

For detailed documentation, please refer to the `docs/` folder in the root repository. 