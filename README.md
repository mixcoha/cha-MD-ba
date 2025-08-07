# Proyecto GROMACS Avanzado

Este es un proyecto avanzado de simulación molecular utilizando GROMACS, diseñado para diferentes tipos de simulaciones: atómicas, coarse grain y proteínas transmembranales.

## Estructura del Proyecto

```
GMX/
├── input/                    # Archivos de entrada organizados por tipo
│   ├── atomic/              # Configuraciones para simulaciones atómicas
│   ├── coarse_grain/        # Configuraciones para simulaciones CG
│   └── transmembrane/       # Configuraciones para proteínas transmembranales
├── output/                  # Archivos de salida organizados por tipo
├── scripts/                 # Scripts de automatización y análisis
│   ├── atomic/              # Scripts para simulaciones atómicas
│   ├── cg/                  # Scripts para simulaciones coarse grain
│   └── transmembrane/       # Scripts para proteínas transmembranales
├── analysis/                # Resultados de análisis organizados por tipo
└── docs/                    # Documentación del proyecto
```

## Requisitos

- GROMACS (versión recomendada: 2023.x o superior)
- Python 3.8+ (para scripts de análisis)
- VMD o PyMOL (para visualización)

## Tipos de Simulación Soportados

### 1. Simulaciones Atómicas
- **Uso**: Para estudios detallados de interacciones moleculares
- **Campo de fuerza**: AMBER99SB-ILDN, CHARMM36, OPLS-AA
- **Script**: `scripts/run_simulation.sh`

### 2. Simulaciones de Coarse Grain
- **Uso**: Para procesos de larga escala temporal
- **Campos de fuerza**: MARTINI, ELNEDIN, SIRAH
- **Script**: `scripts/cg/setup_cg_simulation.sh`

### 3. Proteínas Transmembranales
- **Uso**: Para proteínas embebidas en membranas
- **Tipos de membrana**: POPC, DPPC, DOPC, DMPC
- **Script**: `scripts/transmembrane/setup_membrane.sh`

## Uso

### Simulaciones Atómicas
1. Coloca tu archivo PDB en `input/`
2. Ejecuta: `./scripts/run_simulation.sh nombre_proteina`
3. Analiza: `python scripts/analyze_trajectory.py output/md.xtc output/md.gro`

### Simulaciones Coarse Grain
1. Coloca tu archivo PDB en `input/`
2. Ejecuta: `./scripts/cg/setup_cg_simulation.sh proteina.pdb MARTINI`
3. Analiza con scripts específicos de CG

### Proteínas Transmembranales
1. Coloca tu archivo PDB en `input/`
2. Ejecuta: `./scripts/transmembrane/setup_membrane.sh proteina.pdb POPC`
3. Analiza: `python scripts/transmembrane/analyze_membrane.py output/md.xtc output/md.gro`

## Comandos Útiles

```bash
# Verificar versión de GROMACS
gmx --version

# Listar módulos disponibles
gmx help
```

## Notas

- Asegúrate de tener suficiente espacio en disco para las simulaciones
- Documenta todos los parámetros utilizados en cada simulación
- Haz respaldos regulares de tus datos importantes
