# CHA-MD-BA: Pipeline Avanzado de Simulación de Dinámica Molecular

**CHA-MD-BA** es un pipeline completo y avanzado de simulación de dinámica molecular para el análisis de proteínas, que incluye soporte para simulaciones atómicas, coarse grain y proteínas transmembranales.

## 🚀 Características Principales

- **Simulaciones Atómicas**: Estudios detallados de interacciones moleculares
- **Simulaciones Coarse Grain**: Procesos de larga escala temporal
- **Proteínas Transmembranales**: Análisis de proteínas embebidas en membranas
- **Automatización Completa**: Scripts para todo el pipeline de simulación
- **Análisis Avanzado**: Herramientas de análisis específicas para cada tipo de simulación

## 📁 Estructura del Repositorio

Este repositorio contiene múltiples versiones y enfoques:

### 1. **Versión Bash** (`bash_version/`)
- Implementación original en bash
- Scripts de automatización para MD
- Diseñado para uso directo con GROMACS

### 2. **Versión Python** (`python_version/`)
- Nueva implementación en Python
- Interfaz de línea de comandos mejorada
- Características adicionales y documentación completa

### 3. **Proyecto GROMACS Avanzado** (directorio raíz)
- Estructura organizada por tipo de simulación
- Scripts especializados para cada enfoque
- Configuraciones optimizadas

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

## 🛠️ Tipos de Simulación Soportados

### 1. **Simulaciones Atómicas**
- **Uso**: Estudios detallados de interacciones moleculares
- **Campos de fuerza**: AMBER99SB-ILDN, CHARMM36, OPLS-AA
- **Script**: `scripts/run_simulation.sh`

### 2. **Simulaciones de Coarse Grain**
- **Uso**: Procesos de larga escala temporal
- **Campos de fuerza**: MARTINI, ELNEDIN, SIRAH
- **Script**: `scripts/cg/setup_cg_simulation.sh`

### 3. **Proteínas Transmembranales**
- **Uso**: Proteínas embebidas en membranas
- **Tipos de membrana**: POPC, DPPC, DOPC, DMPC
- **Script**: `scripts/transmembrane/setup_membrane.sh`

## 📋 Requisitos

- **GROMACS** (versión recomendada: 2023.x o superior)
- **Python 3.8+** (para scripts de análisis)
- **VMD o PyMOL** (para visualización)
- **NumPy, MDAnalysis, Rich** (para análisis avanzado)

## 🚀 Instalación

### Instalación del Paquete Python
```bash
pip install cha-md-ba
```

### Instalación Manual
```bash
git clone https://github.com/mixcoha/cha-MD-ba.git
cd cha-MD-ba
pip install -r requirements.txt
```

## 📖 Uso

### Simulaciones Atómicas
```bash
# 1. Coloca tu archivo PDB en input/
# 2. Ejecuta la simulación
./scripts/run_simulation.sh nombre_proteina

# 3. Analiza los resultados
python scripts/analyze_trajectory.py output/md.xtc output/md.gro
```

### Simulaciones Coarse Grain
```bash
# 1. Coloca tu archivo PDB en input/
# 2. Configura la simulación CG
./scripts/cg/setup_cg_simulation.sh proteina.pdb MARTINI

# 3. Analiza con scripts específicos de CG
```

### Proteínas Transmembranales
```bash
# 1. Coloca tu archivo PDB en input/
# 2. Configura la simulación de membrana
./scripts/transmembrane/setup_membrane.sh proteina.pdb POPC

# 3. Analiza la membrana
python scripts/transmembrane/analyze_membrane.py output/md.xtc output/md.gro
```

### Interfaz de Línea de Comandos (Python)
```bash
cha-md-ba prepare --input protein.pdb
cha-md-ba minimize
cha-md-ba nvt --force-constant 1000
cha-md-ba npt
cha-md-ba production
```

### API de Python
```python
from cha_md_ba import prepare, minimize, nvt, npt

# Preparar el sistema
preparator = prepare.SystemPreparator("protein.pdb")
preparator.prepare()

# Minimizar el sistema
minimizer = minimize.Minimizer()
minimizer.minimize()

# Equilibración NVT
nvt_equilibrator = nvt.NVTEquilibrator()
nvt_equilibrator.equilibrate(force_constant=1000)

# Equilibración NPT
npt_equilibrator = npt.NPTEquilibrator()
npt_equilibrator.equilibrate()
```

## 📚 Documentación

Documentación completa disponible en:
- **Versión Bash**: `bash_version/docs/`
- **Versión Python**: `python_version/docs/`
- **Guías de Uso**: `guides/`
- **API Reference**: `docs/api/`

## 🔧 Comandos Útiles

```bash
# Verificar versión de GROMACS
gmx --version

# Listar módulos disponibles
gmx help

# Verificar instalación del proyecto
./scripts/check_gmx.sh
```

## 🤝 Contribuir

¡Las contribuciones son bienvenidas! Por favor lee nuestras [Guías de Contribución](CONTRIBUTING.md) antes de enviar pull requests.

## 📄 Licencia

Este proyecto está licenciado bajo la Licencia MIT - ver el archivo [LICENSE](LICENSE) para detalles.

## 👨‍💻 Autor

**Edgar Mixcoha**

## 🙏 Agradecimientos

- Equipo de desarrollo de GROMACS
- Desarrolladores de MDAnalysis
- Todos los contribuyentes a este proyecto

## 📞 Contacto

Para preguntas y soporte, por favor abre un issue en el repositorio de GitHub.

## ⚠️ Notas Importantes

- Asegúrate de tener suficiente espacio en disco para las simulaciones
- Documenta todos los parámetros utilizados en cada simulación
- Haz respaldos regulares de tus datos importantes
- Verifica la compatibilidad de los campos de fuerza antes de usar
