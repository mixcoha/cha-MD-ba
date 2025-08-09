# Guía de Integración CHA-MD-BA

## 🎯 Integración Completada

Los scripts más recientes de `python_version/cha_md_ba/` han sido integrados exitosamente al proyecto principal, manteniendo intactos los archivos MDP optimizados de los scripts bash.

## 📁 Estructura Integrada

```
GMX/
├── cha_md_ba/                    # Módulo Python principal (ACTUALIZADO)
│   ├── prepare.py               # Preparación avanzada de sistemas
│   ├── minimize.py              # Minimización con soporte GPU
│   ├── nvt.py                   # Equilibración NVT avanzada
│   ├── npt.py                   # Equilibración NPT avanzada
│   ├── analysis.py              # Análisis avanzado con MDAnalysis
│   ├── preprocess.py            # Preprocesamiento de trayectorias
│   └── cli.py                   # Interfaz de línea de comandos
├── bash_version/                # Scripts bash optimizados (PRESERVADOS)
│   ├── analisis.sh             # Análisis con MDP optimizados
│   ├── nvt2.sh                 # NVT con restricciones graduales
│   ├── npt.sh                  # NPT optimizado
│   └── minimiza_atomico.sh     # Minimización atómica
├── scripts/                     # Scripts de integración (NUEVOS)
│   ├── cha_md_ba_wrapper.py    # Wrapper que combina bash + Python
│   └── run_advanced_simulation.py # Script Python puro
└── input/, output/, analysis/   # Estructura original
```

## 🚀 Uso de los Scripts Integrados

### Opción 1: Wrapper Híbrido (Recomendado)

El wrapper combina lo mejor de ambos mundos: **scripts bash con MDP optimizados + análisis Python avanzado**.

```bash
# Pipeline completo
python scripts/cha_md_ba_wrapper.py full_pipeline mi_proteina

# Comandos individuales
python scripts/cha_md_ba_wrapper.py minimize mi_proteina
python scripts/cha_md_ba_wrapper.py nvt mi_proteina
python scripts/cha_md_ba_wrapper.py npt mi_proteina
python scripts/cha_md_ba_wrapper.py analyze_bash mi_proteina
python scripts/cha_md_ba_wrapper.py analyze_advanced mi_proteina
```

### Opción 2: Scripts Python Puros

Para nuevos desarrollos y funcionalidades avanzadas:

```bash
# Preparación de sistema
python scripts/run_advanced_simulation.py prepare proteina.pdb output/

# Minimización
python scripts/run_advanced_simulation.py minimize system.gro topol.top output/

# Equilibración
python scripts/run_advanced_simulation.py nvt em.gro topol.top output/
python scripts/run_advanced_simulation.py npt nvt.gro topol.top output/

# Análisis
python scripts/run_advanced_simulation.py analyze md.xtc md.tpr output/
```

### Opción 3: Scripts Bash Originales

Para máximo control y configuraciones específicas:

```bash
# Usar directamente los scripts bash
./bash_version/minimiza_atomico.sh mi_proteina
./bash_version/nvt2.sh mi_proteina
./bash_version/npt.sh mi_proteina
./bash_version/analisis.sh mi_proteina
```

## 💡 Ventajas de la Integración

### Scripts Bash (Preservados)
- ✅ **Archivos MDP optimizados** - Mantenidos exactamente como estaban
- ✅ **Configuraciones probadas** - Scripts validados en producción
- ✅ **Compatibilidad total** - Funciona exactamente igual que antes

### Funcionalidad Python (Nueva)
- ✅ **Análisis avanzado** - MDAnalysis, matplotlib, rich console
- ✅ **Soporte GPU** - Configuración automática de GPUs
- ✅ **Interfaces modernas** - Click CLI, progress bars
- ✅ **Preprocesamiento** - Limpieza automática de trayectorias

### Wrapper Híbrido (Nuevo)
- ✅ **Lo mejor de ambos** - MDP optimizados + análisis moderno
- ✅ **Pipeline automatizado** - Ejecución secuencial completa
- ✅ **Manejo de errores** - Detección y reporte de problemas
- ✅ **Output colorido** - Rich console para mejor UX

## 🔧 Configuración y Dependencias

### Instalar Dependencias
```bash
pip install -r requirements.txt
```

### Verificar Instalación
```bash
./scripts/check_gmx.sh
python -c "import cha_md_ba; print('✅ CHA-MD-BA integrado correctamente')"
```

## 📊 Análisis Disponible

### Análisis Bash (Optimizado)
- RMSD con estadísticas y desviación estándar
- Radio de giro con gnuplot optimizado
- Gráficas postscript de alta calidad
- Configuraciones de plotting específicas

### Análisis Python (Avanzado)
- RMSD con MDAnalysis
- Radio de giro avanzado
- Análisis de estructura secundaria
- Visualización moderna con matplotlib
- Preprocesamiento de trayectorias

## 🎯 Casos de Uso Recomendados

### Para Producción
**Usar wrapper híbrido**: Combina la confiabilidad de los MDP optimizados con análisis moderno.

```bash
python scripts/cha_md_ba_wrapper.py full_pipeline mi_sistema
```

### Para Desarrollo
**Usar scripts Python**: Permite personalización y nuevas funcionalidades.

```bash
python scripts/run_advanced_simulation.py prepare proteina.pdb output/ --gpu-ids "01"
```

### Para Casos Específicos
**Usar scripts bash directamente**: Máximo control sobre parámetros.

```bash
./bash_version/nvt2.sh mi_sistema  # Restricciones graduales específicas
```

## 🔄 Migración

Los scripts existentes siguen funcionando exactamente igual. La integración añade funcionalidad sin romper compatibilidad.

## 📞 Soporte

- Scripts bash: Funcionalidad original preservada
- Scripts Python: Nuevas funcionalidades con documentación integrada
- Wrapper: Combina ambos enfoques de manera transparente
