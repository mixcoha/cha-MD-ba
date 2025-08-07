# Configuración del Proyecto GROMACS

## Parámetros de Simulación

### Parámetros de Integración
- **dt**: 0.002 ps (2 fs)
- **nsteps**: 50000 (para 100 ps de simulación)
- **nstxout**: 1000 (guardar coordenadas cada 2 ps)
- **nstvout**: 1000 (guardar velocidades cada 2 ps)
- **nstenergy**: 1000 (guardar energía cada 2 ps)
- **nstlog**: 1000 (guardar log cada 2 ps)

### Parámetros de Control de Temperatura
- **tcoupl**: v-rescale
- **tau_t**: 0.1 ps
- **ref_t**: 300 K

### Parámetros de Control de Presión
- **pcoupl**: Parrinello-Rahman
- **tau_p**: 2.0 ps
- **ref_p**: 1.0 bar
- **compressibility**: 4.5e-5 bar^-1

### Parámetros de Restricciones
- **constraints**: h-bonds
- **constraint_algorithm**: LINCS
- **lincs_iter**: 1
- **lincs_order**: 4

### Parámetros de Vecinos
- **cutoff-scheme**: Verlet
- **rlist**: 1.0 nm
- **rvdw**: 1.0 nm
- **rcoulomb**: 1.0 nm

## Fuerzas de Campo

### Campos de Fuerza Recomendados
- **Proteínas**: AMBER99SB-ILDN
- **Lípidos**: CHARMM36
- **Carbohidratos**: GLYCAM06
- **ADN/ARN**: AMBER99SB-ILDN

### Agua
- **TIP3P**: Para AMBER
- **TIP4P**: Para CHARMM
- **SPC/E**: Para OPLS

## Configuración de Análisis

### Análisis de Estructura
- **RMSD**: Backbone y Cα
- **RMSF**: Residuos individuales
- **Rg**: Radio de giro
- **SASA**: Área superficial accesible al solvente

### Análisis de Dinámica
- **MSD**: Coeficiente de difusión
- **VACF**: Función de autocorrelación de velocidades
- **PCA**: Análisis de componentes principales

## Notas de Configuración

- Ajusta los parámetros según tu sistema específico
- Documenta cualquier modificación a estos parámetros
- Verifica la compatibilidad de los campos de fuerza
- Considera el tiempo de simulación necesario para tu análisis
