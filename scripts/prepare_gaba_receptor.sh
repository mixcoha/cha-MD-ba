#!/bin/bash

# Script para preparar el receptor GABA para minimización usando CHA-MD-BA
# Este script usa las configuraciones MDP optimizadas existentes

PROTEIN_NAME="gaba_receptor"
PDB_FILE="gaba_receptor_protonated_pH7.4.pdb"

echo "🧬 Preparando receptor GABA para simulación MD"
echo "================================================"

# Verificar que existe el archivo PDB
if [ ! -f "$PDB_FILE" ]; then
    echo "❌ Error: No se encontró el archivo $PDB_FILE"
    echo "Por favor, asegúrate de que el archivo esté en el directorio actual"
    exit 1
fi

# Crear estructura de directorios
echo "📁 Creando estructura de directorios..."
mkdir -p ${PROTEIN_NAME}/{setup,em,nvt,npt,md,analysis}

# Copiar archivo PDB al directorio de setup
cp "$PDB_FILE" "${PROTEIN_NAME}/setup/"

echo "✅ Estructura de directorios creada:"
echo "   📂 ${PROTEIN_NAME}/"
echo "   ├── 📂 setup/     # Preparación inicial"
echo "   ├── 📂 em/        # Minimización de energía"
echo "   ├── 📂 nvt/       # Equilibración NVT"
echo "   ├── 📂 npt/       # Equilibración NPT"
echo "   ├── 📂 md/        # Simulación de producción"
echo "   └── 📂 analysis/  # Análisis de resultados"

echo ""
echo "🔧 Comandos sugeridos para continuar (requiere GROMACS):"
echo ""
echo "1. Generar topología:"
echo "   cd ${PROTEIN_NAME}/setup"
echo "   gmx pdb2gmx -f $PDB_FILE -o processed.gro -water tip3p -ff amber99sb-ildn"
echo ""
echo "2. Usar scripts CHA-MD-BA optimizados:"
echo "   ./bash_version/minimiza_atomico.sh $PROTEIN_NAME"
echo "   ./bash_version/nvt2.sh $PROTEIN_NAME"
echo "   ./bash_version/npt.sh $PROTEIN_NAME"
echo ""
echo "3. Análisis con scripts optimizados:"
echo "   ./bash_version/analisis.sh $PROTEIN_NAME"
echo ""

# Crear archivo de información del sistema
cat > "${PROTEIN_NAME}/system_info.txt" << EOF
Sistema: Receptor GABA
Archivo original: $PDB_FILE
Protonación: pH 7.4 (PROPKA 3.5.1)
Generado con: MODELLER 10.7
Fecha de preparación: $(date)

Átomos totales: $(grep "^ATOM" "$PDB_FILE" | wc -l)
Residuos únicos: $(grep "^ATOM" "$PDB_FILE" | awk '{print $4}' | sort | uniq | wc -l)

Campo de fuerza recomendado: AMBER99SB-ILDN
Modelo de agua recomendado: TIP3P
Tipo de caja recomendado: Dodecahedron

Notas:
- Sistema protonado a pH fisiológico
- Usar scripts CHA-MD-BA con MDP optimizados
- Monitorear estabilidad durante equilibración
EOF

echo "📋 Información del sistema guardada en: ${PROTEIN_NAME}/system_info.txt"
echo ""
echo "🎯 Preparación completada. El sistema está listo para minimización."
