#!/bin/bash

# Script de ejemplo para ejecutar una simulación básica de GROMACS
# Uso: ./run_simulation.sh [nombre_del_sistema]

set -e  # Salir si hay algún error

# Verificar argumentos
if [ $# -eq 0 ]; then
    echo "Uso: $0 [nombre_del_sistema]"
    echo "Ejemplo: $0 protein_water"
    exit 1
fi

SYSTEM_NAME=$1
INPUT_DIR="input"
OUTPUT_DIR="output"

echo "=== Iniciando simulación para: $SYSTEM_NAME ==="

# Verificar que existe el archivo de estructura
if [ ! -f "$INPUT_DIR/${SYSTEM_NAME}.pdb" ]; then
    echo "Error: No se encontró $INPUT_DIR/${SYSTEM_NAME}.pdb"
    echo "Por favor, coloca tu archivo PDB en el directorio input/"
    exit 1
fi

# Crear directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

echo "1. Generando topología..."
gmx pdb2gmx -f "$INPUT_DIR/${SYSTEM_NAME}.pdb" \
    -o "$OUTPUT_DIR/${SYSTEM_NAME}_processed.gro" \
    -water tip3p \
    -ff amber99sb-ildn \
    -ignh

echo "2. Definiendo caja de simulación..."
gmx editconf -f "$OUTPUT_DIR/${SYSTEM_NAME}_processed.gro" \
    -o "$OUTPUT_DIR/${SYSTEM_NAME}_box.gro" \
    -c -d 1.0 -bt cubic

echo "3. Solvatando el sistema..."
gmx solvate -cp "$OUTPUT_DIR/${SYSTEM_NAME}_box.gro" \
    -cs spc216.gro \
    -o "$OUTPUT_DIR/${SYSTEM_NAME}_solv.gro" \
    -p topol.top

echo "4. Añadiendo iones..."
gmx grompp -f "$INPUT_DIR/ions.mdp" \
    -c "$OUTPUT_DIR/${SYSTEM_NAME}_solv.gro" \
    -p topol.top \
    -o "$OUTPUT_DIR/ions.tpr"

echo "5. Neutralizando el sistema..."
echo "13" | gmx genion -s "$OUTPUT_DIR/ions.tpr" \
    -o "$OUTPUT_DIR/${SYSTEM_NAME}_ions.gro" \
    -p topol.top \
    -pname NA -nname CL -neutral

echo "6. Minimización de energía..."
gmx grompp -f "$INPUT_DIR/minim.mdp" \
    -c "$OUTPUT_DIR/${SYSTEM_NAME}_ions.gro" \
    -p topol.top \
    -o "$OUTPUT_DIR/em.tpr"

gmx mdrun -v -deffnm "$OUTPUT_DIR/em"

echo "7. Equilibración NVT..."
gmx grompp -f "$INPUT_DIR/nvt.mdp" \
    -c "$OUTPUT_DIR/em.gro" \
    -r "$OUTPUT_DIR/em.gro" \
    -p topol.top \
    -o "$OUTPUT_DIR/nvt.tpr"

gmx mdrun -v -deffnm "$OUTPUT_DIR/nvt"

echo "8. Equilibración NPT..."
gmx grompp -f "$INPUT_DIR/npt.mdp" \
    -c "$OUTPUT_DIR/nvt.gro" \
    -r "$OUTPUT_DIR/nvt.gro" \
    -t "$OUTPUT_DIR/nvt.cpt" \
    -p topol.top \
    -o "$OUTPUT_DIR/npt.tpr"

gmx mdrun -v -deffnm "$OUTPUT_DIR/npt"

echo "9. Simulación de producción..."
gmx grompp -f "$INPUT_DIR/md.mdp" \
    -c "$OUTPUT_DIR/npt.gro" \
    -r "$OUTPUT_DIR/npt.gro" \
    -t "$OUTPUT_DIR/npt.cpt" \
    -p topol.top \
    -o "$OUTPUT_DIR/md.tpr"

gmx mdrun -v -deffnm "$OUTPUT_DIR/md"

echo "=== Simulación completada ==="
echo "Los archivos de salida están en: $OUTPUT_DIR/"
echo "Archivo de trayectoria: $OUTPUT_DIR/md.xtc"
echo "Archivo de estructura: $OUTPUT_DIR/md.gro"
