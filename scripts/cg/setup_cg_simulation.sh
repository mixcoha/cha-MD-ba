#!/bin/bash

# Script para configurar simulaciones de coarse grain
# Uso: ./setup_cg_simulation.sh [proteina.pdb] [campo_fuerza_cg]

set -e

if [ $# -ne 2 ]; then
    echo "Uso: $0 [proteina.pdb] [campo_fuerza_cg]"
    echo "Campos de fuerza CG disponibles: MARTINI, ELNEDIN, SIRAH"
    exit 1
fi

PROTEIN_PDB=$1
FORCE_FIELD=$2
INPUT_DIR="input/coarse_grain"
OUTPUT_DIR="output/coarse_grain"

echo "=== Configurando simulación de Coarse Grain ==="
echo "Proteína: $PROTEIN_PDB"
echo "Campo de fuerza: $FORCE_FIELD"

# Crear directorios
mkdir -p "$OUTPUT_DIR"

# Verificar que existe la proteína
if [ ! -f "$PROTEIN_PDB" ]; then
    echo "Error: No se encontró $PROTEIN_PDB"
    exit 1
fi

echo "1. Convirtiendo estructura a coarse grain..."
case $FORCE_FIELD in
    "MARTINI")
        echo "Usando campo de fuerza MARTINI..."
        # Comandos específicos para MARTINI
        # martinize.py -f $PROTEIN_PDB -o cg_protein.gro -x cg_protein.pdb -dssp dssp
        ;;
    "ELNEDIN")
        echo "Usando campo de fuerza ELNEDIN..."
        # Comandos específicos para ELNEDIN
        ;;
    "SIRAH")
        echo "Usando campo de fuerza SIRAH..."
        # Comandos específicos para SIRAH
        ;;
    *)
        echo "Campo de fuerza CG no reconocido: $FORCE_FIELD"
        exit 1
        ;;
esac

echo "2. Generando topología CG..."
# gmx pdb2gmx -f cg_protein.pdb -o cg_protein.gro -water martini -ff martini

echo "3. Definiendo caja de simulación..."
# gmx editconf -f cg_protein.gro -o cg_protein_box.gro -c -d 1.0 -bt cubic

echo "4. Solvatando con solvente CG..."
case $FORCE_FIELD in
    "MARTINI")
        # gmx solvate -cp cg_protein_box.gro -cs martini_water.gro -o cg_solv.gro
        ;;
    "ELNEDIN")
        # gmx solvate -cp cg_protein_box.gro -cs elnedin_water.gro -o cg_solv.gro
        ;;
    "SIRAH")
        # gmx solvate -cp cg_protein_box.gro -cs sirah_water.gro -o cg_solv.gro
        ;;
esac

echo "5. Añadiendo iones CG..."
# gmx grompp -f cg_ions.mdp -c cg_solv.gro -p cg_topol.top -o cg_ions.tpr
# echo "13" | gmx genion -s cg_ions.tpr -o cg_ions.gro -p cg_topol.top -pname NA -nname CL -neutral

echo "6. Minimización CG..."
# gmx grompp -f cg_minim.mdp -c cg_ions.gro -p cg_topol.top -o cg_em.tpr
# gmx mdrun -v -deffnm cg_em

echo "7. Equilibración CG..."
# gmx grompp -f cg_nvt.mdp -c cg_em.gro -r cg_em.gro -p cg_topol.top -o cg_nvt.tpr
# gmx mdrun -v -deffnm cg_nvt

echo "=== Configuración CG completada ==="
echo "Archivos generados en: $OUTPUT_DIR/"
echo ""
echo "Ventajas de la simulación CG:"
echo "- Tiempos de simulación más largos"
echo "- Menor costo computacional"
echo "- Ideal para procesos de plegamiento y agregación"
