#!/bin/bash

# Script para configurar simulaciones de proteínas transmembranales
# Uso: ./setup_membrane.sh [proteina.pdb] [tipo_membrana]

set -e

if [ $# -ne 2 ]; then
    echo "Uso: $0 [proteina.pdb] [tipo_membrana]"
    echo "Tipos de membrana disponibles: POPC, DPPC, DOPC, DMPC"
    exit 1
fi

PROTEIN_PDB=$1
MEMBRANE_TYPE=$2
INPUT_DIR="input/transmembrane"
OUTPUT_DIR="output/transmembrane"

echo "=== Configurando simulación de proteína transmembranal ==="
echo "Proteína: $PROTEIN_PDB"
echo "Membrana: $MEMBRANE_TYPE"

# Crear directorios
mkdir -p "$OUTPUT_DIR"

# Verificar que existe la proteína
if [ ! -f "$PROTEIN_PDB" ]; then
    echo "Error: No se encontró $PROTEIN_PDB"
    exit 1
fi

echo "1. Orientando la proteína en la membrana..."
# Aquí irían comandos específicos para orientar la proteína
# Por ejemplo, usando orientconf o herramientas similares

echo "2. Generando topología de la proteína..."
gmx pdb2gmx -f "$PROTEIN_PDB" \
    -o "$OUTPUT_DIR/protein_processed.gro" \
    -water tip3p \
    -ff amber99sb-ildn \
    -ignh

echo "3. Creando membrana de $MEMBRANE_TYPE..."
# Comando para generar membrana (depende del tipo)
case $MEMBRANE_TYPE in
    "POPC")
        echo "Generando membrana POPC..."
        # gmx membrane -f popc.gro -n 128 -o membrane.gro
        ;;
    "DPPC")
        echo "Generando membrana DPPC..."
        # gmx membrane -f dppc.gro -n 128 -o membrane.gro
        ;;
    "DOPC")
        echo "Generando membrana DOPC..."
        # gmx membrane -f dopc.gro -n 128 -o membrane.gro
        ;;
    "DMPC")
        echo "Generando membrana DMPC..."
        # gmx membrane -f dmpc.gro -n 128 -o membrane.gro
        ;;
    *)
        echo "Tipo de membrana no reconocido: $MEMBRANE_TYPE"
        exit 1
        ;;
esac

echo "4. Insertando proteína en la membrana..."
# Comando para insertar proteína en membrana
# gmx insert-molecules -f membrane.gro -ci protein.gro -o system.gro

echo "5. Solvatando el sistema..."
# gmx solvate -cp system.gro -cs spc216.gro -o solvated.gro

echo "6. Añadiendo iones..."
# gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
# echo "13" | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral

echo "=== Configuración completada ==="
echo "Archivos generados en: $OUTPUT_DIR/"
echo ""
echo "Próximos pasos:"
echo "1. Revisar la orientación de la proteína en la membrana"
echo "2. Ajustar parámetros de simulación si es necesario"
echo "3. Ejecutar minimización y equilibración"
