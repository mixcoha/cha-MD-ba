#!/bin/bash

echo "=== Verificación de GROMACS ==="
echo ""

# Verificar si GROMACS está instalado
GMX_CMD=""
if command -v gmx &> /dev/null; then
    GMX_CMD="gmx"
elif command -v /usr/local/gromacs/bin/gmx_mpi &> /dev/null; then
    GMX_CMD="/usr/local/gromacs/bin/gmx_mpi"
elif command -v gmx_mpi &> /dev/null; then
    GMX_CMD="gmx_mpi"
fi

if [ -n "$GMX_CMD" ]; then
    echo "✓ GROMACS está instalado"
    echo "Comando: $GMX_CMD"
    echo "Versión: $($GMX_CMD --version | head -n 1)"
else
    echo "✗ GROMACS no está instalado o no está en el PATH"
    echo "Por favor, instala GROMACS antes de continuar"
    exit 1
fi

echo ""

# Verificar módulos disponibles
echo "=== Módulos de GROMACS disponibles ==="
gmx help 2>/dev/null | grep -E "^[[:space:]]*[a-zA-Z]" | head -20

echo ""
echo "=== Verificación del entorno ==="

# Verificar variables de entorno importantes
if [ -n "$GMXDATA" ]; then
    echo "✓ GMXDATA está configurado: $GMXDATA"
else
    echo "⚠ GMXDATA no está configurado"
fi

if [ -n "$GMXBIN" ]; then
    echo "✓ GMXBIN está configurado: $GMXBIN"
else
    echo "⚠ GMXBIN no está configurado"
fi

echo ""
echo "=== Verificación completada ==="
