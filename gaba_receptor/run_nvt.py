#!/usr/bin/env python3
"""
Script para equilibración NVT del receptor GABA usando CHA-MD-BA
"""

import os
import sys
import subprocess
from pathlib import Path

# Agregar el directorio donde está el módulo nvt.py
sys.path.insert(0, '/Users/mixcoha/GMX/cha_md_ba')

from nvt import NVTEquilibrator

def main():
    print("🌡️ Iniciando equilibración NVT del receptor GABA")
    
    # Configurar rutas - usar el archivo minimizado como entrada
    input_gro = "em/minim.gro"
    topol_top = "topol.top"
    output_dir = "nvt"
    
    # Verificar que los archivos existen
    if not os.path.exists(input_gro):
        print(f"❌ Error: No se encuentra {input_gro}")
        return 1
    
    if not os.path.exists(topol_top):
        print(f"❌ Error: No se encuentra {topol_top}")
        return 1
    
    print(f"📁 Archivos de entrada:")
    print(f"  • Estructura minimizada: {input_gro}")
    print(f"  • Topología: {topol_top}")
    print(f"  • Directorio de salida: {output_dir}")
    
    # Crear equilibrador NVT
    nvt_equilibrator = NVTEquilibrator(
        input_gro=input_gro,
        topol_top=topol_top
    )
    
    try:
        # Ejecutar equilibración NVT
        print("🚀 Ejecutando equilibración NVT...")
        results = nvt_equilibrator.equilibrate(
            output_dir=output_dir,
            restraint_constant=1000,  # constante de restricción inicial
            simulation_time_ns=0.1   # 100 ps de equilibración
        )
        
        print("✅ Equilibración NVT completada exitosamente!")
        print("📊 Archivos generados:")
        for key, path in results.items():
            if path and os.path.exists(path):
                print(f"  • {key}: {path}")
        
        print("💡 Próximos pasos:")
        print("  1. Revisar la temperatura: gmx energy -f nvt/nvt.edr -o temperature.xvg")
        print("  2. Proceder con equilibración NPT")
        
    except Exception as e:
        print(f"❌ Error durante la equilibración NVT: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
