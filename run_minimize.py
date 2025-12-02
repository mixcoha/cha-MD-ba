#!/usr/bin/env python3
"""
Script directo para minimizar el receptor GABA usando el módulo minimize.py
"""

import os
import sys
import subprocess
from pathlib import Path

# Agregar el directorio donde está el módulo minimize.py
sys.path.insert(0, '/Users/mixcoha/GMX/cha_md_ba')

from minimize import EnergyMinimizer

def main():
    print("🧬 Iniciando minimización del receptor GABA")
    
    # Configurar rutas
    input_gro = "processed.gro"
    topol_top = "topol.top"
    output_dir = "em"
    gmx_cmd = "/usr/local/gromacs/bin/gmx_mpi"
    
    # Verificar que los archivos existen
    if not os.path.exists(input_gro):
        print(f"❌ Error: No se encuentra {input_gro}")
        return 1
    
    if not os.path.exists(topol_top):
        print(f"❌ Error: No se encuentra {topol_top}")
        return 1
    
    print(f"📁 Archivos de entrada:")
    print(f"  • Coordenadas: {input_gro}")
    print(f"  • Topología: {topol_top}")
    print(f"  • Directorio de salida: {output_dir}")
    
    # Crear minimizador
    minimizer = EnergyMinimizer(
        input_gro=input_gro,
        topol_top=topol_top,
        gmx=gmx_cmd
    )
    
    try:
        # Ejecutar minimización
        print("⚡ Ejecutando minimización de energía...")
        results = minimizer.minimize(output_dir=output_dir)
        
        print("✅ Minimización completada exitosamente!")
        print("📊 Archivos generados:")
        for key, path in results.items():
            print(f"  • {key}: {path}")
            
        print("💡 Próximos pasos sugeridos:")
        print("  1. Revisar el archivo de log para verificar convergencia")
        print("  2. Analizar la energía con: gmx energy -f em/ener.edr")
        print("  3. Proceder con equilibración NVT")
        
    except Exception as e:
        print(f"❌ Error durante la minimización: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
