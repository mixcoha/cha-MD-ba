import os
from pathlib import Path
from cha_md_ba.nvt import NVTEquilibrator

# Configuración de la simulación
input_gro = "simulations/1tim/2_minimization/minimized.gro"
topol_top = "simulations/1tim/1_preparation/topol.top"
output_dir = "simulations/1tim/3_nvt"

# Crear equilibrador NVT
equilibrator = NVTEquilibrator(
    input_gro=input_gro,
    topol_top=topol_top
)

# Ejecutar equilibración
results = equilibrator.equilibrate(
    output_dir=output_dir,
    gpu_ids="0"  # Usar GPU 0
)

# Verificar resultados
print("\n🎉 Equilibración NVT completa. Archivos generados:")
for fc, files in results.items():
    print(f"\nConstante de fuerza: {fc}")
    for key, path in files.items():
        print(f"{key}: {path}") 