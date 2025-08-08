import os
from pathlib import Path
from cha_md_ba.prepare import download_pdb, MDSystemPreparator

# Paso 1: Descargar el archivo PDB (si aún no lo tienes)
pdb_id = "1tim"
pdb_path = f"data/{pdb_id}.pdb"

# Asegúrate de que el directorio 'data/' existe
os.makedirs("data", exist_ok=True)

# Descargar desde RCSB
if not os.path.exists(pdb_path):
    success = download_pdb(pdb_id, pdb_path)
    if not success:
        raise RuntimeError(f"No se pudo descargar {pdb_id} desde RCSB.")
else:
    print(f"Archivo {pdb_path} ya existe.")

# Paso 2: Preparar el sistema con GROMACS
preparator = MDSystemPreparator(
    pdb_path=pdb_path,
    forcefield="amber99sb-ildn",
    water_model="tip3p"
)

# Lanzar preparación
result = preparator.prepare_system(
    output_dir="simulations",
    box_type="dodecahedron",
    box_size=1.0,  # Tamaño de la caja en nm
    ions=True,     # Agregar iones para neutralizar
    minimize=True, # Realizar minimización de energía
    gpu_ids="0",   # Usar GPU 0 para la minimización
    replace_group="Backbone"  # Reemplazar backbone con iones
)

# Verificar resultados
print("\n🎉 Preparación completa. Archivos generados:")
print(f"Directorio base: {result['base_dir']}")
print("\nDirectorios creados:")
print(f"1. Preparación: {result['preparation_dir']}")
print(f"2. Minimización: {result['minimization_dir']}")
print(f"3. NVT: {result['nvt_dir']}")
print(f"4. NPT: {result['npt_dir']}")

print("\nArchivos generados:")
for key, path in result.items():
    if key not in ['base_dir', 'preparation_dir', 'minimization_dir', 'nvt_dir', 'npt_dir']:
        print(f"{key}: {path}")
