#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script para visualizar los resultados del docking.
Este script utiliza NGLView para visualizar las estructuras en un navegador web.
"""

import os
import sys
import argparse
import glob
from pathlib import Path
import subprocess

# Configuración de rutas
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RECEPTORS_DIR = BASE_DIR / "receptors"
LIGANDS_DIR = BASE_DIR / "ligands"
RESULTS_DIR = BASE_DIR / "results"

def parse_arguments():
    """Parsear argumentos de línea de comandos."""
    parser = argparse.ArgumentParser(description="Visualizar resultados de docking con AutoDock-Vina")
    
    parser.add_argument("--receptor", type=str, required=True,
                        help="Nombre del archivo del receptor (PDBQT) en el directorio receptors/")
    
    parser.add_argument("--result", type=str, required=True,
                        help="Nombre del archivo de resultado (PDBQT) en el directorio results/")
    
    parser.add_argument("--output", type=str, default="visualization.html",
                        help="Nombre del archivo HTML de salida")
    
    return parser.parse_args()

def check_files(args):
    """Verificar que los archivos necesarios existen."""
    receptor_path = RECEPTORS_DIR / args.receptor
    result_path = RESULTS_DIR / args.result
    
    if not receptor_path.exists():
        print(f"Error: El archivo del receptor {receptor_path} no existe.")
        sys.exit(1)
    
    if not result_path.exists():
        print(f"Error: El archivo de resultado {result_path} no existe.")
        sys.exit(1)
    
    return receptor_path, result_path

def create_visualization_script(receptor_path, result_path, output_path):
    """Crear un script Python para visualizar los resultados con NGLView."""
    script_content = f"""
import nglview as nv
import MDAnalysis as mda

# Cargar estructuras
receptor = mda.Universe("{receptor_path}")
ligand = mda.Universe("{result_path}")

# Crear vista
view = nv.show_mdanalysis(receptor)
view.add_component(ligand, selection="all")

# Configurar vista
view.clear()
view.add_cartoon(selection="protein", color="lightblue")
view.add_ball_and_stick(selection="ligand", color="red")

# Guardar visualización
view.save("{output_path}")
print(f"Visualización guardada en {output_path}")
"""
    
    script_path = RESULTS_DIR / "temp_visualization.py"
    with open(script_path, "w") as f:
        f.write(script_content)
    
    return script_path

def run_visualization(script_path):
    """Ejecutar el script de visualización."""
    try:
        result = subprocess.run(["python", str(script_path)], check=True, capture_output=True, text=True)
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar el script de visualización: {e}")
        print(f"Salida de error: {e.stderr}")
        return False
    finally:
        # Eliminar script temporal
        if os.path.exists(script_path):
            os.remove(script_path)

def main():
    """Función principal."""
    # Parsear argumentos
    args = parse_arguments()
    
    # Verificar archivos
    receptor_path, result_path = check_files(args)
    
    # Crear script de visualización
    script_path = create_visualization_script(receptor_path, result_path, args.output)
    
    # Ejecutar visualización
    success = run_visualization(script_path)
    
    if success:
        print(f"Visualización creada exitosamente: {args.output}")
        print(f"Puedes abrir el archivo en tu navegador web para ver los resultados.")
    else:
        print("La visualización falló.")

if __name__ == "__main__":
    main() 