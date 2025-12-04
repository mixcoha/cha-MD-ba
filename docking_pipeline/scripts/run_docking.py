#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script para automatizar el proceso de docking con AutoDock-Vina.
Este script maneja la preparación de archivos, ejecución de AutoDock-Vina y análisis de resultados.
"""

import os
import sys
import subprocess
import argparse
import glob
import shutil
from pathlib import Path

# Configuración de rutas
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RECEPTORS_DIR = BASE_DIR / "receptors"
LIGANDS_DIR = BASE_DIR / "ligands"
RESULTS_DIR = BASE_DIR / "results"
VINA_EXECUTABLE = Path("/Users/mixcoha/cha-MD-ba/AutoDock-Vina/build/src/main/vina")

def parse_arguments():
    """Parsear argumentos de línea de comandos."""
    parser = argparse.ArgumentParser(description="Pipeline de docking con AutoDock-Vina")
    
    parser.add_argument("--receptor", type=str, required=True,
                        help="Nombre del archivo del receptor (PDBQT) en el directorio receptors/")
    
    parser.add_argument("--ligand", type=str, required=True,
                        help="Nombre del archivo del ligando (PDBQT) en el directorio ligands/")
    
    parser.add_argument("--center_x", type=float, required=True,
                        help="Coordenada X del centro de la caja de búsqueda (Angstrom)")
    
    parser.add_argument("--center_y", type=float, required=True,
                        help="Coordenada Y del centro de la caja de búsqueda (Angstrom)")
    
    parser.add_argument("--center_z", type=float, required=True,
                        help="Coordenada Z del centro de la caja de búsqueda (Angstrom)")
    
    parser.add_argument("--size_x", type=float, default=20.0,
                        help="Tamaño de la caja de búsqueda en dimensión X (Angstrom)")
    
    parser.add_argument("--size_y", type=float, default=20.0,
                        help="Tamaño de la caja de búsqueda en dimensión Y (Angstrom)")
    
    parser.add_argument("--size_z", type=float, default=20.0,
                        help="Tamaño de la caja de búsqueda en dimensión Z (Angstrom)")
    
    parser.add_argument("--exhaustiveness", type=int, default=8,
                        help="Exhaustividad de la búsqueda global (1+)")
    
    parser.add_argument("--num_modes", type=int, default=9,
                        help="Número máximo de modos de unión a generar")
    
    parser.add_argument("--energy_range", type=float, default=3.0,
                        help="Diferencia máxima de energía entre el mejor modo y el peor (kcal/mol)")
    
    parser.add_argument("--cpu", type=int, default=0,
                        help="Número de CPUs a utilizar (0 = automático)")
    
    parser.add_argument("--seed", type=int, default=0,
                        help="Semilla aleatoria explícita (0 = aleatoria)")
    
    parser.add_argument("--scoring", type=str, default="vina",
                        choices=["vina", "vinardo", "ad4"],
                        help="Función de puntuación a utilizar")
    
    return parser.parse_args()

def check_files(args):
    """Verificar que los archivos necesarios existen."""
    receptor_path = RECEPTORS_DIR / args.receptor
    ligand_path = LIGANDS_DIR / args.ligand
    
    if not receptor_path.exists():
        print(f"Error: El archivo del receptor {receptor_path} no existe.")
        sys.exit(1)
    
    if not ligand_path.exists():
        print(f"Error: El archivo del ligando {ligand_path} no existe.")
        sys.exit(1)
    
    return receptor_path, ligand_path

def prepare_output_dir(args):
    """Preparar directorio de salida para los resultados."""
    # Crear nombre de directorio basado en el receptor y ligando
    receptor_name = os.path.splitext(args.receptor)[0]
    ligand_name = os.path.splitext(args.ligand)[0]
    output_dir = RESULTS_DIR / f"{receptor_name}_{ligand_name}"
    
    # Crear directorio si no existe
    os.makedirs(output_dir, exist_ok=True)
    
    return output_dir

def run_docking(args, receptor_path, ligand_path, output_dir):
    """Ejecutar AutoDock-Vina con los parámetros especificados."""
    # Construir comando
    cmd = [
        str(VINA_EXECUTABLE),
        "--receptor", str(receptor_path),
        "--ligand", str(ligand_path),
        "--out", str(output_dir / "result.pdbqt"),
        "--log", str(output_dir / "docking.log"),
        "--center_x", str(args.center_x),
        "--center_y", str(args.center_y),
        "--center_z", str(args.center_z),
        "--size_x", str(args.size_x),
        "--size_y", str(args.size_y),
        "--size_z", str(args.size_z),
        "--exhaustiveness", str(args.exhaustiveness),
        "--num_modes", str(args.num_modes),
        "--energy_range", str(args.energy_range),
        "--cpu", str(args.cpu),
        "--seed", str(args.seed),
        "--scoring", args.scoring
    ]
    
    # Ejecutar comando
    print(f"Ejecutando: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        
        # Guardar salida estándar y error en archivos
        with open(output_dir / "stdout.log", "w") as f:
            f.write(result.stdout)
        
        with open(output_dir / "stderr.log", "w") as f:
            f.write(result.stderr)
            
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar AutoDock-Vina: {e}")
        print(f"Salida de error: {e.stderr}")
        return False

def analyze_results(output_dir):
    """Analizar los resultados del docking."""
    log_file = output_dir / "docking.log"
    
    if not log_file.exists():
        print("No se encontró el archivo de log para análisis.")
        return
    
    # Extraer información relevante del log
    try:
        with open(log_file, "r") as f:
            log_content = f.read()
        
        # Buscar líneas con información de energía
        energy_lines = [line for line in log_content.split("\n") if "-----" in line and "kcal/mol" in line]
        
        if energy_lines:
            print("\nResultados del docking:")
            print("=" * 50)
            for line in energy_lines:
                print(line)
        else:
            print("No se encontraron resultados de energía en el log.")
    except Exception as e:
        print(f"Error al analizar resultados: {e}")

def main():
    """Función principal."""
    # Parsear argumentos
    args = parse_arguments()
    
    # Verificar archivos
    receptor_path, ligand_path = check_files(args)
    
    # Preparar directorio de salida
    output_dir = prepare_output_dir(args)
    
    # Ejecutar docking
    success = run_docking(args, receptor_path, ligand_path, output_dir)
    
    # Analizar resultados si el docking fue exitoso
    if success:
        analyze_results(output_dir)
        print(f"\nResultados guardados en: {output_dir}")
    else:
        print("El proceso de docking falló.")

if __name__ == "__main__":
    main() 