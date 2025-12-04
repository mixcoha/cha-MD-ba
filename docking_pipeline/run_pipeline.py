#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script principal para ejecutar el pipeline completo de docking con AutoDock-Vina.
Este script integra la preparación de archivos, el docking y la visualización de resultados.
"""

import os
import sys
import argparse
import subprocess
import glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import logging

# Configuración de rutas
BASE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = BASE_DIR / "scripts"
RECEPTORS_DIR = BASE_DIR / "receptors"
LIGANDS_DIR = BASE_DIR / "ligands"
RESULTS_DIR = BASE_DIR / "results"

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Asegurar que los directorios existen
os.makedirs(RECEPTORS_DIR, exist_ok=True)
os.makedirs(LIGANDS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

def parse_arguments():
    """Parsear argumentos de línea de comandos."""
    parser = argparse.ArgumentParser(description="Pipeline completo de docking con AutoDock-Vina")
    
    # Opciones de preparación
    parser.add_argument("--prepare_receptor", type=str, help="Archivo PDB del receptor a preparar")
    parser.add_argument("--prepare_ligand", type=str, help="Archivo PDB del ligando a preparar")
    parser.add_argument("--prepare_all", action="store_true", help="Preparar todos los archivos PDB en los directorios")
    
    # Opciones de docking
    parser.add_argument("--receptor", type=str, help="Nombre del archivo del receptor (PDBQT) para docking")
    parser.add_argument("--ligand", type=str, help="Nombre del archivo del ligando (PDBQT) para docking")
    parser.add_argument("--center_x", type=float, help="Coordenada X del centro de la caja de búsqueda (Angstrom)")
    parser.add_argument("--center_y", type=float, help="Coordenada Y del centro de la caja de búsqueda (Angstrom)")
    parser.add_argument("--center_z", type=float, help="Coordenada Z del centro de la caja de búsqueda (Angstrom)")
    parser.add_argument("--size_x", type=float, default=20.0, help="Tamaño de la caja de búsqueda en dimensión X (Angstrom)")
    parser.add_argument("--size_y", type=float, default=20.0, help="Tamaño de la caja de búsqueda en dimensión Y (Angstrom)")
    parser.add_argument("--size_z", type=float, default=20.0, help="Tamaño de la caja de búsqueda en dimensión Z (Angstrom)")
    
    # Opciones de visualización
    parser.add_argument("--visualize", action="store_true", help="Visualizar los resultados del docking")
    parser.add_argument("--result_file", type=str, help="Archivo de resultado para visualizar")
    
    # Opciones de bases de datos
    parser.add_argument("--database", type=str, choices=["zinc15", "zinc20", "pubchem", "chembl", "drugbank", "enamine"],
                      help="Base de datos de la cual descargar compuestos")
    parser.add_argument("--query", type=str, help="Consulta de búsqueda para la base de datos")
    parser.add_argument("--max_compounds", type=int, default=1000,
                      help="Número máximo de compuestos a descargar de la base de datos")
    
    # Opciones de procesamiento paralelo
    parser.add_argument("--threads", type=int, default=4,
                      help="Número de hilos para procesamiento paralelo")
    
    return parser.parse_args()

def run_prepare_pdbqt(args):
    """Ejecutar el script de preparación de archivos PDBQT."""
    cmd = [sys.executable, str(SCRIPTS_DIR / "prepare_pdbqt.py")]
    
    if args.prepare_receptor:
        cmd.extend(["--receptor", args.prepare_receptor])
    
    if args.prepare_ligand:
        cmd.extend(["--ligand", args.prepare_ligand])
    
    if args.prepare_all:
        cmd.append("--all_receptors")
        cmd.append("--all_ligands")
    
    logger.info("Ejecutando preparación de archivos PDBQT...")
    try:
        subprocess.run(cmd, check=True)
        logger.info("Preparación de archivos completada.")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error en la preparación de archivos: {e}")
        return False

def run_docking(args, receptor=None, ligand=None):
    """Ejecutar el script de docking."""
    # Verificar que se proporcionaron los argumentos necesarios
    if not (receptor and ligand and args.center_x and args.center_y and args.center_z):
        logger.error("Error: Para ejecutar el docking, debes proporcionar receptor, ligando y las coordenadas del centro.")
        return False
    
    cmd = [sys.executable, str(SCRIPTS_DIR / "run_docking.py")]
    cmd.extend(["--receptor", receptor])
    cmd.extend(["--ligand", ligand])
    cmd.extend(["--center_x", str(args.center_x)])
    cmd.extend(["--center_y", str(args.center_y)])
    cmd.extend(["--center_z", str(args.center_z)])
    cmd.extend(["--size_x", str(args.size_x)])
    cmd.extend(["--size_y", str(args.size_y)])
    cmd.extend(["--size_z", str(args.size_z)])
    
    logger.info(f"Ejecutando docking para {ligand}...")
    try:
        subprocess.run(cmd, check=True)
        logger.info(f"Docking completado para {ligand}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error en el docking para {ligand}: {e}")
        return False

def run_visualization(args):
    """Ejecutar el script de visualización."""
    # Verificar que se proporcionaron los argumentos necesarios
    if not (args.receptor and args.result_file):
        logger.error("Error: Para visualizar los resultados, debes proporcionar --receptor y --result_file.")
        return False
    
    cmd = [sys.executable, str(SCRIPTS_DIR / "visualize_results.py")]
    cmd.extend(["--receptor", args.receptor])
    cmd.extend(["--result", args.result_file])
    
    logger.info("Generando visualización...")
    try:
        subprocess.run(cmd, check=True)
        logger.info("Visualización generada.")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error en la visualización: {e}")
        return False

def download_compounds(args):
    """Descargar compuestos de la base de datos especificada."""
    if not (args.database and args.query):
        logger.error("Error: Para descargar compuestos, debes proporcionar --database y --query.")
        return False
    
    cmd = [sys.executable, str(SCRIPTS_DIR / "download_compounds.py")]
    cmd.extend(["--database", args.database])
    cmd.extend(["--query", args.query])
    cmd.extend(["--max_compounds", str(args.max_compounds)])
    cmd.extend(["--threads", str(args.threads)])
    
    logger.info(f"Descargando compuestos de {args.database}...")
    try:
        subprocess.run(cmd, check=True)
        logger.info("Descarga de compuestos completada.")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error en la descarga de compuestos: {e}")
        return False

def run_mass_docking(args):
    """Ejecutar docking masivo con todos los ligandos disponibles."""
    if not args.receptor:
        logger.error("Error: Debes proporcionar un receptor para el docking masivo.")
        return False
    
    # Obtener lista de ligandos
    ligand_files = glob.glob(str(LIGANDS_DIR / "*.pdbqt"))
    if not ligand_files:
        logger.error("No se encontraron ligandos en el directorio.")
        return False
    
    logger.info(f"Iniciando docking masivo con {len(ligand_files)} ligandos...")
    
    # Crear función para procesar un ligando
    def process_ligand(ligand_path):
        ligand_name = os.path.basename(ligand_path)
        return run_docking(args, args.receptor, ligand_name)
    
    # Ejecutar docking en paralelo
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = list(executor.map(process_ligand, ligand_files))
    
    successful = sum(results)
    logger.info(f"Docking masivo completado. {successful} de {len(ligand_files)} ligandos procesados exitosamente.")
    return successful > 0

def main():
    """Función principal."""
    # Parsear argumentos
    args = parse_arguments()
    
    # Verificar que al menos una opción fue proporcionada
    if not (args.prepare_receptor or args.prepare_ligand or args.prepare_all or 
            args.receptor or args.visualize or args.database):
        logger.error("Error: Debes especificar al menos una acción para ejecutar.")
        logger.error("Usa --help para ver las opciones disponibles.")
        sys.exit(1)
    
    # Descargar compuestos si se solicitó
    if args.database:
        download_compounds(args)
    
    # Ejecutar preparación de archivos si se solicitó
    if args.prepare_receptor or args.prepare_ligand or args.prepare_all:
        run_prepare_pdbqt(args)
    
    # Ejecutar docking masivo si hay un receptor y ligandos disponibles
    if args.receptor and os.path.exists(LIGANDS_DIR):
        run_mass_docking(args)
    # Ejecutar docking individual si se proporcionaron los argumentos necesarios
    elif args.receptor and args.ligand and args.center_x and args.center_y and args.center_z:
        run_docking(args, args.receptor, args.ligand)
    
    # Ejecutar visualización si se solicitó
    if args.visualize and args.receptor and args.result_file:
        run_visualization(args)

if __name__ == "__main__":
    main() 