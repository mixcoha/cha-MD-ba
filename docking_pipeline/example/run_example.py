#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script para ejecutar un ejemplo completo del pipeline de docking.
Este ejemplo utiliza la proteína 1hsg (HIV-1 protease) y busca compuestos en ZINC15.
"""

import os
import sys
import yaml
import subprocess
import logging
from pathlib import Path

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Obtener el directorio base del pipeline
PIPELINE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EXAMPLE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))

def load_config():
    """Cargar configuración del archivo YAML."""
    config_path = EXAMPLE_DIR / "config.yaml"
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error al cargar el archivo de configuración: {e}")
        sys.exit(1)

def download_receptor(config):
    """Descargar el archivo PDB del receptor."""
    receptor_name = config['receptor']['name']
    receptor_path = PIPELINE_DIR / "receptors" / f"{receptor_name}.pdb"
    
    if not receptor_path.exists():
        logger.info(f"Descargando receptor {receptor_name}...")
        try:
            cmd = [
                "curl",
                "-o", str(receptor_path),
                f"https://files.rcsb.org/download/{receptor_name}.pdb"
            ]
            subprocess.run(cmd, check=True)
            logger.info(f"Receptor descargado exitosamente: {receptor_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error al descargar el receptor: {e}")
            sys.exit(1)
    
    return receptor_path

def prepare_receptor(receptor_path):
    """Preparar el receptor para docking."""
    logger.info("Preparando receptor...")
    try:
        cmd = [
            sys.executable,
            str(PIPELINE_DIR / "scripts" / "prepare_pdbqt.py"),
            "--receptor", str(receptor_path)
        ]
        subprocess.run(cmd, check=True)
        logger.info("Receptor preparado exitosamente")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error al preparar el receptor: {e}")
        sys.exit(1)

def download_compounds(config):
    """Descargar compuestos de la base de datos."""
    logger.info("Descargando compuestos...")
    try:
        cmd = [
            sys.executable,
            str(PIPELINE_DIR / "scripts" / "download_compounds.py"),
            "--database", config['database']['name'],
            "--query", config['database']['query'],
            "--max_compounds", str(config['database']['max_compounds']),
            "--threads", str(config['docking']['threads'])
        ]
        subprocess.run(cmd, check=True)
        logger.info("Compuestos descargados exitosamente")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error al descargar compuestos: {e}")
        sys.exit(1)

def run_docking(config):
    """Ejecutar el docking."""
    logger.info("Ejecutando docking...")
    try:
        cmd = [
            sys.executable,
            str(PIPELINE_DIR / "run_pipeline.py"),
            "--receptor", f"{config['receptor']['name']}.pdbqt",
            "--center_x", str(config['receptor']['center']['x']),
            "--center_y", str(config['receptor']['center']['y']),
            "--center_z", str(config['receptor']['center']['z']),
            "--size_x", str(config['receptor']['size']['x']),
            "--size_y", str(config['receptor']['size']['y']),
            "--size_z", str(config['receptor']['size']['z']),
            "--threads", str(config['docking']['threads'])
        ]
        subprocess.run(cmd, check=True)
        logger.info("Docking completado exitosamente")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error al ejecutar el docking: {e}")
        sys.exit(1)

def visualize_results(config):
    """Visualizar los resultados del docking."""
    logger.info("Generando visualización...")
    try:
        result_file = PIPELINE_DIR / "results" / "docking_results.pdbqt"
        cmd = [
            sys.executable,
            str(PIPELINE_DIR / "scripts" / "visualize_results.py"),
            "--receptor", f"{config['receptor']['name']}.pdbqt",
            "--result", str(result_file),
            "--output", "visualization.html"
        ]
        subprocess.run(cmd, check=True)
        logger.info("Visualización generada exitosamente")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error al generar la visualización: {e}")
        sys.exit(1)

def main():
    """Función principal."""
    try:
        # Cargar configuración
        config = load_config()
        
        # Descargar y preparar receptor
        receptor_path = download_receptor(config)
        prepare_receptor(receptor_path)
        
        # Descargar compuestos
        download_compounds(config)
        
        # Ejecutar docking
        run_docking(config)
        
        # Visualizar resultados
        visualize_results(config)
        
        logger.info("\n¡Ejemplo completado exitosamente!")
        logger.info("Los resultados se encuentran en:")
        logger.info(f"- Receptor preparado: {PIPELINE_DIR}/receptors/{config['receptor']['name']}.pdbqt")
        logger.info(f"- Compuestos descargados: {PIPELINE_DIR}/ligands/")
        logger.info(f"- Resultados del docking: {PIPELINE_DIR}/results/")
        logger.info(f"- Visualización: {PIPELINE_DIR}/visualization.html")
        
    except Exception as e:
        logger.error(f"Error durante la ejecución del ejemplo: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 