#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script para preparar archivos PDBQT para AutoDock-Vina usando OpenBabel y RDKit.
Este script convierte archivos PDB a formato PDBQT, preparando tanto receptores como ligandos.
"""

import os
import sys
import argparse
import glob
import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuración de rutas
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RECEPTORS_DIR = BASE_DIR / "receptors"
LIGANDS_DIR = BASE_DIR / "ligands"

def parse_arguments():
    """Parsear argumentos de línea de comandos."""
    parser = argparse.ArgumentParser(description="Preparar archivos PDBQT para AutoDock-Vina")
    
    parser.add_argument("--receptor", type=str, help="Archivo PDB del receptor a preparar")
    parser.add_argument("--ligand", type=str, help="Archivo PDB del ligando a preparar")
    parser.add_argument("--all_receptors", action="store_true", help="Preparar todos los receptores en el directorio")
    parser.add_argument("--all_ligands", action="store_true", help="Preparar todos los ligandos en el directorio")
    
    return parser.parse_args()

def prepare_receptor(pdb_path):
    """Preparar un receptor para docking usando OpenBabel."""
    try:
        output_path = pdb_path.with_suffix('.pdbqt')
        
        # Comando de OpenBabel para convertir PDB a PDBQT
        cmd = [
            "obabel",
            str(pdb_path),
            "-O", str(output_path),
            "--addhydrogens",
            "--partialcharge", "gasteiger",
            "--gen3d"
        ]
        
        logger.info(f"Preparando receptor: {pdb_path}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Error al preparar el receptor: {result.stderr}")
            return False
        
        logger.info(f"Receptor preparado exitosamente: {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error al preparar el receptor {pdb_path}: {e}")
        return False

def prepare_ligand(pdb_path):
    """Preparar un ligando para docking usando RDKit."""
    try:
        # Leer el archivo PDB
        mol = Chem.MolFromPDBFile(str(pdb_path))
        if mol is None:
            logger.error(f"No se pudo leer el archivo PDB: {pdb_path}")
            return False
        
        # Añadir hidrógenos
        mol = Chem.AddHs(mol)
        
        # Generar conformeros 3D
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Guardar en formato PDBQT
        output_path = pdb_path.with_suffix('.pdbqt')
        Chem.MolToPDBFile(mol, str(output_path))
        
        logger.info(f"Ligando preparado exitosamente: {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error al preparar el ligando {pdb_path}: {e}")
        return False

def prepare_all_receptors():
    """Preparar todos los receptores en el directorio."""
    pdb_files = glob.glob(str(RECEPTORS_DIR / "*.pdb"))
    if not pdb_files:
        logger.warning("No se encontraron archivos PDB en el directorio de receptores")
        return False
    
    success_count = 0
    for pdb_file in pdb_files:
        if prepare_receptor(Path(pdb_file)):
            success_count += 1
    
    logger.info(f"Preparación de receptores completada. {success_count} de {len(pdb_files)} exitosos.")
    return success_count > 0

def prepare_all_ligands():
    """Preparar todos los ligandos en el directorio."""
    pdb_files = glob.glob(str(LIGANDS_DIR / "*.pdb"))
    if not pdb_files:
        logger.warning("No se encontraron archivos PDB en el directorio de ligandos")
        return False
    
    success_count = 0
    for pdb_file in pdb_files:
        if prepare_ligand(Path(pdb_file)):
            success_count += 1
    
    logger.info(f"Preparación de ligandos completada. {success_count} de {len(pdb_files)} exitosos.")
    return success_count > 0

def main():
    """Función principal."""
    # Parsear argumentos
    args = parse_arguments()
    
    # Verificar que al menos una opción fue proporcionada
    if not (args.receptor or args.ligand or args.all_receptors or args.all_ligands):
        logger.error("Error: Debes especificar al menos una opción de preparación.")
        logger.error("Usa --help para ver las opciones disponibles.")
        sys.exit(1)
    
    # Preparar receptor individual si se especificó
    if args.receptor:
        receptor_path = Path(args.receptor)
        if not receptor_path.exists():
            logger.error(f"El archivo del receptor no existe: {receptor_path}")
            sys.exit(1)
        prepare_receptor(receptor_path)
    
    # Preparar ligando individual si se especificó
    if args.ligand:
        ligand_path = Path(args.ligand)
        if not ligand_path.exists():
            logger.error(f"El archivo del ligando no existe: {ligand_path}")
            sys.exit(1)
        prepare_ligand(ligand_path)
    
    # Preparar todos los receptores si se solicitó
    if args.all_receptors:
        prepare_all_receptors()
    
    # Preparar todos los ligandos si se solicitó
    if args.all_ligands:
        prepare_all_ligands()

if __name__ == "__main__":
    main() 