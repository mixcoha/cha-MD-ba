#!/usr/bin/env python3
"""
Script para centrar la proteína en la caja de simulación
Basado en preprocess.py de CHA-MD-BA y herramientas de GROMACS
"""

import os
import subprocess
from pathlib import Path
import tempfile
from typing import Optional, Dict

class ProteinCenterer:
    """
    Clase para centrar proteínas en la caja de simulación usando GROMACS
    Inspirada en el módulo preprocess.py de CHA-MD-BA
    """
    
    def __init__(self, gmx_path: str = "/usr/local/gromacs/bin/gmx_mpi"):
        self.gmx = gmx_path
        
    def center_protein_in_box(self, input_file: str, output_file: str, 
                             topology_file: Optional[str] = None) -> str:
        """
        Centra la proteína en el centro de la caja de simulación
        
        Args:
            input_file: Archivo de entrada (.gro o .pdb)
            output_file: Archivo de salida
            topology_file: Archivo de topología (.tpr) opcional para mayor precisión
            
        Returns:
            Ruta al archivo de salida
        """
        print(f"🎯 Centrando proteína en la caja...")
        print(f"   Entrada: {input_file}")
        print(f"   Salida: {output_file}")
        
        # Usar trjconv para centrar
        cmd = [
            self.gmx, "trjconv",
            "-f", input_file,
            "-o", output_file,
            "-center",
            "-pbc", "mol"
        ]
        
        if topology_file and os.path.exists(topology_file):
            cmd.extend(["-s", topology_file])
            print(f"   Usando topología: {topology_file}")
        
        try:
            # Ejecutar comando con selección automática:
            # 1 = Protein (para centrar)
            # 0 = System (para output)
            process = subprocess.run(
                cmd, 
                input="1\n0\n", 
                text=True, 
                capture_output=True, 
                check=True
            )
            
            print(f"✅ Proteína centrada exitosamente")
            return output_file
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Error centrando proteína: {e}")
            print(f"Stderr: {e.stderr}")
            raise
    
    def center_and_compact(self, input_file: str, output_file: str,
                          topology_file: Optional[str] = None) -> str:
        """
        Centra la proteína y compacta la estructura
        Equivalente a center_and_align del módulo preprocess.py pero con GROMACS
        
        Args:
            input_file: Archivo de entrada
            output_file: Archivo de salida
            topology_file: Archivo de topología opcional
            
        Returns:
            Ruta al archivo final
        """
        print(f"🔄 Centrando y compactando estructura...")
        
        # Paso 1: Centrar
        temp_centered = f"temp_centered_{Path(input_file).stem}.gro"
        
        cmd_center = [
            self.gmx, "trjconv",
            "-f", input_file,
            "-o", temp_centered,
            "-center",
            "-pbc", "mol"
        ]
        
        if topology_file and os.path.exists(topology_file):
            cmd_center.extend(["-s", topology_file])
        
        try:
            # Centrar con proteína
            subprocess.run(
                cmd_center,
                input="1\n0\n",  # Protein para centro, System para output
                text=True,
                capture_output=True,
                check=True
            )
            
            # Paso 2: Compactar
            cmd_compact = [
                self.gmx, "trjconv",
                "-f", temp_centered,
                "-o", output_file,
                "-ur", "compact",
                "-pbc", "mol"
            ]
            
            if topology_file and os.path.exists(topology_file):
                cmd_compact.extend(["-s", topology_file])
            
            subprocess.run(
                cmd_compact,
                input="0\n",  # System
                text=True,
                capture_output=True,
                check=True
            )
            
            # Limpiar archivo temporal
            if os.path.exists(temp_centered):
                os.remove(temp_centered)
            
            print(f"✅ Estructura centrada y compactada: {output_file}")
            return output_file
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Error en el proceso: {e}")
            # Limpiar archivo temporal
            if os.path.exists(temp_centered):
                os.remove(temp_centered)
            raise
    
    def fix_broken_molecules(self, input_file: str, output_file: str,
                           topology_file: Optional[str] = None) -> str:
        """
        Corrige moléculas partidas por condiciones periódicas de frontera
        
        Args:
            input_file: Archivo de entrada
            output_file: Archivo de salida  
            topology_file: Archivo de topología opcional
            
        Returns:
            Ruta al archivo corregido
        """
        print(f"🔧 Corrigiendo moléculas partidas...")
        
        cmd = [
            self.gmx, "trjconv",
            "-f", input_file,
            "-o", output_file,
            "-pbc", "whole"
        ]
        
        if topology_file and os.path.exists(topology_file):
            cmd.extend(["-s", topology_file])
        
        try:
            subprocess.run(
                cmd,
                input="0\n",  # System
                text=True,
                capture_output=True,
                check=True
            )
            
            print(f"✅ Moléculas corregidas: {output_file}")
            return output_file
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Error corrigiendo moléculas: {e}")
            raise
    
    def complete_centering_pipeline(self, input_file: str, output_file: str,
                                  topology_file: Optional[str] = None) -> Dict[str, str]:
        """
        Pipeline completo de centrado inspirado en preprocess.py
        
        1. Corrige moléculas partidas
        2. Centra la proteína
        3. Compacta la estructura
        
        Args:
            input_file: Archivo de entrada
            output_file: Archivo de salida final
            topology_file: Archivo de topología opcional
            
        Returns:
            Diccionario con rutas de archivos generados
        """
        print("🚀 Iniciando pipeline completo de centrado...")
        
        base_name = Path(output_file).stem
        output_dir = Path(output_file).parent
        
        # Paso 1: Corregir moléculas partidas
        step1_file = output_dir / f"{base_name}_step1_whole.gro"
        self.fix_broken_molecules(input_file, str(step1_file), topology_file)
        
        # Paso 2: Centrar proteína
        step2_file = output_dir / f"{base_name}_step2_centered.gro"
        self.center_protein_in_box(str(step1_file), str(step2_file), topology_file)
        
        # Paso 3: Compactar (archivo final)
        cmd_compact = [
            self.gmx, "trjconv",
            "-f", str(step2_file),
            "-o", output_file,
            "-ur", "compact",
            "-pbc", "mol"
        ]
        
        if topology_file and os.path.exists(topology_file):
            cmd_compact.extend(["-s", topology_file])
        
        try:
            subprocess.run(
                cmd_compact,
                input="0\n",
                text=True,
                capture_output=True,
                check=True
            )
            
            print(f"✅ Pipeline completado: {output_file}")
            
            # Limpiar archivos intermedios
            if step1_file.exists():
                step1_file.unlink()
            if step2_file.exists():
                step2_file.unlink()
            
            return {
                "final": output_file,
                "status": "completed"
            }
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Error en compactación final: {e}")
            raise
    
    def convert_to_pdb_centered(self, gro_file: str, pdb_file: str,
                              topology_file: Optional[str] = None) -> str:
        """
        Convierte archivo .gro centrado a .pdb
        
        Args:
            gro_file: Archivo .gro de entrada
            pdb_file: Archivo .pdb de salida
            topology_file: Archivo de topología opcional
            
        Returns:
            Ruta al archivo PDB generado
        """
        print(f"📄 Convirtiendo a PDB centrado...")
        
        cmd = [
            self.gmx, "editconf",
            "-f", gro_file,
            "-o", pdb_file
        ]
        
        if topology_file and os.path.exists(topology_file):
            cmd.extend(["-s", topology_file])
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"✅ PDB centrado generado: {pdb_file}")
            return pdb_file
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Error generando PDB: {e}")
            raise

def main():
    """Función principal para centrar proteínas"""
    print("🎯 Centrador de Proteínas - Receptor GABA")
    print("Basado en preprocess.py de CHA-MD-BA")
    print("=" * 50)
    
    centerer = ProteinCenterer()
    
    # Archivos a procesar
    files_to_process = [
        {
            'input': 'nvt/nvt.gro',
            'output': 'receptor_gaba_nvt_centered.gro',
            'pdb_output': 'receptor_gaba_nvt_centered.pdb',
            'topology': 'nvt/nvt.tpr',
            'description': 'NVT Equilibrated'
        },
        {
            'input': 'em/minim.gro', 
            'output': 'receptor_gaba_minimized_centered.gro',
            'pdb_output': 'receptor_gaba_minimized_centered.pdb',
            'topology': 'nvt/nvt.tpr',  # Usar la topología más reciente
            'description': 'Energy Minimized'
        },
        {
            'input': 'processed.gro',
            'output': 'receptor_gaba_initial_centered.gro', 
            'pdb_output': 'receptor_gaba_initial_centered.pdb',
            'topology': None,
            'description': 'Initial Structure'
        }
    ]
    
    successful_files = []
    
    for file_config in files_to_process:
        if not os.path.exists(file_config['input']):
            print(f"⚠️ Archivo no encontrado: {file_config['input']}")
            continue
        
        print(f"\n🔄 Procesando: {file_config['description']}")
        
        try:
            # Ejecutar pipeline completo de centrado
            result = centerer.complete_centering_pipeline(
                file_config['input'],
                file_config['output'],
                file_config['topology']
            )
            
            # Convertir a PDB
            centerer.convert_to_pdb_centered(
                file_config['output'],
                file_config['pdb_output'],
                file_config['topology']
            )
            
            successful_files.extend([
                file_config['output'],
                file_config['pdb_output']
            ])
            
        except Exception as e:
            print(f"❌ Error procesando {file_config['description']}: {e}")
            continue
    
    # Resumen final
    if successful_files:
        print(f"\n🎉 Centrado completado exitosamente!")
        print(f"📁 Archivos generados:")
        for file in successful_files:
            if os.path.exists(file):
                size_mb = os.path.getsize(file) / (1024 * 1024)
                print(f"   • {file} ({size_mb:.2f} MB)")
    else:
        print(f"\n❌ No se pudieron procesar archivos")
    
    return 0

if __name__ == "__main__":
    exit(main())
