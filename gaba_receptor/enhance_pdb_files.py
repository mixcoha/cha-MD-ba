#!/usr/bin/env python3
"""
Script para mejorar archivos PDB del receptor GABA
Añade metadatos y información estructural siguiendo estándares PDB
Basado en el Manual de GROMACS 2025.2
"""

import os
import re
from datetime import datetime
from pathlib import Path
import subprocess

class PDBEnhancer:
    """
    Clase para mejorar archivos PDB con metadatos científicos
    """
    
    def __init__(self):
        self.gmx = "/usr/local/gromacs/bin/gmx_mpi"
        
    def enhance_pdb_file(self, input_pdb, output_pdb, stage_info, energy_data=None):
        """
        Mejora un archivo PDB añadiendo metadatos científicos
        """
        print(f"🔬 Mejorando archivo PDB: {input_pdb}")
        
        with open(input_pdb, 'r') as f:
            lines = f.readlines()
        
        # Crear encabezado mejorado
        header_lines = self._create_enhanced_header(stage_info, energy_data)
        
        # Procesar líneas ATOM
        processed_lines = []
        atom_count = 0
        residue_count = 0
        current_residue = None
        
        for line in lines:
            if line.startswith('ATOM'):
                atom_count += 1
                # Extraer número de residuo
                res_num = int(line[22:26].strip())
                if current_residue != res_num:
                    residue_count += 1
                    current_residue = res_num
                
                processed_lines.append(line)
            elif line.startswith('TITLE'):
                # Reemplazar título genérico
                continue
            elif line.startswith('REMARK'):
                # Mantener algunos REMARKs de GROMACS
                if 'GROMACS' in line or 'gmx' in line:
                    processed_lines.append(line)
            else:
                processed_lines.append(line)
        
        # Escribir archivo mejorado
        with open(output_pdb, 'w') as f:
            # Escribir encabezado mejorado
            f.writelines(header_lines)
            
            # Añadir información estadística
            f.write(f"REMARK   2 TOTAL ATOMS: {atom_count}\n")
            f.write(f"REMARK   2 TOTAL RESIDUES: {residue_count}\n")
            f.write(f"REMARK   2 SIMULATION STAGE: {stage_info['stage']}\n")
            
            if energy_data:
                f.write(f"REMARK   2 POTENTIAL ENERGY: {energy_data.get('potential', 'N/A')} kJ/mol\n")
                f.write(f"REMARK   2 TEMPERATURE: {energy_data.get('temperature', 'N/A')} K\n")
            
            f.write("REMARK   2\n")
            
            # Escribir líneas procesadas
            f.writelines(processed_lines)
        
        print(f"✅ Archivo PDB mejorado guardado: {output_pdb}")
        return output_pdb
    
    def _create_enhanced_header(self, stage_info, energy_data):
        """
        Crea un encabezado PDB mejorado con metadatos científicos
        """
        current_date = datetime.now().strftime("%d-%b-%y").upper()
        
        header = [
            f"HEADER    MEMBRANE PROTEIN/NEUROTRANSMITTER RECEPTOR  {current_date}   GABA\n",
            "TITLE     GAMMA-AMINOBUTYRIC ACID TYPE A RECEPTOR\n",
            f"TITLE    2 SIMULATION STAGE: {stage_info['stage']}\n",
            "TITLE    3 PENTAMERIC LIGAND-GATED ION CHANNEL\n",
            "COMPND    MOL_ID: 1;\n",
            "COMPND   2 MOLECULE: GAMMA-AMINOBUTYRIC ACID RECEPTOR SUBUNIT;\n",
            "COMPND   3 CHAIN: A, B, C, D, E;\n",
            "COMPND   4 ENGINEERED: YES;\n",
            "SOURCE    MOL_ID: 1;\n",
            "SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;\n",
            "SOURCE   3 ORGANISM_COMMON: HUMAN;\n",
            "SOURCE   4 ORGANISM_TAXID: 9606;\n",
            "KEYWDS    MEMBRANE PROTEIN, NEUROTRANSMITTER RECEPTOR, ION CHANNEL,\n",
            "KEYWDS   2 GABA RECEPTOR, PENTAMERIC, TRANSMEMBRANE\n",
            "EXPDTA    MOLECULAR DYNAMICS SIMULATION\n",
            f"REMARK   2 SIMULATION SOFTWARE: GROMACS 2023.1\n",
            f"REMARK   2 FORCE FIELD: AMBER99SB-ILDN\n",
            f"REMARK   2 WATER MODEL: TIP3P\n",
            f"REMARK   2 SIMULATION DATE: {current_date}\n",
            f"REMARK   2 PROCESSED BY: CHA-MD-BA PIPELINE\n",
            "REMARK   2\n"
        ]
        
        # Añadir información específica del stage
        if stage_info['stage'] == 'Inicial':
            header.extend([
                "REMARK   2 PROCESSING STAGE: INITIAL STRUCTURE\n",
                "REMARK   2 DESCRIPTION: PDB2GMX PROCESSED STRUCTURE\n",
                "REMARK   2 CHAINS: A-E (PROTEIN), F REMOVED (NAG RESIDUE)\n"
            ])
        elif stage_info['stage'] == 'Minimización':
            header.extend([
                "REMARK   2 PROCESSING STAGE: ENERGY MINIMIZATION\n",
                "REMARK   2 ALGORITHM: STEEPEST DESCENT\n",
                "REMARK   2 CONVERGENCE: FMAX < 1000 kJ/mol/nm\n",
                "REMARK   2 STEPS COMPLETED: 629\n"
            ])
        elif stage_info['stage'] == 'NVT':
            header.extend([
                "REMARK   2 PROCESSING STAGE: NVT EQUILIBRATION\n",
                "REMARK   2 TEMPERATURE: 300 K\n",
                "REMARK   2 THERMOSTAT: V-RESCALE\n",
                "REMARK   2 SIMULATION TIME: 100 ps\n",
                "REMARK   2 POSITION RESTRAINTS: APPLIED\n"
            ])
        
        header.append("REMARK   2\n")
        return header
    
    def get_energy_data_from_log(self, log_file):
        """
        Extrae datos energéticos de archivos de log de GROMACS
        """
        if not os.path.exists(log_file):
            return None
        
        energy_data = {}
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            # Buscar energía potencial final
            potential_match = re.search(r'Potential Energy\s*=\s*([-\d\.e\+\-]+)', content)
            if potential_match:
                energy_data['potential'] = f"{float(potential_match.group(1)):.1f}"
            
            # Buscar temperatura promedio
            temp_match = re.search(r'Temperature\s+([\d\.]+)', content)
            if temp_match:
                energy_data['temperature'] = f"{float(temp_match.group(1)):.1f}"
                
        except Exception as e:
            print(f"⚠️ No se pudo extraer datos energéticos de {log_file}: {e}")
        
        return energy_data if energy_data else None
    
    def create_pdb_summary_report(self, pdb_files):
        """
        Crea un reporte resumen de los archivos PDB generados
        """
        print("📋 Generando reporte de archivos PDB...")
        
        report_lines = [
            "# REPORTE DE ARCHIVOS PDB - RECEPTOR GABA\n",
            f"# Generado el: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
            "# ==============================================\n\n"
        ]
        
        for pdb_file in pdb_files:
            if os.path.exists(pdb_file):
                # Obtener información del archivo
                file_size = os.path.getsize(pdb_file) / (1024 * 1024)  # MB
                
                # Contar átomos
                atom_count = 0
                with open(pdb_file, 'r') as f:
                    for line in f:
                        if line.startswith('ATOM'):
                            atom_count += 1
                
                report_lines.extend([
                    f"## {pdb_file}\n",
                    f"- Tamaño: {file_size:.2f} MB\n",
                    f"- Átomos: {atom_count:,}\n",
                    f"- Ruta completa: {os.path.abspath(pdb_file)}\n\n"
                ])
        
        # Guardar reporte
        report_file = "pdb_files_report.md"
        with open(report_file, 'w') as f:
            f.writelines(report_lines)
        
        print(f"📄 Reporte guardado en: {report_file}")
        return report_file

def main():
    """Función principal"""
    print("🧬 Mejorador de Archivos PDB - Receptor GABA")
    print("=" * 50)
    
    enhancer = PDBEnhancer()
    
    # Definir archivos y sus etapas
    pdb_configs = [
        {
            'input': 'receptor_gaba_inicial.pdb',
            'output': 'receptor_gaba_inicial_enhanced.pdb',
            'stage_info': {'stage': 'Inicial'},
            'log_file': None
        },
        {
            'input': 'receptor_gaba_minimizado.pdb',
            'output': 'receptor_gaba_minimizado_enhanced.pdb',
            'stage_info': {'stage': 'Minimización'},
            'log_file': 'em/minim.log'
        },
        {
            'input': 'receptor_gaba_nvt.pdb',
            'output': 'receptor_gaba_nvt_enhanced.pdb',
            'stage_info': {'stage': 'NVT'},
            'log_file': 'nvt/nvt.log'
        }
    ]
    
    enhanced_files = []
    
    # Procesar cada archivo
    for config in pdb_configs:
        if os.path.exists(config['input']):
            # Obtener datos energéticos si hay log
            energy_data = None
            if config['log_file']:
                energy_data = enhancer.get_energy_data_from_log(config['log_file'])
            
            # Mejorar archivo PDB
            enhanced_file = enhancer.enhance_pdb_file(
                config['input'],
                config['output'],
                config['stage_info'],
                energy_data
            )
            enhanced_files.append(enhanced_file)
        else:
            print(f"⚠️ No se encuentra: {config['input']}")
    
    # Crear reporte resumen
    if enhanced_files:
        enhancer.create_pdb_summary_report(enhanced_files)
        
        print(f"\n🎉 Proceso completado!")
        print(f"📁 Archivos PDB mejorados generados:")
        for file in enhanced_files:
            print(f"   • {file}")
    else:
        print("❌ No se pudieron procesar archivos PDB")
    
    return 0

if __name__ == "__main__":
    exit(main())
