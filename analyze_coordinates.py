#!/usr/bin/env python3
"""
Script avanzado para análisis de coordenadas de archivos GROMACS .gro
Basado en el Manual de GROMACS 2025.2 - Análisis estructural y dinámico
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse
from typing import List, Dict, Tuple
import subprocess

class GromacsCoordinateAnalyzer:
    """
    Analizador de coordenadas GROMACS siguiendo las mejores prácticas del manual
    """
    
    def __init__(self, gmx_path: str = "/usr/local/gromacs/bin/gmx_mpi"):
        self.gmx = gmx_path
        self.output_dir = Path("analysis_output")
        self.output_dir.mkdir(exist_ok=True)
        
    def read_gro_file(self, gro_file: str) -> Dict:
        """
        Lee archivo .gro y extrae coordenadas y metadatos
        Formato GROMACS: residue_number+residue_name+atom_name+atom_number+x+y+z
        """
        with open(gro_file, 'r') as f:
            lines = f.readlines()
        
        # Primera línea: título
        title = lines[0].strip()
        
        # Segunda línea: número de átomos
        n_atoms = int(lines[1].strip())
        
        # Última línea: dimensiones de la caja
        box_line = lines[-1].strip().split()
        box_dims = [float(x) for x in box_line]
        
        # Coordenadas de átomos
        atoms = []
        for i in range(2, 2 + n_atoms):
            line = lines[i]
            # Formato fijo de GROMACS .gro
            res_num = int(line[0:5])
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_num = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            
            atoms.append({
                'res_num': res_num,
                'res_name': res_name,
                'atom_name': atom_name,
                'atom_num': atom_num,
                'x': x, 'y': y, 'z': z
            })
        
        return {
            'title': title,
            'n_atoms': n_atoms,
            'box_dims': box_dims,
            'atoms': atoms
        }
    
    def calculate_rmsd(self, structure_file: str, trajectory_file: str = None, 
                      reference_file: str = None) -> str:
        """
        Calcula RMSD usando gmx rms - Análisis de estabilidad estructural
        Manual GROMACS Section 7.3.1: Análisis de trayectorias
        """
        output_file = self.output_dir / "rmsd.xvg"
        
        if trajectory_file:
            cmd = [
                self.gmx, "rms",
                "-s", structure_file,
                "-f", trajectory_file,
                "-o", str(output_file),
                "-tu", "ps"
            ]
        else:
            # Solo estructura final vs inicial
            cmd = [
                self.gmx, "rms",
                "-s", reference_file if reference_file else structure_file,
                "-f", structure_file,
                "-o", str(output_file),
                "-tu", "ps"
            ]
        
        # Seleccionar backbone para RMSD
        process = subprocess.run(cmd, input=b"4\n4\n", capture_output=True, text=True)
        
        if process.returncode == 0:
            print(f"✅ RMSD calculado exitosamente: {output_file}")
        else:
            print(f"❌ Error calculando RMSD: {process.stderr}")
        
        return str(output_file)
    
    def calculate_radius_of_gyration(self, structure_file: str, trajectory_file: str = None) -> str:
        """
        Calcula radio de giro usando gmx gyrate
        Manual GROMACS: Análisis de compactación estructural
        """
        output_file = self.output_dir / "gyrate.xvg"
        
        if trajectory_file:
            cmd = [
                self.gmx, "gyrate",
                "-s", structure_file,
                "-f", trajectory_file,
                "-o", str(output_file)
            ]
        else:
            cmd = [
                self.gmx, "gyrate",
                "-s", structure_file,
                "-f", structure_file,
                "-o", str(output_file)
            ]
        
        # Seleccionar proteína
        process = subprocess.run(cmd, input=b"1\n", capture_output=True, text=True)
        
        if process.returncode == 0:
            print(f"✅ Radio de giro calculado: {output_file}")
        else:
            print(f"❌ Error calculando radio de giro: {process.stderr}")
        
        return str(output_file)
    
    def analyze_energy(self, edr_file: str) -> str:
        """
        Analiza energías usando gmx energy
        Manual GROMACS Section 7.1: Análisis energético
        """
        output_file = self.output_dir / "energy.xvg"
        
        cmd = [
            self.gmx, "energy",
            "-f", edr_file,
            "-o", str(output_file)
        ]
        
        # Seleccionar energía potencial (10), cinética (11), total (12), temperatura (15)
        process = subprocess.run(cmd, input=b"10 11 12 15\n0\n", capture_output=True, text=True)
        
        if process.returncode == 0:
            print(f"✅ Análisis energético completado: {output_file}")
        else:
            print(f"❌ Error en análisis energético: {process.stderr}")
        
        return str(output_file)
    
    def extract_coordinates_summary(self, gro_files: List[str]) -> pd.DataFrame:
        """
        Extrae resumen de coordenadas de múltiples archivos .gro
        """
        summaries = []
        
        for gro_file in gro_files:
            data = self.read_gro_file(gro_file)
            
            # Calcular centro de masa
            coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in data['atoms']])
            center_of_mass = np.mean(coords, axis=0)
            
            # Dimensiones del sistema
            min_coords = np.min(coords, axis=0)
            max_coords = np.max(coords, axis=0)
            system_size = max_coords - min_coords
            
            summaries.append({
                'file': Path(gro_file).name,
                'n_atoms': data['n_atoms'],
                'box_x': data['box_dims'][0] if data['box_dims'] else 0,
                'box_y': data['box_dims'][1] if len(data['box_dims']) > 1 else 0,
                'box_z': data['box_dims'][2] if len(data['box_dims']) > 2 else 0,
                'center_x': center_of_mass[0],
                'center_y': center_of_mass[1],
                'center_z': center_of_mass[2],
                'size_x': system_size[0],
                'size_y': system_size[1],
                'size_z': system_size[2]
            })
        
        return pd.DataFrame(summaries)
    
    def create_analysis_report(self, gro_files: List[str], edr_file: str = None) -> str:
        """
        Genera reporte completo de análisis en formato PDF/HTML
        """
        print("📊 Generando reporte de análisis estructural...")
        
        # 1. Resumen de coordenadas
        coord_summary = self.extract_coordinates_summary(gro_files)
        
        # 2. Análisis RMSD si hay múltiples estructuras
        if len(gro_files) > 1:
            rmsd_file = self.calculate_rmsd(gro_files[0], reference_file=gro_files[0])
        
        # 3. Radio de giro
        gyrate_file = self.calculate_radius_of_gyration(gro_files[-1])  # Última estructura
        
        # 4. Análisis energético si hay archivo .edr
        if edr_file and Path(edr_file).exists():
            energy_file = self.analyze_energy(edr_file)
        
        # 5. Guardar resumen en CSV
        csv_file = self.output_dir / "coordinate_summary.csv"
        coord_summary.to_csv(csv_file, index=False)
        print(f"📋 Resumen guardado en: {csv_file}")
        
        # 6. Crear visualizaciones
        self.create_plots(coord_summary)
        
        return str(self.output_dir)
    
    def create_plots(self, summary_df: pd.DataFrame):
        """
        Crea gráficos de análisis estructural
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Análisis Estructural del Receptor GABA', fontsize=16)
        
        # 1. Dimensiones del sistema
        ax1 = axes[0, 0]
        ax1.bar(['X', 'Y', 'Z'], [summary_df['size_x'].iloc[-1], 
                                  summary_df['size_y'].iloc[-1], 
                                  summary_df['size_z'].iloc[-1]])
        ax1.set_title('Dimensiones del Sistema (nm)')
        ax1.set_ylabel('Tamaño (nm)')
        
        # 2. Centro de masa
        ax2 = axes[0, 1]
        if len(summary_df) > 1:
            ax2.plot(summary_df.index, summary_df['center_x'], 'r-', label='X')
            ax2.plot(summary_df.index, summary_df['center_y'], 'g-', label='Y')
            ax2.plot(summary_df.index, summary_df['center_z'], 'b-', label='Z')
            ax2.set_title('Evolución del Centro de Masa')
            ax2.set_xlabel('Estructura')
            ax2.set_ylabel('Coordenada (nm)')
            ax2.legend()
        else:
            ax2.bar(['X', 'Y', 'Z'], [summary_df['center_x'].iloc[0],
                                      summary_df['center_y'].iloc[0],
                                      summary_df['center_z'].iloc[0]])
            ax2.set_title('Centro de Masa')
        
        # 3. Dimensiones de caja
        ax3 = axes[1, 0]
        ax3.bar(['Box X', 'Box Y', 'Box Z'], [summary_df['box_x'].iloc[-1],
                                              summary_df['box_y'].iloc[-1],
                                              summary_df['box_z'].iloc[-1]])
        ax3.set_title('Dimensiones de la Caja de Simulación (nm)')
        ax3.set_ylabel('Tamaño (nm)')
        
        # 4. Número de átomos
        ax4 = axes[1, 1]
        ax4.text(0.5, 0.5, f'Átomos: {summary_df["n_atoms"].iloc[-1]:,}', 
                ha='center', va='center', fontsize=20, 
                transform=ax4.transAxes)
        ax4.set_title('Sistema Final')
        ax4.axis('off')
        
        plt.tight_layout()
        plot_file = self.output_dir / "structural_analysis.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📈 Gráficos guardados en: {plot_file}")

def main():
    parser = argparse.ArgumentParser(description='Análisis avanzado de coordenadas GROMACS')
    parser.add_argument('--gro-files', nargs='+', required=True, 
                       help='Archivos .gro a analizar')
    parser.add_argument('--edr-file', help='Archivo .edr para análisis energético')
    parser.add_argument('--output-dir', default='analysis_output', 
                       help='Directorio de salida')
    
    args = parser.parse_args()
    
    analyzer = GromacsCoordinateAnalyzer()
    analyzer.output_dir = Path(args.output_dir)
    
    # Verificar que los archivos existen
    for gro_file in args.gro_files:
        if not Path(gro_file).exists():
            print(f"❌ Error: No se encuentra {gro_file}")
            return 1
    
    # Generar reporte
    report_dir = analyzer.create_analysis_report(args.gro_files, args.edr_file)
    
    print(f"\n🎉 Análisis completado!")
    print(f"📁 Resultados en: {report_dir}")
    
    return 0

if __name__ == "__main__":
    exit(main())
