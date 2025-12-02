#!/usr/bin/env python3
"""
Generador de reporte PDF para estructuras de proteínas GROMACS
Basado en análisis de archivos .gro
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from datetime import datetime

class ProteinStructurePDFGenerator:
    """
    Generador de reportes PDF para estructuras de proteínas
    """
    
    def __init__(self):
        self.protein_data = {}
        self.analysis_results = {}
        
    def read_gro_file(self, gro_file):
        """Lee y procesa archivo .gro"""
        with open(gro_file, 'r') as f:
            lines = f.readlines()
        
        title = lines[0].strip()
        n_atoms = int(lines[1].strip())
        
        # Dimensiones de la caja
        box_line = lines[-1].strip().split()
        box_dims = [float(x) for x in box_line]
        
        # Procesar átomos
        atoms_data = []
        residue_info = {}
        chain_info = {}
        
        for i in range(2, 2 + n_atoms):
            line = lines[i]
            res_num = int(line[0:5])
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_num = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            
            atoms_data.append({
                'res_num': res_num,
                'res_name': res_name,
                'atom_name': atom_name,
                'atom_num': atom_num,
                'x': x, 'y': y, 'z': z
            })
            
            # Contar residuos
            if res_name not in residue_info:
                residue_info[res_name] = 0
            residue_info[res_name] += 1
        
        coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms_data])
        
        return {
            'title': title,
            'n_atoms': n_atoms,
            'box_dims': box_dims,
            'atoms': atoms_data,
            'coordinates': coords,
            'residue_info': residue_info,
            'center_of_mass': np.mean(coords, axis=0),
            'min_coords': np.min(coords, axis=0),
            'max_coords': np.max(coords, axis=0),
            'system_size': np.max(coords, axis=0) - np.min(coords, axis=0)
        }
    
    def analyze_protein_structure(self, gro_files):
        """Analiza múltiples estructuras de proteína"""
        results = []
        
        for i, gro_file in enumerate(gro_files):
            if not Path(gro_file).exists():
                continue
                
            data = self.read_gro_file(gro_file)
            stage_name = self._get_stage_name(gro_file)
            
            # Análisis por cadena (asumiendo que los primeros dígitos del residuo indican la cadena)
            chain_analysis = self._analyze_chains(data['atoms'])
            
            result = {
                'stage': stage_name,
                'file': Path(gro_file).name,
                'n_atoms': data['n_atoms'],
                'n_residues': len(set(atom['res_num'] for atom in data['atoms'])),
                'box_volume': np.prod(data['box_dims'][:3]) if len(data['box_dims']) >= 3 else 0,
                'center_x': data['center_of_mass'][0],
                'center_y': data['center_of_mass'][1],
                'center_z': data['center_of_mass'][2],
                'size_x': data['system_size'][0],
                'size_y': data['system_size'][1],
                'size_z': data['system_size'][2],
                'residue_types': len(data['residue_info']),
                'chain_info': chain_analysis,
                'full_data': data
            }
            
            results.append(result)
        
        return results
    
    def _get_stage_name(self, filename):
        """Determina la etapa de simulación basada en el nombre del archivo"""
        if 'processed' in filename:
            return 'Inicial'
        elif 'minim' in filename:
            return 'Minimización'
        elif 'nvt' in filename:
            return 'Equilibración NVT'
        elif 'npt' in filename:
            return 'Equilibración NPT'
        else:
            return 'Desconocido'
    
    def _analyze_chains(self, atoms):
        """Analiza las cadenas de la proteína"""
        chains = {}
        
        for atom in atoms:
            # Identificar cadena basada en rangos de residuos típicos
            res_num = atom['res_num']
            
            if res_num <= 429:
                chain_id = 'A'
            elif res_num <= 858:
                chain_id = 'B'
            elif res_num <= 1306:
                chain_id = 'C'
            elif res_num <= 1754:
                chain_id = 'D'
            else:
                chain_id = 'E'
            
            if chain_id not in chains:
                chains[chain_id] = {'atoms': 0, 'residues': set()}
            
            chains[chain_id]['atoms'] += 1
            chains[chain_id]['residues'].add(res_num)
        
        # Convertir sets a counts
        for chain in chains:
            chains[chain]['n_residues'] = len(chains[chain]['residues'])
            del chains[chain]['residues']
        
        return chains
    
    def generate_pdf_report(self, gro_files, output_filename='receptor_gaba_reporte.pdf'):
        """Genera el reporte PDF completo"""
        
        print("📄 Generando reporte PDF del receptor GABA...")
        
        # Analizar estructuras
        analysis_results = self.analyze_protein_structure(gro_files)
        
        if not analysis_results:
            print("❌ No se pudieron analizar los archivos .gro")
            return
        
        # Crear PDF
        with PdfPages(output_filename) as pdf:
            
            # Página 1: Portada
            self._create_title_page(pdf, analysis_results)
            
            # Página 2: Resumen ejecutivo
            self._create_summary_page(pdf, analysis_results)
            
            # Página 3: Análisis estructural
            self._create_structural_analysis_page(pdf, analysis_results)
            
            # Página 4: Análisis por cadenas
            self._create_chain_analysis_page(pdf, analysis_results)
            
            # Página 5: Evolución estructural
            if len(analysis_results) > 1:
                self._create_evolution_page(pdf, analysis_results)
            
            # Página 6: Coordenadas detalladas
            self._create_coordinates_page(pdf, analysis_results)
        
        print(f"✅ Reporte PDF generado: {output_filename}")
        return output_filename
    
    def _create_title_page(self, pdf, results):
        """Crea la página de título"""
        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis('off')
        
        # Título principal
        ax.text(0.5, 0.8, 'REPORTE ESTRUCTURAL', ha='center', va='center',
                fontsize=24, fontweight='bold', transform=ax.transAxes)
        
        ax.text(0.5, 0.75, 'Receptor GABA', ha='center', va='center',
                fontsize=20, fontweight='bold', color='blue', transform=ax.transAxes)
        
        # Información del sistema
        final_result = results[-1]
        
        info_text = f"""
Sistema: {final_result['full_data']['title']}
Número de átomos: {final_result['n_atoms']:,}
Número de residuos: {final_result['n_residues']:,}
Cadenas analizadas: {len(final_result['chain_info'])}

Dimensiones del sistema:
• X: {final_result['size_x']:.2f} nm
• Y: {final_result['size_y']:.2f} nm  
• Z: {final_result['size_z']:.2f} nm

Volumen de la caja: {final_result['box_volume']:.2f} nm³
"""
        
        ax.text(0.5, 0.45, info_text, ha='center', va='center',
                fontsize=12, transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.7))
        
        # Fecha y herramientas
        ax.text(0.5, 0.1, f'Generado el: {datetime.now().strftime("%d/%m/%Y %H:%M")}', 
                ha='center', va='center', fontsize=10, transform=ax.transAxes)
        
        ax.text(0.5, 0.05, 'Herramientas: GROMACS + CHA-MD-BA', 
                ha='center', va='center', fontsize=10, style='italic', transform=ax.transAxes)
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_summary_page(self, pdf, results):
        """Crea página de resumen ejecutivo"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8.5, 11))
        fig.suptitle('RESUMEN EJECUTIVO', fontsize=16, fontweight='bold')
        
        # Gráfico 1: Evolución del número de átomos
        stages = [r['stage'] for r in results]
        n_atoms = [r['n_atoms'] for r in results]
        
        ax1.bar(range(len(stages)), n_atoms, color='skyblue')
        ax1.set_xticks(range(len(stages)))
        ax1.set_xticklabels(stages, rotation=45)
        ax1.set_ylabel('Número de átomos')
        ax1.set_title('Consistencia del Sistema')
        
        # Gráfico 2: Dimensiones del sistema
        final_result = results[-1]
        dimensions = [final_result['size_x'], final_result['size_y'], final_result['size_z']]
        
        ax2.bar(['X', 'Y', 'Z'], dimensions, color=['red', 'green', 'blue'])
        ax2.set_ylabel('Tamaño (nm)')
        ax2.set_title('Dimensiones del Sistema Final')
        
        # Gráfico 3: Distribución por cadenas
        chain_info = final_result['chain_info']
        chains = list(chain_info.keys())
        chain_atoms = [chain_info[c]['atoms'] for c in chains]
        
        ax3.pie(chain_atoms, labels=[f'Cadena {c}' for c in chains], autopct='%1.1f%%')
        ax3.set_title('Distribución de Átomos por Cadena')
        
        # Gráfico 4: Tabla resumen
        ax4.axis('off')
        summary_data = []
        for result in results:
            summary_data.append([
                result['stage'],
                f"{result['n_atoms']:,}",
                f"{result['n_residues']:,}",
                f"{result['box_volume']:.1f}"
            ])
        
        table = ax4.table(cellText=summary_data,
                         colLabels=['Etapa', 'Átomos', 'Residuos', 'Vol. Caja'],
                         cellLoc='center',
                         loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
        ax4.set_title('Resumen por Etapas', pad=20)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_structural_analysis_page(self, pdf, results):
        """Crea página de análisis estructural"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8.5, 11))
        fig.suptitle('ANÁLISIS ESTRUCTURAL DETALLADO', fontsize=16, fontweight='bold')
        
        final_result = results[-1]
        coords = final_result['full_data']['coordinates']
        
        # Gráfico 1: Distribución espacial (proyección XY)
        ax1.scatter(coords[:, 0], coords[:, 1], alpha=0.1, s=0.1)
        ax1.set_xlabel('X (nm)')
        ax1.set_ylabel('Y (nm)')
        ax1.set_title('Proyección XY de la Proteína')
        ax1.set_aspect('equal')
        
        # Gráfico 2: Distribución espacial (proyección XZ)
        ax2.scatter(coords[:, 0], coords[:, 2], alpha=0.1, s=0.1)
        ax2.set_xlabel('X (nm)')
        ax2.set_ylabel('Z (nm)')
        ax2.set_title('Proyección XZ de la Proteína')
        ax2.set_aspect('equal')
        
        # Gráfico 3: Distribución de coordenadas
        ax3.hist([coords[:, 0], coords[:, 1], coords[:, 2]], 
                bins=50, alpha=0.7, label=['X', 'Y', 'Z'])
        ax3.set_xlabel('Coordenada (nm)')
        ax3.set_ylabel('Frecuencia')
        ax3.set_title('Distribución de Coordenadas')
        ax3.legend()
        
        # Gráfico 4: Información del centro de masa
        ax4.axis('off')
        center = final_result['full_data']['center_of_mass']
        size = final_result['full_data']['system_size']
        
        info_text = f"""
CENTRO DE MASA:
X: {center[0]:.3f} nm
Y: {center[1]:.3f} nm
Z: {center[2]:.3f} nm

DIMENSIONES:
X: {size[0]:.3f} nm
Y: {size[1]:.3f} nm
Z: {size[2]:.3f} nm

ESTADÍSTICAS:
Coord. mín X: {final_result['full_data']['min_coords'][0]:.3f} nm
Coord. máx X: {final_result['full_data']['max_coords'][0]:.3f} nm

Coord. mín Y: {final_result['full_data']['min_coords'][1]:.3f} nm
Coord. máx Y: {final_result['full_data']['max_coords'][1]:.3f} nm

Coord. mín Z: {final_result['full_data']['min_coords'][2]:.3f} nm
Coord. máx Z: {final_result['full_data']['max_coords'][2]:.3f} nm
"""
        
        ax4.text(0.1, 0.9, info_text, transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow"))
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_chain_analysis_page(self, pdf, results):
        """Crea página de análisis por cadenas"""
        fig, axes = plt.subplots(2, 3, figsize=(8.5, 11))
        fig.suptitle('ANÁLISIS POR CADENAS PROTEICAS', fontsize=16, fontweight='bold')
        
        final_result = results[-1]
        chain_info = final_result['chain_info']
        
        # Información de cada cadena
        chains = sorted(chain_info.keys())
        
        for i, chain in enumerate(chains):
            if i < 5:  # Máximo 5 cadenas
                row = i // 3
                col = i % 3
                ax = axes[row, col]
                
                chain_data = chain_info[chain]
                
                # Crear gráfico de barras para la cadena
                categories = ['Átomos', 'Residuos']
                values = [chain_data['atoms'], chain_data['n_residues']]
                
                bars = ax.bar(categories, values, color=['lightblue', 'lightgreen'])
                ax.set_title(f'Cadena {chain}')
                ax.set_ylabel('Cantidad')
                
                # Añadir valores en las barras
                for bar, value in zip(bars, values):
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                           f'{value:,}', ha='center', va='bottom')
        
        # Ocultar ejes no utilizados
        for i in range(len(chains), 6):
            row = i // 3
            col = i % 3
            if row < 2 and col < 3:
                axes[row, col].axis('off')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_evolution_page(self, pdf, results):
        """Crea página de evolución estructural"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8.5, 11))
        fig.suptitle('EVOLUCIÓN ESTRUCTURAL', fontsize=16, fontweight='bold')
        
        stages = [r['stage'] for r in results]
        
        # Evolución del centro de masa
        centers_x = [r['center_x'] for r in results]
        centers_y = [r['center_y'] for r in results]
        centers_z = [r['center_z'] for r in results]
        
        ax1.plot(range(len(stages)), centers_x, 'r-o', label='X')
        ax1.plot(range(len(stages)), centers_y, 'g-o', label='Y')
        ax1.plot(range(len(stages)), centers_z, 'b-o', label='Z')
        ax1.set_xticks(range(len(stages)))
        ax1.set_xticklabels(stages, rotation=45)
        ax1.set_ylabel('Coordenada (nm)')
        ax1.set_title('Evolución del Centro de Masa')
        ax1.legend()
        
        # Evolución del tamaño del sistema
        sizes_x = [r['size_x'] for r in results]
        sizes_y = [r['size_y'] for r in results]
        sizes_z = [r['size_z'] for r in results]
        
        ax2.plot(range(len(stages)), sizes_x, 'r-s', label='X')
        ax2.plot(range(len(stages)), sizes_y, 'g-s', label='Y')
        ax2.plot(range(len(stages)), sizes_z, 'b-s', label='Z')
        ax2.set_xticks(range(len(stages)))
        ax2.set_xticklabels(stages, rotation=45)
        ax2.set_ylabel('Tamaño (nm)')
        ax2.set_title('Evolución del Tamaño del Sistema')
        ax2.legend()
        
        # Cambios relativos
        if len(results) > 1:
            initial = results[0]
            final = results[-1]
            
            # Desplazamiento del centro de masa
            displacement = np.sqrt(
                (final['center_x'] - initial['center_x'])**2 +
                (final['center_y'] - initial['center_y'])**2 +
                (final['center_z'] - initial['center_z'])**2
            )
            
            ax3.axis('off')
            change_text = f"""
CAMBIOS ESTRUCTURALES:

Desplazamiento del centro de masa:
• Total: {displacement:.3f} nm
• ΔX: {final['center_x'] - initial['center_x']:.3f} nm
• ΔY: {final['center_y'] - initial['center_y']:.3f} nm
• ΔZ: {final['center_z'] - initial['center_z']:.3f} nm

Cambios en dimensiones:
• ΔX: {final['size_x'] - initial['size_x']:.3f} nm
• ΔY: {final['size_y'] - initial['size_y']:.3f} nm
• ΔZ: {final['size_z'] - initial['size_z']:.3f} nm

Cambio en volumen:
• Inicial: {initial['box_volume']:.2f} nm³
• Final: {final['box_volume']:.2f} nm³
• Diferencia: {final['box_volume'] - initial['box_volume']:.2f} nm³
"""
            
            ax3.text(0.1, 0.9, change_text, transform=ax3.transAxes, fontsize=10,
                    verticalalignment='top', fontfamily='monospace',
                    bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcyan"))
        
        # Volumen de la caja
        volumes = [r['box_volume'] for r in results]
        ax4.bar(range(len(stages)), volumes, color='gold')
        ax4.set_xticks(range(len(stages)))
        ax4.set_xticklabels(stages, rotation=45)
        ax4.set_ylabel('Volumen (nm³)')
        ax4.set_title('Volumen de la Caja de Simulación')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_coordinates_page(self, pdf, results):
        """Crea página con tabla de coordenadas"""
        fig, ax = plt.subplots(figsize=(8.5, 11))
        fig.suptitle('TABLA DE COORDENADAS DETALLADAS', fontsize=16, fontweight='bold')
        
        ax.axis('off')
        
        # Crear tabla con información detallada
        table_data = []
        headers = ['Etapa', 'Archivo', 'Átomos', 'Residuos', 'Centro X', 'Centro Y', 'Centro Z', 
                  'Tamaño X', 'Tamaño Y', 'Tamaño Z', 'Volumen']
        
        for result in results:
            table_data.append([
                result['stage'],
                result['file'],
                f"{result['n_atoms']:,}",
                f"{result['n_residues']:,}",
                f"{result['center_x']:.3f}",
                f"{result['center_y']:.3f}",
                f"{result['center_z']:.3f}",
                f"{result['size_x']:.3f}",
                f"{result['size_y']:.3f}",
                f"{result['size_z']:.3f}",
                f"{result['box_volume']:.2f}"
            ])
        
        table = ax.table(cellText=table_data, colLabels=headers, cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.2, 2)
        
        # Colorear encabezados
        for i in range(len(headers)):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

def main():
    """Función principal"""
    print("🧬 Generador de Reporte PDF - Receptor GABA")
    print("=" * 50)
    
    # Archivos .gro a analizar
    gro_files = [
        'processed.gro',
        'em/minim.gro',
        'nvt/nvt.gro'
    ]
    
    # Verificar que existen
    existing_files = [f for f in gro_files if Path(f).exists()]
    
    if not existing_files:
        print("❌ No se encontraron archivos .gro")
        return 1
    
    print(f"📂 Archivos encontrados: {len(existing_files)}")
    for f in existing_files:
        print(f"   • {f}")
    
    # Generar PDF
    generator = ProteinStructurePDFGenerator()
    pdf_file = generator.generate_pdf_report(existing_files)
    
    print(f"\n🎉 ¡Reporte PDF generado exitosamente!")
    print(f"📄 Archivo: {pdf_file}")
    
    return 0

if __name__ == "__main__":
    exit(main())
