#!/usr/bin/env python3
"""
Script de análisis básico para trayectorias de GROMACS
Uso: python analyze_trajectory.py [archivo_trayectoria.xtc] [archivo_estructura.gro]
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, diffusion

# Configurar matplotlib para español
rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

def calcular_rmsd(trajectory_file, structure_file, selection='backbone'):
    """
    Calcula el RMSD de la trayectoria
    """
    print(f"Calculando RMSD para selección: {selection}")
    
    # Cargar la trayectoria
    u = mda.Universe(structure_file, trajectory_file)
    
    # Seleccionar átomos de referencia
    protein = u.select_atoms(selection)
    
    # Calcular RMSD
    R = rms.RMSD(u, protein, select='backbone', ref_frame=0)
    R.run()
    
    return R.rmsd

def calcular_rmsf(trajectory_file, structure_file, selection='backbone'):
    """
    Calcula el RMSF por residuo
    """
    print("Calculando RMSF por residuo...")
    
    # Cargar la trayectoria
    u = mda.Universe(structure_file, trajectory_file)
    
    # Seleccionar átomos
    protein = u.select_atoms(selection)
    
    # Calcular RMSF
    R = rms.RMSF(protein, ref_frame=0)
    R.run()
    
    return R.rmsf, protein.resids

def calcular_radio_giro(trajectory_file, structure_file, selection='protein'):
    """
    Calcula el radio de giro
    """
    print("Calculando radio de giro...")
    
    # Cargar la trayectoria
    u = mda.Universe(structure_file, trajectory_file)
    
    # Seleccionar átomos
    protein = u.select_atoms(selection)
    
    # Calcular radio de giro
    rg = []
    for ts in u.trajectory:
        rg.append(protein.radius_of_gyration())
    
    return np.array(rg)

def graficar_resultados(time, rmsd, rmsf, resids, rg, output_prefix):
    """
    Genera gráficas de los análisis
    """
    print("Generando gráficas...")
    
    # Gráfica de RMSD
    plt.figure(figsize=(10, 6))
    plt.plot(time, rmsd[:, 2], 'b-', linewidth=1)
    plt.xlabel('Tiempo (ps)')
    plt.ylabel('RMSD (nm)')
    plt.title('RMSD del Backbone')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_rmsd.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Gráfica de RMSF
    plt.figure(figsize=(12, 6))
    plt.plot(resids, rmsf, 'r-', linewidth=1)
    plt.xlabel('Número de Residuo')
    plt.ylabel('RMSF (nm)')
    plt.title('RMSF por Residuo')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_rmsf.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Gráfica de Radio de Giro
    plt.figure(figsize=(10, 6))
    plt.plot(time, rg, 'g-', linewidth=1)
    plt.xlabel('Tiempo (ps)')
    plt.ylabel('Radio de Giro (nm)')
    plt.title('Radio de Giro de la Proteína')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_rg.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Gráficas guardadas como: {output_prefix}_*.png")

def main():
    if len(sys.argv) != 3:
        print("Uso: python analyze_trajectory.py [archivo_trayectoria.xtc] [archivo_estructura.gro]")
        sys.exit(1)
    
    trajectory_file = sys.argv[1]
    structure_file = sys.argv[2]
    
    print("=== Análisis de Trayectoria de GROMACS ===")
    print(f"Archivo de trayectoria: {trajectory_file}")
    print(f"Archivo de estructura: {structure_file}")
    print()
    
    try:
        # Cargar la trayectoria para obtener información de tiempo
        u = mda.Universe(structure_file, trajectory_file)
        time = np.array([ts.time for ts in u.trajectory])
        
        # Realizar análisis
        rmsd = calcular_rmsd(trajectory_file, structure_file)
        rmsf, resids = calcular_rmsf(trajectory_file, structure_file)
        rg = calcular_radio_giro(trajectory_file, structure_file)
        
        # Generar gráficas
        output_prefix = trajectory_file.replace('.xtc', '_analysis')
        graficar_resultados(time, rmsd, rmsf, resids, rg, output_prefix)
        
        # Guardar datos numéricos
        np.savetxt(f'{output_prefix}_rmsd.dat', 
                  np.column_stack((time, rmsd[:, 2])), 
                  header='Tiempo(ps) RMSD(nm)', fmt='%.3f')
        
        np.savetxt(f'{output_prefix}_rmsf.dat', 
                  np.column_stack((resids, rmsf)), 
                  header='Residuo RMSF(nm)', fmt='%.3f')
        
        np.savetxt(f'{output_prefix}_rg.dat', 
                  np.column_stack((time, rg)), 
                  header='Tiempo(ps) Rg(nm)', fmt='%.3f')
        
        print("\n=== Análisis Completado ===")
        print(f"RMSD promedio: {np.mean(rmsd[:, 2]):.3f} ± {np.std(rmsd[:, 2]):.3f} nm")
        print(f"Radio de giro promedio: {np.mean(rg):.3f} ± {np.std(rg):.3f} nm")
        print(f"RMSF máximo: {np.max(rmsf):.3f} nm (residuo {resids[np.argmax(rmsf)]})")
        
    except Exception as e:
        print(f"Error durante el análisis: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
