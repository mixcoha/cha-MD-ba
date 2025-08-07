#!/usr/bin/env python3
"""
Script de análisis específico para proteínas transmembranales
Incluye análisis de orientación, interacciones con lípidos, etc.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances, contacts

def analizar_orientacion_proteina(trajectory_file, structure_file):
    """
    Analiza la orientación de la proteína en la membrana
    """
    print("Analizando orientación de la proteína...")
    
    u = mda.Universe(structure_file, trajectory_file)
    
    # Seleccionar átomos de la proteína
    protein = u.select_atoms("protein")
    
    # Calcular centro de masa de la proteína
    com_protein = protein.center_of_mass()
    
    # Calcular orientación usando momentos de inercia
    moments, axes = protein.moments_of_inertia()
    
    return com_protein, moments, axes

def analizar_interacciones_lipidos(trajectory_file, structure_file):
    """
    Analiza las interacciones entre la proteína y los lípidos
    """
    print("Analizando interacciones proteína-lípidos...")
    
    u = mda.Universe(structure_file, trajectory_file)
    
    # Seleccionar átomos
    protein = u.select_atoms("protein")
    lipids = u.select_atoms("resname POPC DPPC DOPC DMPC")
    
    # Calcular distancias mínimas
    distances_min = []
    for ts in u.trajectory:
        dist = distances.distance_array(protein.positions, lipids.positions)
        distances_min.append(np.min(dist))
    
    return np.array(distances_min)

def analizar_espesor_membrana(trajectory_file, structure_file):
    """
    Analiza el espesor de la membrana
    """
    print("Analizando espesor de membrana...")
    
    u = mda.Universe(structure_file, trajectory_file)
    
    # Seleccionar átomos de fósforo de los lípidos
    phosphates = u.select_atoms("name P")
    
    thickness = []
    for ts in u.trajectory:
        # Calcular espesor basado en la distribución de fósforos
        z_positions = phosphates.positions[:, 2]
        thickness.append(np.max(z_positions) - np.min(z_positions))
    
    return np.array(thickness)

def analizar_curvatura_membrana(trajectory_file, structure_file):
    """
    Analiza la curvatura de la membrana
    """
    print("Analizando curvatura de membrana...")
    
    # Implementar análisis de curvatura
    # Esto requiere algoritmos más complejos
    pass

def graficar_analisis_membrana(time, thickness, lipid_distances, output_prefix):
    """
    Genera gráficas específicas para análisis de membrana
    """
    print("Generando gráficas de análisis de membrana...")
    
    # Gráfica de espesor de membrana
    plt.figure(figsize=(10, 6))
    plt.plot(time, thickness, 'b-', linewidth=1)
    plt.xlabel('Tiempo (ps)')
    plt.ylabel('Espesor de Membrana (nm)')
    plt.title('Espesor de Membrana vs Tiempo')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_membrane_thickness.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Gráfica de distancias proteína-lípidos
    plt.figure(figsize=(10, 6))
    plt.plot(time, lipid_distances, 'r-', linewidth=1)
    plt.xlabel('Tiempo (ps)')
    plt.ylabel('Distancia Mínima Proteína-Lípidos (nm)')
    plt.title('Interacciones Proteína-Lípidos')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_protein_lipid_distances.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    if len(sys.argv) != 3:
        print("Uso: python analyze_membrane.py [archivo_trayectoria.xtc] [archivo_estructura.gro]")
        sys.exit(1)
    
    trajectory_file = sys.argv[1]
    structure_file = sys.argv[2]
    
    print("=== Análisis de Proteína Transmembranal ===")
    print(f"Archivo de trayectoria: {trajectory_file}")
    print(f"Archivo de estructura: {structure_file}")
    print()
    
    try:
        # Cargar trayectoria
        u = mda.Universe(structure_file, trajectory_file)
        time = np.array([ts.time for ts in u.trajectory])
        
        # Realizar análisis
        com_protein, moments, axes = analizar_orientacion_proteina(trajectory_file, structure_file)
        lipid_distances = analizar_interacciones_lipidos(trajectory_file, structure_file)
        thickness = analizar_espesor_membrana(trajectory_file, structure_file)
        
        # Generar gráficas
        output_prefix = trajectory_file.replace('.xtc', '_membrane_analysis')
        graficar_analisis_membrana(time, thickness, lipid_distances, output_prefix)
        
        # Guardar datos
        np.savetxt(f'{output_prefix}_thickness.dat', 
                  np.column_stack((time, thickness)), 
                  header='Tiempo(ps) Espesor(nm)', fmt='%.3f')
        
        np.savetxt(f'{output_prefix}_lipid_distances.dat', 
                  np.column_stack((time, lipid_distances)), 
                  header='Tiempo(ps) Distancia(nm)', fmt='%.3f')
        
        print("\n=== Análisis Completado ===")
        print(f"Espesor promedio de membrana: {np.mean(thickness):.3f} ± {np.std(thickness):.3f} nm")
        print(f"Distancia mínima promedio proteína-lípidos: {np.mean(lipid_distances):.3f} ± {np.std(lipid_distances):.3f} nm")
        
    except Exception as e:
        print(f"Error durante el análisis: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
