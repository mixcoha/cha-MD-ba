#!/usr/bin/env python3
"""
Análisis simple de coordenadas GROMACS para extracción de datos
"""

import numpy as np
import pandas as pd
from pathlib import Path

def read_gro_file(gro_file):
    """Lee archivo .gro y extrae información básica"""
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    title = lines[0].strip()
    n_atoms = int(lines[1].strip())
    
    # Dimensiones de la caja (última línea)
    box_line = lines[-1].strip().split()
    box_dims = [float(x) for x in box_line]
    
    # Leer coordenadas
    coordinates = []
    residues = []
    
    for i in range(2, 2 + n_atoms):
        line = lines[i]
        res_num = int(line[0:5])
        res_name = line[5:10].strip()
        atom_name = line[10:15].strip()
        atom_num = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        
        coordinates.append([x, y, z])
        residues.append((res_num, res_name, atom_name))
    
    coords_array = np.array(coordinates)
    
    # Estadísticas básicas
    center_of_mass = np.mean(coords_array, axis=0)
    min_coords = np.min(coords_array, axis=0)
    max_coords = np.max(coords_array, axis=0)
    system_size = max_coords - min_coords
    
    return {
        'title': title,
        'n_atoms': n_atoms,
        'box_dims': box_dims,
        'coordinates': coords_array,
        'residues': residues,
        'center_of_mass': center_of_mass,
        'min_coords': min_coords,
        'max_coords': max_coords,
        'system_size': system_size
    }

def analyze_multiple_gro_files(gro_files):
    """Analiza múltiples archivos .gro y genera reporte"""
    
    results = []
    
    for gro_file in gro_files:
        print(f"📊 Analizando: {gro_file}")
        
        if not Path(gro_file).exists():
            print(f"❌ Error: No se encuentra {gro_file}")
            continue
            
        data = read_gro_file(gro_file)
        
        # Agregar información del archivo
        result = {
            'archivo': Path(gro_file).name,
            'titulo': data['title'],
            'num_atomos': data['n_atoms'],
            'caja_x': data['box_dims'][0] if data['box_dims'] else 0,
            'caja_y': data['box_dims'][1] if len(data['box_dims']) > 1 else 0,
            'caja_z': data['box_dims'][2] if len(data['box_dims']) > 2 else 0,
            'centro_masa_x': data['center_of_mass'][0],
            'centro_masa_y': data['center_of_mass'][1],
            'centro_masa_z': data['center_of_mass'][2],
            'tamano_x': data['system_size'][0],
            'tamano_y': data['system_size'][1],
            'tamano_z': data['system_size'][2],
            'coord_min_x': data['min_coords'][0],
            'coord_min_y': data['min_coords'][1],
            'coord_min_z': data['min_coords'][2],
            'coord_max_x': data['max_coords'][0],
            'coord_max_y': data['max_coords'][1],
            'coord_max_z': data['max_coords'][2]
        }
        
        results.append(result)
    
    # Crear DataFrame y guardar
    df = pd.DataFrame(results)
    
    # Guardar como CSV
    output_file = 'coordinate_analysis.csv'
    df.to_csv(output_file, index=False, float_format='%.3f')
    print(f"📋 Análisis guardado en: {output_file}")
    
    # Mostrar resumen
    print("\n📈 RESUMEN DEL ANÁLISIS:")
    print("=" * 50)
    
    for i, row in df.iterrows():
        print(f"\n🔹 {row['archivo']}:")
        print(f"   • Átomos: {row['num_atomos']:,}")
        print(f"   • Caja: {row['caja_x']:.2f} × {row['caja_y']:.2f} × {row['caja_z']:.2f} nm")
        print(f"   • Centro de masa: ({row['centro_masa_x']:.2f}, {row['centro_masa_y']:.2f}, {row['centro_masa_z']:.2f}) nm")
        print(f"   • Tamaño del sistema: {row['tamano_x']:.2f} × {row['tamano_y']:.2f} × {row['tamano_z']:.2f} nm")
    
    # Comparación entre estructuras si hay múltiples
    if len(df) > 1:
        print("\n📊 COMPARACIÓN ENTRE ESTRUCTURAS:")
        print("=" * 50)
        
        initial = df.iloc[0]
        final = df.iloc[-1]
        
        # Desplazamiento del centro de masa
        delta_x = final['centro_masa_x'] - initial['centro_masa_x']
        delta_y = final['centro_masa_y'] - initial['centro_masa_y']
        delta_z = final['centro_masa_z'] - initial['centro_masa_z']
        displacement = np.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
        
        print(f"🔄 Desplazamiento del centro de masa: {displacement:.3f} nm")
        print(f"   • ΔX: {delta_x:.3f} nm")
        print(f"   • ΔY: {delta_y:.3f} nm") 
        print(f"   • ΔZ: {delta_z:.3f} nm")
        
        # Cambios en el tamaño del sistema
        size_change_x = final['tamano_x'] - initial['tamano_x']
        size_change_y = final['tamano_y'] - initial['tamano_y']
        size_change_z = final['tamano_z'] - initial['tamano_z']
        
        print(f"\n📏 Cambios en el tamaño del sistema:")
        print(f"   • ΔX: {size_change_x:.3f} nm")
        print(f"   • ΔY: {size_change_y:.3f} nm")
        print(f"   • ΔZ: {size_change_z:.3f} nm")
    
    return df

def main():
    print("🧬 Análisis de Coordenadas del Receptor GABA")
    print("=" * 50)
    
    # Archivos a analizar
    gro_files = [
        'processed.gro',      # Estructura inicial procesada
        'em/minim.gro',       # Después de minimización
        'nvt/nvt.gro'         # Después de equilibración NVT
    ]
    
    # Ejecutar análisis
    df = analyze_multiple_gro_files(gro_files)
    
    print(f"\n✅ Análisis completado exitosamente!")
    print(f"📁 Resultados disponibles en: coordinate_analysis.csv")
    
    return df

if __name__ == "__main__":
    results = main()
