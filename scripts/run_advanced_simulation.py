#!/usr/bin/env python3
"""
Script avanzado que integra CHA-MD-BA con el proyecto GROMACS
Usa los scripts más recientes manteniendo la generación de archivos MDP intacta
"""

import sys
import argparse
from pathlib import Path

# Añadir el directorio del proyecto al path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cha_md_ba import (
    MDSystemPreparator, 
    EnergyMinimizer, 
    NVTEquilibrator, 
    NPTEquilibrator,
    MDTrajectoryAnalyzer
)

def main():
    parser = argparse.ArgumentParser(
        description="CHA-MD-BA: Pipeline avanzado de simulación de dinámica molecular"
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Comandos disponibles')
    
    # Comando prepare
    prepare_parser = subparsers.add_parser('prepare', help='Preparar sistema para simulación')
    prepare_parser.add_argument('pdb_file', help='Archivo PDB de entrada')
    prepare_parser.add_argument('output_dir', help='Directorio de salida')
    prepare_parser.add_argument('--forcefield', default='amber99sb-ildn', help='Campo de fuerzas')
    prepare_parser.add_argument('--water-model', default='tip3p', help='Modelo de agua')
    prepare_parser.add_argument('--box-type', default='dodecahedron', help='Tipo de caja')
    prepare_parser.add_argument('--gpu-ids', help='IDs de GPUs (ej. "01")')
    
    # Comando minimize
    minimize_parser = subparsers.add_parser('minimize', help='Minimización de energía')
    minimize_parser.add_argument('gro_file', help='Archivo .gro de entrada')
    minimize_parser.add_argument('top_file', help='Archivo .top de topología')
    minimize_parser.add_argument('output_dir', help='Directorio de salida')
    minimize_parser.add_argument('--mdp-file', help='Archivo .mdp personalizado')
    minimize_parser.add_argument('--gpu-ids', help='IDs de GPUs (ej. "01")')
    
    # Comando nvt
    nvt_parser = subparsers.add_parser('nvt', help='Equilibración NVT')
    nvt_parser.add_argument('gro_file', help='Archivo .gro de entrada')
    nvt_parser.add_argument('top_file', help='Archivo .top de topología')
    nvt_parser.add_argument('output_dir', help='Directorio de salida')
    nvt_parser.add_argument('--mdp-file', help='Archivo .mdp personalizado')
    nvt_parser.add_argument('--gpu-ids', help='IDs de GPUs (ej. "01")')
    
    # Comando npt
    npt_parser = subparsers.add_parser('npt', help='Equilibración NPT')
    npt_parser.add_argument('gro_file', help='Archivo .gro de entrada')
    npt_parser.add_argument('top_file', help='Archivo .top de topología')
    npt_parser.add_argument('output_dir', help='Directorio de salida')
    npt_parser.add_argument('--mdp-file', help='Archivo .mdp personalizado')
    npt_parser.add_argument('--gpu-ids', help='IDs de GPUs (ej. "01")')
    
    # Comando analyze
    analyze_parser = subparsers.add_parser('analyze', help='Análisis de trayectoria')
    analyze_parser.add_argument('trajectory_file', help='Archivo de trayectoria (.xtc)')
    analyze_parser.add_argument('topology_file', help='Archivo de topología (.tpr)')
    analyze_parser.add_argument('output_dir', help='Directorio de salida')
    
    args = parser.parse_args()
    
    if args.command == 'prepare':
        print(f"🚀 Preparando sistema desde {args.pdb_file}")
        preparator = MDSystemPreparator(
            args.pdb_file, 
            forcefield=args.forcefield,
            water_model=args.water_model
        )
        
        files = preparator.prepare_system(
            output_dir=args.output_dir,
            box_type=args.box_type,
            gpu_ids=args.gpu_ids
        )
        
        print("✅ Sistema preparado exitosamente!")
        print(f"📁 Archivos generados en: {files['base_dir']}")
        
    elif args.command == 'minimize':
        print(f"⚡ Minimizando energía...")
        minimizer = EnergyMinimizer(
            args.gro_file,
            args.top_file,
            args.mdp_file
        )
        
        files = minimizer.minimize(args.output_dir, args.gpu_ids)
        print("✅ Minimización completada!")
        print(f"📁 Archivos en: {args.output_dir}")
        
    elif args.command == 'nvt':
        print(f"🌡️ Equilibración NVT...")
        equilibrator = NVTEquilibrator(
            args.gro_file,
            args.top_file,
            args.mdp_file
        )
        
        files = equilibrator.equilibrate(args.output_dir, args.gpu_ids)
        print("✅ Equilibración NVT completada!")
        print(f"📁 Archivos en: {args.output_dir}")
        
    elif args.command == 'npt':
        print(f"💧 Equilibración NPT...")
        equilibrator = NPTEquilibrator(
            args.gro_file,
            args.top_file,
            args.mdp_file
        )
        
        files = equilibrator.equilibrate(args.output_dir, args.gpu_ids)
        print("✅ Equilibración NPT completada!")
        print(f"📁 Archivos en: {args.output_dir}")
        
    elif args.command == 'analyze':
        print(f"📊 Analizando trayectoria...")
        analyzer = MDTrajectoryAnalyzer(
            args.trajectory_file,
            args.topology_file
        )
        
        # Calcular RMSD
        times, rmsd = analyzer.calculate_rmsd()
        print(f"📈 RMSD promedio: {rmsd.mean():.3f} ± {rmsd.std():.3f} nm")
        
        # Calcular radio de giro
        times, rog = analyzer.calculate_rog()
        print(f"🔄 Radio de giro promedio: {rog.mean():.3f} ± {rog.std():.3f} nm")
        
        # Generar gráficas
        analyzer.plot_rmsd(f"{args.output_dir}/rmsd_analysis.png")
        analyzer.plot_rog(f"{args.output_dir}/rog_analysis.png")
        
        print("✅ Análisis completado!")
        print(f"📁 Gráficas guardadas en: {args.output_dir}")
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
