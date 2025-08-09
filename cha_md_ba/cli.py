"""
Interfaz de línea de comandos para CHA-MD-BA
"""

import click
from pathlib import Path
from rich.console import Console
from .prepare import MDSystemPreparator
from .analysis import MDTrajectoryAnalyzer
from .minimize import EnergyMinimizer
from .nvt import NVTEquilibrator
from .npt import NPTEquilibrator

console = Console()

@click.group()
def main():
    """CHA-MD-BA: Herramientas para automatizar simulaciones de dinámica molecular"""
    pass

@main.command()
@click.argument('pdb_path', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
@click.option('--forcefield', default='amber99sb-ildn', help='Campo de fuerzas a utilizar')
@click.option('--water-model', default='tip3p', help='Modelo de agua a utilizar')
@click.option('--box-type', default='dodecahedron', help='Tipo de caja de simulación')
@click.option('--box-size', type=float, help='Tamaño de la caja de simulación')
@click.option('--no-ions', is_flag=True, help='No agregar iones al sistema')
@click.option('--no-minimize', is_flag=True, help='No realizar minimización de energía')
@click.option('--gpu-ids', help='IDs de GPUs a utilizar para minimización (ej. "01" para usar GPUs 0 y 1)')
def prepare(pdb_path, output_dir, forcefield, water_model, box_type, box_size, no_ions, no_minimize, gpu_ids):
    """Prepara un sistema para simulación de dinámica molecular"""
    preparator = MDSystemPreparator(pdb_path, forcefield, water_model)
    
    try:
        files = preparator.prepare_system(
            output_dir=output_dir,
            box_type=box_type,
            box_size=box_size,
            ions=not no_ions,
            minimize=not no_minimize,
            gpu_ids=gpu_ids
        )
        
        console.print("[green]Sistema preparado exitosamente!")
        console.print("\nEstructura de directorios generada:")
        console.print(f"- Directorio base: {files['base_dir']}")
        console.print(f"- Preparación: {files['preparation_dir']}")
        console.print(f"- Minimización: {files['minimization_dir']}")
        console.print(f"- Equilibración NVT: {files['nvt_dir']}")
        console.print(f"- Equilibración NPT: {files['npt_dir']}")
        
        console.print("\nArchivos generados:")
        for name, path in files.items():
            if path and isinstance(path, Path) and path.is_file():
                console.print(f"- {name}: {path}")
                
    except Exception as e:
        console.print(f"[red]Error al preparar el sistema: {str(e)}")
        raise click.Abort()

@main.command()
@click.argument('gro_path', type=click.Path(exists=True))
@click.argument('top_path', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
@click.option('--mdp-file', type=click.Path(exists=True), help='Archivo de parámetros de minimización (.mdp)')
@click.option('--gpu-ids', help='IDs de GPUs a utilizar (ej. "01" para usar GPUs 0 y 1)')
def minimize(gro_path, top_path, output_dir, mdp_file, gpu_ids):
    """Realiza minimización de energía de un sistema"""
    minimizer = EnergyMinimizer(gro_path, top_path, mdp_file)
    
    try:
        files = minimizer.minimize(output_dir, gpu_ids)
        
        console.print("[green]Minimización completada exitosamente!")
        console.print("\nArchivos generados en el directorio de minimización:")
        for name, path in files.items():
            console.print(f"- {name}: {path}")
            
    except Exception as e:
        console.print(f"[red]Error durante la minimización: {str(e)}")
        raise click.Abort()

@main.command()
@click.argument('gro_path', type=click.Path(exists=True))
@click.argument('top_path', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
@click.option('--mdp-file', type=click.Path(exists=True), help='Archivo de parámetros NVT (.mdp)')
@click.option('--gpu-ids', help='IDs de GPUs a utilizar (ej. "01" para usar GPUs 0 y 1)')
def nvt(gro_path, top_path, output_dir, mdp_file, gpu_ids):
    """Realiza equilibración NVT de un sistema"""
    equilibrator = NVTEquilibrator(gro_path, top_path, mdp_file)
    
    try:
        files = equilibrator.equilibrate(output_dir, gpu_ids)
        
        console.print("[green]Equilibración NVT completada exitosamente!")
        console.print("\nArchivos generados en el directorio de equilibración NVT:")
        for name, path in files.items():
            console.print(f"- {name}: {path}")
            
    except Exception as e:
        console.print(f"[red]Error durante la equilibración NVT: {str(e)}")
        raise click.Abort()

@main.command()
@click.argument('gro_path', type=click.Path(exists=True))
@click.argument('top_path', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
@click.option('--mdp-file', type=click.Path(exists=True), help='Archivo de parámetros NPT (.mdp)')
@click.option('--gpu-ids', help='IDs de GPUs a utilizar (ej. "01" para usar GPUs 0 y 1)')
def npt(gro_path, top_path, output_dir, mdp_file, gpu_ids):
    """Realiza equilibración NPT de un sistema"""
    equilibrator = NPTEquilibrator(gro_path, top_path, mdp_file)
    
    try:
        files = equilibrator.equilibrate(output_dir, gpu_ids)
        
        console.print("[green]Equilibración NPT completada exitosamente!")
        console.print("\nArchivos generados en el directorio de equilibración NPT:")
        for name, path in files.items():
            console.print(f"- {name}: {path}")
            
    except Exception as e:
        console.print(f"[red]Error durante la equilibración NPT: {str(e)}")
        raise click.Abort()

@main.command()
@click.argument('gro_path', type=click.Path(exists=True))
@click.argument('top_path', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
@click.option('--num-runs', default=1, help='Número de corridas de producción a ejecutar')
@click.option('--gpu-ids', help='IDs de GPUs a utilizar (ej. "01" para usar GPUs 0 y 1)')
def production(gro_path, top_path, output_dir, num_runs, gpu_ids):
    """Ejecuta corridas de producción NPT"""
    equilibrator = NPTEquilibrator(gro_path, top_path)
    
    try:
        results = equilibrator.run_production(output_dir, num_runs, gpu_ids)
        
        console.print(f"[green]Producción completada exitosamente! ({num_runs} corridas)")
        console.print("\nEstructura de directorios generada:")
        for i, run_files in enumerate(results, 1):
            console.print(f"\nCorrida {i}:")
            console.print(f"- Directorio: {run_files['run_dir']}")
            for name, path in run_files.items():
                if name != 'run_dir' and path.is_file():
                    console.print(f"- {name}: {path}")
            
    except Exception as e:
        console.print(f"[red]Error durante la producción: {str(e)}")
        raise click.Abort()

@main.command()
@click.argument('trajectory_path', type=click.Path(exists=True))
@click.argument('topology_path', type=click.Path(exists=True))
@click.option('--output-dir', type=click.Path(), help='Directorio para guardar resultados')
@click.option('--selection', default='protein', help='Selección de átomos para análisis')
@click.option('--skip-preprocess', is_flag=True, help='Omitir preprocesamiento de la trayectoria')
def analyze(trajectory_path, topology_path, output_dir, selection, skip_preprocess):
    """Analiza una trayectoria de dinámica molecular"""
    analyzer = MDTrajectoryAnalyzer(trajectory_path, topology_path)
    
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
    try:
        # Preprocesar trayectoria si es necesario
        if not skip_preprocess:
            console.print("[cyan]Preprocesando trayectoria...")
            processed_files = analyzer.preprocess_trajectory(output_dir, selection)
            console.print("[green]Trayectoria preprocesada exitosamente!")
            
            # Usar la trayectoria procesada para el análisis
            analyzer = MDTrajectoryAnalyzer(str(processed_files["processed"]), topology_path)
            
        # Calcular y graficar RMSD
        if output_dir:
            analyzer.plot_rmsd(output_dir / "rmsd.png")
        else:
            analyzer.plot_rmsd()
            
        # Calcular y graficar radio de giro
        if output_dir:
            analyzer.plot_rog(output_dir / "rog.png")
        else:
            analyzer.plot_rog()
            
        console.print("[green]Análisis completado exitosamente!")
        
    except Exception as e:
        console.print(f"[red]Error durante el análisis: {str(e)}")
        raise click.Abort()

if __name__ == '__main__':
    main() 