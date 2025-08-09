#!/usr/bin/env python3
"""
Wrapper que integra los scripts bash existentes (con archivos MDP optimizados)
con las nuevas funcionalidades Python de CHA-MD-BA
"""

import sys
import os
import subprocess
import argparse
from pathlib import Path

# Añadir el directorio del proyecto al path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cha_md_ba import MDTrajectoryAnalyzer
from rich.console import Console

console = Console()

class CHA_MD_BA_Wrapper:
    """Wrapper que combina scripts bash optimizados con análisis Python avanzado"""
    
    def __init__(self, project_root):
        self.project_root = Path(project_root)
        self.bash_scripts = self.project_root / "bash_version"
        
    def run_bash_script(self, script_name, args):
        """Ejecuta un script bash con los parámetros dados"""
        script_path = self.bash_scripts / script_name
        
        if not script_path.exists():
            console.print(f"[red]❌ Script no encontrado: {script_path}")
            return False
            
        try:
            console.print(f"[cyan]🚀 Ejecutando {script_name}...")
            cmd = [str(script_path)] + args
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            console.print(f"[green]✅ {script_name} completado exitosamente!")
            return True
        except subprocess.CalledProcessError as e:
            console.print(f"[red]❌ Error ejecutando {script_name}: {e}")
            console.print(f"[red]Salida del error: {e.stderr}")
            return False
    
    def minimize_atomico(self, system_name):
        """Ejecuta minimización atómica usando el script bash optimizado"""
        console.print("[bold blue]🔧 Minimización Atómica[/bold blue]")
        console.print("Usando archivos MDP optimizados del script bash...")
        return self.run_bash_script("minimiza_atomico.sh", [system_name])
    
    def nvt_equilibration(self, system_name):
        """Ejecuta equilibración NVT usando el script bash optimizado"""
        console.print("[bold blue]🌡️ Equilibración NVT[/bold blue]")
        console.print("Usando archivos MDP optimizados con restricciones graduales...")
        return self.run_bash_script("nvt2.sh", [system_name])
    
    def npt_equilibration(self, system_name):
        """Ejecuta equilibración NPT usando el script bash optimizado"""
        console.print("[bold blue]💧 Equilibración NPT[/bold blue]")
        console.print("Usando archivos MDP optimizados...")
        return self.run_bash_script("npt.sh", [system_name])
    
    def advanced_analysis(self, system_name, trajectory_file=None, topology_file=None):
        """Ejecuta análisis avanzado usando Python"""
        console.print("[bold blue]📊 Análisis Avanzado[/bold blue]")
        
        # Si no se especifican archivos, usar los del sistema
        if not trajectory_file:
            trajectory_file = f"{system_name}/npt/cell/cell.xtc"
        if not topology_file:
            topology_file = f"{system_name}/npt/cell/topol.tpr"
        
        try:
            # Verificar que existan los archivos
            if not Path(trajectory_file).exists():
                console.print(f"[red]❌ Archivo de trayectoria no encontrado: {trajectory_file}")
                return False
            if not Path(topology_file).exists():
                console.print(f"[red]❌ Archivo de topología no encontrado: {topology_file}")
                return False
            
            # Crear analizador
            analyzer = MDTrajectoryAnalyzer(trajectory_file, topology_file)
            
            # Crear directorio de análisis
            analysis_dir = Path(system_name) / "analysis_advanced"
            analysis_dir.mkdir(exist_ok=True)
            
            console.print("[cyan]📈 Calculando RMSD...")
            times, rmsd = analyzer.calculate_rmsd()
            console.print(f"[green]RMSD promedio: {rmsd.mean():.3f} ± {rmsd.std():.3f} nm")
            
            console.print("[cyan]🔄 Calculando radio de giro...")
            times, rog = analyzer.calculate_rog()
            console.print(f"[green]Radio de giro promedio: {rog.mean():.3f} ± {rog.std():.3f} nm")
            
            # Generar gráficas
            analyzer.plot_rmsd(str(analysis_dir / "rmsd_advanced.png"))
            analyzer.plot_rog(str(analysis_dir / "rog_advanced.png"))
            
            console.print(f"[green]✅ Análisis avanzado completado!")
            console.print(f"[green]📁 Resultados guardados en: {analysis_dir}")
            
            return True
            
        except Exception as e:
            console.print(f"[red]❌ Error en análisis avanzado: {e}")
            return False
    
    def bash_analysis(self, system_name):
        """Ejecuta análisis usando el script bash optimizado"""
        console.print("[bold blue]📊 Análisis con Scripts Bash Optimizados[/bold blue]")
        console.print("Usando configuraciones de gnuplot optimizadas...")
        return self.run_bash_script("analisis.sh", [system_name])
    
    def full_pipeline(self, system_name):
        """Ejecuta el pipeline completo: minimización, equilibración y análisis"""
        console.print("[bold green]🚀 Pipeline Completo CHA-MD-BA[/bold green]")
        
        # 1. Minimización
        if not self.minimize_atomico(system_name):
            return False
        
        # 2. Equilibración NVT
        if not self.nvt_equilibration(system_name):
            return False
        
        # 3. Equilibración NPT
        if not self.npt_equilibration(system_name):
            return False
        
        # 4. Análisis bash (optimizado)
        if not self.bash_analysis(system_name):
            console.print("[yellow]⚠️ Análisis bash falló, continuando con análisis Python...")
        
        # 5. Análisis avanzado Python
        self.advanced_analysis(system_name)
        
        console.print("[bold green]🎉 Pipeline completo finalizado!")
        return True

def main():
    parser = argparse.ArgumentParser(
        description="CHA-MD-BA Wrapper: Integra scripts bash optimizados con análisis Python avanzado"
    )
    
    parser.add_argument('command', choices=[
        'minimize', 'nvt', 'npt', 'analyze_bash', 'analyze_advanced', 'full_pipeline'
    ], help='Comando a ejecutar')
    
    parser.add_argument('system_name', help='Nombre del sistema (carpeta)')
    parser.add_argument('--trajectory', help='Archivo de trayectoria para análisis')
    parser.add_argument('--topology', help='Archivo de topología para análisis')
    
    args = parser.parse_args()
    
    # Crear wrapper
    wrapper = CHA_MD_BA_Wrapper(Path.cwd())
    
    # Ejecutar comando
    if args.command == 'minimize':
        wrapper.minimize_atomico(args.system_name)
    elif args.command == 'nvt':
        wrapper.nvt_equilibration(args.system_name)
    elif args.command == 'npt':
        wrapper.npt_equilibration(args.system_name)
    elif args.command == 'analyze_bash':
        wrapper.bash_analysis(args.system_name)
    elif args.command == 'analyze_advanced':
        wrapper.advanced_analysis(args.system_name, args.trajectory, args.topology)
    elif args.command == 'full_pipeline':
        wrapper.full_pipeline(args.system_name)

if __name__ == "__main__":
    main()
