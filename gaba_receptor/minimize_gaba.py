#!/usr/bin/env python3
"""
Script para minimizar el receptor GABA usando CHA-MD-BA
"""

import sys
import os
sys.path.append('/Users/mixcoha/GMX')

from cha_md_ba.minimize import EnergyMinimizer
from rich.console import Console

console = Console()

def main():
    console.print("[bold green]🧬 Iniciando minimización del receptor GABA[/bold green]")
    
    # Configurar rutas
    input_gro = "processed.gro"
    topol_top = "topol.top"
    output_dir = "em"
    gmx_cmd = "/usr/local/gromacs/bin/gmx_mpi"
    
    # Verificar que los archivos existen
    if not os.path.exists(input_gro):
        console.print(f"[red]Error: No se encuentra {input_gro}[/red]")
        return 1
    
    if not os.path.exists(topol_top):
        console.print(f"[red]Error: No se encuentra {topol_top}[/red]")
        return 1
    
    console.print(f"[cyan]📁 Archivos de entrada:[/cyan]")
    console.print(f"  • Coordenadas: {input_gro}")
    console.print(f"  • Topología: {topol_top}")
    console.print(f"  • Directorio de salida: {output_dir}")
    
    # Crear minimizador
    minimizer = EnergyMinimizer(
        input_gro=input_gro,
        topol_top=topol_top,
        gmx=gmx_cmd
    )
    
    try:
        # Ejecutar minimización
        console.print("[yellow]⚡ Ejecutando minimización de energía...[/yellow]")
        results = minimizer.minimize(output_dir=output_dir)
        
        console.print("[bold green]✅ Minimización completada exitosamente![/bold green]")
        console.print("[cyan]📊 Archivos generados:[/cyan]")
        for key, path in results.items():
            console.print(f"  • {key}: {path}")
            
        console.print("[yellow]💡 Próximos pasos sugeridos:[/yellow]")
        console.print("  1. Revisar el archivo de log para verificar convergencia")
        console.print("  2. Analizar la energía con: gmx energy -f em/ener.edr")
        console.print("  3. Proceder con equilibración NVT")
        
    except Exception as e:
        console.print(f"[red]❌ Error durante la minimización: {e}[/red]")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
