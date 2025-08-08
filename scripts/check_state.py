from cha_md_ba.nvt import NVTEquilibrator
from rich.console import Console
from rich.table import Table

console = Console()

def main():
    # Directorio base de la simulación
    base_dir = "simulations/1tim"
    
    # Crear equilibrador NVT
    equilibrator = NVTEquilibrator(
        input_gro="simulations/1tim/2_minimization/minimized.gro",
        topol_top="simulations/1tim/1_preparation/topol.top"
    )
    
    # Verificar estado
    phase, current_fc, detalles = equilibrator.check_simulation_state(base_dir)
    
    # Crear tabla para mostrar el estado
    table = Table(title="📊 Estado de la Simulación")
    table.add_column("Fase", style="cyan")
    table.add_column("Estado", style="green")
    table.add_column("Detalles", style="yellow")
    
    # Preparación
    table.add_row(
        "Preparación",
        "✅ Completada" if detalles["preparation_completa"] else "❌ Pendiente",
        "Archivos de topología generados" if detalles["preparation_completa"] else "Faltan archivos de topología"
    )
    
    # Minimización
    table.add_row(
        "Minimización",
        "✅ Completada" if detalles["minimization_completa"] else "❌ Pendiente",
        "Estructura minimizada generada" if detalles["minimization_completa"] else "Faltan archivos de minimización"
    )
    
    # NVT
    nvt_constantes = detalles["nvt_constantes"]
    nvt_status = "✅ Completada" if len(nvt_constantes) == 5 else "🔄 En progreso"
    nvt_details = f"Constantes completadas: {', '.join(map(str, nvt_constantes))}" if nvt_constantes else "Ninguna constante completada"
    table.add_row("NVT", nvt_status, nvt_details)
    
    # NPT
    table.add_row(
        "NPT",
        "✅ Completada" if detalles["npt_completa"] else "❌ Pendiente",
        "Equilibración NPT completada" if detalles["npt_completa"] else "Falta equilibración NPT"
    )
    
    # Producción
    table.add_row(
        "Producción",
        "✅ Completada" if detalles["production_completa"] else "❌ Pendiente",
        "Simulación de producción completada" if detalles["production_completa"] else "Falta simulación de producción"
    )
    
    # Mostrar tabla
    console.print(table)
    
    # Mostrar siguiente paso
    console.print("\n⏭️ [bold cyan]Siguiente paso:[/bold cyan]")
    if phase == "preparation":
        console.print("Ejecutar preparación del sistema")
    elif phase == "minimization":
        console.print("Ejecutar minimización de energía")
    elif phase == "nvt":
        console.print(f"Ejecutar equilibración NVT con constante de fuerza {current_fc}")
    elif phase == "npt":
        console.print("Ejecutar equilibración NPT")
    elif phase == "production":
        console.print("Ejecutar simulación de producción")
        
if __name__ == "__main__":
    main() 