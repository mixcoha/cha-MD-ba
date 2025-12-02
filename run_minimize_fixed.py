#!/usr/bin/env python3
"""
Script directo para minimizar el receptor GABA
"""

import os
import sys
import subprocess
from pathlib import Path

def main():
    print("🧬 Iniciando minimización del receptor GABA")
    
    # Configurar rutas
    input_gro = "processed.gro"
    topol_top = "topol.top"
    output_dir = "em"
    gmx_cmd = "/usr/local/gromacs/bin/gmx_mpi"
    
    # Verificar que los archivos existen
    if not os.path.exists(input_gro):
        print(f"❌ Error: No se encuentra {input_gro}")
        return 1
    
    if not os.path.exists(topol_top):
        print(f"❌ Error: No se encuentra {topol_top}")
        return 1
    
    # Crear directorio de salida
    os.makedirs(output_dir, exist_ok=True)
    
    # Crear archivo MDP usando los parámetros optimizados de CHA-MD-BA
    mdp_content = """; Líneas que comienzan con ';' son considerados comentarios
title               = Energy Minimization

; Parameters describing what to do, when to stop and what to save
integrator          = steep     ; steepest descent minimization
emtol               = 1000.0    ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep              = 0.01      ; Energy step size
nsteps              = 50000     ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist             = 10        ; Frequency to update the neighbor list and long range forces
cutoff-scheme       = Verlet    ; Buffered neighbor searching
coulombtype         = PME       ; Treatment of long range electrostatic interactions
rcoulomb            = 1.0       ; Short-range electrostatic cut-off
rvdw                = 1.0       ; Short-range Van der Waals cut-off
pbc                 = xyz       ; Periodic Boundary Conditions in all 3 dimensions
"""
    
    mdp_file = os.path.join(output_dir, "min.mdp")
    with open(mdp_file, 'w') as f:
        f.write(mdp_content)
    
    print(f"📁 Archivos de entrada:")
    print(f"  • Coordenadas: {input_gro}")
    print(f"  • Topología: {topol_top}")
    print(f"  • Parámetros: {mdp_file}")
    print(f"  • Directorio de salida: {output_dir}")
    
    try:
        # 1. Generar archivo .tpr
        print("⚡ Generando archivo .tpr...")
        tpr_path = os.path.join(output_dir, "topol.tpr")
        cmd = [
            gmx_cmd, "grompp",
            "-f", mdp_file,
            "-c", input_gro,
            "-p", topol_top,
            "-o", tpr_path,
            "-maxwarn", "1"  # Ignorar la advertencia de carga neta
        ]
        subprocess.run(cmd, check=True)
        
        # 2. Ejecutar minimización
        print("🚀 Ejecutando minimización...")
        cmd = [
            gmx_cmd, "mdrun", "-v",
            "-s", tpr_path,
            "-deffnm", os.path.join(output_dir, "minim")
        ]
        subprocess.run(cmd, check=True)
        
        print("✅ Minimización completada exitosamente!")
        print("📊 Archivos generados:")
        print(f"  • minim.gro: {os.path.join(output_dir, 'minim.gro')}")
        print(f"  • minim.edr: {os.path.join(output_dir, 'minim.edr')}")
        print(f"  • minim.log: {os.path.join(output_dir, 'minim.log')}")
        
        print("💡 Próximos pasos:")
        print("  1. Revisar convergencia: gmx energy -f em/minim.edr -o potential.xvg")
        print("  2. Proceder con equilibración NVT")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error durante la minimización: {e}")
        return 1
    except Exception as e:
        print(f"❌ Error inesperado: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
