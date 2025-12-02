#!/usr/bin/env python3
"""
Script directo para equilibración NVT del receptor GABA
"""

import os
import subprocess
from pathlib import Path

def main():
    print("🌡️ Iniciando equilibración NVT del receptor GABA")
    
    # Configurar rutas
    input_gro = "em/minim.gro"
    topol_top = "topol.top"
    output_dir = "nvt"
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
    
    # Crear archivo MDP para NVT usando los parámetros optimizados de CHA-MD-BA
    mdp_content = """; Equilibración NVT para receptor GABA
; Parámetros optimizados de CHA-MD-BA
title = NVT equilibration

define = -DPOSRES      ; Activar restricciones posicionales

; RUN CONTROL PARAMETERS
integrator = md        ; leap-frog integrator
nsteps = 50000         ; 2 * 50000 = 100 ps
dt = 0.002             ; 2 fs

; OUTPUT CONTROL OPTIONS
nstxout = 5000         ; save coordinates every 10.0 ps
nstvout = 5000         ; save velocities every 10.0 ps
nstenergy = 5000       ; save energies every 10.0 ps
nstlog = 5000          ; update log file every 10.0 ps

; NEIGHBORSEARCHING PARAMETERS
cutoff-scheme = Verlet
nstlist = 10           ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb = 1.0         ; short-range electrostatic cutoff (in nm)
rvdw = 1.0             ; short-range van der Waals cutoff (in nm)

; ELECTROSTATICS
coulombtype = PME      ; Particle Mesh Ewald for long-range electrostatics
pme_order = 4          ; cubic interpolation
fourierspacing = 0.16  ; grid spacing for FFT

; TEMPERATURE COUPLING  
tcoupl = V-rescale     ; modified Berendsen thermostat
tc-grps = Protein Non-Protein   ; two coupling groups - more accurate
tau_t = 0.1    0.1     ; time constant, in ps
ref_t = 300    300     ; reference temperature, one for each group, in K

; PRESSURE COUPLING
pcoupl = no            ; no pressure coupling in NVT

; PERIODIC BOUNDARY CONDITIONS
pbc = xyz              ; 3-D PBC

; DISPERSION CORRECTION
DispCorr = EnerPres    ; account for cut-off vdW scheme

; VELOCITY GENERATION
gen_vel = yes          ; assign velocities from Maxwell distribution
gen_temp = 300         ; temperature for Maxwell distribution
gen_seed = -1          ; generate a random seed

; CONSTRAINTS
constraints = h-bonds  ; bonds involving H are constrained
constraint_algorithm = LINCS    ; holonomic constraints 
lincs_iter = 1         ; accuracy of LINCS
lincs_order = 4        ; also related to accuracy"""
    
    mdp_file = os.path.join(output_dir, "nvt.mdp")
    with open(mdp_file, 'w') as f:
        f.write(mdp_content)
    
    print(f"📁 Archivos de entrada:")
    print(f"  • Estructura minimizada: {input_gro}")
    print(f"  • Topología: {topol_top}")
    print(f"  • Parámetros NVT: {mdp_file}")
    print(f"  • Directorio de salida: {output_dir}")
    
    try:
        # 1. Generar archivo .tpr
        print("⚡ Generando archivo .tpr para NVT...")
        tpr_path = os.path.join(output_dir, "nvt.tpr")
        cmd = [
            gmx_cmd, "grompp",
            "-f", mdp_file,
            "-c", input_gro,
            "-r", input_gro,  # archivo de referencia para restricciones
            "-p", topol_top,
            "-o", tpr_path,
            "-maxwarn", "1"  # Ignorar advertencias menores
        ]
        subprocess.run(cmd, check=True)
        
        # 2. Ejecutar equilibración NVT
        print("🚀 Ejecutando equilibración NVT (100 ps)...")
        cmd = [
            gmx_cmd, "mdrun", "-v",
            "-s", tpr_path,
            "-deffnm", os.path.join(output_dir, "nvt")
        ]
        subprocess.run(cmd, check=True)
        
        print("✅ Equilibración NVT completada exitosamente!")
        print("📊 Archivos generados:")
        print(f"  • nvt.gro: {os.path.join(output_dir, 'nvt.gro')}")
        print(f"  • nvt.edr: {os.path.join(output_dir, 'nvt.edr')}")
        print(f"  • nvt.log: {os.path.join(output_dir, 'nvt.log')}")
        
        print("💡 Próximos pasos:")
        print("  1. Revisar temperatura: gmx energy -f nvt/nvt.edr -o temperature.xvg")
        print("  2. Proceder con equilibración NPT")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error durante la equilibración NVT: {e}")
        return 1
    except Exception as e:
        print(f"❌ Error inesperado: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
