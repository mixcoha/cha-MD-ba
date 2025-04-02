"""
Módulo para equilibración NVT de sistemas de dinámica molecular
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Dict
from rich.console import Console
from rich.progress import Progress

console = Console()

class NVTEquilibrator:
    """Clase para realizar equilibración NVT de sistemas moleculares"""
    
    def __init__(self, input_gro: str, topol_top: str, mdp_file: Optional[str] = None):
        """
        Inicializa el equilibrador NVT
        
        Args:
            input_gro: Ruta al archivo de coordenadas (.gro)
            topol_top: Ruta al archivo de topología (.top)
            mdp_file: Ruta al archivo de parámetros NVT (.mdp)
        """
        self.input_gro = Path(input_gro)
        self.topol_top = Path(topol_top)
        self.mdp_file = Path(mdp_file) if mdp_file else None
        
    def create_mdp_file(self, output_path: str, title: str = "NVT Equilibration") -> Path:
        """
        Crea un archivo de parámetros para equilibración NVT
        
        Args:
            output_path: Ruta para guardar el archivo .mdp
            title: Título de la simulación
            
        Returns:
            Ruta al archivo .mdp creado
        """
        output_path = Path(output_path)
        
        mdp_content = f"""; Líneas que comienzan con ';' son considerados comentarios
title               = {title}

; Parameters describing what to do, when to stop and what to save
integrator          = md        ; leap-frog integrator
dt                  = 0.002     ; !!! femto second
nsteps              = 50000     ; 100 ps
nstxout             = 500       ; save coordinates every 0.1 ps
nstvout             = 500       ; save velocities every 0.1 ps
nstenergy           = 500       ; save energies every 0.1 ps
nstlog              = 500       ; update log file every 0.1 ps
continuation        = no        ; first dynamics run
constraint_algorithm = lincs     ; holonomic constraints 
constraints         = h-bonds   ; bonds involving H are constrained
lincs_iter          = 1         ; accuracy of LINCS
lincs_order         = 4         ; also related to accuracy

; Nonbonded settings 
cutoff-scheme       = Verlet    ; Buffered neighbor searching
ns_type             = grid      ; search neighboring grid cells
nstlist             = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb            = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr            = EnerPres  ; account for cut-off vdW scheme

; Electrostatics
coulombtype         = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order           = 4         ; cubic interpolation
fourierspacing      = 0.16      ; grid spacing for FFT

; Temperature coupling
tcoupl              = V-rescale ; modified Berendsen thermostat
tc-grps             = Protein Non-Protein ; two coupling groups - more accurate
tau_t               = 0.1       0.1    ; time constant, in ps
ref_t               = 300       300    ; reference temperature, one for each group, in K

; Pressure coupling
pcoupl              = no        ; no pressure coupling in NVT

; Periodic boundary conditions
pbc                 = xyz       ; 3-D PBC

; Velocity generation
gen_vel             = yes       ; assign velocities from Maxwell distribution
gen_temp            = 300       ; temperature for Maxwell distribution
gen_seed            = -1        ; generate a random seed
"""
        
        with open(output_path, 'w') as f:
            f.write(mdp_content)
            
        return output_path
        
    def equilibrate(self, output_dir: str, gpu_ids: Optional[str] = None) -> Dict[str, Path]:
        """
        Realiza la equilibración NVT del sistema
        
        Args:
            output_dir: Directorio para guardar los resultados
            gpu_ids: IDs de GPUs a utilizar (ej. "01" para usar GPUs 0 y 1)
            
        Returns:
            Diccionario con rutas a los archivos generados
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        with Progress() as progress:
            # 1. Crear archivo .mdp si no existe
            if not self.mdp_file:
                task = progress.add_task("[cyan]Creando archivo de parámetros NVT...", total=1)
                self.mdp_file = self.create_mdp_file(output_dir / "nvt.mdp")
                progress.advance(task)
                
            # 2. Generar archivo .tpr
            task = progress.add_task("[cyan]Generando archivo .tpr...", total=1)
            tpr_path = output_dir / "topol.tpr"
            cmd = [
                "gmx", "grompp",
                "-f", str(self.mdp_file),
                "-c", str(self.input_gro),
                "-r", str(self.input_gro),
                "-p", str(self.topol_top),
                "-o", str(tpr_path)
            ]
            subprocess.run(cmd, check=True)
            progress.advance(task)
            
            # 3. Ejecutar equilibración
            task = progress.add_task("[cyan]Ejecutando equilibración NVT...", total=1)
            cmd = ["gmx", "mdrun", "-v", "-s", str(tpr_path)]
            
            if gpu_ids:
                cmd.extend(["-gpu_id", gpu_ids, "-nb", "gpu_cpu", "-tunepme"])
                
            subprocess.run(cmd, check=True)
            progress.advance(task)
            
            # 4. Procesar la estructura final
            task = progress.add_task("[cyan]Procesando estructura final...", total=1)
            confout_gro = output_dir / "confout.gro"
            
            # Centrar molécula
            cmd = ["gmx", "trjconv", "-f", str(confout_gro), "-s", str(tpr_path),
                  "-pbc", "mol", "-center", "-o", str(output_dir / "tmp.gro")]
            subprocess.run(cmd, input=b"1 0\n", check=True)
            
            # Compactar
            cmd = ["gmx", "trjconv", "-f", str(output_dir / "tmp.gro"), "-ur", "compact",
                  "-pbc", "mol", "-o", str(output_dir / "nvt.gro"), "-s", str(tpr_path)]
            subprocess.run(cmd, input=b"0\n", check=True)
            
            # Limpiar archivos temporales
            for tmp_file in output_dir.glob("tmp*.gro"):
                tmp_file.unlink()
            for backup_file in output_dir.glob("#*"):
                backup_file.unlink()
                
            progress.advance(task)
            
        return {
            "tpr": tpr_path,
            "confout": confout_gro,
            "final": output_dir / "nvt.gro",
            "xtc": output_dir / "traj_comp.xtc",
            "edr": output_dir / "ener.edr",
            "log": output_dir / "md.log"
        } 