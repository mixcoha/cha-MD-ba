"""
Módulo para equilibración NPT y producción de trayectorias de dinámica molecular
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Dict, List
from rich.console import Console
from rich.progress import Progress

console = Console()

class NPTEquilibrator:
    """Clase para realizar equilibración NPT y producción de trayectorias"""
    
    def __init__(self, input_gro: str, topol_top: str, mdp_file: Optional[str] = None):
        """
        Inicializa el equilibrador NPT
        
        Args:
            input_gro: Ruta al archivo de coordenadas (.gro)
            topol_top: Ruta al archivo de topología (.top)
            mdp_file: Ruta al archivo de parámetros NPT (.mdp)
        """
        self.input_gro = Path(input_gro)
        self.topol_top = Path(topol_top)
        self.mdp_file = Path(mdp_file) if mdp_file else None
        
    def create_mdp_file(self, output_path: str, title: str = "NPT Equilibration", 
                       is_production: bool = False) -> Path:
        """
        Crea un archivo de parámetros para equilibración NPT o producción
        
        Args:
            output_path: Ruta para guardar el archivo .mdp
            title: Título de la simulación
            is_production: Si es True, configura parámetros para producción
            
        Returns:
            Ruta al archivo .mdp creado
        """
        output_path = Path(output_path)
        
        mdp_content = f"""; Líneas que comienzan con ';' son considerados comentarios
title               = {title}

; Parameters describing what to do, when to stop and what to save
integrator          = md        ; leap-frog integrator
dt                  = 0.002     ; !!! femto second
nsteps              = {'50000' if not is_production else '250000'} ; {'100 ps' if not is_production else '500 ps'}
nstxout             = {'500' if not is_production else '5000'}    ; save coordinates every {'0.1 ps' if not is_production else '1 ps'}
nstvout             = {'500' if not is_production else '5000'}    ; save velocities every {'0.1 ps' if not is_production else '1 ps'}
nstenergy           = {'500' if not is_production else '5000'}    ; save energies every {'0.1 ps' if not is_production else '1 ps'}
nstlog              = {'500' if not is_production else '5000'}    ; update log file every {'0.1 ps' if not is_production else '1 ps'}
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
pcoupl              = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype          = isotropic             ; uniform scaling of box vectors
tau_p               = 2.0                   ; time constant, in ps
ref_p               = 1.0                   ; reference pressure, in bar
compressibility     = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com                   ; absolute position of the center of mass

; Periodic boundary conditions
pbc                 = xyz       ; 3-D PBC

; Velocity generation
gen_vel             = {'yes' if not is_production else 'no'} ; assign velocities from Maxwell distribution
gen_temp            = 300       ; temperature for Maxwell distribution
gen_seed            = -1        ; generate a random seed
"""
        
        with open(output_path, 'w') as f:
            f.write(mdp_content)
            
        return output_path
        
    def equilibrate(self, output_dir: str, gpu_ids: Optional[str] = None) -> Dict[str, Path]:
        """
        Realiza la equilibración NPT del sistema
        
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
                task = progress.add_task("[cyan]Creando archivo de parámetros NPT...", total=1)
                self.mdp_file = self.create_mdp_file(output_dir / "npt.mdp")
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
            task = progress.add_task("[cyan]Ejecutando equilibración NPT...", total=1)
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
                  "-pbc", "mol", "-o", str(output_dir / "npt.gro"), "-s", str(tpr_path)]
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
            "final": output_dir / "npt.gro",
            "xtc": output_dir / "traj_comp.xtc",
            "edr": output_dir / "ener.edr",
            "log": output_dir / "md.log"
        }
        
    def run_production(self, output_dir: str, num_runs: int = 1, 
                      gpu_ids: Optional[str] = None) -> List[Dict[str, Path]]:
        """
        Ejecuta múltiples corridas de producción NPT
        
        Args:
            output_dir: Directorio base para guardar los resultados
            num_runs: Número de corridas a ejecutar
            gpu_ids: IDs de GPUs a utilizar (ej. "01" para usar GPUs 0 y 1)
            
        Returns:
            Lista de diccionarios con rutas a los archivos generados en cada corrida
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = []
        
        for run_num in range(1, num_runs + 1):
            run_dir = output_dir / f"run_{run_num:03d}"
            run_dir.mkdir(exist_ok=True)
            
            with Progress() as progress:
                # 1. Crear archivo .mdp para producción
                task = progress.add_task(f"[cyan]Preparando corrida {run_num}...", total=1)
                mdp_path = self.create_mdp_file(run_dir / "md.mdp", 
                                              title=f"Production Run {run_num}",
                                              is_production=True)
                progress.advance(task)
                
                # 2. Generar archivo .tpr
                task = progress.add_task("[cyan]Generando archivo .tpr...", total=1)
                tpr_path = run_dir / "topol.tpr"
                
                # Usar la estructura final de la corrida anterior como entrada
                input_gro = self.input_gro if run_num == 1 else results[-1]["final"]
                
                cmd = [
                    "gmx", "grompp",
                    "-f", str(mdp_path),
                    "-c", str(input_gro),
                    "-p", str(self.topol_top),
                    "-o", str(tpr_path)
                ]
                subprocess.run(cmd, check=True)
                progress.advance(task)
                
                # 3. Ejecutar producción
                task = progress.add_task("[cyan]Ejecutando producción...", total=1)
                cmd = ["gmx", "mdrun", "-v", "-s", str(tpr_path)]
                
                if gpu_ids:
                    cmd.extend(["-gpu_id", gpu_ids, "-nb", "gpu_cpu", "-tunepme"])
                    
                subprocess.run(cmd, check=True)
                progress.advance(task)
                
                # 4. Procesar la estructura final
                task = progress.add_task("[cyan]Procesando estructura final...", total=1)
                confout_gro = run_dir / "confout.gro"
                
                # Centrar molécula
                cmd = ["gmx", "trjconv", "-f", str(confout_gro), "-s", str(tpr_path),
                      "-pbc", "mol", "-center", "-o", str(run_dir / "tmp.gro")]
                subprocess.run(cmd, input=b"1 0\n", check=True)
                
                # Compactar
                cmd = ["gmx", "trjconv", "-f", str(run_dir / "tmp.gro"), "-ur", "compact",
                      "-pbc", "mol", "-o", str(run_dir / "final.gro"), "-s", str(tpr_path)]
                subprocess.run(cmd, input=b"0\n", check=True)
                
                # Limpiar archivos temporales
                for tmp_file in run_dir.glob("tmp*.gro"):
                    tmp_file.unlink()
                for backup_file in run_dir.glob("#*"):
                    backup_file.unlink()
                    
                progress.advance(task)
                
            results.append({
                "run_dir": run_dir,
                "tpr": tpr_path,
                "confout": confout_gro,
                "final": run_dir / "final.gro",
                "xtc": run_dir / "traj_comp.xtc",
                "edr": run_dir / "ener.edr",
                "log": run_dir / "md.log"
            })
            
        return results 