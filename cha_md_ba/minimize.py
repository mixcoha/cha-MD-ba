"""
Módulo para minimización de energía de sistemas de dinámica molecular
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Dict
from rich.console import Console
from rich.progress import Progress

console = Console()

class EnergyMinimizer:
    """Clase para realizar minimización de energía de sistemas moleculares"""
    
    def __init__(self, input_gro: str, topol_top: str, mdp_file: Optional[str] = None):
        """
        Inicializa el minimizador de energía
        
        Args:
            input_gro: Ruta al archivo de coordenadas (.gro)
            topol_top: Ruta al archivo de topología (.top)
            mdp_file: Ruta al archivo de parámetros de minimización (.mdp)
        """
        self.input_gro = Path(input_gro)
        self.topol_top = Path(topol_top)
        self.mdp_file = Path(mdp_file) if mdp_file else None
        
    def create_mdp_file(self, output_path: str, title: str = "Energy Minimization") -> Path:
        """
        Crea un archivo de parámetros para minimización
        
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
integrator          = steep      ; Algorithm options i.e. steep = steepest descent minimization 
emtol               = 1000.0     ; Stop minimization when the energy changes by less than emtol kJ/mol.
emstep              = 0.01       ; Energy step size
nsteps              = 1000       ; Maximum number of (minimization) steps to perform
nstenergy           = 10         ; Write energies to disk every nstenergy steps
nstxtcout           = 10         ; Write coordinates to disk every nstxtcout steps
xtc_grps            = Protein    ; Which coordinate group(s) to write to disk
energygrps          = Protein    ; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist             = 1          ; Frequency to update the neighbor list and long range forces
ns_type             = grid       ; Method to determine neighbor list (simple, grid)
rlist               = 1.0        ; Cut-off for making neighbor list (short range forces)
coulombtype         = cut-off    ; Treatment of long range electrostatic interactions
rcoulomb            = 1.0        ; long range electrostatic cut-off
rvdw                = 1.0        ; long range Van der Waals cut-off
constraints         = none       ; Bond types to replace by constraints
pbc                 = xyz        ; Periodic Boundary Conditions (yes/no)
"""
        
        with open(output_path, 'w') as f:
            f.write(mdp_content)
            
        return output_path
        
    def minimize(self, output_dir: str, gpu_ids: Optional[str] = None) -> Dict[str, Path]:
        """
        Realiza la minimización de energía del sistema
        
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
                task = progress.add_task("[cyan]Creando archivo de parámetros...", total=1)
                self.mdp_file = self.create_mdp_file(output_dir / "em.mdp")
                progress.advance(task)
                
            # 2. Generar archivo .tpr
            task = progress.add_task("[cyan]Generando archivo .tpr...", total=1)
            tpr_path = output_dir / "topol.tpr"
            cmd = [
                "gmx", "grompp",
                "-f", str(self.mdp_file),
                "-c", str(self.input_gro),
                "-p", str(self.topol_top),
                "-o", str(tpr_path)
            ]
            subprocess.run(cmd, check=True)
            progress.advance(task)
            
            # 3. Ejecutar minimización
            task = progress.add_task("[cyan]Ejecutando minimización...", total=1)
            cmd = ["gmx", "mdrun", "-v", "-s", str(tpr_path)]
            
            if gpu_ids:
                cmd.extend(["-gpu_id", gpu_ids, "-nb", "gpu_cpu", "-tunepme"])
                
            subprocess.run(cmd, check=True)
            progress.advance(task)
            
            # 4. Centrar y procesar la estructura final
            task = progress.add_task("[cyan]Procesando estructura final...", total=1)
            confout_gro = output_dir / "confout.gro"
            
            # Centrar molécula
            cmd = ["gmx", "trjconv", "-f", str(confout_gro), "-s", str(tpr_path),
                  "-pbc", "mol", "-center", "-o", str(output_dir / "tmp.gro")]
            subprocess.run(cmd, input=b"1 0\n", check=True)
            
            # Compactar
            cmd = ["gmx", "trjconv", "-f", str(output_dir / "tmp.gro"), "-ur", "compact",
                  "-pbc", "mol", "-o", str(output_dir / "tmp2.gro"), "-s", str(tpr_path)]
            subprocess.run(cmd, input=b"0\n", check=True)
            
            # Cluster
            cmd = ["gmx", "trjconv", "-f", str(output_dir / "tmp2.gro"), "-pbc", "cluster",
                  "-o", str(output_dir / "tmp.gro"), "-s", str(tpr_path)]
            subprocess.run(cmd, input=b"1 0\n", check=True)
            
            # Compactar y centrar
            cmd = ["gmx", "trjconv", "-f", str(output_dir / "tmp.gro"), "-ur", "compact",
                  "-pbc", "mol", "-center", "-o", str(output_dir / "tmp2.gro"), "-s", str(tpr_path)]
            subprocess.run(cmd, input=b"1 0\n", check=True)
            
            # Ajustar rotación y traslación
            final_pdb = output_dir / "minimized.pdb"
            cmd = ["gmx", "trjconv", "-f", str(output_dir / "tmp2.gro"), "-fit", "rot+trans",
                  "-o", str(final_pdb), "-s", str(tpr_path)]
            subprocess.run(cmd, input=b"4 0\n", check=True)
            
            # Limpiar archivos temporales
            for tmp_file in output_dir.glob("tmp*.gro"):
                tmp_file.unlink()
            for backup_file in output_dir.glob("#*"):
                backup_file.unlink()
                
            progress.advance(task)
            
        return {
            "tpr": tpr_path,
            "confout": confout_gro,
            "final": final_pdb
        } 