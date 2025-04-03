"""
Módulo para minimización de energía de sistemas de dinámica molecular
"""

import os
import subprocess
import shutil
from pathlib import Path
from typing import Optional, Dict
from rich.console import Console

console = Console()

class EnergyMinimizer:
    """Clase para realizar minimización de energía de sistemas moleculares"""
    
    def __init__(self, input_gro: str, topol_top: str, mdp_file: Optional[str] = None, gmx: str = "gmx_mpi"):
        """
        Inicializa el minimizador
        
        Args:
            input_gro: Ruta al archivo de coordenadas (.gro)
            topol_top: Ruta al archivo de topología (.top)
            mdp_file: Ruta al archivo de parámetros de minimización (.mdp)
            gmx: Comando de GROMACS a utilizar
        """
        self.input_gro = Path(input_gro)
        self.topol_top = Path(topol_top)
        self.mdp_file = Path(mdp_file) if mdp_file else None
        self.gmx = gmx
        
    def create_mdp_file(self, output_path: str) -> Path:
        """
        Crea un archivo de parámetros para minimización
        
        Args:
            output_path: Ruta para guardar el archivo .mdp
            
        Returns:
            Ruta al archivo .mdp creado
        """
        output_path = Path(output_path)
        
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
ns_type             = grid      ; Method to determine neighbor list (simple, grid)
coulombtype         = PME       ; Treatment of long range electrostatic interactions
rcoulomb            = 1.0       ; Short-range electrostatic cut-off
rvdw                = 1.0       ; Short-range Van der Waals cut-off
pbc                 = xyz       ; Periodic Boundary Conditions in all 3 dimensions
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
        
        # 1. Crear archivo .mdp si no existe
        console.log("Creando archivo de parámetros de minimización...")
        if not self.mdp_file:
            self.mdp_file = self.create_mdp_file(output_dir / "min.mdp")
                
        # 2. Generar archivo .tpr
        console.log("Generando archivo .tpr...")
        tpr_path = output_dir / "topol.tpr"
        cmd = [
            self.gmx, "grompp",
            "-f", str(self.mdp_file),
            "-c", str(self.input_gro),
            "-p", str(self.topol_top),
            "-o", str(tpr_path)
        ]
        subprocess.run(cmd, check=True)
            
        # 3. Ejecutar minimización
        console.log("Ejecutando minimización...")
        cmd = [
            "gmx_mpi", "mdrun", "-v",
            "-s", str(tpr_path),
            "-gpu_id", gpu_ids,
            "-tunepme"
        ]
        subprocess.run(cmd, check=True)
            
        # 4. Procesar la estructura final
        console.log("Procesando estructura final...")
        
        # Mover archivos generados al directorio de salida
        for file in ["confout.gro", "ener.edr", "md.log"]:
            if Path(file).exists():
                shutil.move(file, str(output_dir / file))
                
        confout_gro = output_dir / "confout.gro"
        
        # Centrar molécula
        cmd = [self.gmx, "trjconv", "-f", str(confout_gro), "-s", str(tpr_path),
              "-pbc", "mol", "-center", "-o", str(output_dir / "tmp.gro")]
        subprocess.run(cmd, input=b"1 0\n", check=True)
        
        # Compactar
        cmd = [self.gmx, "trjconv", "-f", str(output_dir / "tmp.gro"), "-ur", "compact",
              "-pbc", "mol", "-o", str(output_dir / "minimized.gro"), "-s", str(tpr_path)]
        subprocess.run(cmd, input=b"0\n", check=True)
        
        # Limpiar archivos temporales
        for tmp_file in output_dir.glob("tmp*.gro"):
            tmp_file.unlink()
        for backup_file in output_dir.glob("#*"):
            backup_file.unlink()
            
        return {
            "tpr": tpr_path,
            "confout": confout_gro,
            "final": output_dir / "minimized.gro",
            "edr": output_dir / "ener.edr",
            "log": output_dir / "md.log"
        }