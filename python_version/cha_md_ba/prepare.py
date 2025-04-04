"""
Módulo para preparar simulaciones de dinámica molecular
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any
from rich.console import Console
from rich.progress import Progress
import requests
from .minimize import EnergyMinimizer

def download_pdb(pdb_code: str, output_path: str) -> bool:
    """
    Descarga un archivo PDB desde el RCSB PDB.

    Args:
        pdb_code: Código de la estructura en el RCSB (ej. '1tim')
        output_path: Ruta donde guardar el archivo descargado

    Returns:
        True si la descarga fue exitosa, False en caso contrario
    """
    url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            with open(output_path, "w") as f:
                f.write(r.text)
            return True
        else:
            console.print(f"[red]❌ No se pudo descargar {pdb_code} desde el PDB. Código HTTP {r.status_code}")
            return False
    except Exception as e:
        console.print(f"[red]❌ Error al intentar descargar: {e}")
        return False

console = Console()

class MDSystemPreparator:
    """Clase para preparar sistemas de dinámica molecular"""
    
    def __init__(self, pdb_path: str, forcefield: str = "amber99sb-ildn", water_model: str = "tip3p"):
        """
        Inicializa el preparador de sistemas
        
        Args:
            pdb_path: Ruta al archivo PDB
            forcefield: Campo de fuerzas a utilizar
            water_model: Modelo de agua a utilizar
        """
        self.pdb_path = Path(pdb_path)
        self.forcefield = forcefield
        self.water_model = water_model
        self.system_name = self.pdb_path.stem  # Nombre del sistema basado en el PDB
        self.gmx = "gmx_mpi"  # Comando de GROMACS
        
    def _create_ions_mdp(self, output_path: Path) -> Path:
        """Crea el archivo de parámetros para genion"""
        mdp_content = """; Líneas que comienzan con ';' son considerados comentarios
title               = Ions

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
        
    def prepare_system(self, output_dir: str, box_type: str = "dodecahedron", 
                      box_size: Optional[float] = None, ions: bool = True,
                      minimize: bool = True, gpu_ids: Optional[str] = None,
                      replace_group: str = "SOL") -> Dict[str, Path]:
        """
        Prepara el sistema para simulación
        
        Args:
            output_dir: Directorio de salida
            box_type: Tipo de caja de simulación
            box_size: Tamaño de la caja (opcional)
            ions: Si se deben agregar iones
            minimize: Si se debe realizar minimización de energía
            gpu_ids: IDs de GPUs a utilizar para minimización
            replace_group: Grupo de átomos a reemplazar con iones (por defecto "SOL")
            
        Returns:
            Diccionario con rutas a los archivos generados
        """
        # Crear estructura de directorios
        base_dir = Path(output_dir) / self.system_name
        prep_dir = base_dir / "1_preparation"
        min_dir = base_dir / "2_minimization"
        nvt_dir = base_dir / "3_nvt"
        npt_dir = base_dir / "4_npt"
        
        # Crear directorios
        for dir_path in [prep_dir, min_dir, nvt_dir, npt_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
            
        # Copiar el archivo PDB al directorio de preparación
        prep_pdb = prep_dir / f"{self.system_name}.pdb"
        if not prep_pdb.exists():
            import shutil
            shutil.copy2(self.pdb_path, prep_pdb)
        
        with Progress() as progress:
            # 1. Generar topología
            task = progress.add_task("[cyan]Generando topología...", total=1)
            topol_path = self._generate_topology(prep_dir)
            progress.advance(task)
            
            # 2. Definir caja
            task = progress.add_task("[cyan]Definiendo caja de simulación...", total=1)
            box_path = self._define_box(prep_dir, box_type, box_size)
            progress.advance(task)
            
            # 3. Solvatación
            task = progress.add_task("[cyan]Solvatando sistema...", total=1)
            solv_path = self._solvate_system(prep_dir, box_path)
            progress.advance(task)
            
            # 4. Neutralización
            if ions:
                task = progress.add_task("[cyan]Neutralizando sistema...", total=1)
                ions_path = self._neutralize_system(prep_dir, solv_path, replace_group)
                progress.advance(task)
                
            # 5. Minimización de energía
            if minimize:
                task = progress.add_task("[cyan]Minimizando energía...", total=1)
                minimizer = EnergyMinimizer(
                    input_gro=str(ions_path if ions else solv_path),
                    topol_top=str(topol_path),
                    gmx=self.gmx
                )
                min_files = minimizer.minimize(min_dir, gpu_ids)
                progress.advance(task)
                
        result = {
            "base_dir": base_dir,
            "preparation_dir": prep_dir,
            "minimization_dir": min_dir,
            "nvt_dir": nvt_dir,
            "npt_dir": npt_dir,
            "topology": topol_path,
            "box": box_path,
            "solvated": solv_path,
            "ions": ions_path if ions else None
        }
        
        if minimize:
            result.update({
                "minimized": min_files["final"],
                "min_tpr": min_files["tpr"],
                "min_confout": min_files["confout"]
            })
            
        return result
        
    def _generate_topology(self, output_dir: Path) -> Path:
        """Genera la topología del sistema"""
        output_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = output_dir / f"{self.system_name}.pdb"
        
        # Generar topología usando pdb2gmx
        topol_path = output_dir / "topol.top"
        cmd = [
            self.gmx, "pdb2gmx",
            "-f", str(pdb_path),
            "-o", str(output_dir / f"{self.system_name}.gro"),
            "-p", str(topol_path),
            "-ff", self.forcefield,
            "-water", self.water_model
        ]
        subprocess.run(cmd, check=True)
        
        return topol_path
        
    def _define_box(self, output_dir: Path, box_type: str, box_size: Optional[float]) -> Path:
        """Define la caja de simulación"""
        input_gro = output_dir / f"{self.system_name}.gro"
        output_gro = output_dir / f"{self.system_name}_box.gro"
        
        cmd = [self.gmx, "editconf", "-f", str(input_gro), "-o", str(output_gro)]
        
        if box_type:
            cmd.extend(["-bt", box_type])
        if box_size:
            cmd.extend(["-d", str(box_size)])
            
        subprocess.run(cmd, check=True)
        
        return output_gro
        
    def _solvate_system(self, output_dir: Path, box_path: Path) -> Path:
        """Solvata el sistema"""
        input_gro = box_path
        output_gro = output_dir / f"{self.system_name}_solv.gro"
        topol_path = output_dir / "topol.top"
        
        # Solvatación
        cmd = [
            self.gmx, "solvate",
            "-cp", str(input_gro),
            "-cs", "spc216.gro",
            "-o", str(output_gro),
            "-p", str(topol_path)
        ]
        subprocess.run(cmd, check=True)
        
        return output_gro
        
    def _neutralize_system(self, output_dir: Path, solv_path: Path, replace_group: str = "SOL") -> Path:
        """Neutraliza el sistema
        
        Args:
            output_dir: Directorio de salida
            solv_path: Ruta al archivo GRO solvatado
            replace_group: Grupo de átomos a reemplazar con iones (por defecto "SOL")
        """
        input_gro = solv_path
        output_gro = output_dir / f"{self.system_name}_ions.gro"
        topol_path = output_dir / "topol.top"
        
        # Crear archivo .mdp para genion
        mdp_path = self._create_ions_mdp(output_dir / "ions.mdp")
        
        # Generar archivo .tpr para genion
        tpr_path = output_dir / "ions.tpr"
        cmd = [
            self.gmx, "grompp",
            "-f", str(mdp_path),
            "-c", str(input_gro),
            "-p", str(topol_path),
            "-o", str(tpr_path)
        ]
        subprocess.run(cmd, check=True)
        
        # Agregar iones
        cmd = [
            self.gmx, "genion",
            "-s", str(tpr_path),
            "-o", str(output_gro),
            "-p", str(topol_path),
            "-pname", "NA",
            "-nname", "CL",
            "-neutral"
        ]
        subprocess.run(cmd, input=f"{replace_group}\n".encode(), check=True)
        
        return output_gro