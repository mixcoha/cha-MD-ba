"""
Módulo para preparar simulaciones de dinámica molecular
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any
from rich.console import Console
from rich.progress import Progress
from .minimize import EnergyMinimizer

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
        
    def prepare_system(self, output_dir: str, box_type: str = "dodecahedron", 
                      box_size: Optional[float] = None, ions: bool = True,
                      minimize: bool = True, gpu_ids: Optional[str] = None) -> Dict[str, Path]:
        """
        Prepara el sistema para simulación
        
        Args:
            output_dir: Directorio de salida
            box_type: Tipo de caja de simulación
            box_size: Tamaño de la caja (opcional)
            ions: Si se deben agregar iones
            minimize: Si se debe realizar minimización de energía
            gpu_ids: IDs de GPUs a utilizar para minimización
            
        Returns:
            Diccionario con rutas a los archivos generados
        """
        # Crear estructura de directorios
        system_name = self.pdb_path.stem
        base_dir = Path(output_dir) / system_name
        prep_dir = base_dir / "1_preparation"
        min_dir = base_dir / "2_minimization"
        nvt_dir = base_dir / "3_nvt"
        npt_dir = base_dir / "4_npt"
        
        # Crear directorios
        for dir_path in [prep_dir, min_dir, nvt_dir, npt_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
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
                ions_path = self._neutralize_system(prep_dir, solv_path)
                progress.advance(task)
                
            # 5. Minimización de energía
            if minimize:
                task = progress.add_task("[cyan]Minimizando energía...", total=1)
                minimizer = EnergyMinimizer(
                    input_gro=str(ions_path if ions else solv_path),
                    topol_top=str(topol_path)
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
        output_path = output_dir / "topol.top"
        
        cmd = [
            "gmx", "pdb2gmx",
            "-f", str(self.pdb_path),
            "-o", str(output_dir / "conf.gro"),
            "-p", str(output_path),
            "-ff", self.forcefield,
            "-water", self.water_model
        ]
        
        subprocess.run(cmd, check=True)
        return output_path
        
    def _define_box(self, output_dir: Path, box_type: str, box_size: Optional[float]) -> Path:
        """Define la caja de simulación"""
        output_path = output_dir / "box.gro"
        
        cmd = ["gmx", "editconf", "-f", str(output_dir / "conf.gro"),
               "-o", str(output_path), "-c", "-d", "1.0", "-bt", box_type]
        
        if box_size:
            cmd.extend(["-box", str(box_size), str(box_size), str(box_size)])
            
        subprocess.run(cmd, check=True)
        return output_path
        
    def _solvate_system(self, output_dir: Path, box_path: Path) -> Path:
        """Solvata el sistema"""
        output_path = output_dir / "solv.gro"
        
        cmd = [
            "gmx", "solvate",
            "-cp", str(box_path),
            "-cs", "spc216.gro",
            "-o", str(output_path),
            "-p", str(output_dir / "topol.top")
        ]
        
        subprocess.run(cmd, check=True)
        return output_path
        
    def _neutralize_system(self, output_dir: Path, solv_path: Path) -> Path:
        """Neutraliza el sistema con iones"""
        output_path = output_dir / "ions.gro"
        
        # Primero generar el archivo .tpr
        tpr_path = output_dir / "ions.tpr"
        cmd = [
            "gmx", "grompp",
            "-f", "ions.mdp",
            "-c", str(solv_path),
            "-p", str(output_dir / "topol.top"),
            "-o", str(tpr_path)
        ]
        
        subprocess.run(cmd, check=True)
        
        # Luego agregar iones
        cmd = [
            "gmx", "genion",
            "-s", str(tpr_path),
            "-o", str(output_path),
            "-p", str(output_dir / "topol.top"),
            "-pname", "NA",
            "-nname", "CL",
            "-neutral"
        ]
        
        subprocess.run(cmd, check=True)
        return output_path 