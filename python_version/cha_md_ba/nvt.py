"""
Módulo para equilibración NVT de sistemas de dinámica molecular
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Dict, Tuple
from rich.console import Console
from rich.progress import Progress

console = Console()

class NVTEquilibrator:
    """Clase para realizar equilibración NVT de sistemas moleculares"""
    
    def __init__(
        self,
        input_gro: str,
        topol_top: str,
        mdp_file: Optional[str] = None,
        temperature: float = 300.0,
        nsteps: int = 50000,
        gmx: Optional[str] = None,
    ):
        """
        Inicializa el equilibrador NVT
        
        Args:
            input_gro: Ruta al archivo de coordenadas (.gro)
            topol_top: Ruta al archivo de topología (.top)
            mdp_file: Ruta al archivo de parámetros NVT (.mdp)
            temperature: Temperatura de referencia en Kelvin
            nsteps: Número de pasos de integración
            gmx: Ejecutable de GROMACS
        """
        from .gmx_utils import find_gmx

        self.input_gro = Path(input_gro)
        self.topol_top = Path(topol_top)
        self.mdp_file = Path(mdp_file) if mdp_file else None
        self.temperature = temperature
        self.nsteps = nsteps
        self.gmx = gmx or find_gmx()
        
    def check_simulation_state(self, base_dir: str) -> Tuple[str, Optional[int], Dict[str, bool]]:
        """
        Verifica el estado actual de la simulación
        
        Args:
            base_dir: Directorio base de la simulación
            
        Returns:
            Tupla con (fase_actual, constante_actual, detalles) donde:
            - fase_actual: "preparation", "minimization", "nvt", "npt" o "production"
            - constante_actual: Constante de fuerza actual en NVT o None si no está en NVT
            - detalles: Diccionario con información detallada del estado
        """
        base_path = Path(base_dir)
        detalles = {
            "preparation_completa": False,
            "minimization_completa": False,
            "nvt_constantes": [],
            "npt_completa": False,
            "production_completa": False
        }
        
        # Verificar preparación
        prep_dir = base_path / "1_preparation"
        if prep_dir.exists():
            detalles["preparation_completa"] = (prep_dir / "topol.top").exists() and any(
                prep_dir.glob("*.itp")
            )
            if not detalles["preparation_completa"]:
                return "preparation", None, detalles
                
        # Verificar minimización
        min_dir = base_path / "2_minimization"
        if min_dir.exists():
            detalles["minimization_completa"] = all([
                (min_dir / "minimized.gro").exists(),
                (min_dir / "ener.edr").exists()
            ])
            if not detalles["minimization_completa"]:
                return "minimization", None, detalles
                
        # Verificar NVT
        nvt_dir = base_path / "3_nvt"
        if nvt_dir.exists():
            posre_dir = nvt_dir / "posre_constante"
            if posre_dir.exists():
                force_constants = [1000, 800, 600, 400, 200]
                for fc in force_constants:
                    fc_dir = posre_dir / str(fc)
                    if fc_dir.exists():
                        # Verificar si la simulación está completa
                        if all([
                            (fc_dir / "nvt.gro").exists(),
                            (fc_dir / "ener.edr").exists() or (fc_dir / "nvt.edr").exists(),
                        ]):
                            detalles["nvt_constantes"].append(fc)
                            
                if not detalles["nvt_constantes"]:
                    return "nvt", 1000, detalles
                    
                # Encontrar la última constante completada
                last_completed = max(detalles["nvt_constantes"])
                
                # Verificar si es la última constante
                if last_completed == 200:
                    # Verificar NPT
                    npt_dir = base_path / "4_npt"
                    if npt_dir.exists():
                        detalles["npt_completa"] = all([
                            (npt_dir / "npt.gro").exists(),
                            (npt_dir / "ener.edr").exists()
                        ])
                        if not detalles["npt_completa"]:
                            return "npt", None, detalles
                            
                        # Verificar producción
                        prod_dir = base_path / "5_production"
                        if prod_dir.exists():
                            detalles["production_completa"] = all([
                                (prod_dir / "production.gro").exists(),
                                (prod_dir / "ener.edr").exists()
                            ])
                            if not detalles["production_completa"]:
                                return "production", None, detalles
                                
                        return "production", None, detalles
                        
                    return "npt", None, detalles
                    
                # Encontrar la siguiente constante a ejecutar
                next_fc = force_constants[force_constants.index(last_completed) + 1]
                return "nvt", next_fc, detalles
                
            return "nvt", 1000, detalles
            
        return "nvt", None, detalles
        
    def create_posre_files(self, output_dir: Path, force_constants: list[int]) -> None:
        """
        Crea archivos de restricción de posición a partir de los posre de pdb2gmx.
        
        Args:
            output_dir: Directorio base para los archivos de restricción
            force_constants: Lista de constantes de fuerza a usar
        """
        import shutil

        posre_dir = output_dir / "posre_constante"
        posre_dir.mkdir(parents=True, exist_ok=True)
        sources = list(self.topol_top.parent.glob("posre*.itp"))

        for fc in force_constants:
            fc_dir = posre_dir / str(fc)
            fc_dir.mkdir(exist_ok=True)

            if sources:
                for src in sources:
                    text = src.read_text().replace("1000", str(fc))
                    (fc_dir / src.name).write_text(text)
                for topo_file in self.topol_top.parent.iterdir():
                    if topo_file.suffix in {".top", ".itp"} and not topo_file.name.startswith("posre"):
                        shutil.copy2(topo_file, fc_dir / topo_file.name)
                continue

            (fc_dir / "posre.itp").write_text(
                f"""; Position restraint file
[ position_restraints ]
;  i funct       fcx        fcy        fcz
    1    1         {fc}         {fc}         {fc}
"""
            )
        
    def create_mdp_file(
        self,
        output_path: str,
        title: str = "NVT Equilibration",
        temperature: Optional[float] = None,
        nsteps: Optional[int] = None,
        use_posres: bool = True,
        gen_vel: bool = True,
        continuation: bool = False,
    ) -> Path:
        """
        Crea un archivo de parámetros para equilibración NVT
        
        Args:
            output_path: Ruta para guardar el archivo .mdp
            title: Título de la simulación
            temperature: Temperatura de referencia (K)
            nsteps: Número de pasos
            use_posres: Activar restricciones de posición (-DPOSRES)
            gen_vel: Generar velocidades de Maxwell-Boltzmann
            continuation: Continuar desde un estado previo
            
        Returns:
            Ruta al archivo .mdp creado
        """
        output_path = Path(output_path)
        temperature = self.temperature if temperature is None else temperature
        nsteps = self.nsteps if nsteps is None else nsteps
        define_line = "define              = -DPOSRES  ; restrain protein heavy atoms" if use_posres else "; define              = -DPOSRES"
        gen_vel_flag = "yes" if gen_vel else "no"
        continuation_flag = "yes" if continuation else "no"

        mdp_content = f"""; Líneas que comienzan con ';' son considerados comentarios
title               = {title}
{define_line}

; Parameters describing what to do, when to stop and what to save
integrator          = md        ; leap-frog integrator
dt                  = 0.002     ; 2 fs
nsteps              = {nsteps}     ; {nsteps * 0.002:.1f} ps
nstxout             = 500       ; save coordinates every 1 ps
nstvout             = 500       ; save velocities every 1 ps
nstenergy           = 500       ; save energies every 1 ps
nstlog              = 500       ; update log file every 1 ps
continuation        = {continuation_flag}
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
ref_t               = {temperature:g}       {temperature:g}    ; reference temperature, one for each group, in K

; Pressure coupling
pcoupl              = no        ; no pressure coupling in NVT

; Periodic boundary conditions
pbc                 = xyz       ; 3-D PBC

; Velocity generation
gen_vel             = {gen_vel_flag}       ; assign velocities from Maxwell distribution
gen_temp            = {temperature:g}       ; temperature for Maxwell distribution
gen_seed            = -1        ; generate a random seed
"""
        
        with open(output_path, 'w') as f:
            f.write(mdp_content)
            
        return output_path
        
    def equilibrate(self, output_dir: str, force_constants: list[int] = [1000, 800, 600, 400, 200], gpu_ids: Optional[str] = None) -> Dict[int, Dict[str, str]]:
        """
        Realiza la equilibración NVT del sistema
        
        Args:
            output_dir: Directorio de salida para los archivos
            force_constants: Lista de constantes de fuerza a usar
            gpu_ids: IDs de las GPUs a usar. None ejecuta en CPU.
            
        Returns:
            Diccionario con las rutas de los archivos generados para cada constante
        """
        output_path = Path(output_dir)
        results = {}
        current_gro = self.input_gro
        
        # Crear archivos de restricción
        self.create_posre_files(output_path, force_constants)
        
        for index, fc in enumerate(force_constants):
            console.print(f"\n[bold cyan]Ejecutando equilibración NVT con constante de fuerza {fc}[/bold cyan]")
            
            # Crear directorio para esta constante
            fc_dir = output_path / "posre_constante" / str(fc)
            fc_dir.mkdir(parents=True, exist_ok=True)
            
            # Crear archivo de parámetros NVT
            mdp_file = fc_dir / "nvt.mdp"
            self.create_mdp_file(
                mdp_file,
                title=f"NVT fc={fc} T={self.temperature:g} K",
                use_posres=True,
                gen_vel=(index == 0),
                continuation=(index > 0),
            )

            topology = fc_dir / "topol.top"
            if not topology.exists():
                topology = self.topol_top
            
            # Generar archivo .tpr
            tpr_file = fc_dir / "topol.tpr"
            cmd = [
                self.gmx, "grompp",
                "-f", str(mdp_file),
                "-c", str(current_gro),
                "-r", str(self.input_gro),
                "-p", str(topology),
                "-o", str(tpr_file),
                "-maxwarn", "1",
            ]
            
            with Progress() as progress:
                task = progress.add_task("[cyan]Generando archivo .tpr...", total=100)
                subprocess.run(cmd, check=True)
                progress.update(task, completed=100)
            
            # Ejecutar equilibración
            cmd = [
                self.gmx, "mdrun",
                "-v",
                "-s", str(tpr_file),
                "-deffnm", str(fc_dir / "nvt"),
            ]
            if gpu_ids:
                cmd.extend(["-gpu_id", gpu_ids, "-pme", "gpu"])
            
            with Progress() as progress:
                task = progress.add_task(f"[cyan]Ejecutando equilibración NVT para fc={fc}...", total=100)
                subprocess.run(cmd, check=True)
                progress.update(task, completed=100)
            
            # Procesar estructura final
            cmd = [
                self.gmx, "trjconv",
                "-f", str(fc_dir / "nvt.gro"),
                "-s", str(tpr_file),
                "-pbc", "mol",
                "-center",
                "-o", str(fc_dir / "tmp.gro")
            ]
            
            with Progress() as progress:
                task = progress.add_task(f"[cyan]Procesando estructura final para fc={fc}...", total=100)
                subprocess.run(cmd, input=b"1 0\n", check=True)
                progress.update(task, completed=100)

            current_gro = fc_dir / "tmp.gro"
            
            # Guardar resultados
            results[fc] = {
                "gro": str(current_gro),
                "edr": str(fc_dir / "nvt.edr"),
                "log": str(fc_dir / "nvt.log")
            }
            
        return results 