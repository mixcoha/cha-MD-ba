# Artefactos de la corrida 6M03 (NaCl 0.5 M, 310 K)

Copia portable de la simulación ejecutada en el Cloud Agent. Tras
`git pull` de la rama `cursor/6m03-benchmark-nacl-310k-c08b` estos
archivos aparecen en tu clon local **aquí**, no en `simulations/` ni
`data/` (ambos siguen en `.gitignore`).

Estado de la corrida: preparación + minimización (convergida) +
**NVT 100 ps** con POSRES 1000 kJ mol⁻¹ nm⁻² (terminado 2026-08-31).
NPT y producción no se ejecutaron.

## Inventario

| Ruta local | Contenido |
| --- | --- |
| `pdb/6M03_rcsb.pdb` | PDB crudo descargado de RCSB |
| `pdb/6M03.pdb` | PDB limpio (solo ATOM de proteína, cadena A) |
| `1_preparation/` | `6M03.gro`, `6M03_box.gro`, `6M03_ions.gro`, `topol.top`, `posre.itp`, `ions.mdp`, `ions.tpr`, `mdout.mdp` |
| `2_minimization/` | `min.mdp`, `minimized.gro`, `em.gro`, `em.log`, `em.edr`, `topol.tpr` |
| `3_nvt/posre_constante/1000/` | `nvt.mdp`, `nvt.gro`, `nvt.log`, `nvt.edr`, `topol.tpr`, `topol.top`, `posre.itp` |
| `protocol/` | `nvt.mdp`, `npt.mdp`, `md.mdp` de referencia a 310 K |
| `benchmark_report.json` | Composición del sistema y rutas originales |
| `timing_local.json` | Tiempos de minimización y NVT en la VM |

Para continuar NPT/producción localmente usa
`3_nvt/posre_constante/1000/nvt.gro` + `topol.top` + `posre.itp`.

## Omitido (no cabe en git / no es portable)

GitHub advierte a 50 MB y rechaza archivos ≥ 100 MB. Por eso **no** se
incluyó la trayectoria NVT:

| Archivo original en la VM | Tamaño | Motivo |
| --- | --- | --- |
| `simulations/6M03/3_nvt/posre_constante/1000/nvt.trr` | **158 MB** | Trayectoria NVT; supera el límite práctico de git |
| `simulations/6M03/2_minimization/em.trr` | 801 KB | Trayectoria de minimización (omitida junto con `.trr`) |
| `simulations/6M03/3_nvt/posre_constante/1000/nvt.cpt` | 1.6 MB | Checkpoint de reinicio (no necesario para inspeccionar el estado final) |
| `simulations/6M03/3_nvt/posre_constante/1000/tmp.gro` | 4.5 MB | Temporal de GROMACS, duplicado de `nvt.gro` |
| `simulations/6M03/1_preparation/6M03_solv.gro` | 3.0 MB | Intermedio pre-iones; el sistema con NaCl está en `6M03_ions.gro` |

La trayectoria `nvt.trr` quedó solo en la VM del Cloud Agent
(`/workspace/simulations/6M03/.../nvt.trr`). El estado final de NVT
sí está versionado: `nvt.gro` (~4.5 MB) y `nvt.edr` (energías).
