#!/bin/sh

#  kk.sh
#  
#
#  Created by Edgar Mixcoha on 05/01/18.
#  
# Crear los directorios de trabajo para la primera dinámica de los desplegamietos térmicos a 298 310 350 400 450 500 550 600 K

RUNDIR=$PWD

cd $RUNDIR/$1

mkdir $RUNDIR/$1/npt
cd npt
mkdir R1
cd R1

# Crear las carpetas de las diferentes temperaturas

for temp in 298 310 350 400 450 500 550 600
do

mkdir $temp
cd $temp

# Crear los archivos para la ejecución de la dinámica
# Creación del .mdp

cat > nvt_$1_$temp_$nd.mdp << mdp
title        = NPT thermal denaturation simulation of $1 at $temp K
;define        = -DPOSRES    ; No position restrain the protein
; Run parameters
integrator     = md                    ; leap-frog integrator
nsteps         = 10000000              ; 2 * 500000 = 20 ns
dt             = 0.002                 ; 2 fs
; Output control
nstxout        = 5000                  ; save coordinates every 1.0 ps
nstvout        = 5000                  ; save velocities every 1.0 ps
nstenergy      = 5000                  ; save energies every 1.0 ps
nstlog         = 5000                  ; update log file every 1.0 ps
; Bond parameters
continuation              = no         ; first dynamics run
constraint_algorithm      = lincs      ; holonomic constraints
constraints               = all-bonds  ; all bonds (even heavy atom-H bonds) constrained
lincs_iter                = 1          ; accuracy of LINCS
lincs_order               = 4          ; also related to accuracy
; Neighborsearching
cutoff-scheme  = Verlet
ns_type        = grid                  ; search neighboring grid cells
nstlist        = 10                    ; 20 fs, largely irrelevant with Verlet
rcoulomb       = 1.2                   ; short-range electrostatic cutoff (in nm)
rvdw           = 1.2                   ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype    = PME                   ; Particle Mesh Ewald for long-range electrostatics
pme_order      = 4                     ; cubic interpolation
fourierspacing = 0.16                  ; grid spacing for FFT
; Temperature coupling is on
tcoupl         = V-rescale             ; modified Berendsen thermostat
tc-grps        = Protein Non-Protein   ; two coupling groups - more accurate
tau_t          = 0.1      0.1          ; time constant, in ps
ref_t          = $temp       $temp     ; reference temperature, one for each group, in K
; Pressure coupling is on
Pcoupl        = berendsen              ;Pressure coupling method
Pcoupltype    = isotropic              ;Pressure coupling type
tau_p         = 1.0                    ;Pressure relaxation time constant (ps)
compressibility  = 4.5e-5              ;Compressibility (1/bar)
ref_p         = 1.0                    ;Reference Pressure (bar)
; Periodic boundary conditions
pbc           = xyz                    ; 3-D PBC
; Dispersion correction
DispCorr      = EnerPres               ; account for cut-off vdW scheme
; Velocity generation
gen_vel       = yes                    ; assign velocities from Maxwell distribution
gen_temp      = $temp                  ; temperature for Maxwell distribution
gen_seed      = -1                     ; generate a random seed

mdp

# Creación de la carpeta de trabajo y el .tpr para el inicio de la simulación

module purge
module load gromacs51/5.1.4

ls -lhtr > ls
grep drwxr ls > tmp
cat tmp | sort -k 9 -g > tmp1
mv tmp1 tmp
nd=`awk 'END {print $NF+1}' tmp`
rm ls tmp
mkdir $nd
mv *.* $nd
cd $nd

gmx_mpi grompp -f ../nvt_$1_$temp_$nd.mdp -c $RUNDIR/$1/nvt/200/confout.gro -p $RUNDIR/$1/setup/topol.top

# Creación del archivo .lsr para el envío de las dinámicas a la cola de Miztli

cat > $1_$temp_$nd.lsf << estodo1
### Configuración para el envío del cálculo para las colas MPI, Utilizar como segundo argumento del script el nombre de la cola. Como tercer argumento se introduce el número de procesadores a utilizar
#BSUB -q $2  #### q_htc
#BSUB -n $3  #### 32
#BSUB -J $1_$temp_$nd
#BSUB -R "span[hosts=1]"
#BSUB -oo salida
#BSUB -eo error
#BSUB -u edgarmixcoha@gmail.com

cd $RUNDIR/$1/npt/R1/$temp

module purge
module load gromacs51/5.1.4

export OMP_NUM_THREADS=\$LSB_DJOB_NUMPROC

echo "Cantidad de threads \$OMP_NUM_THREADS"
echo "Ejecucion de mdrun de GROMACS"

mpirun -hostfile /hostlist.dat gmx_mpi mdrun -v -s  >& mdrun.log

#exit 0

estodo1

chmod +x $1_$temp_$nd.lsf

bsub < $1_$temp_$nd.lsf

cd $RUNDIR/$1/npt/R1

done







