#!/bin/sh

#  nvt1.sh
#  
#
#  Created by Edgar Mixcoha on 15/12/17.
#  
# Este script leerá la geometría de un sistema minimizado y realizará 1 ns de una simulación NVT y generará los archivos necesarios para disminuir gradualmente la constante de restricción de los átomos pesados en 200 u. Generará las carpetas con el nombre posre_constante-de-fuerza y los mdps para realizar las dinámicas
# el script se ejecutará: ./nvt.sh nombre-de-la-carpeta

cd $1
mkdir nvt
cp setup/topol.top nvt
cd nvt
i=0

###############################
# Realizar las dinámicas comenzando en posre = 1000
###############################

for const in  1000 800 600 400 200
do
mkdir posre_$const
cd posre_$const
cp ../topol.top .

sed "s/1000/${const}/g" ../../setup/posre.itp > ./posre.itp

j=$(( ${const} - ( 200 * $i ) ))
echo "${const}"

i=$(( $i+1 ))

 if [[ "$j" = 1000 ]]
 then
        echo "Inicia la equilibración NVT con posre=$const"

cat > nvt_${const}.mdp << nvt
title        = NVT simulation de $1 con posre de ${const}
define        = -DPOSRES    ; position restrain the protein
; Run parameters
integrator    = md        ; leap-frog integrator
nsteps        = 250000        ; 2 * 250000 = 0.5 ns
dt            = 0.002        ; 2 fs
; Output control
nstxout        = 50        ; save coordinates every 1.0 ps
nstvout        = 50        ; save velocities every 1.0 ps
nstenergy        = 50        ; save energies every 1.0 ps
nstlog        = 50        ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs        ; holonomic constraints
constraints                = all-bonds    ; all bonds (even heavy atom-H bonds) constrained
lincs_iter                = 1            ; accuracy of LINCS
lincs_order                = 4            ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type            = grid        ; search neighboring grid cells
nstlist            = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb        = 1.0        ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype        = PME    ; Particle Mesh Ewald for long-range electrostatics
pme_order        = 4        ; cubic interpolation
fourierspacing    = 0.16    ; grid spacing for FFT
; Temperature coupling is on
tcoupl        = V-rescale                ; modified Berendsen thermostat
tc-grps        = Protein Non-Protein    ; two coupling groups - more accurate
tau_t        = 0.1      0.1           ; time constant, in ps
ref_t        = 298       298           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl        = no         ; no pressure coupling in NVT
; Periodic boundary conditions
pbc        = xyz            ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres    ; account for cut-off vdW scheme
; Velocity generation
gen_vel        = yes        ; assign velocities from Maxwell distribution
gen_temp    = 298        ; temperature for Maxwell distribution
gen_seed    = -1        ; generate a random seed

nvt
		gmx grompp -f nvt_${const}.mdp -c ../../em/confout.gro -p topol.top

    else [[ "$j" < 1000 ]]

        echo "Se continúa la dinámica con una constante de restricción de $j"

k=$(( ${const} + 200 ))

cat > nvt_$const.mdp << nvt1
title        = NVT simulation de $1 con posre de $j
define        = -DPOSRES    ; position restrain the protein
; Run parameters
integrator    = md        ; leap-frog integrator
nsteps        = 25000        ; 2 * 250000 = 0.5 ns
dt            = 0.002        ; 2 fs
; Output control
nstxout        = 50        ; save coordinates every 1.0 ps
nstvout        = 50        ; save velocities every 1.0 ps
nstenergy        = 50        ; save energies every 1.0 ps
nstlog        = 50        ; update log file every 1.0 ps
; Bond parameters
continuation            = yes        ; first dynamics run
constraint_algorithm    = lincs        ; holonomic constraints
constraints                = all-bonds    ; all bonds (even heavy atom-H bonds) constrained
lincs_iter                = 1            ; accuracy of LINCS
lincs_order                = 4            ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type            = grid        ; search neighboring grid cells
nstlist            = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb        = 1.0        ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype        = PME    ; Particle Mesh Ewald for long-range electrostatics
pme_order        = 4        ; cubic interpolation
fourierspacing    = 0.16    ; grid spacing for FFT
; Temperature coupling is on
tcoupl        = V-rescale                ; modified Berendsen thermostat
tc-grps        = Protein Non-Protein    ; two coupling groups - more accurate
tau_t        = 0.1      0.1           ; time constant, in ps
ref_t        = 298       298           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl        = no         ; no pressure coupling in NVT
; Periodic boundary conditions
pbc        = xyz            ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres    ; account for cut-off vdW scheme

nvt1

gmx grompp -f nvt_${const}.mdp -c ../posre_${k}/confout.gro -p topol.top -t ../posre_${k}/state.cpt

    fi

gmx mdrun -gpu_id 01 -nb gpu_cpu -tunepme -v -s >& mdrun.log


cd ..

done

