#!/bin/bash
#
#Este script genera una caja dodecaedrica con una distancia entre la proteína y la caja de 1.5 nm con el nombre del argumento de la carpeta a trabajar con el argumento $1

cd $1

cat > em.mdp < estodo

; Líneas que comienzan con ';' son considerados comentarios
title		= Minimization de la proteina $1 e

; Parameters describing what to do, when to stop and what to save
integrator	= steep		; Algorithm options i.e. steep = steepest descent minimization 
emtol		    = 10.0		; Stop minimization when the energy changes by less than emtol kJ/mol.
nsteps	   	= 1000		; Maximum number of (minimization) steps to perform
nstenergy 	= 10	  	; Write energies to disk every nstenergy steps
nstxtcout 	= 10		  ; Write coordinates to disk every nstxtcout steps
xtc_grps	  = Protein	; Which coordinate group(s) to write to disk
energygrps	= Protein	; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist	  	= 1		    ; Frequency to update the neighbor list and long range forces
ns_type	 	  = grid    ; Method to determine neighbor list (simple, grid)
rlist	    	= 1.4	  	; Cut-off for making neighbor list (short range forces)
coulombtype	= cut-off	; Treatment of long range electrostatic interactions
rcoulomb	  = 1.4	   	; long range electrostatic cut-off
rvdw		    = 1.4		  ; long range Van der Waals cut-off
constraints	= none		; Bond types to replace by constraints
pbc		      = yes		  ; Periodic Boundary Conditions (yes/no)

estodo

cp $1.pqr $1_1.pdb

mkdir setup
mv em.mdp setup
mv $1_1.pdb setup
cd setup

echo 5 | gmx pdb2gmx -f $1_1.pdb -ignh -water spc -o $1_1.gro

gmx editconf -bt dodecahedron -d 2 -c -f $1_1.gro -o $1_1_box.gro

gmx solvate -cp $1_1_box.gro -cs spc216.gro -p topol.top -o prot_wat.gro

gmx grompp -f em.mdp -c prot_wat.gro -p topol.top -o -v

echo 13 | gmx genion -s topol.tpr -p topol.top -o prot_wat_ion.gro -neutral -nname CL -pname NA

gmx grompp -f em.mdp -c prot_wat_ion.gro -p topol.top  -o -v

rm -rf \#*

cd ..
mkdir em
mv setup/topol.tpr em
cd em

gmx mdrun -gpu_id 01 -nb gpu_cpu -tunepme -v -s

rm -rf \#*

echo 1 0 | gmx trjconv -f confout.gro -s topol.tpr -pbc mol -center -o tmp.gro
 
echo 0 | gmx trjconv -f tmp.gro -ur compact -pbc mol -o tmp2.gro -s topol.tpr

echo 1 0 | gmx trjconv -f tmp2.gro -pbc cluster -o tmp.gro -s topol.tpr

echo 1 0 | gmx trjconv -f tmp.gro -ur compact -pbc mol -center -o tmp2.gro -s topol.tpr

echo 4 0 | gmx trjconv -f tmp2.gro -fit rot+trans -o cell.pdb -s topol.tpr

rm tmp*.gro \#*

cd ../..

