#!/bin/bash
#
#Este script genera una caja dodecaedrica con una distancia entre la proteína y la caja de 1.5 nm con el nombre del argumento de la carpeta a trabajar con el argumento $1

cd $1

cp ../em.mdp .

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

