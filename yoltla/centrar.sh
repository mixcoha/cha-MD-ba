#!/bin/bash

cd $1/dynamics

mkdir cell

gmx trjcat -f */traj_comp.xtc -o trajout.xtc 

mv trajout.xtc cell/

cp 1/tprout.tpr cell/tprout.tpr 

cd cell



echo 1 0 | gmx trjconv -f trajout.xtc -s tprout.tpr -pbc mol -center -o tmp.xtc
 
echo 0 | gmx trjconv -f tmp.xtc -ur compact -pbc mol -o tmp2.xtc -s tprout.tpr

echo 1 0 | gmx trjconv -f tmp2.xtc -pbc cluster -o tmp.xtc -s tprout.tpr

echo 1 0 | gmx trjconv -f tmp.xtc -ur compact -pbc mol -center -o tmp2.xtc -s tprout.tpr

echo 1 0 | gmx trjconv -f tmp2.xtc -fit rot+trans -o tmp3.xtc -s tprout.tpr -tu ns

echo 0 | gmx trjconv -f trajout.xtc -o cell.pdb -dump 0 -s tprout.tpr

echo 1 | gmx trjconv -f trajout.xtc -o protein.pdb -dump 0 -s tprout.tpr -n index.ndx

echo 1 | gmx trjconv -f tmp3.xtc -o protein.xtc -s tprout.tpr -n index.ndx

rm tmp2.xtc tmp.xtc

mv tmp3.xtc cell.xtc 

rm \#*

cd ..

