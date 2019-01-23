#!/bin/bash

cd $1/npt

mkdir cell
cd cell

#gmx trjcat -f ../*/traj_comp.xtc -o trajout.xtc
gmx trjcat -f ../*/traj.trr -o trajout.xtc

cp ../1/topol.tpr topol.tpr

echo 1 0 | gmx trjconv -f trajout.xtc -s topol.tpr -pbc mol -center -o tmp.xtc
 
echo 0 | gmx trjconv -f tmp.xtc -ur compact -pbc mol -o tmp2.xtc -s topol.tpr

echo 1 0 | gmx trjconv -f tmp2.xtc -pbc cluster -o tmp.xtc -s topol.tpr

echo 1 0 | gmx trjconv -f tmp.xtc -ur compact -pbc mol -center -o tmp2.xtc -s topol.tpr

echo 4 0 | gmx trjconv -f tmp2.xtc -fit rot+trans -o tmp3.xtc -s topol.tpr -tu ns

echo 0 | gmx trjconv -f trajout.xtc -o cell.pdb -dump 0 -s topol.tpr

echo 1 | gmx trjconv -f trajout.xtc -o protein.pdb -dump 0 -s topol.tpr

echo 1 | gmx trjconv -f tmp3.xtc -o protein.xtc -s topol.tpr

rm tmp2.xtc tmp.xtc

mv tmp3.xtc cell.xtc 

rm \#*

cd ..

