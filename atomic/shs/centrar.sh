#!/bin/bash

cd $1/em


echo 1 0 | gmx trjconv -f confort.gro -s topol.tpr -pbc mol -center -o tmp.gro
 
echo 0 | gmx trjconv -f tmp.gro -ur compact -pbc mol -o tmp2.gro -s topol.tpr

echo 1 0 | gmx trjconv -f tmp2.gro -pbc cluster -o tmp.gro -s topol.tpr

echo 1 0 | gmx trjconv -f tmp.gro -ur compact -pbc mol -center -o tmp2.gro -s topol.tpr

echo 4 0 | gmx trjconv -f tmp2.gro -fit rot+trans -o tmp3.gro -s topol.tpr -tu ns

echo 0 | gmx trjconv -f trajout.xtc -o cell.pdb -dump 0 -s topol.tpr

echo 1 | gmx trjconv -f trajout.xtc -o protein.pdb -dump 0 -s topol.tpr

echo 1 | gmx trjconv -f tmp3.xtc -o protein.xtc -s topol.tpr

rm tmp2.xtc tmp.xtc

mv tmp3.xtc cell.xtc 

rm \#*

cd ..

