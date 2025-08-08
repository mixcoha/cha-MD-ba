#!/bin/bash

cd $1/dynamics

mkdir cell

gmx trjcat -f */traj_*.xtc -o trajout.xtc 
gmx eneconv -f */*.edr

mv trajout.xtc cell/
mv *.edr cell/

cp 1/topol.tpr cell/tprout.tpr 

cd cell

echo -e "1|21 \n q \n" | gmx make_ndx -f tprout.tpr

echo 1 0 | gmx trjconv -f trajout.xtc -s tprout.tpr -pbc mol -center -o tmp.xtc
 
echo 0 | gmx trjconv -f tmp.xtc -ur compact -pbc mol -o tmp2.xtc -s tprout.tpr

echo 1 0 | gmx trjconv -f tmp2.xtc -pbc cluster -o tmp.xtc -s tprout.tpr

echo 1 0 | gmx trjconv -f trajout.xtc -ur compact -pbc mol -center -o tmp2.xtc -s tprout.tpr

echo 1 0 | gmx trjconv -f tmp2.xtc -fit rot+trans -o cell.xtc -s tprout.tpr -tu ns

echo 0 | gmx trjconv -f cell.xtc -o cell.pdb -b 0 -e 0 -s tprout.tpr -tu ps 

#gmx make_ndx -f ../1/confout.gro -o index.ndx

echo 1 | gmx trjconv -f cell.xtc -o protein.xtc -s tprout.tpr -tu ns 
echo 1 | gmx trjconv -f protein.xtc -o protein.pdb -dump 0 -s tprout.tpr



rm tmp2.xtc tmp.xtc

rm \#*

cd ..

