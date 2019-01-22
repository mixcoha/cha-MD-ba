#!/bin/bash
#
#Este script genera una caja dodecaedrica con una distancia entre la proteína y la caja de 1.5 nm con el nombre del argumento de la carpeta a trabajar con el argumento $1

# Se define el directorio actual como directorio de trabajo
RUNDIR=$PWD

# Se crea una carpeta de trabajo con el nombre de la proteina en el argumento $1
mkdir $RUNDIR/$1
cd $RUNDIR/$1/

cat > em.mdp < estodo

; Líneas que comienzan con ';' son considerados comentarios
title		= Minimization de la proteina $1 en agua

; Parameters describing what to do, when to stop and what to save
integrator	= steep		 ; Algorithm options i.e. steep = steepest descent minimization 
emtol		    = 1000.0	 ; Stop minimization when the energy changes by less than emtol kJ/mol.
emstep     = 0.01    ; Energy step size
nsteps	   	= 1000		  ; Maximum number of (minimization) steps to perform
nstenergy 	= 10	  	  ; Write energies to disk every nstenergy steps
nstxtcout 	= 10		    ; Write coordinates to disk every nstxtcout steps
xtc_grps	  = Protein	; Which coordinate group(s) to write to disk
energygrps	= Protein	; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist	  	 = 1		     ; Frequency to update the neighbor list and long range forces
ns_type	 	  = grid    ; Method to determine neighbor list (simple, grid)
rlist	    	 = 1.0	  	 ; Cut-off for making neighbor list (short range forces)
coulombtype	= cut-off	; Treatment of long range electrostatic interactions
rcoulomb	   = 1.0	   	; long range electrostatic cut-off
rvdw		      = 1.0		   ; long range Van der Waals cut-off
constraints	= none		  ; Bond types to replace by constraints
pbc		       = xyz		   ; Periodic Boundary Conditions (yes/no)

estodo

# Al revisar los archivos PDB podemos encontrar que faltan átomos de algunos residuos por lo que es deseable procesarlos previamente
# en un el servidor PDB2PQR http://nbcr-222.ucsd.edu/pdb2pqr_2.1.1/ cuando se obtenga el resultado, se guardará como PDB en la carpeta
# de trabajo.
cp $RUNDIR/$1.pqr $RUNDIR/$1/$1_1.pdb

# Se crea una nueva carpeta donde se generarán los archivos para realizar las simulaciones.
mkdir $RUNDIR/$1/setup
mv $RUNDIR/$1/em.mdp $RUNDIR/$1/setup
mv $RUNDIR/$1/$1_1.pdb $RUNDIR/$1/setup
cd $RUNDIR/$1/setup

# Se elige el campo de Fuerzas AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006) y aguas tipo SPC
echo 5 | gmx pdb2gmx -f $1_1.pdb -ignh -water spc -o $1_1.gro

# Se genera una caja dodecaédrica con una distancia entre la proteína y la caja de 1.5 nm
gmx editconf -bt dodecahedron -d 1.5 -c -f $1_1.gro -o $1_1_box.gro

# Se solvata la proteína con aguas tipo SPC
gmx solvate -cp $1_1_box.gro -cs spc216.gro -p topol.top -o prot_wat.gro

# Se genera un archivo binario TPR para neutralzar el sistema utilizando el archivo em.mdp generado previamente.
gmx grompp -f em.mdp -c prot_wat.gro -p topol.top -o -v

# Se elige al disolvente como medio contínuo para sustituir moléculas de agua con iones Na+ o Cl- para neutralizar el sistema.
# así mismo se agregarán suficientes iones para obtener una concentración de 0.15 M, equivalente a la fisiológica.
echo 13 | gmx genion -s topol.tpr -p topol.top -o prot_wat_ion.gro -neutral -nname CL -pname NA -conc 0.15

# Procesamos los archivos de la proteina con agua e iones para comenzar la minimización
gmx grompp -f em.mdp -c prot_wat_ion.gro -p topol.top  -o -v

rm -rf \#*

# Generamos una nueva carpeta llama em (Energy Minimzation) donde se escribirán los archivos de la minimización.
# Nos movemos a esa carpeta
cd $RUNDIR/$1/
mkdir $RUNDIR/$1/em
mv $RUNDIR/$1/setup/topol.tpr $RUNDIR/$1/em
cd $RUNDIR/$1/em

# Comenzamos a realizar la minización del sistema utilizando las 2 GPUS de la computadora. Generamos la trayectoria de minimización.
gmx mdrun -gpu_id 01 -nb gpu_cpu -tunepme -v -s

# Se eliminan los archivos respaldados
rm -rf \#*

# Proceso para centrar la proteína en la trayectoria de minimización considerando condiciones periódicas a la frontera
echo 1 0 | gmx trjconv -f confout.gro -s topol.tpr -pbc mol -center -o tmp.gro
 
echo 0 | gmx trjconv -f tmp.gro -ur compact -pbc mol -o tmp2.gro -s topol.tpr

echo 1 0 | gmx trjconv -f tmp2.gro -pbc cluster -o tmp.gro -s topol.tpr

echo 1 0 | gmx trjconv -f tmp.gro -ur compact -pbc mol -center -o tmp2.gro -s topol.tpr

echo 4 0 | gmx trjconv -f tmp2.gro -fit rot+trans -o cell.pdb -s topol.tpr

rm tmp*.gro \#*

cd $RUNDIR/

