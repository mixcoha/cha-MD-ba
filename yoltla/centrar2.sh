#!/bin/bash
# _______________
# PREPARANDO LA EQUILIBRACIÓN
# $1 carpeta del calculo

RUNDIR=$PWD

cd $RUNDIR/$1/dynamics/cell/

mkdir cell

module purge
module load intel/15.6.232/impi/5.0.3.49/gromacs/5.1.4-s

#echo -e "1|21 \n q \n" | gmx make_ndx -f tprout.tpr

gmx eneconv -f $RUNDIR/$1/dynamics/*/*.edr

cat > center_$1.slmr << estodo
#!/bin/bash
#SBATCH --job-name=$1_centrar
#SBATCH -p q12h-160p 
#SBATCH -n 20 
#SBATCH -e $1_error_center_%J
#SBATCH -o $1salida_center_%J

cd $RUNDIR/$1/dynamics/cell

module purge
module load intel/15.6.232/impi/5.0.3.49/gromacs/5.1.4-s

export OMP_NUM_THREADS=5

# Only 20 instances, 4 process and 5 threads per process: recommend for this version.

mpirun -n 4 gmx trjcat -f $RUNDIR/$1/dynamics/*/traj_*.xtc -o trajout.xtc

echo 1 0 | mpirun -n 4 gmx trjconv -f trajout.xtc -ur compact -pbc mol -center -o tmp2.xtc -s tprout.tpr

echo 1 0 | mpirun -n 4 gmx trjconv -f tmp2.xtc -fit rot+trans -o cell.xtc -s tprout.tpr -tu ns

rm tmp2.xtc tmp.xtc

rm \#*

exit 0

estodo

cd $RUNDIR/$1/dynamics/cell/
chmod +x center_$1.slmr

sbatch center_$1.slmr

