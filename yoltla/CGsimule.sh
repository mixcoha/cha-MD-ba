#!/bin/sh

#  simule.sh
#  
#
#  Created by Edgar Mixcoha on 23/01/18.
#  
RUNDIR=$PWD

cd $RUNDIR/$1

# Start Simulation

mkdir dynamics

cd $RUNDIR/

cat > ext_1_$1_en$2_$3p.slmr << estodo
#!/bin/bash
#SBATCH --job-name=$1
#SBATCH -p $2
#SBATCH -n $3
#SBATCH -e ${nd}_$1_error_%J
#SBATCH -o ${nd}_$1salida_%J

cd $RUNDIR/$1/dynamics

module purge
module load intel/15.6.232/impi/5.0.3.49/gromacs/5.1.4-s

export OMP_NUM_THREADS=5

# Only 20 instances, 4 process and 5 threads per process: recommend for this version.
# mpirun -n 4 gmx mdrun -v -deffnm file

# > 20 processes
task=\$(( SLURM_NNODES * 4 ))
srun -n \${task} hostname -s | sort > hostlist.dat
# mpirun -np 4 -hostfile ./hostlist.dat gmx mdrun -v -cpi >& mdrun.log
mpiexec.hydra -ppn 4 -hostfile ./hostlist.dat gmx mdrun -v -cpi -noappend >& mdrun.log

exit 0

estodo

cd $RUNDIR/$1/dynamics

module purge
module load intel/15.6.232/impi/5.0.3.49/gromacs/5.1.4-s

gmx grompp -f $RUNDIR/mdps/dynamic.mdp -p $RUNDIR/$1/setup/system.top -c $RUNDIR/$1/equilibration/eq.gro -n $RUNDIR/$1/setup/index.ndx -o topol.tpr -maxwarn 5

cd $RUNDIR

chmod +x ext_1_$1_en$2_$3p.slmr

sbatch ext_1_$1_en$2_$3p.slmr

