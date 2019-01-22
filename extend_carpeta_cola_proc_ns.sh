#!/bin/bash
# _______________
# PREPARANDO LA EQUILIBRACIÓN
# $1 carpeta del calculo
# $2 cola donde se envia el cálculo
# $3 numero de procesadores
# $4 tiempo para extender (ps)
# ————————

RUNDIR=$PWD

cd $RUNDIR/$1/R1/prod

                ls -lhtr > ls
                grep drwxr ls > tmp
                cat tmp | sort -k 9 -g > tmp1
                mv tmp1 tmp
                nd=`awk 'END {print $NF+1}' tmp`
                rm ls tmp
                mkdir $nd
                mv *.* $nd
                cd $nd

cd $RUNDIR/

cat > ext_${nd}_$1_en$2_$3p_$4ns.slmr << estodo
#!/bin/bash
#SBATCH --job-name=$1
#SBATCH -p $2
#SBATCH -n $3
#SBATCH -e ${nd}_$1_error_%J
#SBATCH -o ${nd}_$1salida_%J


cd $RUNDIR/$1/R1/prod

module purge
module load intel/15.6.232/impi/5.0.3.49/gromacs/5.1.4-s

export OMP_NUM_THREADS=5

# Only 20 instances, 4 process and 5 threads per process: recommend for this version.
# mpirun -n 4 gmx mdrun -v -deffnm file

# > 20 processes
task=$(( SLURM_NNODES * 4 ))
srun -n $task hostname -s | sort > hostlist.dat
mpirun -np 4 -hostfile ./hostlist.dat gmx mdrun -v -cpi >& mdrun.log

exit 0

estodo

cd $RUNDIR/$1/R1/prod/

module purge
module load intel/15.6.232/impi/5.0.3.49/gromacs/5.1.4-s

                gmx convert-tpr -s $RUNDIR/$1/R1/prod/$nd/topol.tpr -f $RUNDIR/$1/R1/prod/$nd/state.cpt -o -extend $4
                mv tprout.tpr topol.tpr

                cd $RUNDIR
		chmod +x ext_${nd}_$1_en$2_$3p_$4ns.slmr

sbatch ext_${nd}_$1_en$2_$3p_$4ns.slmr

