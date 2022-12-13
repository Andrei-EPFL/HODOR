#!/bin/bash
#SBATCH -n 	8	# Number of tasks
#SBATCH -J test     # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1          # number of nodes
#SBATCH -c 8          # number of cpus per tasks
#SBATCH --time=00:30:00
#SBATCH -o ./out.out
#SBATCH -e ./err.err

module swap PrgEnv-intel PrgEnv-gnu/6.0.5
module load gsl/2.5
export LD_LIBRARY_PATH=/global/homes/a/avariu/apps/MultiNest/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=8
source activate /global/cscratch1/sd/avariu/hodfit/
export HDF5_USE_FILE_LOCKING=FALSE

srun -n 8 python main.py --config configs/config.ini
