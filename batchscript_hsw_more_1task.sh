#!/bin/bash
#SBATCH -n 1		# Number of tasks
#SBATCH -J 1296_pk     # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1          # number of nodes
#SBATCH -c 64          # number of cpus per tasks
#SBATCH --time=48:00:00
#SBATCH -o ./output/out.para_1296_1000_0.5_64_48_m_DIF_cov_vdisp_pk_1task.out
#SBATCH -e ./output/err.para_1296_1000_0.5_64_48_m_DIF_cov_vdisp_pk_1task.err

module swap PrgEnv-intel PrgEnv-gnu/6.0.5
module load gsl/2.5
export LD_LIBRARY_PATH=/global/homes/a/avariu/apps/MultiNest/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=64
source activate /global/cscratch1/sd/avariu/envs/
export HDF5_USE_FILE_LOCKING=FALSE

srun -n 1 python main.py --config ./configs/config_RSD_diff_1296_1000_0.5_6pars_pk.ini
