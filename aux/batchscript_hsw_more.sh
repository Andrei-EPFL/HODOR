#!/bin/bash
#SBATCH -n 128		# Number of tasks
#SBATCH -J 1536_pk4     # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 16          # number of nodes
#SBATCH -c 8          # number of cpus per tasks
#SBATCH --time=30:00:00
#SBATCH -o ./output/out.para_1536_1000_0.5_DIFF_cov_vdisp_pk4.out
#SBATCH -e ./output/err.para_1536_1000_0.5_DIFF_cov_vdisp_pk4.err

module swap PrgEnv-intel PrgEnv-gnu/6.0.5
module load gsl/2.5
export LD_LIBRARY_PATH=/global/homes/a/avariu/apps/MultiNest/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
source activate /global/cscratch1/sd/avariu/hodfit/
export HDF5_USE_FILE_LOCKING=FALSE

# srun -n 128 python main.py --config configs/FastPM4SLICS/config_RSD_id_1536_1000_0.5_6pars_cf.ini
srun -n 128 python main.py --config configs/FastPM4SLICS/config_RSD_diff_1536_1000_0.5_6pars_pk_0.02_0.4
