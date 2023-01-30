#!/bin/bash
#SBATCH --image=docker:avariu/hodor_new:latest
#SBATCH -n 8	      # Number of tasks
#SBATCH -J img       # Name of the job
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -N 1         # number of nodes
#SBATCH -c 8          # number of cpus per tasks
#SBATCH --time=00:05:00
#SBATCH -o ./out.img.out
#SBATCH -e ./err.img.err


export HDF5_USE_FILE_LOCKING=FALSE

srun -n 8 -c 8 shifter --env=LD_LIBRARY_PATH=/home/apps/MultiNest/MultiNest_v3.12_CMake/multinest/lib:/home/apps/fftw-3.3.10/lib:$LD_LIBRARY_PATH --env=PYTHONPATH=/home/pypowspec:/home/pyfcfc:$PYTHONPATH python3 main_img.py --config configs/config_RSD_diff_1296_0.02_6pars_pk_0.02_0.4.ini
