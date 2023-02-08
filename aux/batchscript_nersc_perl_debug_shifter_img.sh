#!/bin/bash

#SBATCH --image=docker:avariu/hodor_new:latest
#SBATCH -n 8           # Number of tasks
#SBATCH -J debugH      # Name of the job
#SBATCH --qos debug
#SBATCH -C cpu
#SBATCH -N 1           # number of nodes
#SBATCH -c 8           # number of cpus per tasks
#SBATCH --time=00:10:00
#SBATCH -o /pscratch/sd/a/avariu/out.debug.out
#SBATCH -e /pscratch/sd/a/avariu/err.debug.err

export HDF5_USE_FILE_LOCKING=FALSE

# CONFIGFILE=../configs/config_RSD_diff_1296_0.02_6pars_pk_0.02_0.5_wospikes.ini
# CONFIGFILE=../configs/config_RSD_diff_1296_0.02_6pars_pk_0.02_0.4.ini
CONFIGFILE=../configs/config_RSD_diff_1296_0.02_6pars_pk_0.02_0.3.ini
# CONFIGFILE=../configs/config_RSD_diff_1296_6pars_cf_15_50.ini

srun -n 8 -c 8 shifter --env=LD_LIBRARY_PATH=/home/apps/MultiNest/MultiNest_v3.12_CMake/multinest/lib:/home/apps/fftw-3.3.10/lib --env=PYTHONPATH=/home/pypowspec:/home/pyfcfc --env=PATH=/opt/miniconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin python3 ../main.py --config $CONFIGFILE
