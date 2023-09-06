# Halo Occupation Distribution Optimization Routine

## Table of Contents

-   [Introduction](#introduction)
-   [Dependencies](#Dependencies)
-   [Usage](#Usage)
-   [Configuration parameters](#Configuration-parameters)
-   [Co-Developers](#Co-Developers)
-   [Acknowledgements](#Acknowledgements)

## Introduction

**HOD** **O**ptimization **R**outine (HODOR) is an adaptable python code that fits and HOD model on halo catalogues.
Currently, there are two HOD implemented models: the vanila 5-parameters one and the 9-parameters one from Contreras (1301.3497). Additionaly, the velocity dispersion along the line-of-sight is added and cannot be removed (on ToDo list)
The code has been conceived to by highly adaptable to different data-sets, different clustering (2PCF and P(k)) codes and different HOD models.


If you use this tool in research work that results in publications, please cite the following paper:

> Variu et al. 2023 [arXiv:2307.14197](https://arxiv.org/abs/2307.14197).

## Dependencies

The dependencies of this code are as follows:

-   [Python3](https://www.python.org/)  (>= 3.6)
-   [NumPy](https://numpy.org/)
-   [halotools](https://anaconda.org/conda-forge/halotools)
-   [h5py](https://anaconda.org/conda-forge/h5py)
-   [pyfcfc](https://github.com/dforero0896/pyfcfc)
-   [pypowspec](https://github.com/dforero0896/pypowspec)
-   [MultiNest](https://github.com/farhanferoz/MultiNest)
-   [PyMultiNest](https://github.com/JohannesBuchner/PyMultiNest)
Optional:
-   [iminuit](https://anaconda.org/conda-forge/iminuit)
-   [SciPy](https://scipy.org/)

One can replace pyfcfc or pypowspec by different clustering codes. In [compute_2pcf_old.py](hod_pack/alternatives/compute_2pcf_old.py) file there is different ways to compute the 2PCF, e.g. [corrfunc](https://corrfunc.readthedocs.io/en/master/) or [halotools](https://anaconda.org/conda-forge/halotools). Similarly for P(k), in [compute_pspec_old_1.py](hod_pack/alternatives/compute_pspec_old_1.py) and [compute_pspec_old_2.py](hod_pack/alternatives/compute_pspec_old_2.py).

All packages are installed in a [docker](https://www.docker.com/) image that one can pull from [here](https://hub.docker.com/repository/docker/avariu/hodor_new/general).
For a tutorial on how to use docker images on the NERSC super-computer, one can read [this](https://docs.nersc.gov/development/shifter/how-to-use/). However, this [batch-script](aux/batchscript_nersc_perl_shifter_img.sh) is already using the docker image on NERSC.

## Usage

After downloading the code, one has to fill in the configuration file and run [main.py](main.py):

```bash
python main.py --config <CONFIG FILE>
```

By default, the code will run PyMultiNest. However a short test is also provided in the same main file. As an initial check, if the chi2 computation is working and provides reasonable values, a large chunk of the code should work properly. Then, one can try the code together with a minimizer.

The following file
```bash
python main_cat.py
```
is conceived to create the catalogues based on the best-fitting parameters.

Currently, the 'model' catalogues (e.g. FastPM) in [hod_pack/data.py](#hod_pack/data.py) should be and HDF5 file with the following data-set:
```bash
['halo/PID', 'halo/X', 'halo/Y', 'halo/Z', 'halo/VX', 'halo/VY', 'halo/VZ', 'halo/ID', 'halo/Rvir', 'halo/Rs', 'halo/Mvir'].
```
However, the format of the file is not important, as long as the necessary fields are provided. If another format is required, one should modify accordingly the function 
```bash
compute_halocat()
```
in the [hod_pack/data.py](#hod_pack/data.py) file.

If one wants to use different clustering codes, one should modify accordingly the  [hod_pack/compute_2pcf.py](#hod_pack/compute_2pcf.py) [hod_pack/compute_pspec.py](#hod_pack/compute_pspec.py).

If one wants to implement additional HOD models, one should maintain the same structure as the ones in [hod_pack/hod_models.py](#hod_pack/hod_models.py) and then one should modify accordingly the
```bash
compute_model_instance()
```
function in the [hod_pack/hod_models.py](#hod_pack/hod_models.py) file.

If one wants to use a different optimizer, one can simply add it in the [hod_pack/optimizers.py](#hod_pack/optimizers.py) as the MultiNest example, and remove the unnecessary ones.


...To be developped ...

## Configuration parameters
All configuration parameters should be given in the configuration file. An example configuration is [configs/config.ini](configs/config.ini).



## Co-Developers
[Dr. Cheng Zhao](https://github.com/cheng-zhao/)


## Acknowledgements
I thank [Dr. Shadab Alam](https://github.com/shadaba), [Dr. Chia-Hsun Chuang](https://github.com/chia-hsun-chuang/) and [Daniel Felipe Forero-Sánchez](https://github.com/dforero0896) for their suggestions.
