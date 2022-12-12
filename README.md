# Halo Occupation Distribution Optimization Routine

## Table of Contents

-   [Introduction](#introduction)
-   [Dependencies](#Dependencies)
-   [Usage](#Usage)
-   [Configuration parameters](#Configuration-parameters)
-   [Co-Developers](#Co-Developers)
-   [Acknowledgements](#Acknowledgements)

## Introduction

**HOD** **O**ptimization **R**outine (HODOR).

... To be written ...


If you use this tool in research work that results in publications, please cite the following paper:

> Variu et al., in preparation (2023).

## Dependencies

The dependencies of this code are as follows:

-   [Python3](https://www.python.org/)  (>= 3.6)
-   [NumPy](https://numpy.org/)
-   [halotools](https://anaconda.org/conda-forge/halotools)
-   [h5py](https://anaconda.org/conda-forge/h5py)
-   [pypowspec](https://github.com/dforero0896/pypowspec)
-   [pyfcfc](https://github.com/dforero0896/pyfcfc)
-   [MultiNest](https://github.com/farhanferoz/MultiNest)
-   [PyMultiNest](https://github.com/JohannesBuchner/PyMultiNest)

Optional:
-   [iminuit](https://anaconda.org/conda-forge/iminuit)
-   [SciPy](https://scipy.org/)

One can replace pyfcfc or pypowspec by different clustering codes. In [compute_2pcf_old.py](hod_pack/alternatives/compute_2pcf_old.py) file there is different ways to compute the 2PCF, e.g. [corrfunc](https://corrfunc.readthedocs.io/en/master/) or [halotools](https://anaconda.org/conda-forge/halotools). Similarly for P(k), in [compute_pspec_old_1.py](hod_pack/alternatives/compute_pspec_old_1.py) and [compute_pspec_old_2.py](hod_pack/alternatives/compute_pspec_old_2.py).

## Usage

After downloading the code, one has to fill the configuration file and run [main.py](main.py):

```bash
python main.py --config <CONFIG FILE>
```


## Configuration parameters
All configuration parameters should be given in the configuration file. Examples of configuration files can be found in the ([configs/FastPM4SLICS/](configs/FastPM4SLICS/)) folders.



## Co-Developers
[Dr. Cheng Zhao](https://github.com/cheng-zhao/)

[Dr. Chia-Hsun Chuang](https://github.com/chia-hsun-chuang/)
## Acknowledgements
I thank [Dr. Shadab Alam](https://github.com/shadaba) for his suggestions.
