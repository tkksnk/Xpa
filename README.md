# Programming code for "Applying the Explicit Aggregation Algorithm to Heterogeneous Macro Models"

by Takeki Sunakawa (takeki.sunakawa@gmail.com), October 2018

## How to use

The folder `KT` is for Khan-Thomas models and the folder `KMP` is for Krueger-Mitman-Perri models. The code is intended to be compiled by Intel Fortran with Math Kernel Library. In Linux, the code can be compiled and run by
```
$ make solveKS
$ ./solveKS
```
for the KS algorithm and
```
$ make solveXpa
$ ./solveXpa
```
for the Xpa algorithm. The results in Fortran are saved in .json format and Matlab (>=R2016b) is needed to visualize them. Figures and Tables are displayed by running aggresults.m in each KT and KMP folder.

## Required libraries

You need to install

- json-fortran https://github.com/jacobwilliams/json-fortran

If you use the sparse-grid eigenvalue and eigenvector decomposition to calculate the stationary distribution, you also need

- arpack-ng https://github.com/opencollab/arpack-ng

After cloning the arpack-ng repository, run get_arpack_ng.sh (based on https://gist.github.com/galanakis/4069435) to compile libarpack.a.
