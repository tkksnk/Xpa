## Xpa
Programming code for "Applying the Explicit Aggregation Algorithm to Heterogeneous Macro Models"

October 2018 by Takeki Sunakawa

# How to use

The code is intended to be compiled by Intel Fortran with Math Kernel Library. In Linux, the code can be compiled and run by
```
$ make solveKS
$ ./solveKS
```
for the KS algorithm and
```
$ make solveXpa
$ ./solveXpa
```
for the Xpa algorithm. 

The results are saved in .json format. Matlab (>R2016b) is needed to display the results.

# Prerequests

You need to install

- json-fortran https://github.com/jacobwilliams/json-fortran

and if you use eigenvalue and eigenvector decomposition to calculate the stationary distribution, you also need

- arpack-ng https://github.com/opencollab/arpack-ng

After cloning the arpack-ng repository, run get_arpack_ng.sh (based on https://gist.github.com/galanakis/4069435) to compile libarpack.a.
