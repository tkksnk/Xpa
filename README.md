# Xpa
Programming code for "Applying the Explicit Aggregation Algorithm to Heterogeneous Macro Models"

The code is written for intel fortran (or gfortran). 

For Linux, you can compile the code by
```
$ make solveKS
```
for the KS algorithm and
```
$ make solveXpa
```
for the Xpa algorithm. 

Matlab (>R2016b) is needed to display the results.

You need to install 

- json-fortran https://github.com/jacobwilliams/json-fortran

and if you use eigenvalue and eigenvector decomposition to calculate the stationary distribution,

- arpack-ng https://github.com/opencollab/arpack-ng

To be written.
