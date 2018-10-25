# Xpa
Programming code for "Applying the Explicit Aggregation Algorithm to Heterogeneous Macro Models"
October 2018 by Takeki Sunakawa

The code is intended to be compiled by intel fortran 18.0.3. 

In Linux, the code can be compiled and run by
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

You need to install 

- json-fortran https://github.com/jacobwilliams/json-fortran

and if you use eigenvalue and eigenvector decomposition to calculate the stationary distribution, you also need

- arpack-ng https://github.com/opencollab/arpack-ng

To be completed.
