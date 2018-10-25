# Makefile
# TODO, Oct 25 2018: 1. remove include and -fpp (MATLAB related commands),
# 2. integrate calcirf into solve,
# 3. implement spblas in calcss
FC = ifort
FFLAGS = -O3 -qopenmp -fpp
LAPACK = -mkl
INCLUDE = "/opt/MATLAB/R2018a/extern/include"
JSON=-I/home/takeki/lib/json-fortran/lib/intel/mod -L/home/takeki/lib/json-fortran/lib/intel -ljsonfortran

# FC = gfortran
# FFLAGS = -O3 -fopenmp -ffree-line-length-none -cpp
# LAPACK = -llapack -lblas
# INCLUDE = "/usr/local/MATLAB/R2016b/extern/include"
# JSON=-I/home/takeki/lib/json-fortran/lib/gnu/mod -L/home/takeki/lib/json-fortran/lib/gnu -ljsonfortran

%.o: %.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) -c  $<

solveKS: mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_planner.o mod_inner.o mod_simulation.o solveKS.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) $^ -o solveKS $(LAPACK) $(JSON)

solveXpa: mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_planner.o mod_inner.o mod_outer.o mod_simulation.o solveXpa.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) $^ -o solveXpa $(LAPACK) $(JSON)

calcirf: mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_planner.o mod_inner.o mod_simulation.o calcirf.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) $^ -o calcirf $(LAPACK) $(JSON)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *~
