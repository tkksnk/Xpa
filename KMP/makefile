# Makefile
FC = ifort
FFLAGS = -O3 -qopenmp -fpp
# FFLAGS = -O3 -nocheck -inline-level=2 -xCORE-AVX2 -qopenmp -mkl
# FFLAGS = -g -check all -warn all
LAPACK = -mkl
# ARPACK = -L/home/takeki/lib
INCLUDE = "/opt/MATLAB/R2018a/extern/include"
JSON=-I/home/takeki/lib/json-fortran/lib/intel/mod -L/home/takeki/lib/json-fortran/lib/intel -ljsonfortran

# FC = gfortran
# FFLAGS = -O3 -fopenmp -ffree-line-length-none -cpp
# # FFLAGS = "-fopenmp -ffree-line-length-none -cpp -fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -DMATLAB_MEX_FILE "
# LAPACK = -larpack -llapack -lblas
# INCLUDE = "/opt/MATLAB/R2018a/extern/include"
# JSON=-I/home/takeki/lib/json-fortran/lib/gnu/mod -L/home/takeki/lib/json-fortran/lib/gnu -ljsonfortran

%.o: %.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) -c  $<

solveSS: mkl_spblas.o mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o solveSS.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) $^ -o solveSS /home/takeki/lib/libarpack.a $(LAPACK) $(JSON)

solveXpa: mkl_spblas.o mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_inner.o mod_outer.o solveXpa.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) $^ -o solveXpa /home/takeki/lib/libarpack.a $(LAPACK) $(JSON)

solveKS: mkl_spblas.o mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_inner.o mod_outer.o solveKS.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) $^ -o solveKS /home/takeki/lib/libarpack.a $(LAPACK) $(JSON)

solveReiter: mkl_spblas.o mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_reiter.o solveReiter.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) $^ -o solveReiter /home/takeki/lib/libarpack.a $(LAPACK) $(JSON)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *~
