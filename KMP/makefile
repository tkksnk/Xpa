# Makefile
FC = ifort
FFLAGS = -O3 -qopenmp
# FFLAGS = -O3 -nocheck -inline-level=2 -xCORE-AVX2 -qopenmp -mkl
# FFLAGS = -g -check all -warn all
LAPACK = -mkl
JSON=-I/home/takeki/lib/json-fortran/lib/intel/mod -L/home/takeki/lib/json-fortran/lib/intel -ljsonfortran

%.o: %.f90
	$(FC) $(FFLAGS) -c  $<

solveXpa: mkl_spblas.o mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_inner.o mod_outer.o solveXpa.f90
	$(FC) $(FFLAGS) $^ -o solveXpa /home/takeki/lib/libarpack.a $(LAPACK) $(JSON)

solveKS: mkl_spblas.o mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_inner.o mod_outer.o solveKS.f90
	$(FC) $(FFLAGS) $^ -o solveKS /home/takeki/lib/libarpack.a $(LAPACK) $(JSON)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *~
