# Makefile
FC = ifort
FFLAGS = -O3 -qopenmp
LAPACK = -mkl
JSON=-I/home/takeki/lib/json-fortran/lib/intel/mod -L/home/takeki/lib/json-fortran/lib/intel -ljsonfortran

%.o: %.f90
	$(FC) $(FFLAGS) -c  $<

solveKS: mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_planner.o mod_inner.o mod_simulation.o solveKS.f90
	$(FC) $(FFLAGS) $^ -o solveKS $(LAPACK) $(JSON)

solveXpa: mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_planner.o mod_inner.o mod_outer.o mod_simulation.o solveXpa.f90
	$(FC) $(FFLAGS) $^ -o solveXpa $(LAPACK) $(JSON)

calcirf: mod_parameters.o mod_functions.o mod_spline.o mod_calcss.o mod_planner.o mod_inner.o mod_simulation.o calcirf.f90
	$(FC) $(FFLAGS) $^ -o calcirf $(LAPACK) $(JSON)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *~
