# Makefile
FC = ifort
FFLAGS = -O3 -qopenmp -mkl
# FFLAGS = -O3 -nocheck -inline-level=2 -xCORE-AVX2 -qopenmp -mkl
# FFLAGS = -O3 -nocheck -inline-level=2 -xCORE-AVX2 -mkl
# FFLAGS = -g -check all -warn all
JSON=-I/home/takeki/lib/json-fortran/lib/intel/mod -L/home/takeki/lib/json-fortran/lib/intel -ljsonfortran

%.o: %.f90
	$(FC) $(FFLAGS) -c  $< $(JSON)

solvess: mkl_spblas.o mod_parameters.o mod_randnumber.o mod_utils.o mod_spline_v4.o mod_policies_v4.o mod_calcss_v5.o solvess.o
	$(FC) $(FFLAGS) $^ -o solvess /home/takeki/lib/libarpack.a $(JSON)

solveKS: mkl_spblas.o mod_parameters.o mod_randnumber.o mod_utils.o mod_spline_v4.o mod_policies_v4.o mod_calcss_v5.o mod_inner.o mod_outer.o solveKS.o
	$(FC) $(FFLAGS) $^ -o solveKS /home/takeki/lib/libarpack.a $(JSON)

solveXpa: mkl_spblas.o mod_parameters.o mod_randnumber.o mod_utils.o mod_spline_v4.o mod_policies_v4.o mod_calcss_v5.o mod_inner.o mod_outer.o solveXpa.o
	$(FC) $(FFLAGS) $^ -o solveXpa /home/takeki/lib/libarpack.a $(JSON)

# solveReiter: solab.o mkl_spblas.o mod_parameters.o mod_randnumber.o mod_utils.o mod_spline_v4.o mod_policies_v4.o mod_calcss_v5.o mod_reiter.o solveReiter.o
# 	$(FC) $(FFLAGS) $^ -o solveReiter /home/takeki/lib/libarpack.a $(JSON)

# solve: mod_parameters.o mod_randnumber.o mod_utils.o mod_spline_v4.o mod_policies_v4.o mod_calcss_v4.o mod_inner.o mod_outer.o solve.o
# 	$(FC) $(FFLAGS) $^ -o solve
#
# calcirf: mod_parameters.o mod_randnumber.o mod_utils.o mod_spline_v4.o mod_policies_v4.o mod_calcss_v4.o mod_inner.o mod_outer.o mod_irf.o calcirf.o
# 	$(FC) $(FFLAGS) $^ -o calcirf
#
# calcirf_sim: mod_parameters.o mod_randnumber.o mod_utils.o mod_spline_v3.o mod_policies_v3.o mod_calcss.o mod_inner.o mod_outer_sim.o mod_irf_sim.o calcirf_sim.o
# 	$(FC) $(FFLAGS) $^ -o calcirf_sim

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *~
