FC=gfortran
SRCS=gridrefine.f90 derivs.f90 remove.f90 rdot.f90 dsigma_disk.f90 torque_vel_mod.f90 dsigma_solids.f90 interp1d.f90 timestep.f90 remove_dust.f90 macc_rates.f90 write_2d.f90 sigma_init.f90 calculate_quantities.f90 tstep_dust.f90 calc_idx.f90 calc_embryo.f90
#FCFLAGS=-g -fdefault-real-8 -fbounds-check -O3 -Wall -std=legacy -fbacktrace -ffpe-trap=invalid
FCFLAGS= -O5 -fdefault-real-8 -std=legacy

.PHONY: bdwarf
bdwarf:
	$(FC) $(FCFLAGS) -o bdwarf main.f90 $(SRCS)

