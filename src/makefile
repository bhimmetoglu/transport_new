#
#
F90 = gfortran
FFLAGS = -fopenmp 
FFLAGS2 = -fbounds-check

ef: ef.f90
	${F90} ${FFLAGS} cryst_to_car.f90 lint.f90 ef.f90 -o ef.x

dos: dos.f90
	${F90} ${FFLAGS} smearing_mod.f90 cryst_to_car.f90 lint.f90 vband_ibz.f90 dos.f90 -o dos.x

tauef: tauef.f90
	${F90} ${FFLAGS} smearing_mod.f90 cryst_to_car.f90 lint.f90 vband_ibz.f90 vband_kq.f90 kplusq.f90 reducegrid.f90 invtau_nk.f90 tauef.f90 -o tauef.x

tauBZ_plt: tauBZ_plt.f90
	${F90} smearing_mod.f90 cryst_to_car.f90 lint.f90 vband_ibz.f90 vband_kq.f90 kplusq.f90 invtau_nk.f90 tauBZ_plt.f90 -o tauBZ_plt.x

fermi_int_0: fermi_int_0.f90
	${F90} ${FFLAGS} cryst_to_car.f90 lint.f90 vband_ibz.f90 fermi_int_0.f90 -o fermi_int_0.x

fermi_int_1: fermi_int_1.f90
	${F90} ${FFLAGS} cryst_to_car.f90 reducegrid.f90 lint.f90 vband_ibz.f90 fermi_int_1.f90 -o fermi_int_1.x

fermi_int: fermi_int.f90
	${F90} ${FFLAGS} smearing_mod.f90 cryst_to_car.f90 reducegrid.f90 lint.f90 vband_ibz.f90 vband_kq.f90 kplusq.f90 invtau_nk.f90 invtau_nk_Z.f90 fermi_int.f90 -o fermi_int.x
