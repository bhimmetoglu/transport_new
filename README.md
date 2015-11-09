# Boltzmann transport code
## Main Description
This package contains codes to calculate transport properties using Boltzmann Transport theory. 
There are three distinct features of the code:

1. Calculation of scattering rates based on electron-phonon(LO only) interactions. A simple analytical formula is used.
2. Adaptive smearing width choice for the evaluation of energy conserving Dirac delta functions. See documentaion.
3. Adaptive k-grid around the Fermi level to reduce integration time. See documentation.

The code is parallelized using OpenMP. This means you can only use a single node of a cluster for parallelization. Recommended setting is to use nthreads=number of cores in the node. However, you can use nthreads exceeding the number of cores and may get some performance increase. Perform some tests before doing so (In Stampede, I have seen performance increase up to 48 threads. It is probably due to the multi-threading available in the processors).

## Description of each routine in /src
* ef.f90 : Computes Fermi levels for a given doping level in cm-3. 
* dos.f90 : Computes density of states using the adaptive smearing scheme. Recall that DOS is given by
  DOS(e) = sum_{n, k}\, \delta( e - e_{nk} )
* tauef.f90 : This code computes the scattering rate at the Fermi level:
  tau^{-1}(E_F) = \sum_{n,k}\, tau^{-1}_{nk}\, delta( e_{nk} - E_F)
* tauBZ_plt : Visualization code for plotting the scattering rate in the Z = 0 plane of the Brillouin Zone.
* fermi_int_1 : Calculates transport integrals for constant scattering rate at a given Fermi energy.
* fermi_int_0 : Calculatest ransport integrals for a constant scattering rate at a given range of Fermi levels. Useful for visualization.
* fermi_int : Calculates transport integrals for the scattering rate based on electron-phonon interaction. This is the main code. 

makefile contains options for compilation. It is recommended to change gfortran to ifort (so fopenmp to openmp) if you have access to Intel compilers. To compile the code you want simply execute:

make name_of_the_code (e.g. make fermi_int)

## How to use the code
A band structure calculated on a large number of k-points is necessary. One can either calculate the band structure self consistently on this large grid, or do a self-consistent calculation on a small grid and then perform a non-self-consistent field calculation for the larger grid. Below are the basic steps:

<ol>
<li> A band structure calculation on a (large) grid, with the data saved in the prefix.a2Fsave files in QE (see the data folder). </li>
<li> At this point, one can calculate the transprt integrals using fermi_int_0 for a range of Fermi energies (see the example folder eg_fermi_int_0). </li>
<li> If one needs to compute the Fermi integrals for a given doping value, then Fermi level for that doping level needs to be computed using ef.f90 (see the example folder eg_ef). </li>
<li> At this point, one can either compute the Fermi integrals for constant scattering rate at a given Fermi level (fermi_int_1) or using the scattering rate resulting from electron-LO phonon interaction (fermi_int). </li>
<li> For calculation of electron-LO phonon interaction couping constants, one can use the simple Python script vLO.py. For an example, see example folder eg_vLO. </li>
</ol>

So that's basically it! 

## To do
* Extensive tests for a bunch of materials..
* Implementation of more generic scattering rates. Could be a better LO phonon one, in addition to deformation potential scattering ones..
