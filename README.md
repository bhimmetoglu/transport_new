# Boltzmann transport code
## Main Description
This package contains codes to calculate transport properties using Boltzmann Transport theory. 
There are three distinct features of the code:

1. Calculation of scattering rates based on electron-phonon(LO only) interactions. A simple analytical formula is used.
2. Adaptive smearing width choice for the evaluation of energy conserving Dirac delta functions. See documentaion.
3. Adaptive k-grid around the Fermi level to reduce integration time. See documentation.

The code is parallelized using OpenMP. This means you can only use a single node of a cluster for parallelization. Recommended setting is to use nthreads=number of cores in the node. However, you can use nthreads exceeding the nunmber of codes a little and may get some performance increase. Perform some tests before doing so (In Stampede, I have seen performance increase up to 48 threads.)

## Description of each routine in /src
* ef.f90 : Computes Fermi levels for a given doping level in cm-3. 
* dos.f90 : Computes density of states using the adaptive smearing scheme. Recall that DOS is given by
  DOS(e) = sum_{n, k}\, \delta( e - e_{nk} )
* tauef.f90 : This code computes the scattering rate at the Fermi level:
  tau^{-1}(E_F) = \sum_{n,k}\, tau^{-1}_{nk}\, delta( e_{nk} - E_F)
* tauBZ_plt : Visualization code for plotting the scattering rate in the Z = 0 plane of the Brillouin Zone.
* fermi_int_0 : Calculates transport integrals for constant scattering rate at a given Fermi energy.
* fermi_int_1 : Calculatest ransport integrals for a constant scattering rate at a given range of Fermi levels. Useful for visualization.
* fermi_int : Calculates transport integrals for the scattering rate based on electron-phonon interaction. This is the main code. 

makefile contains options for compilation. It is recommended to change gfortran to ifort (so fopenmp to openmp) if you have access to Intel compilers. To compile the code you want simply execute:

make name_of_the_code (e.g. make fermi_int)

## How to use the code
A band structure calculated on a large number of k-points is necessary. One can either calculate the band structure self consistently on this large grid, or do a self-consistent calculation on a small grid and then perform a non-self-consistent field calculation for the larger grid. Below are the basic steps:

<ol>
<li> Band structure calculation on a (large) grid, and the data is saved in the prefix.a2Fsave files in QE (see the /data folder </li>
<li> At this point, one can calculate the transprt integrals using fermi_int_0 for a range of Fermi energies. </li>
</ol>

## To do
* Extensive tests for a bunch of materials..
* Implementation of more generic scattering rates. Could be a better LO phonon one, in addition to deformation potential scattering ones..
