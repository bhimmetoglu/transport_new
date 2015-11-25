#!/usr/bin/env python

import sys
import numpy as np

def main():
  """ This is a simple script which calculates the LO-phonon coupling constants
      based on Toyozowa's Hamiltonian. """
  
  # Open input file
  f = open("input_LO","r")
  f_inp = f.readlines()
  f.close()

  # Assign input variables

  nmod = int(f_inp[0].split()[1]) # Number of split up LO/TO modes 
  a0 = float(f_inp[1].split()[1])   # lattice constant in au 
  vol = float(f_inp[2].split()[1])   # vol in au^3
  eps = float(f_inp[3].split()[1])  # eps_infty

  # LO and TO modes
  LO = np.array([float(f_inp[4].split()[i]) for i in range(1,nmod+1,1)])
  TO = np.array([float(f_inp[5].split()[i]) for i in range(1,nmod+1,1)])

  # Constants to compute the prefactor
  RytoeV = 13.60569253  # Ryd to eV conversion
  tpi = 2.0*np.pi # 2pi
  au = 0.52917720859e-10 # Bohr radius
  me = 9.10938215e-31 # electron mass
  eV = 1.602176487e-19
  hbar = 1.054571628e-34
  mau = 1.66053886e-27 # amu
  Rytocm_1 = 109737.57990 # Ryd to cm-1 conversion
  cm_1toJ = (RytoeV * eV)/Rytocm_1 # cm-1 to J conversion
  eps0 = 8.854187817e-12 # Vac permittivity

  # Prefactor
  pre = eV**2 * cm_1toJ / (2 * eps0 * au * tpi**2) # This is in J**2
  pre = pre / (eV * RytoeV)**2 # This is in Ryd**2
  pre = pre * a0**2 / vol / eps # Add lattice parameter dependencies, and epsilon_infty

  Cnu = np.zeros(nmod)
  # Compute el-ph constants for each polar LO mode
  for nu in range(nmod):
    num = np.prod(np.ones(nmod) - TO**2 / LO[nu]**2)
    den = 1.0
    for j in range(nmod):
      if j != nu:
        den = den * (1.0 - LO[j]**2/LO[nu]**2)
    
    Cnu[nu] = pre * num/den * LO[nu] 

  # Save Cnu
  np.savetxt('Cnu.txt', Cnu)

  #print pre*eps*2*tpi**2


if __name__ == "__main__":
  main()
