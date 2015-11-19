#!/usr/bin/env python

import sys
import numpy as np

def main():
  
  # Open input file
  f = open("input_LO","r")
  f_inp = f.readlines()
  f.close()

  # Assign input variables

  nmod = int(f_inp[0].split()[1]) # Number of split up LO/TO modes 
  a0 = float(f_inp[1].split()[1])   # lattice constant in au 
  eps = float(f_inp[2].split()[1])  # eps_infty

  # LO and TO modes
  LO = np.array([float(f_inp[3].split()[i]) for i in range(1,nmod+1,1)])
  TO = np.array([float(f_inp[4].split()[i]) for i in range(1,nmod+1,1)])

  # Physical constants
  RytoeV = 13.60569253 # Ryd = 13.60569253 eV
  tpi = 2.0*np.pi 
  au = 0.52917720859e-10 # Bohr radius
  eV = 1.602176487e-19
  hbar = 1.054571628e-34
  eps0 = 8.854187817e-12 # Vacuum permittivity
  Rytocm1 = 109737.57990 # Ryd = 109737.57990 cm-1

  # Prefactor 
  pre = eV / ( 2*tpi**2 * eps0 * au ) / RytoeV # This is in Ryd
  pre = pre / eps / a0
   
  Cnu = np.zeros(nmod)
  # Compute el-ph constants for each polar LO mode
  for nu in range(nmod):
    num = np.prod(np.ones(nmod) - TO**2 / LO[nu]**2)
    den = 1.0
    for j in range(nmod):
      if j != nu:
        den = den * (1.0 - LO[j]**2/LO[nu]**2)
    
    Cnu[nu] = pre * num/den * LO[nu] / Rytocm1


  # Save Cnu
  np.savetxt('Cnu.txt', Cnu)

  #print pre*eps*2*tpi**2


if __name__ == "__main__":
  main()
