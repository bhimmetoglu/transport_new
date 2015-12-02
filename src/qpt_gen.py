#!/usr/bin/env python

import sys
import numpy as np

def main():
  """ This script prints out the k-point grid in Cartesian and Crystal
      coordinates from an scf or nscf run of pwscf. 
      ---
      Usage: qpt_gen.py prefix.(n)scf.out  """
  
  f_scf = sys.argv[1] # Read scf.out file
  f = open(f_scf, 'r')
  f_pw = f.readlines()
  f.close()

  i = 0
  j = 0
  # Assign variables from f_pw
  for line in f_pw:
    if "points=" in line:
      nkp = int(line.split()[4])
    elif " cryst. coord." in line:
      ikpt = i + 1 # Store the line number where k-point coordinates start
    elif " cart. coord." in line:
      ikpt_c = j + 1 # Store the line number where k-point coordinates start

    i += 1 # increase i by 1 at each line in f_pw
    j += 1 

  # Store kpoints and weights
  kpoint = []
  kpoint_c = []
  for ik in range(nkp):
      ktext = f_pw[ikpt + ik][20:56].split()
      wktext = f_pw[ikpt + ik][66:75]
      kpoint.append([float(ktext[0]), float(ktext[1]), float(ktext[2]), float(wktext)])
      #
      ktext = f_pw[ikpt_c + ik][20:56].split()
      wktext = f_pw[ikpt_c + ik][66:75]
      kpoint_c.append([float(ktext[0]), float(ktext[1]), float(ktext[2]), float(wktext)])

  # Store in a numpy array and save
  kp = np.array(kpoint)
  kp_c = np.array(kpoint_c)
  np.savetxt("path_cryst.out", kp, fmt='%1.6f')
  np.savetxt("path_cart.out", kp_c, fmt='%1.6f')

if __name__ == "__main__":
  main()

