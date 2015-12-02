#!/usr/bin/env python
# This is a collection of functions that compute the electron-phonon
# interaction using the long-range interaction of polarization (due to 
# longitudinal phonons) with the electronic charge density (treated 
# classically).
#
# Requires output from phonon calculations and matdyn.f90 in Phonon code
# of Quantum Espresso 
#
# Written by Burak Himmetoglu (burakhmmtgl@gmail.com)

import sys
import numpy as np
from numpy.linalg import inv

def read_input():
  """ Read input from path.out, matdyn.modes, prefix.ph.out, and input.dat 
      path.out : The q-points where el-ph and polarization is computed.
                 It has to be given in Cartesian coordinates.
      matdyn.modes : Contains phonon frequencies and eigenvectors computed
                     at the points given at path.out. 
      prefix.ph.out : ph.x output at the Gamma point. Must contain 
                      calculation of E-fields and Z*.
      input.dat : Contains some info on the crystal structure. """
      
  argv = sys.argv

  # Read file names from sd input
  f_dy  = argv[1]  # matdyn.modes
  f_pat  = argv[2] # path.out (should be in crystal coords)
  f_ph  = argv[3]  # ph.x output (Gamma point)

  # Read input card
  f_inp = open("input.dat",'r')
  l1 = f_inp.readline()
  l2 = f_inp.readline()
  l3 = f_inp.readline().split()
  #l4 = f_inp.readline().split()
  f_inp.close()

  # Open files

  f = open(f_dy,'r')   # matdyn.modes 
  f_dyn = f.readlines()
  f.close()

  f = open(f_pat,'r')   # path.out
  f_path = f.readlines()
  f.close()

  f = open(f_ph,'r')   # ph.x output
  f_zs = f.readlines()
  f.close()

  # Assign values to a0, nat, M, nqp
  a0, vol  = float(l1.split()[0]), float(l1.split()[1])
  nat = int(l2) 
  mass = np.zeros(nat)
  for iat in range(nat):
    mass[iat] = float(l3[iat])

  # Read Z* tensor from f_zs
  i = 0
  iz = 0
  zstart = []
  for line in f_zs:
    if "(d P / du)" in line:
      iz = i + 3
    if "Px" in line:
      zstart.append(i)

    i += 1

  # Read the dielectric tensor from f_zs
  i = 0
  ie = 0
  for line in f_zs:
    if "Dielectric constant in cartesian axis" in line:
      ie = i + 2
      break

    i += 1

  # Assign Z* values
  zs = np.zeros((nat,3,3)) # initialize Z*

  for iat in range(nat):
    for ic in range(3):
      ztext = f_zs[zstart[iat]+ic][19:56].split()
      for jc in range(3):
        zs[iat][ic][jc] = float(ztext[jc])

  # Assing the dielectric tensor
  eps = np.zeros((3,3))

  for ic in range(3):
    epstext = f_zs[ie+ic][16:66].split()
    for jc in range(3):
      eps[ic][jc] = float(epstext[jc])

  # We need the inverse of epsilon
  epsinv = inv(eps)

  # Number of modes and q-points
  nmodes = 3 * nat
  nqpt = int(f_path[0].split()[0])

  # Read the q-points
  q = np.zeros((nqpt,4)) # 4th dimension is lenght for q-points on a line, weights for q-points on a grid 
  for iq in range(1,nqpt+1):
    q[iq-1,] = np.array([float(f_path[iq].split()[0]),float(f_path[iq].split()[1]), \
               float(f_path[iq].split()[2]),float(f_path[iq].split()[3])])

  # Read the eigenvalues(om) and eigenvectors(eig) 
  # Initiate 
  om = np.zeros((nmodes,nqpt))
  eig = np.zeros((nmodes,nqpt,nat,3), dtype=complex)   

  # Get the starting lines for each q-pt
  i = 0
  i_q = []
  for line in f_dyn:
    if "q =" in line:
      i_q.append(i+2)
    i += 1

  #Assign values to om and eig
  for iq in range(nqpt):
    for imod in range(nmodes):
      omtext = f_dyn[i_q[iq]+imod*(nat+1)][43:55]
      om[imod][iq] = float(omtext)
      for iat in range(nat):
        etext = f_dyn[i_q[iq]+imod*(nat+1)+iat+1][2:72].split()
        for ic in range(3):
          eig.real[imod][iq][iat][ic]=float(etext[2*ic])*np.sqrt(mass[iat])
          eig.imag[imod][iq][iat][ic]=float(etext[2*ic+1])*np.sqrt(mass[iat])

      #Normalize the eigenvectors
      t1 = eig[imod,iq,:,:]
      t_nu = np.sum(np.sum(np.conjugate(t1)*t1,axis=0))
      eig[imod,iq,:,:] = eig[imod,iq,:,:]/np.sqrt(np.abs(t_nu))

    # Check normalization
    delta = np.zeros((nmodes,nmodes), dtype=complex)
    for iat in range(nat):
      for ic in range(3):
        t2 = eig[:,iq,iat,ic]
        delta += np.outer(np.conjugate(t2),t2)

    unit = np.diag(np.diag(np.ones((nmodes,nmodes)))) # Unit vector
    test = np.abs( (delta-unit) )
    if ( np.max(test) > 1e-3):
       print "Non-orthonormal eigenvector at iq=", q[iq,:]

  return om, eig, q, zs, epsinv, mass, a0, vol, nmodes, nqpt, nat

def ord_modes( om, eig, nmodes, nqpt, nat ):
  """ Order the modes according to overlaps of eigenvectors 
      Somehow, does not work properly, so it is not called
      from the main(). Hopefully, there will be a fix.  """

  # Initiate ordered eigenvalues and eigenvectors
  om_new = np.zeros((nmodes,nqpt))
  eig_new = np.zeros((nmodes,nqpt,nat,3), dtype=complex)   

  om_new[:,0] = om[:,0] ; eig_new[:,0,:,:] = eig[:,0,:,:]

  overlap = np.zeros(nmodes,dtype=complex)
  t1 = np.array((nat,3),dtype=complex)
  t2 = np.array((nat,3),dtype=complex)

  # Loop over q-points 
  nn = 1
  for iq in range(0,nqpt-nn):
    for imu in range(nmodes):
      t1 = eig[imu,iq,:,:]
      # Calculate overlaps
      overlap[:] = 0.0
      for inu in range(nmodes):
        t2 = eig[inu,iq+nn,:,:]
        overlap[inu] = np.sum(np.sum(np.conjugate(t2)*t1,axis=0))

      # When overlap has maximum absolute value: inu(max) --> imu
      ind = np.argmax(np.abs(overlap))
      if ( np.abs(om[imu,iq]-om[ind,iq+nn]) < 20.0):
        eig_new[imu,iq+nn,:,:] = eig[ind,iq+nn,:,:]
        om_new[imu,iq+nn] = om[ind,iq+nn] 
      else:
        eig_new[imu,iq+nn,:,:] = eig[imu,iq+nn,:,:]
        om_new[imu,iq+nn] = om[imu,iq+nn] 

  return om_new, eig_new

def polarization(q,zs,eig,nmodes,nqpt,nat):
  """ This function computes the polarization generated by phonon modes.
      It illustrates that at generic k-points, the separation between
      longitudinal and transverse modes gets lost """

  # Induced longitudinal polarizations 
  pol = np.zeros((nmodes,nqpt), dtype=complex)
  l_mod = np.zeros((nmodes,nqpt), dtype=complex)

  # The l_mod is nonzero for longitudinal modes
  for imod in range(nmodes):
    for ic in range(3):
      for iat in range(nat):
        l_mod[imod,:] += q[:,ic]*eig[imod,:,iat,ic]

  l_mask = np.abs(l_mod) > 1e-6 # Region in the imod,iq space where l_mod is nonzero
  
  # Now compute the polarization using Z*
  for imod in range(nmodes):
    for ic in range(3):
      for iat in range(nat):
        for jc in range(3):
          pol[imod,:] += q[:,ic]*zs[iat][ic][jc]*eig[imod,:,iat,jc]

  # Use l_mask to mask out the longitudinal components
  #pol = pol #* l_mask

  return pol  

def el_ph(om,eig,q,zs,mass,epsinv,nmodes,nqpt,nat):
  """ Electron-phonon interaction matrix, based on a 
      continuum interaction of longitudinal polarization 
      field with electrons. """

  # Initiate
  g = np.zeros((nqpt,nmodes),dtype=complex)

  # 1/q^2
  q2 = np.sum(q[:,1:4]**2,axis=1)
  inv_q2 = 1.0 / (q2+1e-10)

  for imod in range(nmodes):
    for ia in range(3):
      for ib in range(3):
        for ic in range(3):
          for iat in range(nat):
            g[:,imod] += inv_q2[:]*q[:,ia]*epsinv[ia,ib]*zs[iat,ib,ic]*eig[imod,:,iat,ic] \
                            / np.sqrt(2.0*mass[iat]*np.abs(om[imod,:])+1e-10)

  return g

# Main prog  

def main():
  """ Usage: mode_sep.py matdyn.modes path.out prefix.ph.out """  

  # Get omega, eigenvectors, q-pts (and their dimensions)
  om, eig, q, zs, epsinv, mass, a0, vol, nmodes, nqpt, nat = read_input()  

  # Get ordered omega and eigenvectors (some issues..)
  #om_new, eig_new = ord_modes(om,eig,nmodes,nqpt,nat)

  # Get polarization (not needed, just for curiosity)
  #pol = polarization(q,zs,eig,nmodes,nqpt,nat)

  # Get el_ph
  g = el_ph(om,eig,q,zs,mass,epsinv,nmodes,nqpt,nat)

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

  pre = eV**2 * hbar / (eps0 * au**2 * tpi * np.sqrt(mau * cm_1toJ)) # Prefactor in J (i.e. SI)
  pre = pre / eV / RytoeV * (a0 / vol) # prefactor in Ryd

  g *= pre

  # # # SAVING # # #
  ##
  # Phonon frequencies
  om_plt = np.transpose(om)
  np.savetxt("freq.out",om_plt,fmt='%1.4e')

  # El-ph
  g_plt = g.real**2 + g.imag**2 
  np.savetxt("g2.out",g_plt,fmt='%1.4e')

  # Saving both real and complex parts
  #g_plt = np.zeros((nqpt,2*nmodes))
  #g_plt[:,0:nmodes] = g.real[:,0:nmodes]
  #g_plt[:,nmodes:2*nmodes] = g.imag[:,0:nmodes]
  #np.savetxt("g.out",g_plt,fmt='%1.4e')


if __name__ == "__main__":
  main()
