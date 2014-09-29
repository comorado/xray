#!/usr/bin/env python
import os, sys
import subprocess
import lebedev.grids as grids

COMPTON_PATH = "/home/clevac/old/feff9/bin/MPI/compton"

skel=""" run compton module?
           1
pqmax, npq
   5.000000            1000
ns, nphi, nz, nzp
  32  32  32 192
smax, phimax, zmax, zpmax
      {smax}      6.28319      {zmax}     {zpmax}
jpq? rhozzp? force_recalc_jzzp?
 T F T
window_type (0=Step, 1=Hann), window_cutoff
           1  0.0000000E+00
temperature (in eV)
      0.00000
set_chemical_potential? chemical_potential(eV)
 F  0.0000000E+00
rho_xy? rho_yz? rho_xz? rho_vol? rho_line?
 F F F F F
qhat_x qhat_y qhat_z
  {qhat} 
"""

def write_compton_input(qhat, scale=1.0):
  smax = 2.3 * scale
  zmax = 2.3 * scale
  zpmax = 25 * scale

  qhatstr = ' '.join(['%.10f'%f for f in qhat])
  inpstr = skel.format(smax=smax, zmax=zmax, zpmax=zpmax, qhat=qhatstr)
  with open("compton.inp", "w") as f:
    f.write(inpstr)

def run_compton(nprocs, nodefile, tag):

#  ret = subprocess.call(["echo","-n",">", "compton.{tag}.dat".format(tag=tag)]) 
#  if ret != 0:
#    print("Error copying. Aborting.")
#    return False


  cmd = ["mpirun", "-np", str(nprocs), "-bynode", "-machinefile", nodefile, COMPTON_PATH]
  ret = subprocess.call(cmd)

  if ret != 0:
    print("Error. Aborting.")
    return False

  ret = subprocess.call(["cp", "compton.dat", "compton.{tag}.dat".format(tag=tag)]) 

  if ret != 0:
    print("Error copying. Aborting.")
    return False

  return True
def copy(tag, qhat, npts):
  res = subprocess.call(["ls","-l","compton.{tag}.dat".format(tag=tag)])
  if (res == 0):
       return True
  for j in grids.grid_sizes:
    if j < npts:
      grid2 = grids.grids[j]

      for n,(x2,y2,z2,w2) in enumerate(grid2):
        if qhat == [x2,y2,z2]:
          tag2 = "L{j}.{n}".format(j=j, n=n)
          res = subprocess.call(["ls","-l","compton.{tag2}.dat".format(tag2=tag2)])
          num_lines = file_len("compton.{tag2}.dat".format(tag2=tag2))

          #Check if the file exist and number of lines in file is greater then 1000
          if (res == 0 and num_lines > 1000):
            subprocess.call(["cp", "compton.{tag2}.dat".format(tag2=tag2), "compton.{tag}.dat".format(tag=tag)]) 
            with open("compton.{tag}.dat".format(tag=tag), "a") as myfile:
              myfile.write("# Original file compton.{tag2}.dat\n".format(tag2=tag2))
            return True
          elif (res == 0):
            return True

  return False

#checks number of lines in the filename provided
def file_len(fname):
    i = 0
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
 
nprocs = sys.argv[1]
nodefile = sys.argv[2]
totpts = 0
for npts in grids.grid_sizes:
  if totpts > 2000: break
  grid = grids.grids[npts]
  for i,(x,y,z,w) in enumerate(grid):
      if (npts > 4 and npts < 50):
        qhat = [x,y,z]
        tag = "L{npts}.{i}".format(npts=npts, i=i)
        if (copy(tag, qhat, npts) == False):
          write_compton_input(qhat)
          subprocess.call(["cp", "compton.inp", "compton.{tag}.inp".format(tag=tag)]) 
          with open("compton.{tag}.dat".format(tag=tag), "w") as f:
            f.write("")

          ret = run_compton(nprocs, nodefile, tag)
          if not ret:
            exit(1)
          totpts += 1

