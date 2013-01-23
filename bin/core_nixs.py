#!/usr/bin/env python2
# calculate core contribution of S(q,w) in IA
import sys
import argparse
import numpy as np
from xray import core_nixs, feff

# 1s, 2s, 2p1/2, 2p3/2, 3s, 3p1/2, 3p3/2, 3d...
ORBITAL = [0, 0, 1, 1, 0, 1, 1, 2, 2]

def setup_argparser():
  parser = argparse.ArgumentParser(description='Calculate core-shell contribution to NIXS in IA', formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('wf_file', help='FEFF atomic wavefunctions file', default='atomic_wf.dat')
  parser.add_argument('-q', type=float, default=10.0, help='momentum transfer in A^-1')
  parser.add_argument('-w', type=str, help='energy transfer range start:end[:step]', default='0:1500:5')
  parser.add_argument('-s', '--shells', type=str, help='comma separated list of shells to include (e.g. \'0,1\')', default=None)
  return parser

def parse_energy_transfer(wstr):
  pieces = [int(s) for s in wstr.split(':')]

  if len(pieces) == 3:
    w = np.arange(pieces[0], pieces[1]+pieces[2], pieces[2])
  elif len(pieces) == 2:
    w = np.arange(pieces[0], pieces[1]+pieces[2], 5)
  else:
    raise ArgumentError("Unknown energy transfer range: %s" % wstr)

  return w

def parse_shells(shellstr):
  shells = [int(s) for s in shellstr.split(',')]
  return shells

def main():
  parser = setup_argparser()
  args = parser.parse_args()

  wf = feff.load_atomic_wf(args.wf_file)
  w = parse_energy_transfer(args.w)
  q = args.q

  shellstr = args.shells
  if shellstr is None:
    shells = np.arange(len(wf))
    shellstr = '(All)'
  else:
    shells = parse_shells(args.shells)

  print "# w: ", args.w
  print "# q: ", args.q
  print "# shells: ", shellstr
  print "# Ls: ", ','.join([str(ORBITAL[i]) for i in shells])
  print "# w     S"

  #r2 = np.exp(-8.8 + 0.0125 * np.arange(1004))
  n = 4
  r2 = np.exp(-8.8 + (0.05/n) * np.arange(251*n))
  wf = [[r2, np.interp(r2, x[0], np.real(x[1])) + 1j * np.interp(r2, x[0], np.imag(x[1]))] for x in wf]
  S = [core_nixs.sqw(wf[i][0], wf[i][1],q,w, l=ORBITAL[i]) for i in shells]

  np.savetxt(sys.stdout, np.vstack([w]+S).T)
  

if __name__ == "__main__":
  main()

