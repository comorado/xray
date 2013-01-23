#!/usr/bin/env python
from argparse import ArgumentParser
from matplotlib import pyplot as pt
import numpy as np
import ldos
from xray import fermi

COLORS = [
  "#0000aa",
  "#ba0000",
  "#009600",
  "#ff4c00",
  "#6800ba",
  "#0094ff",
  "#e57400",
  "#d000ff",
  "#008662",
  "#ff8661"
]

parser = ArgumentParser(description="Plot FEFF LDOS file")

parser.add_argument('-t', '--total', help='plot total DOS (otherwise, individual lDOS)', action='store_true', default=False)
parser.add_argument('-f', '--filled', help='plot filled DOS', action='store_true', default=False)
parser.add_argument('-T', '--temperature', help='temperature in eV', type=float, default=0)
parser.add_argument('-m', '--mu', help='Chemical potential (or Fermi energy for T=0)', type=float, default=None)
parser.add_argument('ldos_file', help='File containing ldos data', default='ldos00.dat')

args = parser.parse_args()

l = ldos.LDOS(args.ldos_file)

x = l.energy

mu = args.mu
if mu is None:
  mu = l.mu

if args.total:
  y = np.sum(l.dos,0)
  pt.plot(x, y, color=COLORS[0])
  if args.filled:
    yf = y * fermi.f(x, args.temperature, mu)
    pt.plot(x,yf, color=COLORS[0], dashes=(2,1))
    pt.fill_between(x,yf, alpha=0.5, color=COLORS[0])
else:
  for i,y in enumerate(l.dos):
    pt.plot(x,y, color=COLORS[i])
    if args.filled:
      yf = y * fermi.f(x, args.temperature, mu)
      pt.plot(x,yf, color=COLORS[i], dashes=(2,1))
      pt.fill_between(x,yf, alpha=0.5, color=COLORS[i])


pt.show()
