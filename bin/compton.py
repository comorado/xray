#!/usr/bin/env python2
"""
Compton Profile Convertor

Usage:
  {cmd} jpq2sqw <filename> <theta> <energy> [options]
  {cmd} sqw2jpq <filename> <theta> <energy> [options]
  {cmd} jpq2rhop <filename> [options]
  {cmd} xrts2sqw <filename> <energy> [options]
  {cmd} -h

Options:
  -h --help               Show this screen
  -w --omega=<omega>      Energy transfer "start,end[,step]" (e.g., "0,1001,0.5")
  -s --scattered          <energy> is scattered photon energy instead of
                          incident photon energy
  -n --num-electrons=N    Number of electrons to normalize J(pq) to.
  -a --atomic             Input file is atomic units
  -o --output=<filename>  Filename to save output to [default: stdout]
  -S --symmetrize         Symmetrize CP
  -c --column=<column>    Column in file containing data (only in xrts2sqw)
  -d --debug              Show debugging output
"""
from docopt import docopt
import sys, os
import numpy as np
from schema import Schema, Or, And, Use
from xray import compton, analysis, const

is_positive = lambda x: x > 0
open_output = lambda o: sys.stdout if o == 'stdout' else open(o, 'w')

def parse_omega(omega):
  if omega == None:
    w = np.linspace(0,1500,1501)
  else:
    tokens = omega.split(',')
    wstep = 1.0
    if len(tokens) == 3:
      wstep = float(tokens.pop())

    if len(tokens) < 2:
      raise ValueError("--omega must be of the form 'start,end[,step]'")

    w1,w2 = [float(t) for t in tokens]
    w = np.arange(w1,w2+wstep/2,wstep)

    return w



opts_schema = Schema({
  '<filename>': Schema(os.path.exists, error='Input file does not exist.'),
  '<theta>': Or(None, And(Use(float), lambda x: 0 <= x <= 180,
                 error='Invalid theta value')),
  '<energy>': Or(None, And(Use(float), is_positive,
                  error='Invalid energy')),
  '--omega': Use(parse_omega, error='Invalid omega'),
  '--num-electrons': Or(None,
                        And(Use(float), is_positive)
                       ), 
  '--scattered': bool,
  '--atomic': bool,
  '--symmetrize': bool,
  '--column': Or(None, And(Use(int))),
  '--output': Use(open_output),
  '--debug': bool,
  '--help': bool,
  'jpq2sqw': bool,
  'sqw2jpq': bool,
  'jpq2rhop': bool,
  'xrts2sqw': bool,
})


DEBUGGING = False
def debug(msg):
  if DEBUGGING:
    sys.stderr.write(str(msg))
    sys.stderr.write("\n")

def jpq2sqw(args):
    filename = args['<filename>']
    theta = args['<theta>']*np.pi/180
    energy = args['<energy>']
    scattered = args['--scattered']
    num_electrons = args['--num-electrons']
    atomic_units = args['--atomic']
    out = args['--output']

    cp = compton.ComptonProfile.from_file(filename, num_electrons=num_electrons, atomic_units=atomic_units)

    w = np.linspace(0,1500,1501)

    energy_key = 'E2' if scattered else 'E1'
    kwargs = {energy_key: energy}
    sqw = cp.to_sqw_theta(w, theta, **kwargs)

    header = ("# S(q,w) generated from J(p_q)\n"
              "# Source: {filename:s}\n"
              "# theta: {theta:.2f}\n"
              "# {which_energy:s}: {energy:.2f}\n"
              "#\n# w  S  q\n").format(
                  filename=os.path.abspath(filename),
                  theta=args['<theta>'],
                  which_energy=energy_key,
                  energy=energy)
    out.write(header)
    sqw.save(out)

def sqw2jpq(args):

  filename = args['<filename>']
  theta = args['<theta>']
  energy = args['<energy>']
  scattered = args['--scattered']
  num_electrons = args['--num-electrons']
  atomic_units = args['--atomic']
  out = args['--output']

  sqw = analysis.Curve.from_file(filename)

  if np.abs(np.log(energy) - np.log(sqw.x).mean()) < 1:
    sys.stderr.write("Input file appears to be in terms of scanned energy instead of energy transfer. Treating it as such.\n")
    if scattered: # scattered energy is fixed, sqw.x is incident
      sqw.x = (sqw.x - energy)[::-1]
    else:
      sqw.x = (energy - sqw.x)[::-1]
    sqw.y = sqw.y[::-1]

  if scattered:
    E1, E2 = energy + sqw.x, energy
  else:
    E1, E2 = energy, energy - sqw.x

  q = compton.momentum_transfer(E1, E2, theta*np.pi/180)
  cp = compton.sqw_to_jpq(sqw.y, q, sqw.x)

  if num_electrons:
    # XXX check that tails are present. otherwise warn...
    cp = cp.normalize_integral(num_electrons)

  if atomic_units:
    cp.x *= const.BOHR
    cp.y /= const.BOHR

  cp.save(out)


def jpq2rhop(args):
    filename = args['<filename>']
    num_electrons = args['--num-electrons']
    atomic_units = args['--atomic']
    symmetrize = args['--symmetrize']
    out = args['--output']

    cp = compton.ComptonProfile.from_file(filename, num_electrons=num_electrons, atomic_units=atomic_units)

    if symmetrize:
      cp = cp.symmetrize()

    rhop = cp.to_rhop()

    header = ("# rho(p) generated from J(p_q)\n"
              "# Source: {filename:s}\n"
              "#\n# p       rho(p)\n").format(
                  filename=os.path.abspath(filename),
                  )
    out.write(header)
    rhop.save(out)

def xrts2sqw(args):
    filename = args['<filename>']
    w1 = args['<energy>']
    column = args['--column'] or 1
    out = args['--output']


    data = analysis.Curve.from_file(filename, y=column)

    # switch to energy tranfer
    w2 = data.x
    w = (w1 - w2)[::-1]
    S = (data.y * (w1/w2)**2)[::-1] 

    sqw = analysis.Curve(w, S)

    header = ("# S(q,w) generated from XRTS output\n"
              "# Source: {filename:s}\n"
              "#\n# w       S\n").format(
                  filename=os.path.abspath(filename),
                  )
    out.write(header)
    sqw.save(out)

if __name__ == "__main__":
  usage = __doc__.format(cmd=os.path.basename(sys.argv[0]))
  args = docopt(usage, version="Compton 0.1)")

  if args['--debug']: DEBUGGING = True
  debug(args)

  try:
    args = opts_schema.validate(args)
  except Exception as e:
    #sys.stdout.write(usage)
    #sys.stdout.write("\n")
    sys.exit(e.code)
    #sys.stderr.write(str(e) + "\n")

  cmds = ['jpq2sqw', 'sqw2jpq', 'jpq2rhop', 'xrts2sqw']

  for cmd in cmds:
    if args[cmd]:
      locals()[cmd](args)
      exit()
