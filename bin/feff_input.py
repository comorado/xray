#!/usr/bin/env python2
"""
FEFF Input File Modifier

Usage:
    {cmd} scale <filename> <scale> [--scf] [--fms] [--ldos] [options]
    {cmd} rotate <filename> <qx> <qy> <qz> [options]
    {cmd} (-h | --help)

Generic Options:
    <filename>              FEFF input file
    -h --help               Show this screen
    -o --output=<output>    Output filename

Command Specific Options:

  scale:
    <scale>                 Scale factor
    -s --scf                Scale SCF radius
    -F --fms                Scale FMS radius
    -l --ldos               Scale LDOS parameters

  rotate:
    <qx>, <qy>, <qz>        Vector for new z-axis to point along

"""

from docopt import docopt
from xray import feff
import sys, os

if __name__ == "__main__":
  arg0 = os.path.basename(sys.argv[0])
  args = docopt(__doc__.format(cmd=arg0), version="FEFF Input File Modifier 0.1)")

  #print(args)

  if args['scale']:
    f = feff.InputFile(args['<filename>'])
    scale = float(args['<scale>'])

    f.scale_atoms_uniformly(scale)

    def modify_line(line):
      if ((args['--scf'] and line.find('SCF') >= 0) or
          (args['--fms'] and line.find('FMS') >= 0)):
        components = line.strip().split()
        card = components[0]
        val = float(components[1])*scale
        rest = ' '.join(components[2:])
        return ' {card:s}       {val:.3f} {rest:s}\n'.format(**locals())
      elif args['--ldos'] and line.find('LDOS') >= 0:
        components = line.strip().split()
        card = components[0]
        start = float(components[1]) / scale**2
        end = float(components[2]) / scale**2
        rest = ' '.join(components[3:])
        return ' {card:s}       {start:.3f} {end:.3f} {rest:s}\n'.format(**locals())
      else:
        return line

    f.pre_atoms = [modify_line(l) for l in f.pre_atoms]

    filename = args['--output'] or sys.stdout
    f.save(filename)

  elif args['rotate']:
    infile = args['<filename>']
    outfile = args['--output'] or sys.stdout

    qx = float(args['<qx>'])
    qy = float(args['<qy>'])
    qz = float(args['<qz>'])

    f = feff.InputFile(infile)
    f.rotate_atoms_to([qx, qy, qz])
    f.save(outfile)

