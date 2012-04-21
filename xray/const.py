from math import pi

HBAR  = 6.58211899e-16 # in eV s
HBARC = 1973.269631 # eV * angstroms
HC    = 2 * pi * HBARC 
MC2   = 5.10998903e5 # electron mass in eV
C     = 2.99792458e18 # speed of light (angstroms / s)
Me    = MC2 / C**2
#BOHR  = 0.52917720859 # Bohr radius of H in angstroms (codata?)
BOHR  = 0.52917721092 # Bohr radius of H in angstroms (codata 1.27.2012)
HARTREE = 27.2116 # Hartree in eV
R0    = 2.81794029e-5 # Thomson scattering length in angstroms

AVOGADRO = 6.0221415e23

BOLTZMANN = 8.6173324e-5 # eV/K

DEGREE = pi/180 # in radians


def rad(deg):
  """
  Convert from degrees to radians
  """
  return deg * DEGREE

def deg(rad):
  """
  Convert from radians to degrees
  """
  return rad / DEGREE
