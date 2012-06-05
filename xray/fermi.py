"""
Fermi Gas calculations 
"""

import numpy as np
from scipy.integrate import quad
from scipy.integrate.quadpack import Inf
from .special import sph_jn
import const

pi = np.pi
T_MIN = 1e-8

SPIN_DEGENERACY = 2

# effective Fermi gas parameters
elt_info = {
    'Li': {
      'N': 1,
      'a': 3.51,
      'V': 21.62,
      'kF': 1.11,
    },
    'Be': {
      'N': 2, # number of valence electrons per atom
      'a': 2.291, # lattice param in A
      'c': 3.581, # lattice param in A
      'V': 8.1374, # atomic volume in A^3
      'kF': 1.938, # fermi momentum in A^-1
    },
    'Na': {
      'N': 1,
      'a': 4.2906,
      'V': 39.4934,
      'kF': 0.9084,
    },
    'Mg': {
      'N': 2,
      'a': 3.2094,
      'c': 5.2108,
      'V': 23.241,
      'kF': 1.366,
    },
}

CONVERSION = {
    'N': 1,
    'a': 1/const.BOHR,
    'b': 1/const.BOHR,
    'c': 1/const.BOHR,
    'V': 1/const.BOHR**3,
    'kF': const.BOHR,
}

# convert to atomic units
elt_info_au = {k:{kk:CONVERSION[kk]*elt_info[k][kk] for kk in elt_info[k]} for k in elt_info}


def f(E, T, mu):
  """
  Fermi distribution

  Parameters:
    E: energy
    T: temperature 
    mu: chemical potential

  All values should be given in energy units (that is, `T` is really k_B * T).

  If T < 1e-8, it is set to 0.

  """

  # for small values of T, use T=0 distribution
  if (T < T_MIN):
    return (E <= mu) * 1.0

  return 1.0/(np.exp((E-mu)/T)+1)

def fermi_momentum(density):
  """
  Fermi momentum for a given electronic density
  """
  return (6 * pi**2 * density / SPIN_DEGENERACY)**(1/3.)

def fermi_energy(density, m=1.0):
  """
  Fermi energy for a given electronic density

  Parameters:
    density: electronic density (in electrons / bohr^-3)
    m: effective mass (in units of electron mass)

  Returns:
    Fermi energy in Hartree
  """
  return fermi_momentum(density)**2 / (2 * m)

def rhop(p, pF=1.0, T=0, mu=None):
  """
  Occupied momentum density for a Fermi gas

  Parameters:
    p: momenta to evaluate density at
    pF: fermi momentum
  """

  rho = SPIN_DEGENERACY / (2*pi)**3

  if mu is None:
    mu = pF**2/2

  if T < T_MIN:
    return rho * (pF >= p)
  else:
    return rho * f(p**2/2, T, mu)

def rhoe(energy, V=1.0, m=1.0):
  """
  DoS for a Fermi Gas

  Parameters:
    energy: energy values to evaluate rhoE at
    V: unit-cell volume (defaults to 1.0)
    m: effective mass (in units of electron mass)
  """

  return m**1.5 * SPIN_DEGENERACY / (np.sqrt(2)*pi**2) * np.sqrt(energy) * V

def rholk(k, l, V):
  """
  Local l-projected momentum density for Fermi Gas

  Parameters:
    k: momentum values to calculate density at
    l: angular momentum
    V: Unit cell volume 
  """

  R = (3/(4*pi)*V)**(1/3.)

  single_k = np.isscalar(k)
  k = np.atleast_1d(k)

  ret = np.array([quad(lambda u: u**2 * sph_jn(l,u)**2, 0, ki*R)[0] for ki in k])
  i = k>0
  ret[i] *= (4/pi)*(2*l+1)/k[i]**3
  ret[k==0] = V/pi**2 if l == 0 else 0

  if single_k: ret = ret[0]

  return ret

def rhole(energy, l, V):
  """
  Local l-projected DoS for Fermi Gas

  Parameters:
    energy: energy values to calculate rhole at
    l: angular momentum
    V: Unit cell volume 
  """

  k = np.sqrt(2*energy)
  return rholk(k, l, V) * k

def mu_T(rhoe_func, n, T, args=(), disp=False):
  """
  Chemical potental at temperature T

  Parameters:
    rhoe_func: density of states function
      must take energy as its first parameter
    n: electronic density (electrons / bohr^3)
    T: temperature (in energy units)
    args: any additional arguments to pass on to rhoe_func
  """

  def occupied(energy, mu):
    return rhoe_func(energy, *args) * f(energy, T, mu)

  def nmu (mu):
    return quad(occupied, 0, Inf, mu)[0]

  from scipy.optimize import fmin

  ret = fmin(lambda mu: np.abs(n-nmu(mu)), fermi_energy(n), full_output=True, disp=disp)

  mu = ret[0][0]

  # XXX handle non-convergence
  #warnflag = ret[4]

  return mu

