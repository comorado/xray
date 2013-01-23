import numpy as np
from .const import HBARC, BOHR, HARTREE

def sph_fourier(r,f,p,l=0):
  """
  Fourier transform radial function f_l(r) to F_l(p)

  TODO: modify this to determine p based on extent of f(r)
  """
  from utils import sph_jn
  jn = sph_jn(l,p[:,np.newaxis]*r[np.newaxis,:])
  return np.sqrt(2/np.pi)*np.trapz(r*f[np.newaxis,:]*jn, r)

def momentum_density(u):
  """
  Radial momentum distribution

  FIXME: this assumes 1s wavefunction

  Parameters:
    u: radial momentum-space wavefunction
  """
  return abs(u)**2 / (4*np.pi)

def pwffa(p, u, q, w, B):
  """
  Calculate S(q,w) in plane-wave form-factor approximation
  (Schumacher & Smend 1975)

  u(p): momentum space wavefunction in a.u.
  q, w: momentum-, energy transfer in a.u (Hartree)
  B: binding energy in Hartree
  """
  rhop = momentum_density(u)
  pmin = abs(np.sqrt(2*(w-B))-q);
  pmax = np.sqrt(2*(w-B))+q
  pmin[w<B] = 0
  pmax[w<B] = 0
  indices = ((p1<=p)&(p<=p2) for p1,p2 in zip(pmin,pmax))
  S = np.array([np.trapz(p[i]*rhop[i], p[i]) for i in indices])

  return (2*np.pi/q) * S

def impulse(p, u, q, w):
  """
  Calculate S(q,w) in impulse approximation

  u(p): momentum space wavefunction in a.u.
  q, w: momentum-, energy transfer in a.u (Hartree)
  """

  rhop = momentum_density(u)
  pmin = abs(w/q-q/2)
  indices = ((p1<=p) for p1 in pmin)
  S = np.array([np.trapz(p[i]*rhop[i], p[i]) for i in indices])

  return (2*np.pi/q) * S

def hydrogenic1s(r, Z):
  """ Impulse approximation expression for hydrogenic 1s shell """
  return 2 * r * Z**1.5 * np.exp(-Z*r)

def momentum_transfer(E1, E2, theta):
  """
  E1,E2: incident and emitted energy in eV
  theta: scattering angle in radians
  returns q in inv. angstroms
  """
  return np.sqrt(E1**2 + E2**2 - 2 * E1 * E2 * np.cos(theta))/HBARC

def S_hydrogenic_1s(w1, w2, theta, Z, B=0, type='IA'):
  """
  w1,w2 - incident, scattered energies in eV
  theta - scattering angle in degrees
  Z - effective nuclear charge

  returns: S(q,w) for 1 electron in inv. eV
  """
  r = np.linspace(0,20,1001)
  p = np.linspace(0,20,1001)
  f = hydrogenic1s(r,Z)
  u = sph_fourier(r,f, p)

  w = (w1 - w2) / HARTREE
  q = momentum_transfer(w1,w2,theta*np.pi/180) * BOHR

  if type == 'IA':
    S = impulse(p, u, q, w)
  elif type == 'PWFFA':
    S = pwffa(p, u, q, w, B/HARTREE)

  return S / HARTREE

def sthetaw(r, psi, w1, w2, theta, B=0, type='IA'):
  """
  Parameters:
    r: radial coordinate (a.u.)
    psi: radial wave function (a.u.)
    w1: incident enery (eV)
    w2: incident enery (eV)
    theta: scattering angle (degrees)
    B: binding energy (eV) (only used for PWFFA)
    type: 'IA' for impulse approximation or 'PWFFA' for plane-wave form-factor
  """

  w = (w1 - w2)
  q = momentum_transfer(w1,w2,theta*np.pi/180)
  return sqw(r,psi,q,w,B,type)

def sqw(r,psi, q, w, B=0, type='IA', l=0):
  """
  Parameters:
    r: radial coordinate (a.u.)
    psi: radial wave function (a.u.)
    q: momentum transfer (inv. A)
    w: enery transfer (eV)
    B: binding energy (eV) (only used for PWFFA)
    type: 'IA' for impulse approximation or 'PWFFA' for plane-wave form-factor
    l: orbital quantum number for shell
  """
  p = np.linspace(0,40,2001) # XXX be careful about this upper bound
  u = sph_fourier(r,psi, p, l)

  if type == 'IA':
    S = impulse(p, u, q*BOHR, w/HARTREE)
  elif type == 'PWFFA':
    S = pwffa(p, u, q*BOHR, w/HARTREE, B/HARTREE)

  return S / HARTREE


def eisenberger1s(w, q, Z, B=None):
  """
  Parameters:
    w: energy transfer in Hartree
    q: momentum transfer in inv. Bohr
    Z: effective nuclear charge
    B: binding energy (leave as None for correct physics!)

  Returns:
    S(q,w) for a hydrogenic system in inv. Hartree

  This combines equations (6), (8), (9) and (10) of Eisenberger and Plazman PRA 2, 2 (1970) for a single hydrogenic electron.

  (Note: bound states are not included in the sum over final states. should they be?)
  """

  if B is None:
    B = Z**2 / 2.0
  p = np.sqrt(2*(w-B))

  a = 1.0 / Z
  pa = p / Z
  qa = q / Z

  S = (np.pi**2 * 8**3 * a**2  * (1 - np.exp(-2*np.pi/pa))**(-1) *
       np.exp(-2/pa * np.arctan2(2 * pa , (1 + qa**2 - pa**2))) *
       ((q/Z)**4 + (q/Z)**2 * (1 + (p/Z)**2) / 3.0) *
       (((q/Z)**2 + 1 - (p/Z)**2)**2 + 4 * (p/Z)**2)**(-3) * 
       4*np.pi / 8 / np.pi**3  # factors from integration
       ) 
  S[np.isnan(S)] = 0
  return S
