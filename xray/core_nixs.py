import numpy as np
from .const import HBARC, BOHR, HARTREE
from .special import sph_jn
from .utils import trapezoid
from scipy.misc import factorial

def sph_fourier(r,f,p,l=0):
  """
  Fourier transform radial function f_l(r) to F_l(p)

  TODO: modify this to determine p based on extent of f(r)
  """
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
  #indices = ((p1<=p)&(p<=p2) for p1,p2 in zip(pmin,pmax))
  #S = np.array([np.trapz(p[i]*rhop[i], p[i]) for i in indices])
  S = np.array([ trapezoid(p, p*rhop, p1,p2) for p1,p2 in zip(pmin,pmax)])

  return (2*np.pi/q) * S

def impulse(p, u, q, w):
  """
  Calculate S(q,w) in impulse approximation

  u(p): momentum space wavefunction in a.u.
  q, w: momentum-, energy transfer in a.u (Hartree)
  """

  rhop = momentum_density(u)
  pmin = abs(w/q-q/2)
  #indices = ((p1<=p) for p1 in pmin)
  #S = np.array([np.trapz(p[i]*rhop[i], p[i]) for i in indices])
  S = np.array([ trapezoid(p, p*rhop, pm) for pm in pmin])

  return (2*np.pi/q) * S

def hydrogenic1s(r, Z):
  """ Impulse approximation expression for hydrogenic 1s shell """
  return 2 * r * Z**1.5 * np.exp(-Z*r)

def hydrogenic_radial(r, Z, n=1,l=0):
  """Radial wavefunction (times r) for a hydrogenic system"""
  if n == 1 and l == 0:
    return 2 * r * Z**1.5 * np.exp(-Z*r)
  elif n == 2 and l == 0:
    return (np.sqrt(2)/4) * r * Z**1.5 * (2 - r*Z) * np.exp(-Z*r/2)
  elif n == 2 and l == 1:
    return (np.sqrt(6)/12) * r * Z**1.5 * (r*Z) * np.exp(-Z*r/2)
  else:
    raise Exception("Unimplemented")


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

def hydrogenic_bound_bound(n,q,Z):
  """
  Eisenberger & Platzman eq. (32)
  Parameters:
    n: principal quantum number
    q: momentum transfer in inv. Bohr
    Z: effective nuclear charge
  """
  a = 1.0 / Z
  qa = q*a
  n = float(n)
    
  A = (2**8/3.) * qa**2 / n**3 * (3 * qa**2 + (n**2-1)/n**2)
  B = ((n-1)**2/n**2 + qa**2)**(n-3)
  C = (((n+1)/n)**2 + qa**2)**(n+3)
  S = A * B/C

  return S

def schumacher_phi(n,l,y):
    if n == 1:
        if l == 0:
            return 1.0 / (3.0 * y**3) # was this a typo in Schumacher (which has y^2)?
    elif n == 2:
        if l == 0:
            return 4.0 * (1.0 / (3.0 * y**3) -  1.0 / y**4 + 4.0 / (5.0 * y**5))
        elif l == 1:
            return 1.0 / (4 * y**4) - 1.0 / (5 * y**5)
    raise NotImplementedError

def schumacher_IA_jpq(pq,n,l,Z=1):
    Za = Z
    return 2**(4*l+3)/np.pi * factorial(n-l-1)/factorial(n+l) * n**2 * factorial(l)**2/Za * schumacher_phi(n,l,1+pq**2*n**2/Za**2)

def schumacher_IA_sqw(q,w,n,l,Z=1):
    pq = w/q - q/2
    return schumacher_IA_jpq(pq, n, l, Z) / q



def xrts_IA(w1, w2, theta, n=1, l=0, Z=3.81, B=112,  scale_rk=True, scale_sinc=True):
  """
  Reproduce truncated IA calculation from XRTS code.

  w1 in eV
  w2 in eV
  theta in degrees

  If `scale_rk` is True, then a factor 1-|f|^2 is included
  If `scale_sinc` is True than factor of (1+q^2/2 / w1)^-3 is included
  """

  w = w1 - w2
  q = momentum_transfer(w1, w2, theta*np.pi/180)*BOHR

  S = schumacher_IA_sqw(q,w/HARTREE,n,l,Z)/HARTREE
  S[w<B] = 0
  if scale_rk:
    fi2 = hydrogenic_bound_bound(n,q,Z)
    S *= 1 - fi2 # verified
  if scale_sinc:
    sinc = (1 + (q**2/2 * HARTREE / w1))**-3
    S *= sinc # verified
    S *= (w2/w1)**2 # verified

  return S

def xrts_PWFFA(w1, w2, theta, n=1, l=0, Z=3.81, B=112,  scale_rk=True, scale_sinc=True):
  if n != 1 or l != 0:
    raise Exception("Unimplemented. Only 1s supported.")

  S = S_hydrogenic_1s(w1, w2, theta, Z, B, type='PWFFA')
  q = momentum_transfer(w1, w2, theta*np.pi/180)*BOHR
  if scale_rk:
    fi2 = hydrogenic_bound_bound(n,q,Z)
    S *= 1 - fi2
  if scale_sinc:
    sinc = (1 + (q**2/2 * HARTREE / w1))**-3
    S *= sinc
    S *= (w2/w1)**2

  return S

def xrts_scale_factor(w1, q, Z, n=1):
  fi2 = hydrogenic_bound_bound(n,q,Z)
  rk = 1 - fi2 # r_k factor

  recoil = (1 + (q**2/2 * HARTREE / w1))**-3 # Breit-Dirac recoil

  return rk * recoil
