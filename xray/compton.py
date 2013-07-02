import numpy as np
from . import const
from .analysis import Curve

def isotropic_profile(rho_p_func, pq, args):
  """
  Calculates compton profile for an isotropic momentum density
  """

  from scipy.integrate import quad
  from scipy.integrate.quadpack import Inf

  integrand = lambda p, *args: rho_p_func(p, *args) * p * 2 * np.pi

  # Integration to Inf appears to be numerically unstable for some pq
  # So, we use \int_pq^\inf = \int_0^inf - \int_0^pq instead.
  #
  #data = np.array([quad(integrand, pqi, Inf, args) for pqi in pq])

  J0 = quad(integrand, 0, Inf, args, epsabs=1e-3, epsrel=1e-3)[0]
  dJ = np.array([quad(integrand, 0, pqi, args, epsabs=1e-3, epsrel=1e-3)[0] for pqi in pq])

  J = J0 - dJ
  return J

def fermi_profile(pq, N, V, T=0, mu=None):
  from .fermi import fermi_momentum, mu_T, rhop

  if T < 1e-5:
    pf = fermi_momentum(N/float(V))
    return V / (2 * np.pi)**2 * (pf**2 - pq**2) * (pf >= np.abs(pq))
  else:
    if mu is None:
      mu = mu_T(N,V,T)
    return isotropic_profile(rhop, pq, (N,V,T,mu))

def compton_shift(E, theta):
    """
    Calculate the Compton shift at a fixed scattering angle

    Parameters:
      E: fixed energy in eV
      theta: angle in degrees"

    Returns:
      Compton shift in eV

    For fixed incident energy experiment (e.g. laser-shock XRTS), the peak will be located at the scattered energy E - Ec.
    For fixed scattered energy (e.g. LERIX), the peack will be located at incident energy E + Ec.
    """

    from xray.const import HARTREE
    c = 137.035999074 # inverse fine structure from CODATA

    E = E / HARTREE # convert to a.u.
    theta = theta * np.pi / 180 # convert to radians
    A = 1
    B = -2*E*(1-np.cos(theta)) - 2 *c**2
    C = 2*E**2*(1-np.cos(theta))
    return (-B - np.sqrt(B**2 - 4 * A * C)) / 2 * HARTREE

def momentum_transfer(E1, E2, theta):
  """
  Square of momentum transfer in energy units for inelastically scattering from `E1` to `E2` by angle `theta`.

  Parameters:
    E1: incident  photon momentum (eV)
    E2: scattered photon momentum (eV)
    theta: scattering angle (radians)

  Returns:
    q: momentum transfer in inverse Angstroms
  """
  return np.sqrt(E1 * E1 + E2 * E2 - 2 * E1 * E2 * np.cos(theta)) / const.HBARC


def projected_momentum(w, q):
  """
  Magnitude of electron momentum in direction of momentum transfer in Impulse Approximation.

  Parameters:
    w: energy transfer in eV
    q: momentum transfer in inv. A
  Returns:
    p_q: momentum projected onto direction of momentum transfer in inverse Angstroms
  """
  return w * const.MC2 / q / const.HBARC**2 - q / 2.0



def sqw_to_jpq(S, q, w):
  """
  Parameters:
    S: dynamic structure factor in 1/eV
    q: momentum transfer in 1/A
    w: energy transfer in eV

  Returns:
    Curve(pq, J)

    J: Compton profile in A
    pq: projected electron momentum in 1/A
  """

  pq = projected_momentum(w, q)
  J = q * S * const.HBARC**2 / const.MC2

  return Curve(pq, J)

def jpq_to_sqw(J, pq, w, q):
  """
  Parameters:
    J: Compton profile in Angstroms
    pq: projected momentum in 1/A
    w: energy loss in eV
    q: momentum transfer in inverse Angstroms

  Returns:
    Curve containing S(w; q)
  """

  pq_interp = projected_momentum(w,q)
  J_interp = np.interp(pq_interp, pq, J)
  S = J_interp / q * const.MC2 / const.HBARC**2

  return S

class ComptonProfile(Curve):

  def to_sqw(self, w, q):
    """
    Parameters:
      w: energy transfer in eV
      q: momentum transfer in 1/A

    Returns:
      Curve with (x=w,y=S(q, w))
    """
    S = jpq_to_sqw(self.y, self.x, w, q)
    sqw = Curve(w.copy(), S, q=q)
    return sqw

  def to_sqw_theta(self, w, theta, E1=None, E2=None):
    """
    Parameters:
      w: energy transfer in eV
      theta: angle in radians
      E1,E2: fixed incident and scattered energy (in eV). Only one of these
             may be specified (the other is determined by the condition
             w = E1 - E2)

    Returns:
      Curve containing S(w)
    """

    assert not (E1 is None and E2 is None), "One of E1 and E2 must be given"
    assert E1 is None or E2 is None, "Only one of E1 and E2 may be given"

    if E1 is None: E1 = E2 + w
    else: E2 = E1 - w

    q = momentum_transfer(E1, E2, theta)

    sqw = self.to_sqw(w,q)
    if E1 is not None: sqw.meta['E1'] = E1
    else: sqw.meta['E2'] = E2
    sqw.meta['theta'] = theta

    return sqw

  @classmethod
  def from_file(cls, filename, num_electrons=None, atomic_units=False,
                **kwargs):
    """
    Parameters:
      filename: name of file to load
      num_electrons: number of electrons to normalize curve to
      atomic_units: if True, convert from atomic units to Angstroms
      kwargs: see Curve.from_file for additional arguments
    """

    c = super(ComptonProfile, cls).from_file(filename, **kwargs)
    if atomic_units:
      c.x /= const.BOHR
      c.y *= const.BOHR

    if (np.all(c.x >= 0)):
      c = c.extend_symmetric()

    if num_electrons is not None:
      c = c.normalize_integral(num_electrons)

    return c

  def to_rhop(self, use_neg=False):
      der = self.differentiate()

      if not use_neg:
        i = der.x > 0
        der.x = der.x[i]
        der.y = (der.y[i] / (-2.0 * np.pi * der.x))
      else:
        i = der.x < 0
        der.x = -der.x[i]
        der.y = (der.y[i] / (2.0 * np.pi * der.x))
        der.x = der.x[::-1]
        der.y = der.y[::-1]

      return der

