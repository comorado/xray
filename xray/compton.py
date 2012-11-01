import numpy as np

def isotropic_profile(rho_p_func, pq, *args):
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

  J0 = quad(integrand, 0, Inf, args)[0]
  dJ = np.array([quad(integrand, 0, pqi, args)[0] for pqi in pq])

  J = J0 - dJ
  return J

def fermi_profile(pq, N, V):
  from .fermi import fermi_momentum 

  pf = fermi_momentum(N/float(V))
  return V / (2 * np.pi)**2 * (pf**2 - pq**2) * (pf >= np.abs(pq))



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
    B = 2*E*(1-np.cos(theta)) - 2*c**2
    C = 2*E**2*(1-np.cos(theta))
    return (-B - np.sqrt(B**2 - 4 * A * C)) / 2 * HARTREE
