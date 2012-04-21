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
  return V / (2 * np.pi)**2 * (pf**2 - pq**2) * (pf >= pq)
