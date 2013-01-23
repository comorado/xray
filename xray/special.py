"""
Special functions
"""
import numpy as np

def sph_jn(v, x):
  """
  Spherical Bessel Function j_v(x)
  """
  from scipy.special import jn

  single_x = np.isscalar(x)
  x = np.atleast_1d(x)

  y = jn(v + 0.5, x)
  i = (x != 0)
  y[i] *= np.sqrt(np.pi/(2*x[i]))
  y[x==0] = 1 if v == 0 else 0

  if single_x:
    y = y[0]

  return y

def laguerre_generator(alpha, x):
  """
  Calculate Associated Laguerre polynomials, L^alpha_n(x) using a recurrence relation

  Parameters:
    alpha: the upper index
    x: values to evaluate the polynomial at

  Returns:
    g: A generator of Associated Laguerre polynomials

  The (n+1)th call to g.next() will return L^alpha_n(x)
  """
  l1 = 0
  l = np.ones_like(x)

  i = 1.0
  while True:
    yield l
    l0 = l1
    l1 = l
    l = ((2 * i - 1 + alpha - x) * l1 - (i - 1 + alpha) * l0) / i
#    print "[",i,"]", l1, l0, "(", (2 * i - 1 + alpha - x), (i - 1 + alpha), ") => ", l
    i += 1.0

def laguerre_list(nmax, alpha, x):
  """
  Calculate Associated Laguerre polynomials, L^alpha_n(x) using a recurrence relation for n in [0,nmax]
  """

  g = laguerre_generator(alpha, x)
  return [g.next() for i in range(nmax+1)]

def laguerre(n, alpha, x):
  """
  Calculate Associated Laguerre polynomials, L^alpha_n(x) using a recurrence relation for n in [1,nmax]
  """

  g = laguerre_generator(alpha, x)
  for i in range(n+1):
    l = g.next()

  return l

def hydrogenic_radial(n,l,r, Z=1):
  """
  r assumed to be in units of Bohr
  """
  from scipy.misc import factorial
  from scipy.special import genlaguerre

  # Messiah's L's are different from ours, so this norm is incorrect
  #N = (2.0 / n**2) * np.sqrt(factorial(n-l-1) / factorial(n+l)**3)

  # correct norm:
  N = (2.0/n**2) * np.sqrt(factorial(n-l-1) / factorial(n+l))
  x = 2.0 * r * Z / n

  #L = laguerre(n-l-1, 2*l+1, x)
  L = genlaguerre(n-l-1, 2*l+1)(x)
  F = x**l * np.exp(-x/2.0) * L

  R = Z**(1.5) *  F * N

  return R
