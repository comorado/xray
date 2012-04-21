"""
Special functions
"""
import numpy as np

def sph_jn(v, x):
  """
  Spherical Bessel Function j_v(x)
  """
  from scipy.special import jn

  single_x = False
  if not hasattr(x, '__getitem__'):
    x = [x]
    single_x = True

  x = np.array(x)

  y = jn(v + 0.5, x)
  i = (x != 0)
  y[i] *= np.sqrt(np.pi/(2*x[i]))
  y[x==0] = 1 if v == 0 else 0

  if single_x:
    y = y[0]

  return y

