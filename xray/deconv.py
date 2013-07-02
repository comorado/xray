import numpy as np

class ConvergenceFailed(Exception):
  """
  An exception indicating that deconvolution failed to converge.
  The final result is stored in the `result` parameter.
  """
  def __init__(self, result=None):
    super(ConvergenceFailed, self).__init__()
    self.result = result


def richardson_lucy_step(data, point_spread, prior, noise_filter=None):
  cur = prior * np.convolve(data/np.convolve(prior, point_spread, 'same'),
                            point_spread, 'same')
  if noise_filter is not None:
    cur = np.convolve(cur, noise_filter, 'same')

  return cur

def richardson_lucy(data, point_spread, noise_filter=None, abserr=1e-5, maxiter=1e4):
  prior = data
  cur = None

  delta = abserr*2
  iter = 0
  while delta > abserr and iter < maxiter:

    cur = richardson_lucy_step(data, point_spread, prior, noise_filter)
    delta = max(abs(prior - cur))
    prior = cur
    iter += 1

  if iter >= maxiter:
    raise ConvergenceFailed(cur)

  return cur

def unit_gaussian(x, sigma):
  y = np.exp(-(x/sigma)**2 / 2.0)
  return y / y.sum()

def unit_gaussian_grid(dx, sigma):
  xmax = dx * np.floor(4*sigma / dx)
  x = np.arange(-xmax, xmax+dx, dx)
  return unit_gaussian(x, sigma)

def deconvolve_gaussian(x, y, sigma, nsigma=None, abserr=1e-5, maxiter=1e4):
  dx = x[1] - x[0]

  # generate normalized gaussian for deconvolution
  gy = unit_gaussian_grid(dx, sigma)

  ny = None
  if nsigma is not None:
    ny = unit_gaussian_grid(dx, nsigma)

  result = richardson_lucy(y, gy, ny, abserr, maxiter)
  if result.shape != y.shape:
    # XXX trim result to y.shape
    pass

  return result
