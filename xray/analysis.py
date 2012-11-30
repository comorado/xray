# This file uses the text encoding: utf-8
import sys
import operator
import numpy as np

class Curve(object):
  def __init__(self, x, y, s=None, **kwargs):
    self.x = np.asarray(x)
    self.y = np.asarray(y)

    if s is not None:
      s = np.asarray(s)
    self.s = s

    self.meta = kwargs

  def __getattr__(self, name):
    try:
      return self.meta[name]
    except KeyError:
      raise AttributeError, "'{cls}' object has no attribute '{name}'".format(cls=self.__class__, name=name), sys.exc_info()[2]

  @classmethod
  def from_file(cls, filename, x=0, y=1, s=None, **kwargs):
    """
    Load curve from file

    Parameters:
      filename - name of file to load
      x: column for x values
      y: column for y values
      s: column for s values

    Additionally, any named parameters from numpy.loadtxt may be provided. However, the value
    of the `unpack` parameter is ignored.

    The x,y,s column values can be either integers or a callable that is
    passed the full data as a transposed numpy array. See second example
    below.

    Examples:

    >>> from StringIO import StringIO
    >>> s = StringIO(\"\"\"
    ...       0 0 1 2
    ...       1 0 2 1
    ...       2 1 3 2
    ...       3 1 4 1
    ...     \"\"\")
    >>> c = Curve.from_file(s)
    >>> c.x
    array([ 0.,  1.,  2.,  3.])
    >>> c.y
    array([ 0.,  0.,  1.,  1.])

    >>> s.seek(0) # rewind to beginning of data
    >>> c = Curve.from_file(s, 0, lambda d: d[1:].sum(0))
    >>> c.y # contains sums of all columns except first
    array([ 3.,  3.,  6.,  6.])

    """

    kwargs['unpack'] = True

    data = np.loadtxt(filename, **kwargs)

    x = x(data) if callable(x) else data[x]
    y = y(data) if callable(y) else data[y]

    if s is not None:
      s = s(data) if callable(s) else data[s]

    return cls(x,y,s, filename=filename)

  def copy(self):
    return self.__class__(self.x.copy(), self.y.copy(), self.s.copy() if self.s else None, **self.meta)

  def plot(self, **kwargs):
    from matplotlib import pyplot as pt

    axes = kwargs.pop('axes', None)
    if axes is None: axes = pt.gca()
    errorbars = kwargs.pop('errorbars', False)

    if errorbars:
      return axes.errorbar(self.x, self.y, self.s, **kwargs)
    else:
      return axes.plot(self.x, self.y, **kwargs)

  def save(self, filename):
    if self.s:
      data = np.c_[self.x, self.y, self.s]
    else:
      data = np.c_[self.x, self.y]

    np.savetxt(filename, data)

  def moment(self, n=0, normed=False):
    """
    Calculate the nth moment

    Parameters:
      n: which moment to calculate
      normed: whether to normalize by 0th moment

    Return:
      If normed is False (default):  ∫dx y x^n
      If normed is True: ∫dx y x^n / ∫dx y 

    """
    ret = np.trapz(self.y * self.x**n, self.x)
    if normed: ret /= np.trapz(self.y, self.x)

    return ret

def _ranges_to_index(x, ranges):
  """
  """
  if ranges:
    return reduce(np.logical_or, [(x1 <= x) & (x <= x2) for x1,x2 in ranges])
  else:
    return np.ones(x.shape, dtype=bool)

def processor(func):
  """
  Function decorator that adds an identically named method to the Curve class.

  `func` must take a Curve as its first parameter, and should return a Curve

  TODO: add ability to customize name
  """

  setattr(Curve, func.__name__, func)
  return func

@processor
def interpolate(curve, new_x, left=0, right=0):
  """
  Linearly interpolate `curve` onto new grid.
  """
  new_y = np.interp(new_x, curve.x, curve.y, left, right)

  new_s = None
  if curve.s is not None:
    new_s = np.interp(new_x, curve.x, curve.s, left, right)

  return curve.__class__(new_x, new_y, new_s)

@processor
def interpolate_spline(curve, new_x, order=3, smooth=None, der=0):
  from scipy.interpolate import splrep, splev
  w = None
  if curve.s is not None:
    w = 1.0/curve.s

  tck = splrep(curve.x, curve.y, w, k=order, s=smooth)
  new_y = splev(new_x, tck, der, ext=1)

  new_s = None
  if curve.s is not None:
    new_s = np.interp(new_x, curve.x, curve.s, 0, 0)

  return curve.__class__(new_x, new_y, new_s)

@processor
def extend_symmetric(curve):
  """
    Include curve flipped about x==0

    Note: this the x-values must be all positive or all negative
  """
  if all(curve.x >= 0):
    start = 0
    if curve.x[0] == 0:
      start = 1
    new_x = np.append(-curve.x[start:][::-1], curve.x)
    new_y = np.append(curve.y[start:][::-1], curve.y)
    new_s = None
    if curve.s is not None:
      new_s = np.append(curve.s[start:][::-1], curve.s)
  elif all(curve.x <= 0):
    new_s = None
    if curve.x[-1] == 0:
      new_x = np.append(curve.x, -curve.x[:-1][::-1])
      new_y = np.append(curve.y, curve.y[:-1][::-1])
      if curve.s is not None:
        new_s = np.append(curve.s, curve.s[:-1][::-1])
    else:
      new_x = np.append(curve.x, -curve.x[:][::-1])
      new_y = np.append(curve.y, curve.y[:][::-1])
      if curve.s is not None:
        new_s = np.append(curve.s, curve.s[:][::-1])
  else:
    raise ValueError("Data to be symmetrized must have x-values all positive or all negative.")

  return curve.__class__(new_x, new_y, new_s)

@processor
def broaden(curve, sigma):
  """
  Broaden curve by convolving with a gaussian with parameter `sigma`.

  Preconditions:
    curve.x must be evenly spaced, and should be dense relative to sigma for results to be accurate
  """
  d = np.diff(curve.x)
  #assert np.allclose(d, d[0])

  dx = d[0]
  gx = np.arange(-4*sigma, 4*sigma+dx, dx)
  gy = np.exp(-gx**2 / (2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)
  new_y = np.convolve(curve.y,gy,'same') / gy.sum()

  if len(gx) > len(curve.x):
    new_y = np.interp(curve.x,gx,new_y)

  # XXX How should uncertainty be handled (stripped for now)?
  return curve.__class__(curve.x, new_y)

@processor
def broaden_fwhm(curve, fwhm):
  """
  Broaden curve by convolving with a gaussian with full-width half-max `fwhm`.

  Preconditions:
    curve.x must be evenly spaced, and should be dense relative to sigma for results to be accurate
  """

  sigma = fwhm / (2 * np.sqrt(2*np.log(2)))
  return curve.broaden(sigma)


@processor
def normalize_integral(curve, norm=1.0, xmin=None, xmax=None):
  """
  Normalize curve to set integral to a given value.
  """
  from scipy.interpolate import splrep, splint

  if xmin is None: xmin = curve.x[0]
  if xmax is None: xmax = curve.x[-1]

  tck = splrep(curve.x, curve.y)
  n = splint(xmin, xmax, tck)

  new_y = curve.y * norm / n
  new_s = None
  if curve.s is not None:
    new_s = curve.s * norm / n

  return curve.__class__(curve.x, new_y, new_s)


@processor
def normalize_integral_trapz(curve, norm=1.0):
  """
  Normalize curve to set integral to a given value.
  """
  n = np.trapz(curve.y, curve.x)

  new_y = curve.y * norm / n
  new_s = None
  if curve.s is not None:
    new_s = curve.s * norm / n

  return curve.__class__(curve.x, new_y, new_s)

@processor
def background_subtract(curve, fit_ranges, degree):
  """

  """

  index = _ranges_to_index(curve.x, fit_ranges)
  x = curve.x[index]
  y = curve.y[index]

  p = np.polyfit(x, y, degree)
  bg = np.polyval(p, curve.x)

  # XXX handle uncertainty...

  return curve.__class__(curve.x, curve.y - bg, curve.s, background_subtract=curve.__class__(curve.x, bg))

@processor
def fit_gaussian(curve, fit_ranges=None, guess=(1.0,0.0,1.0)):
  index = _ranges_to_index(curve.x, fit_ranges)
  x = curve.x[index]
  y = curve.y[index]

  def model_gaussian(params, x):
    A,x0,sigma = params
    return A * np.exp(-(x-x0)**2/2/sigma**2)

  def error_gaussian(params, *args):
    x,y = args
    return model_gaussian(params, x) - y

  from scipy.optimize import leastsq
  fit, ier = leastsq(error_gaussian, guess, (x,y))

  # XXX check for failure of fit
  new_y = model_gaussian(fit, curve.x)
  return curve.__class__(curve.x, new_y, curve.s, gaussian_parameters=fit)

@processor
def fit_lorentzian(curve, fit_ranges=None, guess=(1.0,0.0,1.0)):
  index = _ranges_to_index(curve.x, fit_ranges)
  x = curve.x[index]
  y = curve.y[index]

  def model_lorentzian(params, x):
    A,x0,gamma = params
    return A * 1.0 / ((x-x0)**2 + (gamma/2.0)**2)

  def error_lorentzian(params, *args):
    x,y = args
    return model_lorentzian(params, x) - y

  from scipy.optimize import leastsq
  fit, ier = leastsq(error_lorentzian, guess, (x,y))

  # XXX check for failure of fit
  new_y = model_lorentzian(fit, curve.x)
  return curve.__class__(curve.x, new_y, curve.s, lorentzian_parameters=fit)


def voigt(x,amp,pos,fwhm,shape):
  from scipy import special
  """
  Voigt profile

  V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
  z = (x+i*gam)/(sig*sqrt(2))

  http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
  """

  tmp = 1/special.wofz(np.zeros((len(x))) \
        +1j*np.sqrt(np.log(2.0))*shape).real
  tmp = tmp*amp* \
        special.wofz(2*np.sqrt(np.log(2.0))*(x-pos)/fwhm+1j* \
        np.sqrt(np.log(2.0))*shape).real
  return tmp

@processor
def fit_voigt(curve, fit_ranges=None, guess=(1.0,0.0,1.0, 0.1)):
  index = _ranges_to_index(curve.x, fit_ranges)
  x = curve.x[index]
  y = curve.y[index]

  def error_voigt(params, *args):
    x,y = args
    return voigt(x, *params) - y

  from scipy.optimize import leastsq
  fit, ier = leastsq(error_voigt, guess, (x,y))

  # XXX check for failure of fit
  new_y = voigt(curve.x, *fit)
  return curve.__class__(curve.x, new_y, curve.s, voigt_parameters=fit)



@processor
def fit_voigt_bg(curve, fit_ranges=None, guess=(1.0,0.0,1.0, 0.1, 0, 0)):
  """
  guess: tuple of initial fit parameters (A, x0, fwhm, shape, bg0, bg1)
  """
    
  index = _ranges_to_index(curve.x, fit_ranges)
  x = curve.x[index]
  y = curve.y[index]

  def voigt_bg(params, x):
    A,x0,fwhm,shape,bg0,bg1 = params
    v = voigt(x, A,x0,fwhm,shape)
    bg = bg0 + x*bg1
    return v + bg

  def error(params, *args):
    x,y = args
    return voigt_bg(params, x) - y

  from scipy.optimize import leastsq
  fit, ier = leastsq(error, guess, (x,y))

  # XXX check for failure of fit
  new_y = voigt_bg(fit, curve.x)
  return curve.__class__(curve.x, new_y, curve.s, voigt_parameters=fit)

def _savitzky_golay(y, window_size, order, deriv=0):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [0] http://www.scipy.org/Cookbook/SavitzkyGolay
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m, y, mode='valid')

@processor
def savitzky_golay(curve, window_size, order, deriv=0):
  new_y = _savitzky_golay(curve.y, window_size, order, deriv)
  return curve.__class__(curve.x, new_y, curve.s)



def _spline_integral(x, y, xmin=None, xmax=None):
  """
  Calculate spline integral of curve from xmin to xmax
  """
  from scipy.interpolate import splrep, splint

  if xmin is None: xmin = x[0]
  if xmax is None: xmax = x[-1]

  tck = splrep(x, y)
  return splint(xmin, xmax, tck)

@processor
def spline_integral(curve, xmin=None, xmax=None):
  return _spline_integral(curve.x, curve.y, xmin, xmax)

@processor
def center_of_mass(curve, xmin=None, xmax=None):
  return _spline_integral(curve.x, curve.y*curve.x, xmin, xmax) / _spline_integral(curve.x, curve.y, xmin, xmax)


@processor
def fit_polynomial(curve, fit_ranges, degree):
  """

  """

  index = _ranges_to_index(curve.x, fit_ranges)
  x = curve.x[index]
  y = curve.y[index]

  p = np.polyfit(x, y, degree)
  bg = np.polyval(p, curve.x)

  return curve.__class__(curve.x, bg, params=p)

@processor
def subtract(curve, curve2, interpolate=None, **kwargs):
  """
  Subtract a curve2 from curve 

  Parameters:
    curve - Curve to subtract from
    curve2 - Curve to subtract
    interpolate - None (default), 'linear', or 'spline'
    kwargs - keyword args to pass on to the interpolation method
             see `interpolate()` and `interpolate_spline()`

  If curve2.x differs from curve.x, then curve2 is interpolated onto curve.x using the method specified by `interpolate`. If  `interpolate` is None, a ValueError is raised instead.
  """

  if not np.all(curve.x == curve2.x):
    if interpolate == 'linear' :
      curve2 = curve2.interpolate(curve.x, **kwargs)
    elif interpolate == 'spline':
      curve2 = curve2.interpolate_spline(curve.x, **kwargs)
    else:
      raise ValueError("curve and curve2 must have same x values, or an interpolation method must be specificied")

  new_x = curve.x
  new_y = curve.y - curve2.y

  if curve.s is not None and curve2.s is not None:
    new_s = np.sqrt(curve.s**2 +curve2.s**2)
  elif curve.s is not None:
    new_s == curve.s
  elif curve2.s is not None:
    new_s == curve2.s

  return curve.__class__(new_x, new_y, new_s)

@processor
def combine(curve, curve2, op, interpolate=None, **kwargs):
  """
  Combine y values of `curve` and `curve2` using `operator`.

  Parameters:
    curve - left side argument
    curve2 - right side argument
    op - operator used to combine y values
    interpolate - None (default), 'linear', or 'spline'
    kwargs - keyword args to pass on to the interpolation method
             see `interpolate()` and `interpolate_spline()`

  If curve2.x differs from curve.x, then curve2 is interpolated onto curve.x using the method specified by `interpolate`. If  `interpolate` is None, a ValueError is raised instead.
  """

  if (hasattr(curve2, 'x')):
    if not np.all(curve.x == curve2.x):
      if interpolate == 'linear' :
        curve2 = curve2.interpolate(curve.x, **kwargs)
      elif interpolate == 'spline':
        curve2 = curve2.interpolate_spline(curve.x, **kwargs)
      else:
        raise ValueError("curve and curve2 must have same x values, or an interpolation method must be specificied")
    y2 = curve2.y
    s2 = curve2.s
  else:
    y2 = curve2
    s2 = None

  new_x = curve.x
  new_y = op(curve.y, y2)

  if curve.s is not None and s2 is not None:
    # XXX this isn't correct for all operations...
    new_s = np.sqrt(curve.s**2 +s2**2)
  elif curve.s is not None:
    new_s == curve.s
  elif s2 is not None:
    new_s == s2
  else:
    new_s = None

  return curve.__class__(new_x, new_y, new_s)

# standard operations
@processor
def __add__(curve, curve2, interpolate='linear', **kwargs):
  return combine(curve, curve2, operator.add, interpolate, **kwargs)

@processor
def __sub__(curve, curve2, interpolate='linear', **kwargs):
  return combine(curve, curve2, operator.sub, interpolate, **kwargs)

@processor
def __mul__(curve, curve2, interpolate='linear', **kwargs):
  return combine(curve, curve2, operator.mul, interpolate, **kwargs)

@processor
def __div__(curve, curve2, interpolate='linear', **kwargs):
  return combine(curve, curve2, operator.div, interpolate, **kwargs)

@processor
def fsum(curve, q, N):
  from .const import HBARC, MC2
  norm = N * (HBARC**2) * q**2 / MC2
  ret = curve.copy()
  ret.y *=  norm / np.trapz(ret.x*ret.y, ret.x)
  return ret


##################
# Additional API #
##################

def plot_curves(curves, styles=None):
  if styles is None:
    styles = [{}]
  elif isinstance(styles, dict):
    styles = [styles]
  elif not isinstance(styles, list):
    raise ValueError("Invalid styles")

  styles *= len(curves) / len(styles) + 1
  for c, s in zip(curves, styles):
    c.plot(**s)

