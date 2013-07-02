import warnings
import numpy as np

def unsqueeze(arr, dim, axes=None):
  """
  Expands an array to dim dimensions, placing the original axes along the dimensions specified by axes.

  If axes is None, the old axes are placed at the end of the list. 

  For example:

    arr = [1,2,3]
  """
  arr = np.asarray(arr)

  assert(dim >= arr.ndim)

  if (dim == arr.ndim):
    return arr

  if axes is None:
    axes = np.arange(len(arr.shape)) + dim - len(arr.shape)

  assert(len(axes) == len(arr.shape))

  shape = [1] * dim
  for frm, to in zip(arr.shape, axes):
    shape[to] = frm

  return np.reshape(arr, shape)

def indices_bounding(arr, val):
  """
  Finds indices of the elements of an ordered array that are
  just below and just above a given value.

  i.e., finds i,j s.t. arr[i] <= val <= arr[j] and j-1 <= 1

  Parameters:
      arr - array to search (must be in ascending order)
      val - value to search for

  Returns:
    (i,j) such that arr[i] is the closest element below (or equal to) val and arr[j] is the closest element above (or equal to) val

    If val < arr[0], then (None,0) is returned
    If val > arr[-1], then (len(arr)-1,0) is returned

    Finally, if val is in arr, then (i,i) is returned where i is the index of val in arr.

  This uses a binary search, and thus runs in O(log(N)) time.
  """

  i1 = 0
  i2 = len(arr)-1

  if (val == arr[0]):  return 0,0
  if (val <  arr[0]):  return None,0
  if (val == arr[-1]): return i2,i2
  if (val >  arr[-1]): return i2,None

  while (i2 - i1 > 1):
    i = i1 + (i2-i1)/2
    if (arr[i] >= val): i2 = i
    if (arr[i] <= val): i1 = i

  return i1,i2


def trapezoid(x, y, x0=None, x1=None, axis=0):
  """
  Integrates y(x) between x0 and x1 using the trapezoid rule, including endpoint corrections.

  Parameters:
    x: array of x values
    y: array of y values (same length as x)

    x0: lower bound of integration range (or None for x[0])
    x1: upper bound of integration range (or None for x[-1])
  """

  x = np.asarray(x)
  y = np.asarray(y)

  if x0 is None: x0 = x[0] 
  if x1 is None: x1 = x[-1] 

  i0,j0 = indices_bounding(x, x0)
  i1,j1 = indices_bounding(x, x1)

  if i0 is None:
    warnings.warn("Lower bound is below range of data. Truncating integral to data range.")
    i0 = 0
    x0 = x[i0]

  if j1 is None:
    warnings.warn("Upper bound is above range of data. Truncating integral to data range.")
    j1 = len(x) - 1
    x1 = x[j1]

  if i1 is None:
    warnings.warn("Lower bound is above data range. Nothing to integrate.")
    return 0.0

  if j0 is None:
    warnings.warn("Upper bound is below data range. Nothing to integrate.")
    return 0.0

  # trapezoid
  shape = [y.shape[0]] + [1] * len(y.shape)
  shape[0] = y.shape[0]

  dx = [ (x1 if i == i1 else x[i+1]) -
         (x0 if i == j0 else x[i-1])
         for i in xrange(j0,i1+1) ]
  dx = unsqueeze(dx, y.ndim, [axis]) # match dimensions of y

  #yy = y[j0:i1+1]
  yy = np.take(y, np.arange(j0,i1+1), axis)
  ret = np.sum( yy * dx / 2.0, axis)

  # a helper to simplify grabbing slicing along one axis
  # this is the same as y[:,...,ii,...,:] where the ii is in the slot specified by 'axis'
  ty = lambda ii: np.squeeze(np.take(y, [ii], axis))

  # left endpoint correction
  if (j0 != i0):
    f = (x0 - x[i0]) / float(x[j0] - x[i0])
    y0 = (1 - f) * ty(i0) + f * ty(j0) 
    corr = y0 * (x[j0] - x0) / 2.0
    ret += corr

  # right endpoint correction
  if (j1 != i1):
    f = (x1 - x[i1]) / float(x[j1] - x[i1])
    y1 = (1 - f) * ty(i1) + f * ty(j1)
    corr = y1 * (x1 - x[i1]) / 2
    ret += corr

  return ret

def trapezoidn(xs, y, x0s=None, x1s=None, axes=None):
  if axes is None:
    axes = range(y.ndim)

  assert(len(xs) == len(axes))

  """
  if x0s is None:
    x0s = [x[0] for x in xs] 
  else:
    x0s = [x[0] if x0 is None else x0 for x0,x in zip(x0s, xs)]

  if x1s is None:
    x1s = [x[-1] for x in xs] 
  else:
    x1s = [x[-1] if x1 is None else x1 for x1,x in zip(x1s, xs)]
  """

  if x0s is None: x0s = [None] * len(xs)
  if x1s is None: x1s = [None] * len(xs)

  assert(len(x0s) == len(x1s) == len(xs)) 

  tmp = y
  for i,(x, x0, x1, axis,) in enumerate(zip(xs, x0s, x1s, axes)):
    tmp = trapezoid(x, tmp, x0, x1, axis - i)

  return tmp 


def interpolate(x0s, x, y, axis=0):
  """
  Interpolate y(x) onto points in x0s using piecewise linear interpolation

  If y is multidimensional then `axis` specifies the axis along which x varies
  """
  y0s = [interpolate1(x0, x, y, axis) for x0 in x0s]
  return np.concat([unsqueeze(y0, y.ndim, axis) for y0 in y0s], axis)

def interpolate1(x0, x, y, axis=0):
  """
  Interpolate y(x) onto the point x0 using piecewise linear interpolation

  If y is multidimensional then `axis` specifies the axis along which x varies
  """
  i,j = indices_bounding(x, x0)
  if i is None or j is None:
    raise ValueError("x0 is outside of range of x")

  s = [slice(None)] * y.ndim
  s[axis] = i
  yi = y[tuple(s)]

  if i == j: return yi

  s[axis] = j
  yj = y[tuple(s)]

  f = (x0 - x[i]) / (x[j] - x[i])
  y0 = yi * (1-f) + yj * f

  return  y0
