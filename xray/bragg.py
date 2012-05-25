import numpy as np
from const import HC, pi, deg, rad

UNKNOWN    = 0
DIAMOND    = 1
ROCKSALT   = 2
ZINCBLENDE = 3
WURTZITE   = 4
GAP        = 5

crystal_info = {
    'Diamond': (DIAMOND,    3.56679),
    'Ge':      (DIAMOND,    5.65735),
    'Si':      (DIAMOND,    5.4309),
    'NaCl':    (ROCKSALT,   5.63978),
    'LiF':     (ROCKSALT,   4.0263),
    'AlAs':    (ZINCBLENDE, 5.6611),
    'GaAs':    (ZINCBLENDE, 5.6533),
    'GaP':     (GAP, 5.4512),
    'InAs':    (ZINCBLENDE, 6.0584),
    'InSb':    (ZINCBLENDE, 6.4782),
    'SiC':     (ZINCBLENDE, 4.348),
    'InP':     (UNKNOWN,    5.8688),
    'CsF':     (UNKNOWN,    6.008),
    'KCl':     (UNKNOWN,    6.29294),
    'CsCl':    (UNKNOWN,    7.02),
    }

lattice_constants = {k:crystal_info[k] for k in crystal_info}
#
# FILTERING
#

def all_even(arr, axis=None):
  return ((arr % 2) == 0).all(axis)

def all_odd(arr, axis=None):
  return ((arr % 2) == 1).all(axis)

def sum_not_congruent(arr, modulus, remainder, axis=None):
  return np.logical_not(arr.sum(axis) % modulus == remainder)

def sum_congruent(arr, modulus, remainder, axis=None):
  return arr.sum(axis) % modulus == remainder

def filter_hkl_diamond(hkl):
  """
  Selection Rules:
    h,k,l either all odd or all even and (h+k+l) = 4n
  """
  index = np.logical_or(
      all_odd(hkl,1),
      np.logical_and(
        all_even(hkl,1),
        sum_congruent(hkl,4,0,axis=1)
        )
      )

  return hkl[index]

def filter_hkl_rocksalt(hkl):
  index = np.logical_or(all_even(hkl,1), all_odd(hkl,1))
  return hkl[index]

def filter_hkl_gap(hkl):
  """
  GaP is only available in a few cuts: 100, 110 and 111
  XXX finish this up...
  """
  hkl = filter_hkl_rocksalt(hkl)
  index = np.array([np.all((row == row[0])|(row == 0)) for row in hkl])
  return hkl[index]
  

def filter_hkl_dummy(hkl):
  return hkl

filters = [
    filter_hkl_dummy,
    filter_hkl_diamond,
    filter_hkl_rocksalt,
    filter_hkl_rocksalt,
    filter_hkl_dummy,
    filter_hkl_gap,
    ]

#
# hkl matrix operations
#

def generate_hkl_matrix(maximum):
  return np.vstack([np.array([h,k,l]) 
    for h in range(0,maximum+1)
    for k in range(0, h+1)
    for l in range(0 if h+k > 0 else 1, k+1)
    ])

#
# angle calculations
#

def angle_between(v1, v2):
  """
  Calculate the angle between two vectors of arbitrary dimension
  """

  return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def bragg_angle(energy, d, m=None, degrees=False):
  """
  Calculate the Bragg angle for reflection of `energy` off of the plane specified by `m` with lattice constant `d`

  Parameters
  ----------
    energy: the energy of incident light
    d: the lattice constant in Angstroms
    m: a 3 component vector containt he Miller indices (h,k,l)
    degrees: if True, return result in degrees, otherwise radians

  Returns
  -------
    The angle in degrees if `degrees` is True, otherwise radians.
    If the Bragg condition is never satsified for the given parameters, np.nan will be returned 
  """
 
  if m is not None:
    d /= np.linalg.norm(m)

  sintheta = HC / (2 * d * energy)

  oldsettings = np.seterr(all='ignore')
  ret = np.arcsin(sintheta)
  np.seterr(**oldsettings)

  if degrees:
    ret = deg(ret)

  return ret

def bragg_energy(theta, d, m=None, degrees=False):
  """
  Energy scattered at Bragg angle `theta` from the `m` planes of a crystal with lattice spacing `d` 

  Parameters
  ----------
    theta - Bragg angle (units depend on value of `degrees`
    d - lattice spacing in Angstroms
    m - 3-vector of Miller indices (h,k,l)
    degrees - if True, theta is in degrees, otherwise radians

  Returns
  -------
    energy in eV
  """
  if m is not None:
    d /= np.linalg.norm(m)

  if degrees:
    theta = rad(theta)

  return HC/(2*d) / np.sin(theta)

def bragg_denergy(theta, d, m=None, degrees=False):
  """
  -dE/dtheta at Bragg angle `theta` from the `m` planes of a crystal with lattice spacing `d` 

  Parameters
  ----------
    theta - Bragg angle (units depend on value of `degrees`
    d - lattice spacing in Angstroms
    m - 3-vector of Miller indices (h,k,l)
    degrees - if True, theta is in degrees, otherwise radians

  Returns
  -------
    energy in eV
  """
  if m is not None:
    d /= np.linalg.norm(m)

  if degrees:
    theta = rad(theta)

  return HC/(2*d) * np.cos(theta) / np.sin(theta)**2


#
# entry points
#

def find_reflections(xtal_names, energy, min_angle, max_index=10):
  reflections = []

  hkl_full = generate_hkl_matrix(max_index)

  theta_min = min_angle * pi / 180.
  theta_max = pi / 2
  
  for xtal_name in xtal_names:
    xtal_type, lattice_const = crystal_info[xtal_name]
    hkl = filters[xtal_type](hkl_full)
    hkl_norm = np.apply_along_axis(np.linalg.norm, 1, hkl)
    unique_norm = np.unique(hkl_norm)
    for n in unique_norm:
      d = lattice_const / n
      theta = bragg_angle(energy, d)
      if theta_min <= theta <= theta_max:
        index = np.where(hkl_norm == n)[0]
        planes = [list(hkl[i]) for i in index]
        reflections.append((xtal_name, theta*180.0/np.pi, planes))

  return reflections
  
def find_alternate_reflection0(xtal_name, m0, energy, min_angle, max_index=10):
  #hkl_full = np.vstack([ [h,k,l]
  #  for h in range(0,max_index)
  #  for k in range(0,max_index)
  #  for l in range(0,max_index)
  #  ])
  hkl_full = generate_hkl_matrix(max_index)

  xtal_type, lattice_const = crystal_info[xtal_name]
  hkl = filters[xtal_type](hkl_full)
  #hkl_norm = np.linalg.norm(hkl_full)

  theta_min = min_angle * pi / 180.
  theta_max = pi / 2

  reflections = []

  for m in hkl:
    phi = np.abs(angle_between(m, m0))
    theta = bragg_angle(energy, lattice_const, m)

    r = 1.0*m[2] / m[1]
    if r < 1.0/3.0:
      theta2 = theta + phi
    else:
      theta2 = theta - phi

    if theta_min <= theta2 <= theta_max:
      reflections.append((theta2 * 180.0/np.pi, list(m)))

  return reflections


def find_alternate_reflections(xtal_name, m0, E1, E2, max_angle, max_index=10, include_forbidden=False):
  #hkl_full = np.vstack([ [h,k,l]
  #  for h in range(0,max_index)
  #  for k in range(0,max_index)
  #  for l in range(0,max_index)
  #  ])
  hkl_full = generate_hkl_matrix(max_index)

  xtal_type, lattice_const = crystal_info[xtal_name]
  if include_forbidden:
    hkl = hkl_full
  else:
    hkl = filters[xtal_type](hkl_full)

  #hkl_norm = np.linalg.norm(hkl_full)

  theta1 = bragg_angle(E1, lattice_const, m0)*180.0/pi
  theta2 = bragg_angle(E2, lattice_const, m0)*180.0/pi

  reflections = []

  for m in hkl:
    phi = np.abs(angle_between(m, m0))*180.0/pi
    if phi < max_angle:
      energies = [
          bragg_energy((theta1+phi)*pi/180.0, lattice_const, m),
          bragg_energy((theta2+phi)*pi/180.0, lattice_const, m),
          ]
      reflections.append((phi, list(m), energies))
  

  return reflections

def find_alternate_energy_ranges(theta1, theta2, max_index=5, pixels=195, xtals=None):
  hkl_full = np.vstack([ [h,k,l]
    for h in range(0,max_index)
    for k in range(0,h+1)
    for l in range(0,k+1)
    ])

  ranges = []

  if xtals is None:
    xtals = crystal_info

  for xtal_name in xtals:
    xtal_type, d = crystal_info[xtal_name]
    hkl = filters[xtal_type](hkl_full)

    for m in hkl:
      emax = bragg_energy(theta1, d, m)
      emin = bragg_energy(theta2, d, m)
      res =  abs(emax-emin)/pixels

      ranges.append([xtal_name, m, emin, emax, res])

  return ranges

if __name__ == "__main__":
  import sys

  if len(sys.argv) < 2:
    print "Usage: python bragg.py energy [min angle]\n"
    exit()

  energy = float(sys.argv[1])

  min_angle = 70
  if len(sys.argv) > 2:
    min_angle = float(sys.argv[2])

  reflections = find_reflections(['Si', 'Ge'], energy, min_angle)

  for xtal, angle, cuts in reflections:
    print xtal, angle, cuts
   

