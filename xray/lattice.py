from numpy import sqrt, sin, cos, pi

# lattice structures
ISOMETRIC     = 0
TETRAGONAL    = 1
HEXAGONAL     = 2
TRIGONAL      = 3
ORTHORHOMBIC  = 4
MONOCLINIC    = 5
TRICLINIC     = 6
STRUCTURE_MAX = 7

structure_names = [
  "isometric",
  "tetragonal",
  "hexagonal",
  "trigonal",
  "orthorhombic",
  "monoclinic",
  "triclinic",
  ]

def volume(structure, a, b=None, c=None, alpha=None, beta=None, gamma=None):
  """
  Calculate volume of unit cell
  """
  if not 0 <= structure < STRUCTURE_MAX:
    raise ValueError("Unknown structure\n")

  errmsg = "%s are required to calculate the %s lattice unit-cell volume "
  error = lambda req: ValueError(errmsg % (req, structure_names[structure]))

  if structure == ISOMETRIC:
    return a**3
  elif structure == TETRAGONAL:
    if c is None: raise error("a and c")
    return a*2 * c
  elif (structure == HEXAGONAL or
        structure == TRIGONAL):
    if c is None: raise error("a and c")
    return a**2 * c * sin(pi/3)
  elif structure == ORTHORHOMBIC:
    if c is None: raise error("a, b, and c")
    return a*b*c
  elif structure == MONOCLINIC:
    if c is None: raise error("a, b, c, and beta")
    return a*b*c * sin(beta)
  elif structure == TRICLINIC:
    if c is None: raise error("a, b, c, alpha, beta, and gamma")
    return a*b*c*(1-cos(alpha)**2-cos(b)**2-cos(gamma)**2) + 2*sqrt(cos(a) + cos(beta) + cos(gamma))


