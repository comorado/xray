"""
NIST Database of x-ray absorption edges and emission lines from 0 to 10keV

  >>> import nistdb
  >>> fe = nistdb.elements()['Fe']
  >>> fe.edges['K']
  AbsorptionEdge('Fe', 'K', 7110.747000, 0.020000, TYPE_DIRECT)
  >>> print '\n'.join([str(l) for l in nistdb.lines() if 5020 < l.energy < 5040])
  Pr L3M5 5033.79 eV (La1)
  Xe L2N4 5036.95 eV (Lg1)
  Xe L2N5 5038.78 eV

  >>> nistdb.closest_line(7060)
  [EmissionLine('Fe', 'KM2', 7058.175000, 0.016000, TYPE_DIRECT)]
"""
import os

DB_FILE = os.path.join(os.path.dirname(__file__), 'data', 'NIST_xray_all_elements_0_to_30keV.txt')

"""
TYPE_* indicates the source for the emission line / absorption edge data
"""
TYPE_UNKNOWN = 0
TYPE_THEORY = 1
TYPE_DIRECT = 2
TYPE_VAPOR = 3
TYPE_COMBINED = 4

type_reprs = [
  "TYPE_UNKNOWN",
  "TYPE_THEORY", "TYPE_DIRECT",
  "TYPE_VAPOR",
  "TYPE_COMBINED",
  ]
type_names = [
  "Unknown",
  "Theory",
  "Direct",
  "Vapor",
  "Combined",
  ]

# partial mapping between siegbahn and IUPAC notation
siegbahn_map = {
    'Ka1': 'KL3',
    'Ka2': 'KL2',
    'Kb1': 'KM3',
    'Kb2': 'KN3',
    'Kb3': 'KM2',
    'La1': 'L3M5',
    'La2': 'L3M4',
    'Lb1': 'L2M4',
    'Lb2': 'L3N5',
    'Lb3': 'L1M3',
    'Lg1': 'L2N4',
    'Lg2': 'L1N2',
    'Lg3': 'L1N3',
    'Ma1': 'M5N7',
    'Ma2': 'M5N6',
    'Mb': 'M4N6',
    'Mg': 'M3N5',
    }

# inverse of siegbahn_map
iupac_map = {siegbahn_map[k]:k for k in siegbahn_map}

class EmissionLine(object):
  """
  X-ray emission line data

  Properties:
    elt_name:
      name of element (e.g. 'Au')
    trans:
      transition name in IUPAC notation (e.g. 'KL3')
    energy:
      energy of emission line (eV)
    uncertainty:
      uncertainty in energy value
    type:
      source for energy value
  """
  def __init__(self, elt_name, trans, energy, uncertainty=0.0, type=TYPE_THEORY):
    self.elt_name = elt_name
    self.trans = trans
    self.energy = energy
    self.uncertainty = uncertainty
    self.type = type

  def __repr__(self):
    return "%s('%s', '%s', %f, %f, %s)" % (self.__class__.__name__,
               self.elt_name, self.trans, self.energy, self.uncertainty, type_reprs[self.type])

  def __str__(self, include_elt_name=True):

    s = "%s %.2f eV" % (self.trans, self.energy)
    if include_elt_name:
      s = "%s %s" % (self.elt_name, s)

    if self.trans in iupac_map:
      s += " (%s)" % iupac_map[self.trans]
    return s

class AbsorptionEdge(object):
  """
  X-ray absorption edge data

  Properties:
    elt_name:
      name of element (e.g. 'Au')
    edge:
      name of edge (e.g. 'L3')
    energy:
      energy of absorption edge (eV)
    uncertainty:
      uncertainty in energy value
    type:
      source for energy value
  """
  def __init__(self, elt_name, edge, energy, uncertainty=0.0, type=TYPE_THEORY):
    self.elt_name = elt_name
    self.edge = edge
    self.energy = energy
    self.uncertainty = uncertainty
    self.type = type

  def __repr__(self):
    return "%s('%s', '%s', %f, %f, %s)" % (self.__class__.__name__,
               self.elt_name, self.edge, self.energy, self.uncertainty, type_reprs[self.type])

  def __str__(self, include_elt_name=True):
    s = "%s %.2f" % (self.edge, self.energy)
    if include_elt_name:
      s = "%s %s" % (self.elt_name, s)
    return s

class Element(object):
  """
  An element

  Properties:
    name:
      Chemical symbol for element (e.g. 'Ag')
    lines:
      dictionary of EmissionLine objects keyed by IUPAC transition name
    edges:
      dictionary of AbsorptionEdge objects keyed by IUPAC shell name
  """

  def __init__(self, name, lines={}, edges={}):
    from collections import OrderedDict
    self.name = name
    self.lines = OrderedDict()
    self.edges = OrderedDict()

  def siegbahn(self):
    lines = {}
    for sym in siegbahn_map:
      line = self.lines.get(siegbahn_map[sym])
      if line: lines[sym] = line

    return lines

  def __repr__(self):
    return "%s('%s')" % (self.__class__.__name__, self.name)

def load_nist_database(filename):
  f = open(filename)

  by_elt = {}
  lines = []
  edges = []

  for fline in f:
    # skip comments and empty rows
    if not fline.strip() or fline[0] == '#': continue

    elt_name, A, trans, E_th, sigma_th, E_exp, sigma_exp, E_combined, sigma_combined,  E_vapor, sigma_vapor, blend, ref  = fline.split('\t')

    elt = by_elt.get(elt_name)
    if not elt:
      elt = Element(elt_name)
      by_elt[elt_name] = elt

    # prefer direct experimental measurement, then vapor, then theory
    if E_exp:
      E = float(E_exp)
      sigma = float(sigma_exp)
      type = TYPE_DIRECT
    elif E_vapor:
      E = float(E_vapor)
      sigma = float(sigma_vapor)
      type = TYPE_VAPOR
    elif E_th:
      E = float(E_th)
      sigma = float(sigma_th)
      type = TYPE_THEORY

    if 'edge' in trans:
      edge_name = trans[:-5]
      edge = AbsorptionEdge(elt_name, edge_name, E, sigma, type)
      elt.edges[edge_name] = edge
      edges.append(edge)
    else:
      line = EmissionLine(elt_name, trans, E, sigma, type)
      elt.lines[trans] = line
      lines.append(line)

  return by_elt, lines, edges

class Database(object):
  _singleton = None
  def __init__(self, db_file=DB_FILE):
    self.db_file = db_file
    self.by_elt, self.lines, self.edges = load_nist_database(db_file)

  @classmethod
  def instance(cls):
    if cls._singleton is None:
      cls._singleton = cls()
    return cls._singleton

"""
Simplified API below
"""
SHARED_DB = None
def db():
  global SHARED_DB
  if not SHARED_DB:
    SHARED_DB = load_nist_database(DB_FILE)
  return SHARED_DB

def elements():
  return db()[0]

def lines():
  return db()[1]

def edges():
  return db()[2]

def _energies(obj_list):
  return [o.energy for o in obj_list]

def _closest_in_energy(obj_list, energy, before=0, after=0):
  import numpy as np
  energies = np.array(_energies(obj_list))
  closest_index = np.argmin(np.abs(energies - energy))
  lower = max(closest_index - before, 0)
  upper = min(closest_index + after + 1, len(energies))
  return obj_list[lower:upper]

def _in_energy_range(obj_list, energy_min, energy_max):
  return [o for o in obj_list if energy_min < o.energy < energy_max ]

def line_energies():
  return _energies(lines())

def edge_energies():
  return _energies(edges())

def closest_lines(energy, before=0, after=0):
  return _closest_in_energy(lines(), energy, before, after)

def closest_line(energy):
  try:
    return closest_lines(energy)[0]
  except IndexError:
    return None

def lines_in_range(energy_min, energy_max):
  """
  Return a list of emission lines in the supplied energy range

  Parameters:
    energy_min: minimum energy bound
    energy_max: maximum energy bound
  """
  return _in_energy_range(lines(), energy_min, energy_max)

def closest_edges(energy, before=0, after=0):
  return _closest_in_energy(edges(), energy, before, after)

def closest_edge(energy):
  try:
    return closest_edges(energy)[0]
  except IndexError:
    return None

def edges_in_range(energy_min, energy_max):
  """
  Return a list of absortion edges in the supplied energy range

  Parameters:
    energy_min: minimum energy bound
    energy_max: maximum energy bound
  """
  return _in_energy_range(edges(), energy_min, energy_max)

