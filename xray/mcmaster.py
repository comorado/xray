"""
X-ray cross sections form McMaster tables
"""

import mucal
import numpy as np

TYPES = ['total', 'coherent', 'incoherent', 'photoelectric']
VALID_TYPES = set(TYPES + ['all'])


def check_xsec_type(type):
  if not type in VALID_TYPES:
    raise ValueError("`type` must be one of 'total', 'coherent', 'incoherent' or 'photoelectric'")

class Element(object):
  def __init__(self, elt):
    """
    Parameters
      elt: chemical symbol for element
    """

    self.elt = elt

    # extract energy independent parameters
    info = mucal.mucal(elt, 1.0)
    self.energies = info['energies']

    xsec = info['xsec']
    self.atwt = xsec['atwt']
    self.density = xsec['density']
    self.conversion = xsec['conversion']

  def xsec(self, energy, type='total', barns=True):
    """
    Cross-section in barns/atom (or cm^2/g if `barns` is False)
   
    Parameters
    ----------
      energy: photon energy in eV
      type: total, coherent, incoherent or photoelectric
      barns: whether to use barns/atom or cm^2/g units
    """

    check_xsec_type(type)

    info = mucal.mucal(self.elt, energy/1000.0)['xsec']

    if type == 'all':
      xsec = np.array([info[t] for t in TYPES])
    else:
      xsec = info[type]

    if not barns:
      xsec /= self.conversion

    return xsec

  def mu(self, energy, type='total'):
    """
    Inverse absorption length in 1/um
   
    Parameters:
      energy: photon energy in eV
      type: total, coherent, incoherent or photoelectric
    """

    check_xsec_type(type)
    info = mucal.mucal(self.elt, energy/1000.0)['xsec']
    conv = self.density / self.conversion / 1e4
   
    if type == 'all':
      mu = np.array([info[t]*conv for t in TYPES])
    else:
      mu = info[type] * conv

    return mu

e = Element()
