import numpy as np
from itertools import cycle
from scipy.interpolate import splrep, splint

"""
This file contains some old cruft that is duplicated elsewhere.
For now the main purpose is to bring in the code to read in FEFF LDOS files. 
"""

def fermi(E,T,mu): 
  """
  fermi distribution

  Parameters:
    E: energy values at which to calculate distribution
    T: temperature
    mu: chemical potential
  """
  if T == 0:
    return (E <= mu) * 1.0

  return 1.0/(np.exp((E-mu)/T) + 1)

def N(mu, E,T,g):
  """
  Find number of electrons by integrating g(E) * f(E,T,mu) over all E

  Parameters:
    mu: chemical potential
    E: energy grid
    T: temperature
    g: DOS
  """

  tck = splrep(E, g*fermi(E,T,mu))
  return splint(E[0], E[-1], tck)

  #return np.trapz( g * fermi(E,T,mu), E )

def N_constraint(mu, E, T, g, N0):
  return abs(N0 - N(mu, E, T, g))

class LDOS(object):
  def __init__(self, filename=None):
    self.mu = 0
    self.charge_transfer = 0
    self.num_atoms = 0
    self.broadening = 0
    self.counts = []
    self.energy = np.array([])
    self.dos = np.array([[]])

    if filename:
      self.load(filename)

  def load(self, filename):
    """
    Loads ldosXX.dat file from FEFF.

    Returns:
      (mu, energy, ldos)
      mu: Fermi Energy
      energy: energy grid
      ldos: array of DOS for each value of angular momentum l
    """

    # read in fermi level
    with open(filename) as f:
      in_counts = False

      # read in header
      for line in f:
        if in_counts:
          pieces = line[1:].strip().split()
          if len(pieces) == 2:
            try:
              int(pieces[0])
              n = float(pieces[1])
              self.counts.append(n)
              continue
            except ValueError:
              pass
          in_counts = False

        if "Fermi level" in line:
          self.mu = float(line.split()[-1])
        elif "Charge transfer" in line:
          self.charge_transfer = float(line.split()[-1])
        elif "Electron counts" in line:
          in_counts = True
        elif "Number of atoms" in line:
          self.num_atoms = int(line.split()[-1])
        elif "Lorentzian broadening" in line:
          self.broadening = float(line.split()[-2])
        elif "-----" in line:
          break #done with header

    # read rest of data
    data = np.loadtxt(filename).T 
    self.energy = data[0]
    self.dos = data[1:]
    self.total_dos = self.dos.sum(0)


  def mu_of_T(self, T, n0=None):
    """
    Find chemical potential as function of temperature
    """
    #from scipy.optimize import fminbound
    from scipy.optimize import fmin

    if n0 is None:
      n0 = sum(self.counts)
    dos = sum(self.dos, 0)

    #mu, dn, ierr, numfunc  = fminbound(N_constraint, self.energy[0], self.energy[-1], (self.energy, T, dos, n0), full_output=True)
    mu, dn, iter, funccalls, warnflag = fmin(N_constraint, self.energy[len(self.energy)/2], (self.energy, T, dos, n0), full_output=True, disp=False)

    mu = mu[0] # why is this needed now?

    #print "n0,n, mu: ", n0, N(mu, self.energy, T, dos), mu
    if warnflag > 0:
      mu = np.nan

    return mu

  def electron_counts(self, T, mu=None):
    """
    Calculate number of electrons for each angular momentum
    """
    if mu is None:
      mu = self.mu_of_T(T)

    # interpolate onto finer grid for integration
    #E = np.linspace(self.energy[0], self.energy[-1], 1001)
    #dos = [np.interp(E, self.energy, d) for d in self.dos]

    return [N(mu, self.energy, T, d) for d in self.dos]

def plot_ldos(ldos, T, mu):
  from matplotlib import pyplot as pt
  pt.plot(ldos.energy, ldos.dos.sum(0))

  colors = cycle("bgrcmy")
  f = fermi(ldos.energy,T, mu)
  pt.plot(ldos.energy, f, color='k', linestyle='dashed')

  dos = ldos.dos.sum(0)
  pt.plot(ldos.energy, dos , color='gray')
  pt.fill_between(ldos.energy, f*dos, color='gray', linestyle='dotted')

  for dos in ldos.dos[::-1]:
    c = colors.next()
    print c
    pt.plot(ldos.energy, dos, color=c)
    pt.fill_between(ldos.energy, f*dos, linestyle='dotted', color=c)

    
  pt.show()


if __name__ == "__main__":
  import sys

  fn = sys.argv[1]
  ldos = LDOS(fn)

  def calc(T):
    mu = ldos.mu_of_T(T)
    return [T, mu] + ldos.electron_counts(T,mu) + [fermi(-43,T,mu)]

  data = [calc(T) for T in np.arange(0,10,0.1)]

  print "# T  mu  Np"
  for line in data:
    print ' '.join(str(s) for s in line)

    #print fminret
    #print "mu: ", mu
    #print "N_check: ", N(mu, E, T, dos)

  #plot_ldos(T,mu)
