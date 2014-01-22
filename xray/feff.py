import numpy as np
from itertools import cycle
from scipy.interpolate import splrep, splint
from numpy import fft
#from .lib import fortranfile
import fortranfile, random
from operator import itemgetter, attrgetter
from . import const

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
  if T == 0.0:
    return (E <= mu) * 1.0

  return 1.0/(np.exp((E-mu)/T) + 1.0)

def N(mu, E, T, g, n=None):
  """
  Find number of electrons by integrating g(E) * f(E,T,mu) * E**n over all E

  Parameters:
    mu: chemical potential
    E: energy grid
    T: temperature
    g: DOS
    n: momentum
  """
  if n is None:
   tck = splrep(E, g*fermi(E,T,mu))
   return splint(E[0], E[-1], tck)
  else:
   tck = splrep(E, g*fermi(E,T,mu)*E**n)
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

  def __call__(self, E, l=None):
    y = self.total_dos if l is None else self.dos[l]
    return np.interp(E, self.energy, y, left=0, right=0)

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

  def _find_occupation(self, Emin, N):
    """
    Finds Emax such that the integral of total_dos from Emin to Emax equals N
    """
    from .utils import trapezoid
    from scipy.optimize import fminbound


    constraint = lambda Emax: (trapezoid(self.energy,
                                               self.total_dos,
                                               Emin, Emax
                                              ) - N)**2

    Emax, dn, warnflag, funccalls = fminbound(constraint,
                                              Emin,
                                              self.energy[-1],
                                              full_output=True,
                                              disp=0)


    #print dn, warnflag, funccalls
    if warnflag != 0 or dn > 1e-5:
      return None

    return Emax

  def electron_counts(self, T, mu=None, n=None):
    """
    Calculate number of electrons for each angular momentum
    n : momentum to calculate e.g. E**n*fermi(T)*g(E)
    """
    if mu is None:
      mu = self.mu_of_T(T)

    # interpolate onto finer grid for integration
    #E = np.linspace(self.energy[0], self.energy[-1], 1001)
    #dos = [np.interp(E, self.energy, d) for d in self.dos]

    #print "Print T = ",T
    #print "Print mu = ",mu

    if n is None:
      return [N(mu, self.energy, T, d) for d in self.dos]
    else:
      return [N(mu, self.energy, T, d, n) for d in self.dos]


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

class DensityBin(object):
  def __init__(self, filename):
    f = fortranfile.FortranFile(filename)
    self.ndims = f.readInts()[0]
    self.origin = f.readReals('d')

    gridspec = [(f.readReals('d'), f.readInts()[0]) for i in range(self.ndims)]
    self.axes = [axis for axis,n in gridspec]
    self.npts = [n for axis,n in gridspec]
    self.steps = [axis/(n-1) for axis,n in gridspec]


    self.rho = f.readReals('d').reshape(self.npts[::-1])


  def point_at_index(self, idx):
    assert(idx.ndims == self.ndims)
    self.origin + sum([i*s for i,s in zip(idx,self.steps)], 0)

  def idx_iter(self):
    idx_iter(self.npts)

  def is_orthog(self):
    comparisons = [(0,1), (0,2), (1,2)]
    return all([np.dot(self.axes[i], self.axes[j]) < 1e-10
                for i,j in comparisons])

  def orthog_coords(self):
    # XXX this currently assumes that axes are ordered x-y-z...
    assert(self.is_orthog())
    return [x0 + np.linalg.norm(ax) / (npts-1) * np.arange(npts)
            for x0,ax,npts in zip(self.origin,self.axes, self.npts)]

def idx_iter(shape):
  from itertools import product
  return product(*[xrange(n) for n in shape])

class DensityMap(object):
  def __init__(self, filename, float_precision='f'):
    f = fortranfile.FortranFile(filename)
    self.energy = f.readReals('d')
    self.x = f.readReals('d')
    self.y = f.readReals('d')
    self.density = f.readReals(float_precision).reshape(len(self.x), len(self.y), len(self.energy))

  def occupied(self, E_fermi, T=0):
    """
    Integrate over Fermi distribution to obtain occupied density

    For T != 0, `E_fermi` should be the chemical potential!
    """
    from .fermi import f
    from .utils import trapezoid

    if T > 0:
      return np.trapz(self.density[:,:,:] * f(self.energy, T, E_fermi), self.energy)
    else:
      return trapezoid(self.energy, self.density, x1=E_fermi, axis=2)

  def occupied2(self, Emin, Emax):
    """
    Integrate over energies in [Emin,Emax] to obtain density
    """
    pass


def load_atomic_wf(filename):
  """
  Load atomic wavefunction output from FEFF

  returns list of (x,y) tuples (one for each L value)
  """
  wf = np.loadtxt(filename).T
  indices = [wf[0] == i for i in np.unique(wf[0])]
  return [(wf[2][i], wf[3][i] + 1j * wf[4][i]) for i in indices]


#
# FEFF Input code below
#

def rot_axis_angle(a, b):
  """
  Find rotation axis and angle that rotates vector a to vector b
  """

  c = np.cross(a,b)
  norm = np.linalg.norm
  theta = np.arcsin(norm(c) / (norm(a) * norm(b)))

  return (c,theta)

def rotation_matrix(axis, theta):
  R = np.zeros((3,3))

  u = axis / np.linalg.norm(axis)
  c = np.cos(theta)
  d = 1- c
  s = np.sin(theta)

  R[0,0] = c + u[0]**2 * d
  R[0,1] = u[0] * u[1] * d - u[2] * s
  R[0,2] = u[0] * u[2] * d + u[1] * s
  R[1,0] = u[1] * u[0] * d + u[2] * s
  R[1,1] = c + u[1]**2 * d
  R[1,2] = u[1] * u[2] * d - u[0] * s
  R[2,0] = u[2] * u[0] * d - u[1] * s
  R[2,1] = u[2] * u[1] * d + u[0] * s
  R[2,2] = c + u[2]**2 * d

  return R

class Lattice(object):
  def __init__(self, a, b, c, alpha=0, beta=0, gamma=0, x=1, y=1, z=1, n=10, tag="", filename=""):
    self.a = a
    self.b = b
    self.c = c
    self.x = x
    self.y = y
    self.z = z
    self.alpha = alpha*np.pi/180.0
    self.beta = beta*np.pi/180.0
    self.gamma = gamma*np.pi/180.0
    self.n = n;
    self.tag = tag;
    self.filename = filename;

  def partial_coordinates_transform_matrix(self):
    R = np.zeros((3,3))
    sgamma = np.sin(self.gamma)
    cgamma = np.cos(self.gamma)
    sbeta  =  np.sin(self.beta)
    cbeta  =  np.cos(self.beta)
    salpha =  np.sin(self.alpha)
    calpha =  np.cos(self.alpha)

    V = np.abs(1-calpha**2-cbeta**2-cgamma**2+2*calpha*cbeta*cgamma)**(0.5)
 
    R[0,0] = self.a
    R[0,1] = self.b * cgamma
    R[0,2] = self.c * cbeta
    R[1,0] = 0.0
    R[1,1] = self.b * sgamma
    R[1,2] = self.c * (calpha-cbeta*cgamma)/sgamma
    R[2,0] = 0.0
    R[2,1] = 0.0
    R[2,2] = self.c*V/sgamma

    return R

  def lattice(self,structure='hcp', shake=1,d=0.01): 
    self.atoms=[]
    R = self.partial_coordinates_transform_matrix()
    dx,dy,dz = np.dot(R, np.array([self.x,self.y,self.z]))
    x1,y1,z1=(0.0,0.0,0.0)
    for i in range(-self.n,self.n+1):
     for j in range(-self.n,self.n+1):
      for k in range(-self.n,self.n+1):
       if np.abs(i)+ np.abs(j)+ np.abs(k) <= self.n:
        xl,yl,zl = np.dot(R, np.array([i,j,k]))
        if shake==1:
          x1,y1,z1=self.random(d)

        self.atoms.append(Atom(
      	  x = xl+x1,
      	  y = yl+y1,
          z = zl+z1,
      	  pot = 1,
      	  tag = self.tag,
      	  r = ((xl+x1)**2+(yl+y1)**2+(zl+z1)**2)**(0.5),
      	  n = np.abs(i) + np.abs(j) + np.abs(k)))
        #if structure=='hcp':  
        if shake==1:
          x1,y1,z1=self.random(d)
        self.atoms.append(Atom(
          x = xl+dx+x1,
          y = yl+dy+y1,
          z = zl+dz+z1,
          pot = 1,
          tag = self.tag,
          r =((xl+dx+x1)**2 +
            (yl+dy+y1)**2 +
            (zl+dz+z1)**2)**(0.5),
          n = np.abs(i) + np.abs(j) + np.abs(k)))

    self.atoms.sort(key=attrgetter('r'))
    f = open(self.filename + '.xyz', 'w')
    g = open(self.filename + '.inp', 'w')
    i = 0
    for a in self.atoms:
      a.n = i
      g.write(str(a)+"\n")
      f.write(a.xyz()+"\n")
      i+=1
  def random(self,d):
    return (self.x*random.uniform(-d,d),self.y*random.uniform(-d,d),self.z*random.uniform(-d,d))
    

  def volume(self,V=0):
    sgamma = np.sin(self.gamma)
    cgamma = np.cos(self.gamma)
    sbeta  =  np.sin(self.beta)
    cbeta  =  np.cos(self.beta)
    salpha =  np.sin(self.alpha)
    calpha =  np.cos(self.alpha)
    self.V = self.a*self.b*self.c*sgamma*self.c*(1-cbeta*(((calpha/(cbeta)-cgamma)*1.0/sgamma)**2+1.0)**(0.5))**(.5)#/(const.BOHR)**3
    V=self.V
    return self.V

class Atom(object):
  def __init__(self, x, y, z, pot=0, tag='', r=0.0, n=0):
    self.x = x
    self.y = y
    self.z = z
    self.pot = pot
    self.tag = tag
    self.r = r
    self.n = n

  def __str__(self):
    return " %10.5f %10.5f %10.5f  %3d  %10s %10.5f % 5d" % (self.x,self.y,self.z,self.pot,self.tag,self.r,self.n)

  def xyz(self):
    return " %2s %10.5f %10.5f  %10.5f" % (self.tag,self.x,self.y,self.z)

  def __repr__(self):
    return "xray.feff.Atom(%10.5f, %10.5f, %10.5f)" % (self.x, self.y, self.z)

class InputFile(object):
  def __init__(self, filename):
    if filename:
      self.load(filename)

  def load(self, filename):
    self.filename = filename

    self.pre_atoms = []
    self.post_atoms = []
    self.atoms = []
    self.cur = self.pre_atoms

    with open(filename) as f:
      in_atoms = False

      for line in f:
        if not in_atoms:
          self.cur.append(line)
          if line.strip().startswith('ATOMS'):
            in_atoms = True

        else:
          sline = line.strip()
          if sline.startswith('END'):
            self.cur = self.post_atoms
            self.cur.append(line)
            in_atoms = False
          elif sline[0] == '*':
            self.cur.append(line)
          else:
            pieces = line.split()
            self.atoms.append(Atom(
                x = float(pieces[0]),
                y = float(pieces[1]),
                z = float(pieces[2]),
                pot = int(pieces[3]),
                tag = pieces[4],
                r = float(pieces[5]),
                n = int(pieces[6])
                ))

  def scale_atoms_uniformly(self, scale):
    self.scale_atoms(scale, scale, scale)

  def scale_atoms(self, scale_x, scale_y, scale_z):
    for a in self.atoms:
      a.x *= scale_x
      a.y *= scale_y
      a.z *= scale_z
      a.r = np.sqrt(a.x**2 + a.y**2 + a.z**2)

  def rotate_atoms_to(self, zhat, axis="z"):
    """
    Rotates atoms so that z axis now points along `zhat`
    """
    if axis =="y":
       vec = (np.array([0.,1.,0.])) 
    elif axis == "x":
       vec = (np.array([1.,0.,0.])) 
    else:
       vec = (np.array([0.,0.,1.])) 

    ax,ang = rot_axis_angle(vec, zhat)
    R = rotation_matrix(ax, -ang)
    self.rotate_atoms(R)

  def rotate_atoms(self, rotation_matrix):
    for a in self.atoms:
      a.x,a.y,a.z = np.dot(rotation_matrix, np.array([a.x,a.y,a.z]))

  def truncate_cylinder(self, radius):
    # filter out atoms further from z axis than radius
    self.atoms = [a for a in self.atoms if a.x**2 + a.y**2 <= radius**2]
    # renumber atom list
    for i,a in enumerate(self.atoms):
      a.n = i

  def save(self, filename):
    with to_filehandle(filename, "w") as (f, fname):
      for line in self.pre_atoms:
        f.write(line)

      for a in self.atoms:
        f.write(str(a) + "\n")

      for line in self.post_atoms:
        f.write(line)

def vasp_atoms(feff_input, bounds):
  xmin,xmax = bounds[0]
  ymin,ymax = bounds[1]
  zmin,zmax = bounds[2]

  dx = xmax - xmin
  dy = ymax - ymin
  dz = zmax - zmin

  # limit to atoms within bounds
  atoms_in_box = [a for a in feff_input.atoms if
                    xmin <= a.x <= xmax and
                    ymin <= a.y <= ymax and
                    zmin <= a.z <= zmax]

  # scale to box coordinates
  vatoms = []
  for a in atoms_in_box:
    x = (a.x - xmin) / dx
    y = (a.y - ymin) / dy
    z = (a.z - zmin) / dz
    vatoms.append((x,y,z))


  print "Comment"
  print "    1.000"
  print "    %10.5f %10.5f %10.5f" % (dx,0,0)
  print "    %10.5f %10.5f %10.5f" % (0,dy,0)
  print "    %10.5f %10.5f %10.5f" % (0,0,dz)
  print "Selective dynamics"
  print "Direct"
  for va in vatoms:
    print "    %10.5f %10.5f %10.5f F F F" % va

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


# XXX the following two functions should go elsewhere...
def writable(fname):
  return hasattr(fname, 'write')

from contextlib import contextmanager
@contextmanager
def to_filehandle(fname_or_fh, flag="r"):
  """
  Either open a new file or use existing filehandle for IO
  """
  try:
    if writable(fname_or_fh):
      fh = fname_or_fh
      own_fh = False
      fname_or_fh = None
    else:
      fh = open(fname_or_fh, flag)
      own_fh = True
    yield fh, fname_or_fh
  finally:
    if own_fh:
      fh.close()

def load_pot(fname):
  info = {}
  with open(fname) as f:

    # read in first row which contains dimensions
    vars = f.readline().strip().split()
    vals = f.readline().split()
    for k,v in zip(vars,vals):
      info[k.strip(',')] = int(v)

    info['title'] = []
    for i in range(info['ntitle']):
      info['title'].append(f.readline().strip())

    for i in range(13):
      key, val = f.readline().strip().split()
      info[key] = float(val)

    # read in everything else
    key = None
    for line in f:
      line = line.strip()
      if not (line[0].isdigit() or line[0] == '-'):
        if key is not None:
          info[key] = np.array(info[key])
        key = line
        info[key] = []
      else:
        info[key] += [float(s) for s in line.split()]
    info[key] = np.array(info[key])

    return info

def luecks_grid(x0=-8.8, dx=0.05, n=251):
    return np.exp(x0 + dx * np.arange(n))

class PotDat(object):
  def __init__(self, fname, nr=251):
    pot = load_pot(fname)

    # convert from dict to attributes
    for k in pot:
      self.__setattr__(k, pot[k])

    self.lx = 30

    # fix up dimensions
    self.dgc = self.dgc.reshape(self.nph+1, -1, nr)
    self.dpc = self.dpc.reshape(self.nph+1, -1, nr)
    self.adgc = self.adgc.reshape(self.nph+1, -1, 10)
    self.adpc = self.adpc.reshape(self.nph+1, -1, 10)
    self.edens = self.edens.reshape(self.nph+1, nr)
    self.vclap = self.vclap.reshape(self.nph+1, nr)
    self.vtot = self.vtot.reshape(self.nph+1, nr)
    self.edenvl = self.edenvl.reshape(self.nph+1, nr)
    self.vvalgs = self.vvalgs.reshape(self.nph+1, nr)
    self.dmag = self.dmag.reshape(self.nph+1, nr)
    self.xnval = self.xnval.reshape(self.nph+1, -1)
    self.iorb = self.iorb.reshape(self.nph+1, 8) # XXX check this
    self.xnmues = self.xnmues.reshape(self.nph+1, -1)

    self.lx = self.xnmues.shape[1]
    self.r = luecks_grid()

  def density(self, ipot):
    """
    Calculate *atomic* density for unique potential ipot

    Returns:
      radial density for each shell
    """

    if ipot > self.nph:
      raise ValueError("ipot must be less than self.nph")

    return [d/(4*np.pi*self.r**2) for d in self.dgc[ipot]**2 + self.dpc[ipot]**2 if not (d == 0.0).all()]

  def atomic_formfactors(self, ipot):
    x = np.linspace(0,self.r[-1], self.r[-1]*101)
    ys = (np.interp(x, self.r, d) for d in self.density(ipot))
    tmp = [sph_fft(x, y) for y in ys]

    q = tmp[0][0]
    fs = [t[1] for t in tmp]
    return  q, fs

    """
    q = fft.fftshift(fft.fftfreq(len(x), dx)) * 2 * np.pi
    fs = [-np.imag(fft.fftshift(fft.fft(np.interp(x, self.r, d/self.r)*dx)))/
          for d in self.density(ipot)]
    return q, fs
    """
  

def sph_fft(x, y):
  dx = x[1] - x[0]
  q = fft.fftshift(fft.fftfreq(len(x), dx)) * 2 * np.pi
  f = -(4*np.pi)*fft.fftshift(np.imag(fft.fft(y*dx * x)))
  f /= q
  f[q == 0] = np.trapz(4*np.pi*y*x**2, x)
  return q, f

def load_binmatrix(fname):
  f = fortranfile.FortranFile(fname)
  dims = f.readInts()
  shape = list(dims)[::-1] + [2]
  data = f.readReals('d').reshape(shape)
  return data[...,0] + 1j * data[...,1]
