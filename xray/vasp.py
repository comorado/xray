import numpy as np
import math
import random
from operator import itemgetter, attrgetter

"""
This file contains tool box for working with XDATCAR file 
"""
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

def readHeader(f):
      """
      Reads the header information form the XDATCAR file
      Leaves file handler at the position past the header
      """
      nextLine = True
      scale   = 1.0
      natoms  = 0
      xvector = 0.0
      yvector = 0.0
      zvector = 0.0
      totalLines = 0
      atomLabel = []
      l = 0
      while (nextLine):
        line = f.readline()
        # Read label
        if totalLines == 1:
            pieces = line.split()
            scale = float(pieces[l])
        if totalLines == 2:
            pieces = line.split()
            xvector = (float(pieces[l]),float(pieces[l+1]),float(pieces[l+2]))
        if totalLines == 3:
            pieces = line.split()
            yvector = (float(pieces[l]),float(pieces[l+1]),float(pieces[l+2]))
        if totalLines == 4:
            pieces = line.split()
            zvector = (float(pieces[l]),float(pieces[l+1]),float(pieces[l+2]))
        if totalLines == 5:
            pieces = line.split()
            for atom in pieces:
              atomLabel.append(atom)
        if totalLines == 6:
            pieces = line.split()
            natoms = float(pieces[l])
        if (line.strip() != '' and totalLines > 5):
            nextLine = False
            line = f.readline()

        totalLines += 1
      return (natoms, xvector, yvector, zvector, atomLabel)

def countFrames(f):
      """
      Helper method counts number of frames in the file
      Starts counting number of frames where the file handler was left
      """
      numberOfAtoms = readHeader(f)[0] 
      nextLine = True
      blankLines = 0
      totalLines = 0
      while(nextLine):
         line = f.readline()
         totalLines += 1
         if line.strip() == '' or line.find('Direct') != -1:
            blankLines += 1
            frameCount = float(totalLines - blankLines) / numberOfAtoms  
            assert (frameCount.is_integer()), "Not consistent number of entries in frame "+ str(blankLines)

         if (line == ""):
            nextLine = False
      f.seek(0) 
      return int(frameCount)

def readFrame(f):
      """
      Helper method, reads next available frame
      """
      l = 0
      atoms = []
      line = f.readline()

      while(line.find('Direct') == -1 and line.strip() != ''):
            if (line.find('Direct') != -1):
                print line.find('Direct') != -1

            pieces = line.split()
            atoms.append(Atom(
              x = float(pieces[l  ]),
              y = float(pieces[l+1]),
              z = float(pieces[l+2]),
              r = 0))
   
            line = f.readline()
      return atoms

def skipFrame(f):
      """
      Skips next available frame
      """
      line = f.readline()
            
      while(line.find('Direct') == -1 and line.strip() != ''):
            line = f.readline()
            if (line.find('Direct') != -1):
                print line.find('Direct') != -1
            line = f.readline()

def readFrameNum(f, frame):
      """
      Returns array of atoms at the frame
      """
      i = 1
      while (i < frame):
         skipFrame(f) 
         i+=1
    
      return readFrame(f) 

def vaspFrameToFeff(f, frame, cutOff, w=None):
      """
      Returns feff style input for the atom position, starting from random atom position
      in the frame
      """

      outputAtoms = []
      numberOfAtoms, xvector, yvector, zvector, atomLabel  = readHeader(f) 
      atoms = readFrameNum(f,frame)

      random.seed()
      randInt = random.randrange(len(atoms))
      centerAtom = atoms[randInt]

      # Assume one kind of atoms in the XDATCAR input
      for nextAtom in atoms:
        for ix in range(-1,2):
           dist = 0.0
           xtran =  xvector[0] * (ix + nextAtom.x - centerAtom.x)
           dist += (xtran)**2
           #if (dist**(0.5) < cutOff):
           for iy in range(-1,2):
             ytran =  yvector[1] * (iy + nextAtom.y - centerAtom.y)
          
             #if ((dist+(ytran)**2)**(0.5) < cutOff):
             dist += (ytran)**2
             for iz in range(-1,2):
               ztran =  zvector[2] * (iz + nextAtom.z - centerAtom.z) 
               #print "Hello",ix,iy,iz 
               dist = xtran**2 + ytran**2 + ztran**2
               if dist**(0.5) < cutOff:
                   #dist += ztran**2
                   outputAtoms.append(Atom(
                     x = xtran,
                     y = ytran,
                     z = ztran,
                     pot = 1,
                     tag = atomLabel[0],
                     r = dist**(0.5),
                     n = 0))
      outputAtoms.sort(key=attrgetter('r'))
      i = 0

      if w is None:
        print "* ",'{0}'.format(frame)," frame used to generate POT input"
      else:
         print >> w,"* ",'{0}'.format(frame)," frame used to generate POT input"

      for atom in outputAtoms:
         if (i > 13):
            pot = 13
         else:
            pot = i

         if w is None:
           print ' {0: .5f} {1: .5f} {2: .5f} {3:3d} {4} {5: .5f} {6}'.format(atom.x, atom.y, atom.z, pot, atom.tag, atom.r, i)
         else:
           print >> w, ' {0: .5f} {1: .5f} {2: .5f} {3:3d} {4} {5: .5f} {6}'.format(atom.x, atom.y, atom.z, pot, atom.tag, atom.r, i)
         i += 1
      if (w is None):
        print 'END'
      else:
        print >> w, 'END'


def comV(f, frameMax = None):
    """
    Calculates center of mass velocity by going through frames
    """
    if frameMax == None:
       nframe = countFrames(f)
    else:
       nframe = frameMax
    
    numberOfAtoms, xvector, yvector, zvector  = readHeader(f) 
    #print numberOfAtoms, xvector, yvector, zvector 
    previousFrame = readFrame(f)	

    #print len(previousFrame)
    for i in xrange(2, nframe-1):
        currentFrame = readFrame(f)
        #print i, "\t".join(map(str,comVFrame(f, currentFrame, previousFrame,numberOfAtoms, xvector, yvector, zvector)))
        xv, yv, zv = comVFrame(f, currentFrame, previousFrame,numberOfAtoms, xvector, yvector, zvector)
        print " %04d\t%10.5f\t%10.5f\t%10.5f\t%10.5f" % (i, xv, yv, zv,(xv**2+yv**2+zv**2)**(0.5))
        previousFrame = currentFrame

def comVFrame(f, currentFrame, previousFrame, numberOfAtoms, xvector, yvector, zvector):
   """
   Helper method for calculating velocity of COM between two adjecent frames
   Checks if there is a translation by the size of unit cell
   """
   comVelocity = []
   dvx = 0
   dvy = 0
   dvz = 0
   tol = 0.1
   for i in xrange(0, int(numberOfAtoms)):
        dx = (currentFrame[i].x-previousFrame[i].x) * xvector[0]
        dy = (currentFrame[i].y-previousFrame[i].y) * yvector[1]
        dz = (currentFrame[i].z-previousFrame[i].z) * zvector[2]

        if abs(dx) + tol >= xvector[0]:
           dx -= math.copysign(xvector[0],dx)
        if abs(dy) + tol >= yvector[1]:
           dy -= math.copysign(yvector[1],dy)
        if abs(dz) + tol >= zvector[2]:
           dz -= math.copysign(zvector[2],dz)
        
        dvx += dx
        dvy += dy
        dvz += dz

   return (dvx/numberOfAtoms,dvy/numberOfAtoms,dvz/numberOfAtoms)

def pairDistance(atoms, numberOfAtoms, xvector, yvector, zvector, cutoff=None, xyz=None):
    """
    Helper method for calculating pair distance in provided array atoms coordinates
    Uses cutoff distance to optimize calculation, if not specified sets to 4 angs 
    """
    #print cutoff
    if cutoff == None:
       cutoff = min([xvector[0],yvector[1],zvector[2]])
       if cutoff > 5:
          cutoff = 4

    dist = []
    m = 0
    d = 0
    if xyz==None:
      for atom1 in atoms:
         for atom2 in atoms:
            for ix in range(-1,2):
              dx = (atom1.x - (ix * xvector[0] + atom2.x))**2
              if (dx**(0.5) < cutoff):
                for iy in range(-1,2):
                  dy = dx + (atom1.y - (iy * yvector[1] + atom2.y))**2
                  if (dy**(0.5) < cutoff):
                    for iz in range(-1,2):
                      d =(dy +(atom1.z - (iz * zvector[2] + atom2.z))**2)**(0.5)
  
                      if d != 0 and d < cutoff:
                        dist.append(d)
                        m += 1

    else:
       for atom1 in atoms:
         for atom2 in atoms:
                  d = 0.0
                  d += (atom1.x - atom2.x)**2
                  d += (atom1.y - atom2.y)**2
                  d += (atom1.z - atom2.z)**2

                  d = d **(0.5)
  
                  if d != 0 and d < cutoff:
                    dist.append(d)
                    m += 1
       
    return dist
  
