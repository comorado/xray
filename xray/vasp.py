import numpy as np
import math

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
      l = 0
      while (nextLine):
        line = f.readline()
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
        if totalLines == 6:
            pieces = line.split()
            natoms = float(pieces[l])
        if (line.strip() != '' and totalLines > 5):
            nextLine = False
            line = f.readline()

        totalLines += 1
      return (natoms, xvector, yvector, zvector)

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
  
