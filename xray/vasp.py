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
      """
      #f = open(filename,"r")
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
        totalLines += 1
        if (totalLines > 6):
            nextLine = False
      # Go to the top of the file
      f.seek(0)
      return (natoms, xvector, yvector, zvector)

def countFrames(f):
      """
      Helper method counts number of frames in the file
      """
      numberOfAtoms = readHeader(f)[0] 
      nextLine = True
      blankLines = 0
      totalLines = 0
      headerLines = 8
      while(nextLine):
         line = f.readline()
         totalLines += 1
         if line.strip() == '' and totalLines > headerLines:
            blankLines += 1
            frameCount = float(totalLines - headerLines - blankLines) / numberOfAtoms  
            assert (frameCount.is_integer()), "Not consistent number of entries in frame "+ str(blankLines)

         if (line == ""):
            nextLine = False
      f.seek(0) 
      return int(frameCount)

def readFrame(f, frameNumber = 1, fileName = None):
      """
      Helper method, reads specified frame from the file handler,
      returns array of objects of atoms
      """
      nextLine = True
      blankLines = 0
      totalLines = 0
      headerLines = 8
      l = 0
      atoms = []
      if (fileName == None):
         numberOfAtoms, xvector, yvector, zvector  = readHeader(f) 
   
         while(nextLine):
            line = f.readline()
            totalLines += 1
   
            if (line.strip() == '' and totalLines > headerLines):
               blankLines += 1
   
               while (blankLines == frameNumber):
                  line = f.readline()
   
                  if line.strip() != '':
                     pieces = line.split()
                     atoms.append(Atom(
                          x = float(pieces[l  ]) * xvector[0],
                          y = float(pieces[l+1]) * yvector[1],
                          z = float(pieces[l+2]) * zvector[2],
                          r = 0))
                  else:
                     blankLines += 1        
                     nextLine = False
   
         f.seek(0) 
         assert (len(atoms) == numberOfAtoms), "Number of atoms is wrong " + str(len(atoms))
         assert (len(atoms) != 0), "Frame number is out of bounds " + str(len(atoms))
         
         return atoms

      else: 
         fileHandler = open(fileName,'r') 
         numberOfAtoms, xvector, yvector, zvector  = readHeader(fileHandler) 
         while(nextLine):
            line = f.readline()
            totalLines += 1
   
            if (line.strip() == ''):
               blankLines += 1
   
               while (blankLines == 1):
                  line = f.readline()
   
                  if (line.strip() != ''):
                     pieces = line.split()
                     atoms.append(Atom(
                          x = float(pieces[l  ]) * xvector[0],
                          y = float(pieces[l+1]) * yvector[1],
                          z = float(pieces[l+2]) * zvector[2],
                          r = 0))
                  else:
                     blankLines += 1        
   
         assert (len(atoms) == numberOfAtoms), "Number of atoms is wrong " + str(len(atoms))
         assert (len(atoms) != 0), "Frame number is out of bounds " + str(len(atoms))
         
         return atoms

def comV(f, frameMax = None):
    """
    Calculates center of mass velocity by going through frames
    """
    previousFrame = readFrame(f)	
    if frameMax == None:
       nframe = countFrames(f)
    else:
       nframe = frameMax
    
    for i in xrange(2, nframe-1):
        currentFrame = readFrame(f, i)
        print i, "\t".join(map(str,comVFrame(f, currentFrame, previousFrame)))
        previousFrame = currentFrame

def comVFrame(f, currentFrame, previousFrame):
   """
   Helper method for calculating velocity of COM between two adjecent frames
   Checks if there is a translation by the size of unit cell
   """
   comVelocity = []
   numberOfAtoms, xvector, yvector, zvector  = readHeader(f) 
   #print readHeader(f)
   dvx = 0
   dvy = 0
   dvz = 0
   for i in xrange(0, int(numberOfAtoms)):
        dx = currentFrame[i].x-previousFrame[i].x
        dy = currentFrame[i].y-previousFrame[i].y
        dz = currentFrame[i].z-previousFrame[i].z

        if abs(dx)+0.1 >= xvector[0]:
           dx -= math.copysign(xvector[0],dx)
        if abs(dy)+0.1 >= yvector[1]:
           dy -= math.copysign(yvector[1],dy)
        if abs(dz)+0.3 >= zvector[2]:
           dz -= math.copysign(zvector[2],dz)
        
        dvx += dx
        dvy += dy
        dvz += dz

   return (dvx/numberOfAtoms,dvy/numberOfAtoms,dvz/numberOfAtoms)
   
