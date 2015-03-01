import numpy as np

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

def readFrame(f, frameNumber):
      nextLine = True
      blankLines = 0
      totalLines = 0
      headerLines = 8
      while(nextLine):
         line = f.readline()
         totalLines += 1
         if line.strip() == '' and totalLines > headerLines and blankLines == frameNumber:
            blankLines += 1
            #while ():

            atoms.append(AtomXYZ(
                 x = float(pieces[l]) * self.xvector[0],
                 y = float(pieces[l+1]) * self.yvector[1],
                 z = float(pieces[l+2]) * self.zvector[2],
                 r = 0))

         if (line == ""):
            nextLine = False

      f.seek(0) 
      return atoms

