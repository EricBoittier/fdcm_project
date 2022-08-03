#!/usr/bin/env python3

# attempt to define local DCM frames automatically to convert a charge model in the
# global axis to a charge model with local axis definitions used in CHARMM's DCM
# module

import sys
import math
import subprocess

def usage():
  print ("Usage: python3 get_frames.py <pdbfile> <resname>")

if len(sys.argv) < 3:
  usage()
  exit()

# vdw radii taken from https://physlab.lums.edu.pk/images/f/f6/Franck_ref2.pdf
# H & Noble gases taken from Table 12 of https://pubs.acs.org/doi/10.1021/jp8111556

pdbfile = sys.argv[1]
resname = sys.argv[2]
types=[]
x=[]
y=[]
z=[]

# read coordinates from PDB file
try:
  with open(pdbfile, 'r') as fin:
    for line in fin:
      a=line.split()
      if len(a) > 4:
        if a[0].lower() == "atom":
          tt = line[12:16]
          tt2 = ''.join([i for i in tt if (not i.isdigit() and not i == ' ')])
          types.append(tt2.lower())
          x.append(float(line[30:38]))
          y.append(float(line[38:46]))
          z.append(float(line[46:54]))
except IOError as e:
  print("Could not open PDB file %s" % pdbfile)
  print( "File error ({0}): {1}".format(e.errno, e.strerror))
  self.terminate_code()

natm=len(types)

# now find bonded neighbors based on vdw radii
sdfFile=open('babel.sdf','w')
errFile=open('babel.err','w')
subprocess.run(['obabel','-ipdb',pdbfile,'-osdf'],
        stdout=sdfFile,stderr=errFile)
sdfFile.close()
errFile.close()

# read bonded pairs from newly created SDF file
nl=0
skip=0
nbond=0
bonds=[]
try:
  with open('babel.sdf', 'r') as sdf:
    for line in sdf:
      nl=nl+1
      a=line.split()
      if len(a) == 11:
        if a[0] != str(natm):
          raise Exception("No. atoms in babel.sdf does not match "+pdbFile)
          self.terminate_code()
        else:
          nbond = int(a[1])
          skip=nl
      if nl > skip+natm and nl <= skip+natm+nbond and nbond > 0:
          bonds.append([int(a[0]),int(a[1])])
except IOError as e:
  print("Could not open Babel SDF file babel.sdf, is obabel working?")
  print( "File error ({0}): {1}".format(e.errno, e.strerror))
  self.terminate_code()

if nbond == 0:
  raise Exception("Error reading bonded pairs from babel.sdf")
  self.terminate_code()

stoich=[]
neighbor=[]
for i in range(0,natm):
  stoich.append(0)
  neighbor.append([])
  for j in range(0,nbond):
    if bonds[j][0] == i+1 or bonds[j][1] == i+1:
      stoich[i]=stoich[i]+1
      if bonds[j][0] == i+1:
        neighbor[i].append(bonds[j][1])
      else:
        neighbor[i].append(bonds[j][0])

# Check that all atoms are involved in at least one bond:
inmol=[]
inframe=[]
for i in range(0,natm):
  if stoich[i] == 0:
    raise Exception("No bonded neighbors detected for atom "+str(i+1))
    self.terminate_code()
  inmol.append(0)
  inframe.append(0)
#  print("atom "+str(i+1)+" stoichiometry: "+str(stoich[i]))
#  print("neighbors: "+str(neighbor[i]))

# Check that all atoms are joined together:
def check_connectivity(atm):
  inmol[atm]=1
  for i in range(0,stoich[atm]):
    if inmol[neighbor[atm][i]-1]==0:
      check_connectivity(neighbor[atm][i]-1)

inmol[0]=1
check_connectivity(0)

for i in range(0,natm):
  if inmol[i] == 0:
    raise Exception("Discontinuity detected, atom "+str(i+1)+" is not connected to molecule")
    self.terminate_code()

# Now start creating frames:
# First the terminal atoms:
frames=[]
for i in range(0,natm):
  tframe=[]
  at2=-1
  at3=-1
  if stoich[i] == 1 and inframe[i] == 0:
    tframe.append(i+1)
    inframe[i]=1
    tframe.append(neighbor[i][0])
    at2=neighbor[i][0]-1
    inframe[at2]=1
    # try to add an atom that's not yet in any frame:
    for j in range(0,stoich[at2]):
      if inframe[neighbor[at2][j]-1] == 0:
        at3=neighbor[at2][j]-1
        tframe.append(at3+1)
        inframe[at3]=1
        break
    # otherwise add an atom that's already in a frame
    if len(tframe) != 3:
      for j in range(0,stoich[at2]):
        if neighbor[at2][j]-1 != i:
          at3=neighbor[at2][j]-1
          tframe.append(at3+1)
          break
    # if we still haven't completed the frame then something went wrong...
    if len(tframe) != 3:
      raise Exception("Could not complete frame "+str(tframe))
      self.terminate_code()
    frames.append(tframe)

# Now add any non-terminal atoms that are left:
for i in range(0,natm):
  tframe=[]
  at2=-1
  at3=-1
  if inframe[i] == 0:
    # try to add an atom that's not yet in any frame:
    tframe.append(i+1)
    inframe[i]=1
    for j in range(0,stoich[i]):
      if inframe[neighbor[i][j]-1] == 0:
        at2=neighbor[i][j]-1
        tframe.append(at2+1)
        inframe[at2]=1
        break
    # else just add any bonded neighbor:
    if at2 == -1:
      at2=neighbor[i][0]-1
      tframe.append(at2+1)
    # try to add a 3rd atom that's not yet in any frame:
    for j in range(0,stoich[at2]):
      if inframe[neighbor[at2][j]-1] == 0:
        at3=neighbor[at2][j]-1
        tframe.append(at3+1)
        inframe[at3]=1
        break
    # otherwise add an atom that's already in a frame
    if len(tframe) != 3:
      for j in range(0,stoich[i]):
        if neighbor[i][j]-1 != at2:
          at3=neighbor[i][j]-1
          tframe.append(at3+1)
          break
    # if we still haven't completed the frame then something went wrong...
    if len(tframe) != 3:
      raise Exception("Could not complete frame "+str(tframe))
      self.terminate_code()
    frames.append(tframe)

fr=open('frames.txt','w')
fr.write(resname+'\n')
for i in range(0,len(frames)):
  fr.write("%-4i %-4i %-4i\n" % (frames[i][0],frames[i][1],frames[i][2]))
fr.close()

print("\nfile frames.txt has been written\n")
