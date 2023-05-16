#!/usr/bin/env python
# encoding: utf-8

#from __future__ import unicode_literals
#from __future__ import print_function
#from __future__ import division
from StringIO import StringIO
import math
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from osgeo import ogr,osr # GDAL library
from osgeo.gdalconst import * # GDAL constants
import gdal
from gdalconst import *


#To change
meshFile='./bin/test_case_deb3.msh'
meshFileOut='./bin/test_case_deb_bathy.msh'

 
#End to change



def bathy(x,y):
   xmin1,xmax1=0  , 300 
   xmin2,xmax2=0  , 300 
   ymin1,ymax1=50 , 100 
   ymin2,ymax2=0  , 50

   x=float(x)
   y=float(y)

   if ((xmin1<=x) and (xmax1>=x) and (ymin1<=y) and (ymax1>=y)):
      z=-(1./300)*y+2.0

   if ((xmin2<=x) and (xmax2>=x) and (ymin2<=y) and (ymax2>=y)):
      z=-(2./300)*y+2.0

   else :
      z =0.0

   return z

def deleteSpaceTab(line):
   line=line.replace('\n','')
   line=line.replace('\t','')
   line=line.replace(' ','')
   line=str(line)
   return line
x,y=[],[]
#
# Read of the file
#
fileM = open(meshFile,'r')
fileMOut = open(meshFileOut,'w')

# Get number of nodes
line=''
while (line!='$Nodes'):
   line=fileM.readline()
   fileMOut.write(line)
   line=deleteSpaceTab(line)

line=fileM.readline()
fileMOut.write(line)
line=deleteSpaceTab(line)
numberOfNode=int(line)
print 'Number of nodes : ' + str(numberOfNode)

# Get nodes
Toread=''
for i in range(numberOfNode):
   line=fileM.readline()
   lineSplit=line.split(' ')
   x=lineSplit[1]
   y=lineSplit[2]
   lineSplit[3]=str(bathy(x,y))
   line=' '.join(lineSplit)+'\n'
   fileMOut.write(line)


# Get number of element
line=fileM.readline() #'$EndNodes'
fileMOut.write(line)
line=fileM.readline() #'$Elements'
fileMOut.write(line)
line=fileM.readline()
fileMOut.write(line)
line=deleteSpaceTab(line)
numberOfElement=int(line)
print 'Number of elements : ' + str(numberOfElement)

# Get Elements
for i in range(numberOfElement):
   line=fileM.readline()
   fileMOut.write(line)

line=fileM.readline()
fileMOut.write(line)

