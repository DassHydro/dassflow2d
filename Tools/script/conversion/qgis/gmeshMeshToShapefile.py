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
#meshFile='./bin/test_case_deb.msh'
#outputDirectory='./bin/mesh'
meshFile='./garonne_jack.geo'
outputDirectory='.'
TraceNode=True
TraceEdge=True
EPGScode=2154
#End to change

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

# Get number of nodes
line=''
while (line!='$Nodes'):
   line=fileM.readline()
   line=deleteSpaceTab(line)

line=fileM.readline()
line=deleteSpaceTab(line)
numberOfNode=int(line)
print 'Number of nodes : ' + str(numberOfNode)

# Get nodes
Toread=''
for i in range(numberOfNode):
   line=fileM.readline()
   Toread=Toread+str(line)

c=StringIO(Toread) 
nbr,x,y,z=np.loadtxt(c,unpack=True)


# Get number of element
line=fileM.readline() #'$EndNodes'
line=fileM.readline() #'$Elements'
line=fileM.readline()
line=deleteSpaceTab(line)
numberOfElement=int(line)
print 'Number of elements : ' + str(numberOfElement)

# Get Elements
Toread=''
for i in range(numberOfElement):
   line=fileM.readline()
   lineSplit=line.split(' ')
   if str(lineSplit[1])!='2':
      pass
   else :
      Toread=Toread+str(line)

c=StringIO(Toread) 
nbr,a,b,c,d,pts1,pts2,pts3=np.loadtxt(c,unpack=True)

numberOfElement=len(a)


print 'Number of elements : ' + str(numberOfElement)

#Creation of output directory shapefile
if os.path.isdir(outputDirectory)==False:
	os.mkdir(outputDirectory)

# Creation of new shapefile
driver = ogr.GetDriverByName('ESRI Shapefile')
fshpout=outputDirectory
	
#Delete of old shapefile if exist
if os.path.exists(fshpout):
	driver.DeleteDataSource(fshpout)
dataout = driver.CreateDataSource(fshpout)
proj = osr.SpatialReference()

proj.ImportFromEPSG(EPGScode) # 4326 = EPSG code for lon/lat WGS84 projection


#Trace node
if TraceNode==True :
   print 'Creation node'
   layerout = dataout.CreateLayer('nodes',proj,geom_type=ogr.wkbPoint)  
   fieldDef1 = ogr.FieldDefn('z_b', ogr.OFTReal)
   layerout.CreateField(fieldDef1)
   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)

   for i in range(0,numberOfNode):
      #- Create Geometry Point with pixel coordinates
      pixel_point = ogr.Geometry(ogr.wkbPoint)
      pixel_point.AddPoint(x[i],y[i])
      #- Add the geometry to the feature
      feature_out.SetGeometry(pixel_point)
      #- Set feature attributes
      feature_out.SetField('z_b', float(z[i]))
      layerout.CreateFeature(feature_out)
      #- Delete point geometry
      pixel_point.Destroy()

#Trace Edge of cell
if TraceEdge==True :
   print 'Creation edge'
   
   layerout = dataout.CreateLayer('cells',proj,geom_type=ogr.wkbLineString) 
   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)
   for i in range(numberOfElement):
      line = ogr.Geometry(ogr.wkbLineString)
      line.AddPoint(x[int(pts1[i]-1)],y[int(pts1[i]-1)])
      line.AddPoint(x[int(pts2[i]-1)],y[int(pts2[i]-1)])
      line.AddPoint(x[int(pts3[i]-1)],y[int(pts3[i]-1)])
      line.AddPoint(x[int(pts1[i]-1)],y[int(pts1[i]-1)])
      feature_out.SetGeometry(line)
      layerout.CreateFeature(feature_out)
      line.Destroy()

feature_out.Destroy()
dataout.Destroy()



