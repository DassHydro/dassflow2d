#!/usr/bin/env python
# encoding: utf-8

#from __future__ import unicode_literals
#from __future__ import print_function
#from __future__ import division

import math
import os, sys
os.chdir("/home/livillenave/Documents/distant/svn/dassflow-2d/trunk/tools/conversion/qgis/")
import numpy as np
#from netCDF4 import Dataset
import matplotlib.pyplot as plt
from osgeo import ogr,osr # GDAL library
from osgeo.gdalconst import * # GDAL constants
import gdal
from gdalconst import *
from function import *


# Create from a DassFlow mesh a shapefile file of the mesh ( node and edge)

#To change
pathInput='/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/input.txt'
meshFile='/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/automaticaly_generated_mesh.txt'
TraceNode=False
TraceEdge=True
EPGScode=2154

outputDirectoryQgis="/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"
outputDirectory=outputDirectoryQgis+'mesh'

#End to change

#
#meshFile=readFromInputFile(pathInput,['mesh_name'])
#meshFile='./bin/'+str(meshFile[0])

print( 'Working :')
print( 'Argument 1 : Path of the mesh file (default : '+str(meshFile) +')')
print( 'Argument 2 : Trace node 1-> True , 0 -> False (default : '+str(TraceNode) +')')
print('Argument 3 : Trace edge 1-> True , 0 -> False (default : '+str(TraceEdge) +')\n')



# Get mesh file path
if len(sys.argv)==2:
   argument1=sys.argv[1]
   try:
      meshFile=str(argument1)
   except:
      pass

# Get mesh file path + node bool + edge bool
if len(sys.argv)==4:
   argument1=sys.argv[1]
   try:
      meshFile=str(argument1)
   except:
      pass

   try:
      TraceNode=int(argument1)
      if TraceNode==1:
         TraceNode=True
      else: 
         TraceNode=False
   except:
      pass

   try:
      TraceEdge=int(argument1)
      if TraceEdge==1:
         TraceEdge=True
      else: 
         TraceEdge=False
   except:
      pass




x,y=[],[]
#Get number of nodes, element and size:
try :
   fileM = open(meshFile,'r')
except :
   print( 'File '+str(meshFile) + ' does not exist\n')

line=fileM.readline()
line=fileM.readline().split()
numberOfNode=int(line[0])
numberOfElement=int(line[1])
print( 'number of nodes :   ' + str(numberOfNode))
print( 'number of element : ' + str(numberOfElement))
Node=fileM.readlines()


#Creation of output directory Qgis
if os.path.isdir(outputDirectoryQgis)==False:
	os.mkdir(outputDirectoryQgis)

#Creation of output directory shapefile
if os.path.isdir(outputDirectory)==False:
	os.mkdir(outputDirectory)


#Delete of old shapefile if exist
if TraceNode==True:
   if os.path.isfile(outputDirectory+'/nodes.shp'):
      filepath=outputDirectory+'/nodes'
      os.remove(filepath+'.dbf')
      os.remove(filepath+'.prj')
      os.remove(filepath+'.shp')
      os.remove(filepath+'.shx')

if TraceEdge==True:
   if os.path.isfile(outputDirectory+'/cells.shp'):
      filepath=outputDirectory+'/cells'
      os.remove(filepath+'.dbf')
      os.remove(filepath+'.prj')
      os.remove(filepath+'.shp')
      os.remove(filepath+'.shx')

	
# Creation of new shapefile
driver = ogr.GetDriverByName('ESRI Shapefile')
dataout = driver.CreateDataSource(outputDirectory)
proj = osr.SpatialReference()

proj.ImportFromEPSG(EPGScode) # 4326 = EPSG code for lon/lat WGS84 projection


#Trace node
if TraceNode==True :
   print( 'Creation node')
   layerout = dataout.CreateLayer('nodes',proj,geom_type=ogr.wkbPoint)  
   fieldDef1 = ogr.FieldDefn('z_b', ogr.OFTReal)
   layerout.CreateField(fieldDef1)
   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)

   for i in range(0,numberOfNode):
      #- Create Geometry Point with pixel coordinates
      pixel_point = ogr.Geometry(ogr.wkbPoint)
      x.append(float(Node[i+1].split()[1]))
      y.append(float(Node[i+1].split()[2]))
      pixel_point.AddPoint(x[i],y[i])
      #- Add the geometry to the feature
      feature_out.SetGeometry(pixel_point)
      #- Set feature attributes
      feature_out.SetField('z_b', float(Node[i+1].split()[3]))
      layerout.CreateFeature(feature_out)
      #- Delete point geometry
      pixel_point.Destroy()

#Trace Edge of cell
if TraceEdge==True :
   print( 'Creation edge')
   if TraceNode==False:
      for i in range(0,numberOfNode):
         x.append(float(Node[i+1].split()[1]))
         y.append(float(Node[i+1].split()[2]))
   
   layerout = dataout.CreateLayer('cells',proj,geom_type=ogr.wkbLineString) 
   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)
   for i in range(numberOfNode+1,numberOfNode+1+numberOfElement):
      line = ogr.Geometry(ogr.wkbLineString)
      for j in range(1,5):
         line.AddPoint(x[int(Node[i+1].split()[j])-1],y[int(Node[i+1].split()[j])-1])
      line.AddPoint(x[int(Node[i+1].split()[1])-1],y[int(Node[i+1].split()[1])-1])
      feature_out.SetGeometry(line)
      layerout.CreateFeature(feature_out)
      #line.__del__

feature_out.Destroy()
dataout.Destroy()

