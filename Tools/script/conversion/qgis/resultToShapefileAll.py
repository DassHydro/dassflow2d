#!/usr/bin/env python
# encoding: utf-8

#from __future__ import unicode_literals
#from __future__ import print_function
#from __future__ import division
import math
import os, sys
os.chdir("/home/livillenave/Documents/distant/svn/dassflow-2d/trunk/tools/conversion/qgis/")
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from osgeo import ogr,osr # GDAL library
from osgeo.gdalconst import * # GDAL constants
import gdal
from gdalconst import *
from function import *


bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"

outputDirectoryQgis=f'{bin_dir}/mesh/'
outputDirectory=outputDirectoryQgis+'results'

#Creation of output directory Qgis
if os.path.isdir(outputDirectoryQgis)==False:
	os.mkdir(outputDirectoryQgis)

#Creation of output directory shapefile
if os.path.isdir(outputDirectory)==False:
	os.mkdir(outputDirectory)

rep=f'{bin_dir}/res/'


print('Working :')
print('Argument 1 : Path of the result directory (default : '+str(rep) +')\n')



# Get mesh file path
if len(sys.argv)==2:
   argument1=sys.argv[1]
   try:
      rep=str(argument1)
   except:
      pass





FileListResult=os.listdir(rep)
FileListResultDat=[]


for i in FileListResult:
   if i.find('.dat')!=-1:
      FileListResultDat.append(i)
   else :
      pass

ListeNumero=[]


for resultFile in FileListResultDat :
   try :
      numero=getNumeroResultFile(rep+resultFile)
      numero=float(numero)
      ListeNumero.append(numero)
   except:
      pass


FileListResultDat=[]
ListeNumero.sort()

if os.path.isfile(rep+'/result_initial.dat'):
   FileListResultDat.append(f'{bin_dir}/res/result_initial.dat')

for i in ListeNumero:
   FileListResultDat.append(f'{bin_dir}/res/result_'+"{:12.6E}".format(i)+'.dat')


if os.path.isfile(rep+'/result_final.dat'):
   FileListResultDat.append(f'{bin_dir}/res/result_final.dat')

itera=0
for nameFile in FileListResultDat:
   itera=itera+1
   resultFile=str(nameFile)
   EPGScode=2154
   #End to change


   try :
      print('Reading of '+str(resultFile))
      i, x,y,bathy,h,zs,Manning,u,v = np.loadtxt(resultFile,unpack=True)
   except :
      #print('Reprocessing of '+str(resultFile))
      #reprocessing(resultFile)
      i, x,y,bathy,h,zs,Manning,u,v=0,0,0,0,0,0,0,0,0

   numero=getNumeroResultFile(resultFile)



   # Creation of new shapefile
   driver = ogr.GetDriverByName('ESRI Shapefile')


   

   #Delete of old shapefile if exist
   if os.path.isfile(outputDirectory+'/result_'+str(numero)+'.shp'):
      filepath=outputDirectory+'/result_'+str(numero)
      os.remove(filepath+'.dbf')
      os.remove(filepath+'.prj')
      os.remove(filepath+'.shp')
      os.remove(filepath+'.shx')
   

   dataout = driver.CreateDataSource(outputDirectory)
   proj = osr.SpatialReference() 
   proj.ImportFromEPSG(EPGScode)


   layerout = dataout.CreateLayer('result_'+str(numero),proj,geom_type=ogr.wkbPoint)  
   fieldDef1 = ogr.FieldDefn('bathy', ogr.OFTReal)
   fieldDef2 = ogr.FieldDefn('h', ogr.OFTReal)
   fieldDef3 = ogr.FieldDefn('zs', ogr.OFTReal)
   fieldDef4 = ogr.FieldDefn('Manning', ogr.OFTReal)
   fieldDef5 = ogr.FieldDefn('u', ogr.OFTReal)
   fieldDef6 = ogr.FieldDefn('v', ogr.OFTReal)

   layerout.CreateField(fieldDef1)
   layerout.CreateField(fieldDef2)
   layerout.CreateField(fieldDef3)
   layerout.CreateField(fieldDef4)
   layerout.CreateField(fieldDef5)
   layerout.CreateField(fieldDef6)


   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)

   for i in range(0,len(x)):
      #- Create Geometry Point with pixel coordinates
      pixel_point = ogr.Geometry(ogr.wkbPoint)
      pixel_point.AddPoint(x[i],y[i])
      #- Add the geometry to the feature
      feature_out.SetGeometry(pixel_point)
      #- Set feature attributes
      feature_out.SetField('bathy'  , bathy[i])
      feature_out.SetField('h'      , h[i])
      feature_out.SetField('zs'     , zs[i])
      feature_out.SetField('Manning', Manning[i])
      feature_out.SetField('u'      , u[i])
      feature_out.SetField('v'      , v[i])
      layerout.CreateFeature(feature_out)
      #- Delete point geometry
      #pixel_point.Destroy()


   feature_out.Destroy()
   dataout.Destroy()
