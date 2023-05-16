#!/usr/bin/env python
# encoding: utf-8

#from __future__ import unicode_literals
#from __future__ import print_function
#from __future__ import division
import math
import os, sys
import numpy as np

# Description : Cree un fichier nameFile de cellules d'observation ('station_grp_xxx') a partir d'une zone definie par des bornes inf et sup et d'un fichier de resultat ( pour recuperer le centre des cellules)

# TO change 
nameFile='station_grp_0001.txt'
resultFile='./bin/res/result_initial.dat'
xlimInf = 450
xlimSup = 550
ylimInf = -300
ylimSup = 300
outputDirectory='./bin/'
#

try :
   x,y,bathy,h,zs,Manning,u,v=np.loadtxt(resultFile,unpack=True)
except :
   print 're traitement'
   fileInput=open(resultFile,'r')
   fileOutput=open(resultFile+'_post','w')
   Lines=fileInput.readlines()
   for i in range(len(Lines)):
      Line=Lines[i]
      Line=Line.replace('-',' -')
      Line=Line.replace('E -','E-')
      fileOutput.write(Line)
   fileOutput.close()
   x,y,bathy,h,zs,Manning,u,v=np.loadtxt(resultFile+'_post',unpack=True)

StationX=[]
StationY=[]
for i in range(len(x)):
   if ((xlimInf<x[i]) and (xlimSup>x[i]) and (ylimInf<y[i]) and (ylimSup>y[i])):
      StationX.append(x[i])
      StationY.append(y[i])


print 'Station : '+ str(len(StationX)) + ' points.' 

file=open(outputDirectory+str(nameFile), 'w')
file.write('points\n')
for j in range(len(StationX)) :
   file.write(str(StationX[j]) +  '  '+str(StationY[j])+'\n')
file.close()
   
