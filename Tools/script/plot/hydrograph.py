#!/usr/bin/env python
# encoding: utf-8

#from __future__ import unicode_literals
#from __future__ import print_function
#from __future__ import division
from io import StringIO
import math
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from osgeo import ogr,osr # GDAL library
from osgeo.gdalconst import * # GDAL constants
import gdal
from gdalconst import *
from function import *

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl


bin_dir =  "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"

pathPlot=f'{bin_dir}/plot/'
pathOutFile=pathPlot+f'{bin_dir}/hydrograph/'
pathTargetHydrograph=f'{bin_dir}/hydrograph_target.txt'
pathFirstGuess=f'{bin_dir}/hydrograph_prior.txt'

i=0


##########
#To change
##########
legendTime        = 'hours' #Plot legend or not (examples: 'days' 'seconds' 's')
legend            = True
plotMin           = True
plotObs           = False
plotIntermediaire = False
maxplotInter      = 5

labelTarget    ='$Q^{real}$'
labeFirstGuess ='$Q^{FG}$'

titlefigMin    = 'Inflow discharge identification' 
titlefigInter  = 'Inflow discharge identification (with intermediary iterations)'

##########
#End to change
##########

pathobstxt=f'{bin_dir}/obs.txt'
pathparamobs=f'{bin_dir}/param_obs.txt'

factorTime=factorTimeFromLegend(legendTime)

Qplot=[]

#Get off set param_obs.txt file
l,offset,x1,y1,x2,y2=np.loadtxt(pathparamobs, unpack=True,skiprows=1)
if type(offset)!=type(np.array([])):
   nbrObs=int(1)
   offset=[offset]
   print(1)
else :
   nbrObs=len(offset)

#Get observation time in obs.txt file
ObsTime=[]
obsfile=open(pathobstxt,'r')
ObsFileLine=obsfile.readlines()
j=0
Toread=''
numberLineToRead=0
for i in ObsFileLine:
   print(i)
   if numberLineToRead>0 :
      tempi=i.replace('\n','')
      tempi=tempi.replace('\t','')
      tempi=tempi.replace(' ','')
      if tempi!='':
         Toread=Toread+str(i)
         numberLineToRead=numberLineToRead-1
   if ('stations_with_grp' in i):
      i=i.replace('stations_with_grp','')
      i=i.replace('\n','')
      i=i.replace('\t','')
      i=i.replace(' ','')
      numberLineToRead=int(i)
         
print(Tore
c=StringIO(Toread) 
nbr,Tobs,weight=np.loadtxt(c,unpack=True)

if type(Tobs)!=type(np.array([])):
   Tobs=[Tobs]



#Get Q(t) from target hydrograph
t,q=np.loadtxt(pathTargetHydrograph,skiprows=8,unpack=True)
tend=t[-1]

#For each swot band, get observations times on computation time domain.
DifferentTime=[]
for i in range(nbrObs):
   tempArray=[]
   t=float(offset[i])
   tempArray.append(t)
   while t < tend:
      t=t+Tobs[i]
      tempArray.append(t)
      if (t in DifferentTime)==False :
         DifferentTime.append(t)
   del tempArray[-1]
   ObsTime.append(tempArray)

   



#Creation of outputfile 
if (os.path.isdir(pathPlot)==False):
   os.mkdir(pathPlot)

if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)


#
# First plot : hydrograph.png
#

t,q=np.loadtxt(pathTargetHydrograph,skiprows=8,unpack=True)
plt.plot(t*factorTime,q,'b',linewidth=1,label=labelTarget)
Qplot.append(q)

t,q=np.loadtxt(pathFirstGuess,skiprows=8,unpack=True)
plt.plot(t*factorTime,q,'k--',linewidth=1,label=labeFirstGuess)
Qplot.append(q)


if plotMin:
   i=0
   stop=True
   plt.plot([-2,-1],[0,1])
   plt.plot([-2,-1],[0,1])
   while stop :
      a="{0:03d}".format(i+1)
      try :
         t,q=np.loadtxt('./bin/min/hydrograph_001.'+str(a),unpack=True)
         
         if plotIntermediaire:
            if i<(maxplotInter-1):
               plt.plot(t*factorTime,q,label='$Q^{ident}$ (ite: $'+str(i+1)+'$)')
               Qplot.append(q)
         i=i+1
         
             
      except :
         stop=False

   i=i-1
   a="{0:03d}".format(i+1)
   t,q=np.loadtxt('./bin/min/hydrograph_001.'+str(a),unpack=True)

   
   if plotIntermediaire==False :
      plt.plot(t*factorTime,q,color='g',label='$Q^{ident}$ (ite: $'+str(i+1)+'$)')
   else :
      plt.plot(t*factorTime,q,color='g',label='$Q^{ident}$ last ite\n (ite: $'+str(i+1)+'$)')
   Qplot.append(q)
      
   plt.ylabel('$Q_{in}$ $(m^3.s^{-1})$')
   plt.xlabel('$Time$ $(hours)$')
   plt.legend()
  


if plotObs:
   for i in range(nbrObs):
      tempArray=ObsTime[i]
      zero=list(np.zeros(len(tempArray))+Ymin)
      plt.plot(np.array(tempArray)*factorTime,zero,'x',linewidth=15,markersize=12)



if legend :
   plt.legend()
   

if plotMin:
   titlefig = titlefigMin
elif plotIntermediaire:
   titlefig = titlefigInter
else :
   titlefig = '' 

plt.title(titlefig)
plt.ylabel('$Q_{in}$ $(m^3.s^{-1})$')
plt.xlabel('$Time$ $('+str(legendTime)+ ')$')

[Ymin,Ymax]=computeYminYmax(Qplot,0.1)
plt.ylim(Ymin,Ymax)
plt.xlim(t[0],t[-1]*factorTime)
pathFile=pathOutFile+'plotHydrograph.png'
plt.savefig(pathFile)


#
# Second plot : hydrographObs.png
#
plt.clf()
t,q=np.loadtxt(pathTargetHydrograph,skiprows=8,unpack=True)
plt.plot(t*factorTime,q,'b',linewidth=1,label=labelTarget)
t,q=np.loadtxt(pathFirstGuess,skiprows=8,unpack=True)
plt.plot(t*factorTime,q,'k--',linewidth=1,label=labeFirstGuess)


if plotObs:
   for i in range(nbrObs):
      tempArray=ObsTime[i]
      zero=list(np.zeros(len(tempArray))+Ymin)
      plt.plot(np.array(tempArray)*factorTime,zero,'*')

   for i in range(nbrObs):
      tempArray=ObsTime[i]
      tempArray=np.array(tempArray)*factorTime
      for j in range(len(tempArray)):
         plt.axvline(x=tempArray[j],linewidth=0.5)

d=pathOutFile+'plotHydrographObs.png'
if legend :
   plt.legend()

if plotMin:
   titlefig = '$Q_{in}$ target and $Q_{in}$ first guess' 
else : 
   titlefig = 'Evolution of $Q_{in}$ target and first guess\n' + str(len(DifferentTime)) +' observation times, ' + str(len(ObsTime)) + ' SWOT biefs'


plt.title(titlefig)
plt.ylabel('$Q_{in}$ $(m^3.s^{-1})$')
plt.xlabel('$Time$ $( '+str(legendTime)+ ')$')
plt.savefig(d)
os.system('eog '+str(pathFile) + '&' ) 
