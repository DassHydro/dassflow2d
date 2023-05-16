#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys,os
from function import *

# Plot H and w observed in time at each SWOT bief


##########
#To change
##########
nbrStation=10
stations=range(nbrStation)
legendTime        = 'hours' 


##########
#End to change
##########


pathPlot='./bin/plot/'
pathOutFile=pathPlot+'./obs/'
pathInput='./bin/obs/'
factorTime=factorTimeFromLegend(legendTime)

#Creation of outputfile 
if (os.path.isdir(pathPlot)==False):
   os.mkdir(pathPlot)

if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)


print 'Working :'
print 'Argument 1: Number of station (default: '+str(nbrStation) +')'
print 'Argument 2: Legend time (default: '+str(legendTime) +')'
print 'Argument 3: Directory files (default: '+str(pathInput) +')\n'


# Get mesh file path
if len(sys.argv)==4:
   argument1=sys.argv[1]
   try:
      stations=range(int(argument1))
   except:
      pass

   argument1=sys.argv[2]
   try:
      legendTime=str(argument1)
   except:
      pass

   argument1=sys.argv[3]
   try:
      pathInput=str(argument1)
   except:
      pass


H=[]
fig = plt.figure()
ax = plt.subplot(111)


for i in stations:



   a="{0:04d}".format(i+1)
   try:
      data=np.loadtxt(str(pathInput)+'obs_station_'+str(a)+'.plt',skiprows=1)
   except:
      data=np.loadtxt(str(pathInput)+'obs_station_'+str(a)+'.dat',skiprows=1)  

   t=data[:,0]
   h=data[:,1]




   ax.plot(t*factorTime,h,label='Bief $'+str(int(a))+'$')
   titlefig='$h^{obs}$ observed in time at each SWOT bief'
   plt.title(titlefig)
   plt.ylabel('$h^{obs}$ $(m)$')
   plt.xlabel('$Time$ $( '+str(legendTime)+ ')$')
   
   H.append(h)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
[Ymin,Ymax]=computeYminYmax(H,0.1)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(Ymin,Ymax)
plt.xlim(t[0]*factorTime,t[-1]*factorTime)
pathFile=pathOutFile+'obsxh.png'
plt.savefig(pathFile)


#
plt.clf()


W=[]

fig = plt.figure()
ax = plt.subplot(111)


for i in stations:



   a="{0:04d}".format(i+1)
   try:
      data=np.loadtxt(str(pathInput)+'obs_station_'+str(a)+'.plt',skiprows=1)
   except:
      data=np.loadtxt(str(pathInput)+'obs_station_'+str(a)+'.dat',skiprows=1)  
   t=data[:,0]
   w=data[:,4]

   titlefig='$w^{obs}$ observed in time at each SWOT bief'
   plt.title(titlefig)
   plt.ylabel('$w^{obs}$ $(m)$')
   plt.xlabel('$Time$ $( '+str(legendTime)+ ')$')
   plt.plot(t*factorTime,w,label='Bief $'+str(int(a))+'$')
   plt.legend()
   W.append(w)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
[Ymin,Ymax]=computeYminYmax(W,0.1)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(Ymin,Ymax)
plt.xlim(t[0]*factorTime,t[-1]*factorTime)
pathFile=pathOutFile+'obsxw.png'
plt.savefig(pathFile)

os.system('eog '+str(pathFile) + '&' ) 


