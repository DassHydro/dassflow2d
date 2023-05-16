#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from function import *

# Plot  w observed in time at each SWOT bief


##########
#To change
##########

legendTime        = 'hours' 
stations=range(10)

##########
#End to change
##########


pathPlot='./bin/plot/'
pathOutFile=pathPlot+'./obs/'
factorTime=factorTimeFromLegend(legendTime)

#Creation of outputfile 
if (os.path.isdir(pathPlot)==False):
   os.mkdir(pathPlot)

if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)

W=[]
for i in stations:



   a="{0:04d}".format(i+1)
   data=np.loadtxt('./bin/obs/obs_station_'+str(a)+'.plt',skiprows=1)

   t=data[:,0]
   w=data[:,4]

   titlefig='$w^{obs}$ observed in time at each SWOT bief'
   plt.title(titlefig)
   plt.ylabel('$w^{obs}$ $(m)$')
   plt.xlabel('$Time$ $( '+str(legendTime)+ ')$')
   plt.plot(t*factorTime,w,label='Bief $'+str(int(a))+'$')
   plt.legend()
   W.append(w)

[Ymin,Ymax]=computeYminYmax(W,0.1)
plt.ylim(Ymin,Ymax)
plt.xlim(t[0]*factorTime,t[-1]*factorTime)
pathFile=pathOutFile+'obsxw.png'
plt.savefig(pathFile)

os.system('eog '+str(pathFile) + '&' ) 
