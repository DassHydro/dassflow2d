#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys,os
from function import *
from matplotlib.pyplot import cm 
from matplotlib.legend_handler import HandlerLineCollection
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D


class HandlerDashedLines(HandlerLineCollection):
   def create_artists(self, legend, orig_handle,
                      xdescent, ydescent, width, height, fontsize, trans):
       # figure out how many lines there are
       numlines = len(orig_handle.get_segments())
       xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                            width, height, fontsize)
       leglines = []
       # divide the vertical space where the lines will go
       # into equal parts based on the number of lines
       ydata = ((height) / (numlines + 1)) * np.ones(xdata.shape, float)
       # for each line, create the line at the proper location
       # and set the dash pattern
       for i in range(numlines):
           legline = Line2D(xdata, ydata * (numlines - i) - ydescent)
           self.update_prop(legline, orig_handle, legend)
           # set color, dash pattern, and linewidth to that
           # of the lines in linecollection
           try:
               color = orig_handle.get_colors()[i]
           except IndexError:
               color = orig_handle.get_colors()[0]
           try:
               dashes = orig_handle.get_dashes()[i]
           except IndexError:
               dashes = orig_handle.get_dashes()[0]
           try:
               lw = orig_handle.get_linewidths()[i]
           except IndexError:
               lw = orig_handle.get_linewidths()[0]
           if dashes[0] != None:
               legline.set_dashes(dashes[1])
           legline.set_color(color)
           legline.set_transform(trans)
           legline.set_linewidth(lw)
           leglines.append(legline)
       return leglines


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
pathInput1='./bin/obs/'
pathInput2='./bin/res/'

factorTime=factorTimeFromLegend(legendTime)

#Creation of outputfile 
if (os.path.isdir(pathPlot)==False):
   os.mkdir(pathPlot)

if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)


print 'Working :'
print 'Argument 1: Number of station (default: '+str(nbrStation) +')'
print 'Argument 2: Legend time (default: '      +str(legendTime) +')'
print 'Argument 3: Directory files 1 (default: '+str(pathInput1) +')'
print 'Argument 4: Directory files 2 (default: '+str(pathInput2) +')\n'


# Get argument
if len(sys.argv)==4:
   argument1=sys.argv[1]
   try:
      nbrStation=int(argument1)
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

color=cm.brg(np.linspace(0,1,nbrStation))
Txt,Line=[],[]
for i in stations:



   a="{0:04d}".format(i+1)
   #Data 1:
   try:
      data=np.loadtxt(str(pathInput1)+'obs_station_'+str(a)+'.plt',skiprows=1)
   except:
      data=np.loadtxt(str(pathInput1)+'obs_station_'+str(a)+'.dat',skiprows=1)  

   t1=data[:,0]
   h1=data[:,1]

   #Data 2:
   try:
      data=np.loadtxt(str(pathInput2)+'obs_station_'+str(a)+'.plt',skiprows=1)
   except:
      data=np.loadtxt(str(pathInput2)+'obs_station_'+str(a)+'.dat',skiprows=1)  

   t2=data[:,0]
   h2=data[:,1]



   ax.plot(t1*factorTime,h1,color=color[i])
   ax.plot(t2*factorTime,h2,label='Bief $'+str(int(a))+'$',linestyle='dashed',color=color[i])

   line = [[(0, 0)]]
   l = LineCollection(2 * line, linestyles = ['solid', 'dashed'], colors = [color[i],color[i]])
   txt='Bief $'+str(int(a))+'$'
   Line.append(l)
   Txt.append(txt)

   titlefig='$H^{obs}$ observed in time at each SWOT bief\n Groups 1: solid line, Groups 2: dashed line'
   plt.title(titlefig)
   plt.ylabel('$h^{obs}$ $(m)$')
   plt.xlabel('$Time$ $( '+str(legendTime)+ ')$')


   
   H.append(h1)
   H.append(h2)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
[Ymin,Ymax]=computeYminYmax(H,0.1)
ax.legend(Line,Txt,handler_map = {type(l) : HandlerDashedLines()}, handlelength = 2.5,loc='center left', bbox_to_anchor=(1, 0.5))


plt.ylim(Ymin,Ymax)

plt.xlim(t1[0]*factorTime,max([t1[-1],t2[-1]])*factorTime)
pathFile=pathOutFile+'obsxh_compare.png'
plt.savefig(pathFile)


#
plt.clf()


fig = plt.figure()
ax = plt.subplot(111)


W=[]
for i in stations:
   a="{0:04d}".format(i+1)

   #Data 1:
   try:
      data=np.loadtxt(str(pathInput1)+'obs_station_'+str(a)+'.plt',skiprows=1)
   except:
      data=np.loadtxt(str(pathInput1)+'obs_station_'+str(a)+'.dat',skiprows=1)  

   t1=data[:,0]
   w1=data[:,4]

   #Data 2:
   try:
      data=np.loadtxt(str(pathInput2)+'obs_station_'+str(a)+'.plt',skiprows=1)
   except:
      data=np.loadtxt(str(pathInput2)+'obs_station_'+str(a)+'.dat',skiprows=1)  

   t2=data[:,0]
   w2=data[:,4]

   titlefig='$w^{obs}$ observed in time at each SWOT bief\n Groups 1: solid line, Groups 2: dashed line'
   plt.title(titlefig)
   plt.ylabel('$w^{obs}$ $(m)$')
   plt.xlabel('$Time$ $( '+str(legendTime)+ ')$')
   ax.plot(t1*factorTime,w1,color=color[i])
   ax.plot(t2*factorTime,w2,label='Bief $'+str(int(a))+'$',linestyle='dashed',color=color[i])
   plt.legend()
   W.append(w1)
   W.append(w2)


box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
[Ymin,Ymax]=computeYminYmax(W,0.1)
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.legend(Line,Txt,handler_map = {type(l) : HandlerDashedLines()}, handlelength = 2.5,loc='center left', bbox_to_anchor=(1, 0.5))


plt.ylim(Ymin,Ymax)
plt.xlim(t1[0]*factorTime,max([t1[-1],t2[-1]])*factorTime)
pathFile=pathOutFile+'obsxw_compare.png'
plt.savefig(pathFile)





os.system('eog '+str(pathFile) + '&' ) 


