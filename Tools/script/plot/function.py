import numpy as np
from io import StringIO
import os 

def Q_from_SF(t,a0,n,A,B,T):
   Q=float(a0)/2
   for i in range(len(A)):
      Q=Q+A[i]*np.cos((n[i]*2*t*np.pi)/T)+B[i]*np.sin((n[i]*2*t*np.pi)/T)
   return Q


def ReadLineClean(line,split):
   line=line.split(split)
   sizeline=len(line)
   nbrDel=0
   for i in range(sizeline):
      if line[i-nbrDel]=='':
         del line[i-nbrDel]
         nbrDel=nbrDel+1
   return line


def readCoefFromOriFile(pathfile):
   A,B,n=[],[],[]
   minHydro_FS_file=open(pathfile,'r')
   Lines=minHydro_FS_file.readlines()
   a0=float(Lines[4])
   #for i in range(len(Lines)-5):
   #   line=Lines[i+5]
   #   print line
   #   line=ReadLineClean(line,'\t')
   #   line=ReadLineClean(line,' ')
   #   n.append(int(float(line[0])))
   #   A.append(float(line[1]))
   #   B.append(float(line[2]))
   ToRead=''
   for i in range(len(Lines)-5):
      ToRead=ToRead+Lines[i+5]
   print(ToRead)
   ToRead= StringIO(ToRead)

   n,A,B=np.loadtxt(ToRead,unpack=True)

   n=list(n)
   A=list(A)
   B=list(B)
   return [a0,n,A,B]


def readFromInputFile(path,parameters):
   input_file=open(path,'r')
   Lines=input_file.readlines()

   #Delete Comment
   for i in range(len(Lines)):
      DeleteComment=Lines[i].find('!')
      if DeleteComment!=-1:
         Lines[i]=Lines[i][0:DeleteComment]

   #Find parameters
   Todelete=[' ','\n','=',',','\t','\'']
   parametersIdentificationName=[]
   for i in range(len(parameters)):
      parametersIdentificationName.append(str(parameters[i]))
   parametersIdentification=[]

   
   for i in parametersIdentificationName:
      parameters=[x for x in Lines if i in x]
      if len(parameters)==1:
         parameters=parameters[0]
         for j in Todelete :
            parameters=parameters.replace(j,'')

         parameters=parameters.replace(i,'')

      else :
         parameters=int(0)
      parametersIdentification.append((parameters))

   return parametersIdentification

def readCoefFromMinFile(pathfile):
   A,B,n=[],[],[]
   minHydro_FS_file=open(pathfile,'r')
   Lines=minHydro_FS_file.readlines()
   a0=float(Lines[0])
   for i in range(len(Lines)-1):
      line=Lines[i+1]
      line=ReadLineClean(line,' ')
      n.append(int(line[0]))
      A.append(float(line[1]))
      B.append(float(line[2]))
   return [a0,n,A,B]

def computeYminYmax(data,deltacoeff):

   try:
      lenData=len(data)
   except:
      lenData=0

   try :
      lenData1=len(data[0])
      Ymax=max(data[0][1])
      Ymin=min(data[0][1])
      for i in range(lenData1):
         ymax=max(data[0][i])
         ymin=min(data[0][i])
         Ymin=min(Ymin,ymin)    
         Ymax=max(Ymax,ymax)
   except :
      Ymax=max(data[0])
      Ymin=min(data[0])

   for i in range(1,lenData):
      temp=data[i]

      try :
         ymax=max(temp)
         ymin=min(temp)
      except:
         ymax=temp
         ymin=temp

      try : 
         ymax=ymax[0]
         ymin=ymin[0]
      except :
         pass

      Ymin=min(Ymin,ymin)
      Ymax=max(Ymax,ymax)

   deltay=Ymax-Ymin
   Ymin=Ymin-deltacoeff*deltay
   Ymax=Ymax+deltacoeff*deltay
   return [Ymin,Ymax]


def adaptedYlim(data):
   try:
      lenData=len(data)
   except:
      lenData=0

   for i in range(lenData):
      temp=data[i]
      try :
         
         ymax=max(temp)
         ymin=min(temp)
      except:
         ymax=temp
         ymin=temp

      try : 
         ymax=ymax[0]
         ymin=ymin[0]
      except :
         pass

      if i==0:
         Ymin=ymin
         Ymax=ymax
      else :      
         Ymin=min(Ymin,ymin)
         Ymax=max(Ymax,ymax)

   deltay=Ymax-Ymin
   Ymin=Ymin-0.1*deltay
   Ymax=Ymax+0.1*deltay
   return [Ymin,Ymax]



def factorTimeFromLegend(legendTime):
   if legendTime=='years':
      factorTime=1./(3600*24*365)
   elif legendTime=='days':
      factorTime=1./(3600*24)
   elif legendTime=='hours' :
      factorTime=1./3600
   elif legendTime=='minutes' :
      factorTime=1./60
   elif legendTime=='seconds':
      factorTime=1.
   elif legendTime=='s':
      factorTime=1.
   else:
      factorTime=1.
   return factorTime


def factorSpace(legendTime):
   if legendTime=='km' :
      factorTime=1./1000
   elif legendTime=='m':
      factorTime=1.
   else:
      factorTime=1.
   return factorTime




def getObsTime(pathobs,pathparam,T):

   tobs_offSet = []


   #Read obs.txt
   obsfile=open(pathobs,'r')
   obsFileLine=obsfile.readlines()
   numberLine=len(obsFileLine)

   #Read stations
   i=0
   line=obsFileLine[i]
   while (('stations' in line)==False):
      line=obsFileLine[i]
      i=i+1

   line=line.replace('\n','')
   line=line.replace('\t','')
   line=line.replace('stations','')
   line=line.replace(' ','')

   numberOfStation=int(line)


   if numberOfStation!=0:
      toread=''
      for j in range(i+1,i+1+numberOfStation):
         toread=toread+obsFileLine[j]


      c=StringIO(toread) 
      x,y,tobs_period,weight=np.loadtxt(c,unpack=True)

   #Read stations_with_grp
   i=0
   line=obsFileLine[i]
   while (('stations_with_grp' in line)==False):
      line=obsFileLine[i]
      i=i+1

   line=line.replace('\n','')
   line=line.replace('\t','')
   line=line.replace('stations_with_grp','')
   line=line.replace(' ','')

   numberOfStation=int(line)

   if numberOfStation!=0:
      toread=''
      for j in range(i,i+numberOfStation):
         toread=toread+obsFileLine[j]

      c=StringIO(toread) 
      x,tobs_period,weight=np.loadtxt(c,unpack=True)


   #
   # read param_obs.txt
   #
   L,t_obs_off_set=np.loadtxt(pathparam,skiprows=1,unpack=True)

   numberGroups=[]
   listeCouple=[]
   groups=[]
   for i in range(len(tobs_period)):
      couple=[tobs_period[i],t_obs_off_set[i]]
      if couple in listeCouple:
         index=listeCouple.index(couple)
         temp=groups[index]
         temp.append(i)
         groups[index]=temp
         numberGroups[index]=numberGroups[index]+1
      else :
         listeCouple.append(couple)
         numberGroups.append(1)
         groups.append([i])

   obsTimeT=[]

   #For each couple, generate ObsTimeT
   for i in range(len(listeCouple)):
      dt      = listeCouple[i][0]
      t       = listeCouple[i][1]
      t_temp=[t]
      while t+dt<=T:
         t=t+dt
         t_temp.append(t)
   
      obsTimeT.append(t_temp)
 
   
   return [obsTimeT,numberGroups,groups]



def getObsTime_from_obs_file(pathobsdirectory):
   Tobs,numberGroups=[],[]
   i=1
   a="{0:04d}".format(i)
   pathobsfile=pathobsdirectory+'obs_station_'+str(a)+'.dat'
   print(pathobsfile)
   while os.path.exists(pathobsfile) :
      t,a,b,c,d=np.loadtxt(pathobsfile,unpack=True)
      t=list(t)
      if t in Tobs:
         index=Tobs.index(t)
         numberGroups[index]=numberGroups[index]+1
      else :
         Tobs.append(t)
         numberGroups.append(1)
      i=i+1
      a="{0:04d}".format(i)
      pathobsfile=pathobsdirectory+'obs_station_'+str(a)+'.dat'

   return [Tobs,numberGroups]
