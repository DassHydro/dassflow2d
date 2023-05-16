import numpy as np

# From groups of station locate by points, this script create a groups of station locate by indexes.


#Modifiable
nomDossierInput='./inputs/stations'
nomDossierOuput='./outputs/stations'
nomFichierSortie='./inputs/result_initial.dat'
nombreStation=22
#Finmodifiable

#Reccuperation des noms des stations d'entree

nameStationGrpPoints=[]
nameStationsGrpIndex=[]
for i in range(1,nombreStation+1) :
   a="{0:04d}".format(i) 
   nameStationGrpPoints.append(nomDossierInput+'/station_grp_'+a+'.txt')
   nameStationsGrpIndex.append(nomDossierOuput+'/station_grp_'+a+'.txt')

#print nameStationGrpPoints

#Reccuperation des coordonnees des centres des cellules
X,Y,bathy,h,zs,Manning,u,v=np.loadtxt(nomFichierSortie,unpack=True)


# Lecture des stations d'entrees
#Pour chaque station ...
i=0
numFichier=0
index=[]
for fichierStationPoints in nameStationGrpPoints :
   #Reccuperation des l'ensemble des points
   print fichierStationPoints
   data=np.loadtxt(fichierStationPoints,skiprows=1)
   x=data[:,0]
   y=data[:,1]
   index=[]
   #Pour chaque points 
   for i in range(len(x)):
      #print np.where(X==x[i])[0],np.where(Y==y[i])[0]
      index.append(np.intersect1d(np.where(X==x[i])[0],np.where(Y==y[i])[0])+1)
      #print
      #i=i+1 
      if len(index[-1]) !=1 :
         print 'pbl'
      #if i%100==0 :
      #   print i
   #print 1
   obsfile=open(nameStationsGrpIndex[numFichier],'w+') # Create group station
   obsfile.write('indexes\n')
   for j in range(0,len(index)) : #For each point ...
      ind=index[j]
      obsfile.write(str(ind[0])+'\n') # ... write coor
   numFichier=numFichier+1
   obsfile.close()


