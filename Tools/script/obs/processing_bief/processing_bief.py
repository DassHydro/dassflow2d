import ogr, os, sys
import numpy as np


#To change
SWOTBandDir= '../SWOT_like_data/orbit891km_Garonne_SWOTswathbands1km/orbit891km_Garonne_SWOTswathbands1km' #Path of output directory of SWOT_like_data SWOT bands
NodesGrpDir= '../SWOT_like_data/orbit891km_Garonne_SWOTswathbands1km/orbit891km_Garonne_SWOTswathbands1km/groups' #Path of output directory of SWOT_like_data groups of nodes
centerline=  './inputs/centerline.shp' #Path of center line
outputDir='../SWOT_like_data/ouputs'
#End to change

SWOTBandFiles=[]           #SWOT Band files
numberOf1IntersectionT=0   #Number of swot band which intersects only one time center line
numberOfSWOTBandIntersect=0#Number of swot band which intersects at least center line

outputDirPts=outputDir+'/points'
outputDirIndex=outputDir+'/indexes'


if os.path.isdir(outputDir)==False:
	os.mkdir(outputDir)

if os.path.isdir(outputDirPts)==False:
	os.mkdir(outputDirPts)

if os.path.isdir(outputDirIndex)==False:
	os.mkdir(outputDirIndex)

lenght,leftCoor,rigthCoor,Bandnumber,OutDistance,ObsnodeFile,timesList=[],[],[],[],[],[],[]

driver = ogr.GetDriverByName('ESRI Shapefile')

listFile=os.listdir(SWOTBandDir) #List of swot band

for fileB in listFile: #Selection of useful swot band
   if (fileB.find('s_l93.shp')!=-1):
      SWOTBandFiles.append('/'+fileB) 


dsSWOTPolyline = driver.Open(centerline) #Load center line



for SWOTBandFile in SWOTBandFiles:
   times=SWOTBandFile[fileB.find('km_time')+8:fileB.find('s_l93.shp')-8]
   times=str(np.round(float(times)))
   print 'Time :' + times
   dsSWOTBand     = driver.Open(SWOTBandDir+SWOTBandFile)
   SWOTBand       = dsSWOTBand.GetLayerByIndex(0)
   SWOTPolyline   = dsSWOTPolyline.GetLayerByIndex(0)

   c    =SWOTPolyline.GetFeature(0)
   Line  =c.GetGeometryRef()

   NumberOfBand=SWOTBand.GetFeatureCount()
   for i in range(0,NumberOfBand): #For each under band of swot band

      a     =SWOTBand.GetFeature(i)
      Band  =a.GetGeometryRef()
      outDis=a.GetField('Out dist')

      
      intersec=Line.Intersection(Band)
      if intersec.GetPointCount()!=0: # show that the under swot band intersects only one time the center line
         numberOf1IntersectionT=numberOf1IntersectionT+1
         numberOfSWOTBandIntersect=numberOfSWOTBandIntersect+1
         numberOfIntersection=1
      else : #Else
         numberOfIntersection=intersec.GetGeometryCount()
         
      

      if (numberOfIntersection==1): #If the under swot band is intersect only one time
         print 'Band ' +str(i) + ', outdis '+str(outDis)+' intersects river ' +str(numberOfIntersection) + ' times.'
         print '\tLength          : '+str(intersec.Length())
         print '\tLeft  Extremity : '+str(intersec.GetPoints()[0])
         print '\tRight Extremity : '+str(intersec.GetPoints()[-1])
         lenght.append(float(intersec.Length()))
         leftCoor.append((intersec.GetPoints()[0]))
         rigthCoor.append((intersec.GetPoints()[-1]))
         Bandnumber.append(int(i))
         OutDistance.append(np.round(float(outDis)))
         ObsnodeFile.append(str('gr_nodes_'+str(times)+'_d'+str(OutDistance[-1])+'.txt'))
         timesList.append(times)

      elif (numberOfIntersection>1) :# If the under swot band is intersect more than one time
         numberOfSWOTBandIntersect=numberOfSWOTBandIntersect+1 
         #print 'Band ' +str(i) + ', outdis '+str(outDis)+' intersects river ' +str(numberOfIntersection) + ' times.'   
         for j in range(0,numberOfIntersection):
            
            geo=intersec.GetGeometryRef(j)
            #print 'File ' + SWOTBandFile
            #print 'Intersection ' +str(j) + ' : '
            #print '\tLength          : '+str(geo.Length())
            #print '\tLeft  Extremity : '+str(geo.GetPoints()[0])
            #print '\tRight Extremity : '+str(geo.GetPoints()[-1])
      else : # Else
         pass
         #print 'Band ' +str(i) + ' does not intersect river '

print 'Number of SWOT Band intersects : ' + str(numberOfSWOTBandIntersect)
print 'Number of SWOT Band with 1 intersection : ' +str(numberOf1IntersectionT)

#Creation of obs file points
i=0
for obsfilepath in ObsnodeFile:  #For each group station ...
   i=i+1
   x,y,z=np.loadtxt(NodesGrpDir+'/'+obsfilepath,unpack=True) # load associate file
   a="{0:04d}".format(i) 
   obsfile=open(outputDirPts+'/station_grp_'+a+'.txt','w+') # Create group station
   obsfile.write('points\n')
   for j in range(0,len(x)) : #For each point ...
      x1="{:12.6f}".format(x[j]) #Adapted format  
      x2="{:12.6f}".format(y[j]) #Adapted format
      obsfile.write(x1+'   '+x2+'\n') # ... write coor


#Modif
#numerotationCell= '../SWOT_like_data/inputs/result.shp' #Path of output directory of SWOT_like_data groups of nodes
#dsCell     = driver.Open(numerotationCell)
#datanodes = driver.Open(finvect, 0)
#layernodes = datanodes.GetLayer()
#nb_feat_layernodes = layernodes.GetFeatureCount()
#(lonmin_nds, lonmax_nds, latmin_nds, latmax_nds) = layernodes.GetExtent()
#featurenode = layernodes.GetNextFeature()
#      nbnodes_interswotb = 0
#      lon_nds = []
#      lat_nds = []
#      swotb_nds = []
#      while featurenode:
#   	fnodegeom = featurenode.GetGeometryRef()
#	  lon_nds = np.append(lon_nds,fnodegeom.GetX())
#	  lat_nds = np.append(lat_nds,fnodegeom.GetY())


#Fin modif


#Creation of obs file indexes
i=0
for obsfilepath in ObsnodeFile:  #For each group station ...
   i=i+1
   x,y,z=np.loadtxt(NodesGrpDir+'/'+obsfilepath,unpack=True) # load associate file
   a="{0:04d}".format(i) 
   obsfile=open(outputDirIndex+'/station_grp_'+a+'.txt','w+') # Create group station
   obsfile.write('indexes\n')
   for j in range(0,len(x)) : #For each point ...
      x1="{:12.6f}".format(x[j]) #Adapted format  
      x2="{:12.6f}".format(y[j]) #Adapted format
      obsfile.write(x1+'   '+x2+'\n') # ... write coor


      
#Creation of obs_param file
obs_paramfile=open(outputDir+'/param_obs.txt','w+')  #Creation of file
obs_paramfile.write(str(numberOf1IntersectionT)+'\n') #Number of group stations
 
for i in range(0,numberOf1IntersectionT): #For each group station ...
   obs_paramfile.write(str(lenght[i])+ '  '+str(timesList[i]) + '  '+str(leftCoor[i][0]) +  '   '+str(leftCoor[i][1]) +  '   '+str(rigthCoor[i][0]) +  '   '+str(rigthCoor[i][1]) + '\n') # ... write descrition 








