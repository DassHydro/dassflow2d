import ogr, os, sys
import numpy as np
from osgeo import ogr,osr # GDAL library
import shutil 

# See readme.txt 

#To change
SWOTBandDir= '../SWOT_like_data/orbit891km_Garonne_SWOTswathbands1km/orbit891km_Garonne_SWOTswathbands1km' #Path of output directory of SWOT_like_data SWOT bands
centerPoints=  './inputs/centerPoints.shp' #Path of center line
outputDir='./ouputs/'
#End to change

SWOTBandFiles=[]           #SWOT Band files
numberOf1IntersectionT=0   #Number of swot band which intersects only one time center line

#
# Function which do groups contiguous from a array none continuguous
# x [2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 22, 25, 26, 28, 51, 52]
# return [[2, 3, 4, 5], [12, 13, 14, 15, 16, 17], [22], [25, 26], [28], [51, 52]]
def group(L):
   l=L[0]
   Lsort=[]
   localLsort=[l]
   for i in range(1,len(L)):
      l=L[i]
      if (l)==(localLsort[-1]+1):
         localLsort.append(l)
      else:
         Lsort.append(localLsort)
         localLsort=[l]            
   Lsort.append(localLsort)
   return Lsort



#
# Remove and creation of ouput directory
# 
outputDirIndex=outputDir+'/indexes'

if os.path.isdir(outputDir)==False:
	os.mkdir(outputDir)

if os.path.isdir(outputDirIndex)==False:
	os.mkdir(outputDirIndex)



#
# Search of SWOT reach 1D
#


lenght,leftCoor,rigthCoor,Bandnumber,OutDistance,ObsnodeFile,timesList=[],[],[],[],[],[],[]

driver = ogr.GetDriverByName('ESRI Shapefile')

listFile=os.listdir(SWOTBandDir) #List of swot band

for fileB in listFile: #Selection of useful swot band (l93 band)
   if (fileB.find('s_l93.shp')!=-1):
      SWOTBandFiles.append('/'+fileB) 


dsCenterPoints = driver.Open(centerPoints) #Load centre cell points

numberOf1IntersectionT,numberOfIntersectionT,numberTotalReach=0,0,0
pointsObserved,pointsObservedT,timesList=[],[],[]


CenterPoint    = dsCenterPoints.GetLayerByIndex(0)
numberCenterPoint=CenterPoint.GetFeatureCount(0)


for SWOTBandFile in SWOTBandFiles: #For each Band
   
   # Get time 
   times=SWOTBandFile[fileB.find('km_time')+8:fileB.find('s_l93.shp')-8]
   times=str(np.round(float(times)))
   print 'Time :' + times

   # Get Band
   dsSWOTBand     = driver.Open(SWOTBandDir+SWOTBandFile)
   SWOTBand       = dsSWOTBand.GetLayerByIndex(0)

   NumberOfBand=SWOTBand.GetFeatureCount() # Number of under band
   
   for i in range(0,NumberOfBand): #For each under band of swot band

      #Get under band
      a     =SWOTBand.GetFeature(i)
      Band  =a.GetGeometryRef()
      outDis=a.GetField('Out dist')


      IndexesPointsInBand=[] # Initialization of points array

      for j in range(numberCenterPoint): # For each center cell point
         
         #Get center cell point j
         c    =CenterPoint.GetFeature(j)
         Points  =c.GetGeometryRef()
        
         Contains=Band.Contains(Points) # True is Points in the under SWOT band Band

         
         if Contains==True:
            IndexesPointsInBand.append(j+1) #append to the points arrays


      if len(IndexesPointsInBand)!=0 : # If there are at last one groups of cell observed
         
         IndexesPointsInBandSort=group(IndexesPointsInBand)  # Do groups from the array 
         numberOfGroups=len(IndexesPointsInBandSort) 
         numberOfIntersectionT=numberOfIntersectionT+1
         numberTotalReach=numberTotalReach+numberOfGroups

         if numberOfGroups==1: # If there are only one groups of cell observed :
            numberOf1IntersectionT=numberOf1IntersectionT+1
            pointsObserved.append(IndexesPointsInBandSort[0]) # Save the list
            timesList.append(times)
            print 'Band ' +str(i) + ', outdis '+str(outDis)+' intersects river ' +str(numberOfGroups) + ' times.'
            
         if numberOfGroups>1: # If there are more of one groups of cell observed
            pass
         pointsObservedT.append(IndexesPointsInBandSort) # Savec all list



#
# Computing percentage cell observed
#
percentagePtsObsArray =[0]*(numberCenterPoint)
percentagePtsObsTArray=[0]*(numberCenterPoint)

# Percentage of cell observed directly usable (one under band -> one reach)
for tab in pointsObservedT:
   for tab2 in tab:
      for i in tab2:
         if percentagePtsObsTArray[i-1]==0:
            percentagePtsObsTArray[i-1]=1
percentagePtsObsT=float(sum(percentagePtsObsTArray ))/numberCenterPoint


# Percentage of cell observed 
for tab in pointsObserved:
   for i in tab:
      if percentagePtsObsArray[i-1]==0:
         percentagePtsObsArray[i-1]=1
percentagePtsObs=float(sum(percentagePtsObsArray ))/numberCenterPoint

#
# Print stats
#
print 'Number of SWOT Band intersects : ' + str(numberOfIntersectionT)
print 'Number of SWOT Band with 1 intersection : ' +str(numberOf1IntersectionT)
print 'Number total of reaches : ' +str(numberTotalReach)
print 'Percentage cell observed directly usable : ' + str(percentagePtsObs*100) +' % '
print 'Percentage cell observed : ' + str(percentagePtsObsT*100)+ ' %'

#     
# Creation of shapefile
#
print 'Creation of shapefile ...' 
#Creation of output directory shapefile
outputDirectoryShapefile=outputDir+'shapefile'

if os.path.exists(outputDirectoryShapefile):
   shutil.rmtree(outputDirectoryShapefile)
   os.mkdir(outputDirectoryShapefile)




for i in range(len(pointsObserved)): 
   EPGScode=2154
   nameFile='obs_'+str(i+1)

   # Creation of new shapefile
   driver = ogr.GetDriverByName('ESRI Shapefile')
   fshpout=outputDirectoryShapefile#+str(nameFile)
	
   dataout = driver.CreateDataSource(fshpout)
   proj = osr.SpatialReference()

   proj.ImportFromEPSG(EPGScode) # 4326 = EPSG code for lon/lat WGS84 projection

   
   #Trace node
   layerout = dataout.CreateLayer(str(nameFile),proj,geom_type=ogr.wkbPoint)  
   fieldDef1 = ogr.FieldDefn('x', ogr.OFTReal)
   fieldDef2 = ogr.FieldDefn('y', ogr.OFTReal)
 
   layerout.CreateField(fieldDef1)
   layerout.CreateField(fieldDef2)
 

   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)

   
   points=pointsObserved[i]
   for j in range(0,len(points)):
      c    =CenterPoint.GetFeature(points[j]-1)
      a    =c.geometry()
      coord=a.GetPoint()
      x=coord[0]
      y=coord[1]
      #- Create Geometry Point with pixel coordinates
      pixel_point = ogr.Geometry(ogr.wkbPoint)
      pixel_point.AddPoint(x,y)
      #- Add the geometry to the feature
      feature_out.SetGeometry(pixel_point)
      #- Set feature attributes
      feature_out.SetField('x'  , x)
      feature_out.SetField('y'  , y)
      layerout.CreateFeature(feature_out)
      #- Delete point geometry
      pixel_point.Destroy()


   feature_out.Destroy()
   dataout.Destroy()

#
# Creation of obs file indexes
#
print 'Creation of station_grp_XXX.txt files ...'
for i in range(len(pointsObserved)):  # For each groups of cell
   a="{0:04d}".format(i+1) 
   obsfile=open(outputDirIndex+'/station_grp_'+a+'.txt','w+') # Create group station
   obsfile.write('indexes\n')
   points=pointsObserved[i]
   for j in range(0,len(points)) : #For each point ...
      pts=int(points[j]) #"{:12.6f}".format(points[j]) #Adapted format  
      obsfile.write(str(pts)+'\n') # ... write coor

#
# Creation of obs_param file
#
print 'Creation of param_obs.txt files ...'
obs_paramfile=open(outputDirIndex+'/param_obs.txt','w+')  #Creation of file
obs_paramfile.write(str(numberOf1IntersectionT)+'\n') #Number of group stations
 
for i in range(0,numberOf1IntersectionT): #For each group station ...
   obs_paramfile.write(str(i+1)+ '  '+str(timesList[i]) + '\n') # ... write 
