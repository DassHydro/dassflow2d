import os
import numpy as np


def ReadDassFlowResults(pathResultFile):
   print('Read DassFlow result  : '+str(pathResultFile))
   try :
      i, x,y,bathy,h,zs,Manning,u,v=np.loadtxt(pathResultFile,unpack=True)
   except :
      pathResultFilePost=pathResultFile+'_post'
      print( 'Error file result parsing ...')
      if os.path.isfile(pathResultFilePost):
         print ('File : '+ str(pathResultFilePost) + ' exists already')
         i, x,y,bathy,h,zs,Manning,u,v=np.loadtxt(pathResultFilePost,unpack=True)
      else:
         print( 'Creation of  '+ str(pathResultFilePost) + ' file')
         fileInput=open(pathResultFile,'r')
         fileOutput=open(pathResultFilePost,'w')
         Lines=fileInput.readlines()
         for i in range(len(Lines)):
            Line=Lines[i]
            Line=Line.replace('-',' -')
            Line=Line.replace('E -','E-')
            fileOutput.write(Line)
         fileOutput.close()
         i, x,y,bathy,h,zs,Manning,u,v=np.loadtxt(pathResultFilePost,unpack=True)

   return [i, x,y,bathy,h,zs,Manning,u,v]


def cleanLine(line):
   line=line.replace('\n','')
   line=line.replace('\t','')
   return line

def cleanLineSplit(line):
   count=line.count('')
   if count!=0:
      for i in range(count):
         line.remove('')
   return line

def ReadDassFlowMesh(pathInputMesh):
   InputMesh=open(pathInputMesh,'r')
   Lines=InputMesh.readlines()

   # Get number of node and number of cell
   line=Lines[1]
   line=cleanLine(line)
   line=line.split(' ')
   cleanLineSplit(line)
   numberOfNode=int(line[0])
   numberOfCell=int(line[1])
   print( 'DassFlow mesh : ' + pathInputMesh)
   print( '\t -'+str(numberOfNode) + ' nodes')
   print( '\t -'+str(numberOfCell) + ' cells')
   print( '')

   # Get nodes
   Xnode,Ynode=[],[]
   gap=  1 + 1 + 1 # Comment 1+nbr node, nbr cell, scale+Comment 2

   print( 'Read nodes ...')
   for i in range(numberOfNode):
      line=Lines[gap+i]
      line=cleanLine(line)
      line=line.split(' ')
      cleanLineSplit(line)
      Xnode.append(float(line[1]))
      Ynode.append(float(line[2]))

   # Get connectivity
   connect1, connect2, connect3, connect4=[],[],[],[]
   bathy=[]
   gap= 3 +  numberOfNode +  1 #  comment 1 and 2 + nbr node and cell+ number of nodes+ comment 3
   print( 'Read connectivity ...')
   for i in range(numberOfCell):
      line=Lines[gap+i]
      line=cleanLine(line)
      line=line.split(' ')
      cleanLineSplit(line)
      connect1.append(int(int(line[1])-1))
      connect2.append(int(int(line[2])-1))
      connect3.append(int(int(line[3])-1))
      connect4.append(int(int(line[4])-1))
      bathy   .append(float(line[6]))

   # Connectivity array
   Connect=[connect1,connect2,connect3,connect4]
   Connect.append(connect1)

   return [numberOfNode,numberOfCell,Xnode,Ynode,Connect,bathy]

def ComputeBarycenter(Xnode,Ynode,Connect):
   numberOfCell=len(Connect[0])
   if Connect[0]==Connect[3]:
      s=3
   else :
      s=4

   xT,yT=[],[]
   print( s)
   for i in range(numberOfCell):
      x,y=0,0
      for j in range(s):
         connect=Connect[j]

         x=x+Xnode[connect[i]]
         y=y+Ynode[connect[i]]

      xT.append(float(x)/s)
      yT.append(float(y)/s)

   return [xT,yT]

def distance2Points(pts1,pts2):
   return np.sqrt(np.power(pts2[0]-pts1[0],2)+np.power(pts2[1]-pts1[1],2))


def castBathyOnNode(Connect,bathy,Xnode,Ynode):

   numberOfNode=len(Xnode)
   numberOfCell=len(Connect[0])

   PointsCells    = [list([]) for _ in range(numberOfNode)]
   PointsCellsDis = [list([]) for _ in range(numberOfNode)]

   x,y=ComputeBarycenter(Xnode,Ynode,Connect)


   if Connect[0]==Connect[3]:
      s=3
   else :
      s=4

   ranges=range(s)
   for i in range(numberOfCell): #For each cell
      for j in ranges: #For each points of the cell i
         
         connect=Connect[j] #Get correct connectivity array
         # Node
         node=connect[i] 
         # X,Y node
         xnode=Xnode[node]
         ynode=Ynode[node]
         # x,y center
         xcenter=x[i]
         ycenter=y[i]
         value=bathy[i]
         distance=distance2Points([xnode,ynode],[xcenter,ycenter]) #Distance between pts and center cell


         #Save value of center to the array PointsCell
         temp=PointsCells[node]
         temp.append(value)
         PointsCells[node]=temp
         
         #Save distance pts-center cell, to the PointsCell array
         temp=PointsCellsDis[node]
         temp.append(distance)
         PointsCellsDis[node]=temp

   
   #
   # Computes value at the node
   #
   ValueNode=[]
   print( 'Compute nodes values ...')
   for i in range(numberOfNode):
      pointsCells    = PointsCells[i]
      pointsCellsDis = PointsCellsDis[i]

      #print pointsCellsDis
      mean=0
      for j in range(len(pointsCells)):
         mean=mean+pointsCellsDis[j]*pointsCells[j]
      mean=mean/sum(pointsCellsDis)
      ValueNode.append(mean)
   return ValueNode




def writeVariableVTK(file,nameVariable,variable):
   
   file.write('SCALARS         '+str(nameVariable)+'         double 1       \n')
   file.write('LOOKUP_TABLE default           \n')
   for i in range(len(variable)):
      file.write('{:15.8E}'.format(variable[i])+'\n')
   #file.write(
#variable



