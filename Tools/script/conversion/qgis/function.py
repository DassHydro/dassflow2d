import os 

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


def reprocessing(resultFile):
   fileInput=open(resultFile,'r')
   fileOutput=open(resultFile+'_temp','w')
   Lines=fileInput.readlines()
   for i in range(len(Lines)):
      Line=Lines[i]
      Line=Line.replace('-',' -')
      Line=Line.replace('E -','E-')
      fileOutput.write(Line)
    
   fileOutput.close()
   fileInput.close()
   os.remove(resultFile)
   os.rename(resultFile+'_temp', resultFile)


def getNumeroResultFile(resultFile):
   numero=resultFile
   numero=numero.split('/')
   numero=numero[-1]
   numero=str(numero[7:-4])
   return numero


def cleanString( string_i):
   string_i=string_i.replace('\n','')
   string_i=string_i.replace('\t','')
   string_i=string_i.replace('  ','')
   return string_i

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
      Xnode.append(float(line[1]))
      Ynode.append(float(line[2]))

   # Get connectivity
   connect1, connect2, connect3=[],[],[]
   gap= 3 +  numberOfNode +  1 #  comment 1 and 2 + nbr node and cell+ number of nodes+ comment 3
   print( 'Read connectivity ...')
   for i in range(numberOfCell):
      line=Lines[gap+i]
      line=cleanLine(line)
      line=line.split(' ')
      connect1.append(int(int(line[1])-1))
      connect2.append(int(int(line[2])-1))
      connect3.append(int(int(line[3])-1))

   # Connectivity array
   Connect=[connect1,connect2,connect3]
   Connect.append(connect1)

   return [numberOfNode,numberOfCell,Xnode,Ynode,Connect]


def ComputeBarycenter(Xnode,Ynode,Connect):
   numberOfCell=len(Connect[0])
   if Connect[0]==Connect[3]:
      s=3
   else :
      s=4

   xT,yT=[],[]
   for i in range(numberOfCell):
      x,y=0,0
      for j in range(s):
         connect=Connect[j]

         x=x+Xnode[connect[i]]
         y=y+Ynode[connect[i]]

      xT.append(float(x)/s)
      yT.append(float(y)/s)

   return [xT,yT]

def cellsObservedToPointsObserved(Cell,Xnode,Ynode,Connect):

   X,Y=[],[]

   x,y=ComputeBarycenter(Xnode,Ynode,Connect)


   # For each cell observed
   for i in range(len(Cell)):

      #Compute the barycenter
      cell=int(Cell[i])-1

   
      # Append x,y of the barycenter
      X.append(x[cell])
      Y.append(y[cell])

   return [X,Y]

