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



