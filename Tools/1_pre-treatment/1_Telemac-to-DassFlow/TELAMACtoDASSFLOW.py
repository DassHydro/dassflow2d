import os
############
# pytel libraires directory
lib_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/PRE-TREATMENT/Telemac-to-DassFlow/pytel"
bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/PRE-TREATMENT/Telemac-to-DassFlow/example_bin"

mesh_name = '/Cartino2D.slf'
cli_name =  '/Cartino2D.cli'

print(os.path.isdir(bin_dir))
print(os.path.isfile(bin_dir + '/Cartino2D.slf'))
print("-----------------")
# work in bin directory
cwd = os.getcwd()
print(cwd)
os.chdir(bin_dir)

cwd = os.getcwd()
print(cwd)
###########
import sys
import numpy as np

sys.path.insert(0, lib_dir) # To adapt (get python parser of telemac)
sys.path.insert(0, lib_dir + "utils" ) # To adapt (get python parser of telemac)
sys.path.insert(0, lib_dir + "data_manip" ) # To adapt (get python parser of telemac)
print(sys.path)

# # source libraries
from utils.files import put_file_content
from data_manip.formats.selafin import Selafin
from data_manip.formats.conlim import Conlim


# 
# # Transform a telemac mesh into a Dassflow mesh
# 


# To change
mesh = Selafin(bin_dir + mesh_name )               # path of telemac mesh file
cond = Conlim (bin_dir + cli_name)                 # path of boundary condition telemac mesh file
DassFlowMesh = open(bin_dir + '/mesh.geo','w')   # path of new dassflow mesh file
# End to change


print('Mesh title : ' + str(mesh.title))
print('Number points   : ' + str(mesh.npoin3))
print('Number elements : ' + str(mesh.nelem3))

x     =np.copy(mesh.meshx)
y     =np.copy(mesh.meshy)
zb_n  =np.copy(mesh.get_values(0)[0])

ele   =np.copy(mesh.ikle3)
zb_c  =[]
DassFlowMesh.write('# Converted mesh : ' +str(mesh.title)+'\n')
DassFlowMesh.write(str(mesh.npoin3)+ ' ' +str(mesh.nelem3)+ ' '+ str(1.0000)+'\n')
DassFlowMesh.write('# Nodes '+'\n')
for i in range(mesh.npoin3):
   DassFlowMesh.write(str(i+1) + ' '+str(x[i])+ ' ' +str(y[i])+ ' '+ str(zb_n[i])+'\n')

DassFlowMesh.write('# Cells '+'\n')

for i in range(mesh.nelem3):
   zb=(1./3)*(zb_n[ele[i][0]]+zb_n[ele[i][1]]+zb_n[ele[i][2]]) 
   DassFlowMesh.write(str(i+1)+ ' ' + str(ele[i][0]+1) + ' ' + str(ele[i][1]+1) + ' ' + str(ele[i][2]+1)+ ' ' + str(ele[i][0]+1)+ ' 1 ' +  str(zb) + '\n')

DassFlowMesh.write('# Boundaries '+'\n')
#OUTPUT
nbrInput=0
BorderPoint=[]
for i in range(cond.nptfr) :
   if ((cond.bor[i][0]==5)and (cond.bor[i][1]==4)) :
      nbrInput=nbrInput+1
      BorderPoint.append(cond.bor[i][11])

print('Number output border point :' +str(len(BorderPoint)))




ElemNodeT=[]
for i in range(len(BorderPoint)):
   ElemNode=np.where(ele==BorderPoint[i]-1)
   for j in range(len(ElemNode[0])):
      ElemNodeT.append(ElemNode[0][j])

ElemNode=[]
for i in  range(len(ElemNodeT)):
   if (ElemNodeT.count(ElemNodeT[i])>1 and ElemNodeT[i] not in ElemNode):
      ElemNode.append(ElemNodeT[i])

print('Number output border cell : ' + str(len(ElemNode)))

DassFlowMesh.write('OUTLET ' +str(len(ElemNode))+ ' 0 \n')




# len(ele[ElemNode[i]])

for i in range(len(ElemNode)):
   connectivityBorder=0
   for j in range(len(BorderPoint)):
      a=np.where(ele[ElemNode[i]]==BorderPoint[j]-1)[0]
      if len(a)!=0:
         connectivityBorder=connectivityBorder+a[0]
   if connectivityBorder==2:
      connectivityBorder=3
   elif connectivityBorder==3:
      connectivityBorder=2
   zb=(1./3)*(zb_n[ele[ElemNode[i]][0]]+
      zb_n[ele[ElemNode[i]][1]]+
      zb_n[ele[ElemNode[i]][2]]) #To recalculate
   DassFlowMesh.write(str(ElemNode[i]+1)+ ' ' +str(connectivityBorder)+' ' +'1 '+ str(zb) +'\n')
    
    
    
    


#INPUT
nbrInput=0
BorderPoint=[]
for i in range(cond.nptfr) :
   if ((cond.bor[i][0]==4)and (cond.bor[i][1]==5)) :
      nbrInput=nbrInput+1
      BorderPoint.append(cond.bor[i][11])

print('Number input border point :' +str(len(BorderPoint)))

ElemNodeT=[]
for i in range(len(BorderPoint)):
   ElemNode=np.where(ele==BorderPoint[i]-1)
   for j in range(len(ElemNode[0])):
      ElemNodeT.append(ElemNode[0][j])

ElemNode=[]
for i in  range(len(ElemNodeT)):
   if (ElemNodeT.count(ElemNodeT[i])>1 and ElemNodeT[i] not in ElemNode):
      ElemNode.append(ElemNodeT[i])


print('Number input border cell : ' + str(len(ElemNode)))

DassFlowMesh.write('INLET ' +str(len(ElemNode))+ ' 0 \n')

for i in range(len(ElemNode)):
   connectivityBorder=0
   for j in range(len(BorderPoint)):
      a=np.where(ele[ElemNode[i]]==BorderPoint[j]-1)[0]
      if len(a)!=0:
         connectivityBorder=connectivityBorder+a[0]
   if connectivityBorder==2:
      connectivityBorder=3
   elif connectivityBorder==3:
      connectivityBorder=2
   zb=(1./3)*(zb_n[ele[ElemNode[i]][0]]+zb_n[ele[ElemNode[i]][1]]+zb_n[ele[ElemNode[i]][2]]) #To recalculate
   DassFlowMesh.write(str(ElemNode[i]+1)+ ' ' +str(connectivityBorder)+ ' ' + '8 '+str(zb) +'\n')


   #A finir elevation fantomes + OUTPUT sur le input         

