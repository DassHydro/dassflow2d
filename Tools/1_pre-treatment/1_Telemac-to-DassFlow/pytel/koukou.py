
#
# Transform a telemac mesh into a Dassflow mesh


#!/usr/bin/env python
import sys
sys.path.append( "/home/lilian.villenave/Documents/save/GIT/hydrological-modeling/2021-Lilian-VILLENAVE/5_DassFlow_prise-en-main/2_scripts/python3" )

#from data_manip.formats.selafin import Selafin # conlim

from utils.files import is_newer, put_file_content

import numpy as np


####
# To change
mesh = SELAFIN('bathy_test1_v2.slf')                   # path of telemac mesh file
cond = CONLIM ('bathy_test1_v2.bc2')                 # path of boundary condition telemac mesh file
DassFlowMesh = open('dassflox_test1.geo','w')   # path of new dassflow mesh file

########## -------
# End to change

print 'Mesh title : ' + str(mesh.TITLE)
print 'Number points   : ' + str(mesh.NPOIN3)
print 'Number elements : ' + str(mesh.NELEM3)

x     =np.copy(mesh.MESHX)
y     =np.copy(mesh.MESHY)
zb_n  =np.copy(mesh.getVALUES(0)[0])

ele   =np.copy(mesh.IKLE)
zb_c  =[]
DassFlowMesh.write('# Converted mesh : ' +str(mesh.TITLE)+'\n')
DassFlowMesh.write(str(mesh.NPOIN3)+ ' ' +str(mesh.NELEM3)+ ' '+ str(1.0000)+'\n')
DassFlowMesh.write('# Nodes '+'\n')

for i in range(mesh.NPOIN3):
   DassFlowMesh.write(str(i+1) + ' '+str(x[i])+ ' ' +str(y[i])+ ' '+ str(zb_n[i])+'\n')

DassFlowMesh.write('# Cells '+'\n')

for i in range(mesh.NELEM3):
   zb=(1./3)*(zb_n[ele[i][0]]+zb_n[ele[i][1]]+zb_n[ele[i][2]])
   DassFlowMesh.write(str(i+1)+ ' ' + str(ele[i][0]+1) + ' ' + str(ele[i][1]+1) + ' ' + str(ele[i][2]+1)+ ' ' + str(ele[i][0]+1)+ ' 1 ' +  str(zb) + '\n')

DassFlowMesh.write('# Boundaries '+'\n')


#OUTPUT
nbrInput=0
BorderPoint=[]
for i in range(cond.NPTFR) :
   if ((cond.BOR[i][0]==5)and (cond.BOR[i][1]==4)) :
      nbrInput=nbrInput+1
      BorderPoint.append(cond.BOR[i][11])

print 'Number output border point :' +str(len(BorderPoint))

ElemNodeT=[]
for i in range(len(BorderPoint)):
   ElemNode=np.where(ele==BorderPoint[i]-1)
   for j in range(len(ElemNode[0])):
      ElemNodeT.append(ElemNode[0][j])

ElemNode=[]
for i in  range(len(ElemNodeT)):
   if (ElemNodeT.count(ElemNodeT[i])>1 and ElemNodeT[i] not in ElemNode):
      ElemNode.append(ElemNodeT[i])

print 'Number input border cell : ' + str(len(ElemNode))

DassFlowMesh.write('OUTLET ' +str(len(ElemNode))+ ' 0 \n')


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
   DassFlowMesh.write(str(ElemNode[i]+1)+ ' ' +str(connectivityBorder)+' ' +'1 '+ str(zb) +'\n')




#INPUT
nbrInput=0
BorderPoint=[]
for i in range(cond.NPTFR) :
   if ((cond.BOR[i][0]==4)and (cond.BOR[i][1]==5)) :
      nbrInput=nbrInput+1
      BorderPoint.append(cond.BOR[i][11])

print 'Number input border point :' +str(len(BorderPoint))

ElemNodeT=[]
for i in range(len(BorderPoint)):
   ElemNode=np.where(ele==BorderPoint[i]-1)
   for j in range(len(ElemNode[0])):
      ElemNodeT.append(ElemNode[0][j])

ElemNode=[]
for i in  range(len(ElemNodeT)):
   if (ElemNodeT.count(ElemNodeT[i])>1 and ElemNodeT[i] not in ElemNode):
      ElemNode.append(ElemNodeT[i])

print 'Number input border cell : ' + str(len(ElemNode))

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
