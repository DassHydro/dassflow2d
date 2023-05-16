import os
os.chdir("/home/livillenave/Documents/distant/svn/dassflow-2d/trunk/tools/conversion/dat_to_vtk/")
from functions import *
import numpy as np

# Transform a .vtk file into a .dat file


bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"
# TO change
meshFileName   =f'{bin_dir}/automaticaly_generated_mesh.txt'
inputFileName  =f'{bin_dir}/res/result_final.dat'
outputFileName  =f'{bin_dir}/res/result_final.vtk'
# End to change

#
# Load result file
#
print( '')
print( '------------------------------------------------')
print( '        Read DassFlow results file              ')
print( '------------------------------------------------')
[i, x,y,bathy,h,zs,Manning,u,v]=ReadDassFlowResults(inputFileName)
#x,y,bathy,h,zs,Manning,u,v=va[0],va[1],va[2],va[3],va[4]


#
# Read DassFlow mesh
#
print( '')
print( '------------------------------------------------')
print( '        Read DassFlow mesh              ')
print( '------------------------------------------------')
values=ReadDassFlowMesh(meshFileName)
numberOfNode,numberOfCell,Xnode,Ynode,Connect,bathy2=values[0],values[1],values[2],values[3],values[4],values[5]
connect1,connect2,connect3,connect4=Connect[0],Connect[1],Connect[2],Connect[3]


Znode=castBathyOnNode(Connect,bathy,Xnode,Ynode)


#
# write VTK
#
outputFile=open(outputFileName,'w+')
outputFile.write('# vtk DataFile Version 3.0     \n')
outputFile.write('DassFlow Output File           \n')
outputFile.write('ASCII          \n')
outputFile.write('DATASET UNSTRUCTURED_GRID      \n')

# Write points
outputFile.write('POINTS          '+str('{:16d}'.format(numberOfNode))+' double        \n')
for i in range(numberOfNode):
   outputFile.write(str('{:15.8E}'.format(Xnode[i]))+' '+str('{:15.8E}'.format(Ynode[i]))+' '+str('{:15.8E}'.format(Znode[i]))+'\n')

# Write cells
outputFile.write(' CELLS         '+str('{:16d}'.format(numberOfCell))+str('{:16d}'.format(numberOfCell*5))+'\n')

# Quadrangle
if connect1[0]!=connect4[0]:
   for i in range(numberOfCell):
      outputFile.write(str('{:15d}'.format(4))+str('{:16d}'.format(connect1[i]))+str('{:16d}'.format(connect2[i]))+str('{:16d}'.format(connect3[i]))+str('{:16d}'.format(connect4[i]))+'\n')

# Triangle
else:
   for i in range(numberOfCell):
      outputFile.write(str('{:15d}'.format(3))+str('{:16d}'.format(connect1[i]))+str('{:16d}'.format(connect2[i]))+str('{:16d}'.format(connect3[i]))+'\n')

outputFile.write('CELL_TYPES     '+str('{:16d}'.format(numberOfCell))+'\n')
for i in range(numberOfCell):
   outputFile.write(str('{:15d}'.format(9))+'\n')

outputFile.write('CELL_DATA                   '+str(numberOfCell)+'\n')
writeVariableVTK(outputFile,'u      ',u)
writeVariableVTK(outputFile,'v      ',v)
writeVariableVTK(outputFile,'bathy  ',bathy)
outputFile.close()
outputFile=open(outputFileName,'a')
writeVariableVTK(outputFile,'h      ',h)
writeVariableVTK(outputFile,'zs     ',zs)
writeVariableVTK(outputFile,'n',Manning)
outputFile.close()
#writeVariableVTK(+str('{:16d}'.format(connect1[i]))
   #print( 'ok'
#str("{:13.7E}".format(h)
#I16
#outputFile.write('mesh%nn
#outputFile.write(10,rec=10,fmt='(A15,A1)') ' double        ', char(10)
