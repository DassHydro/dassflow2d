#!/usr/bin/env python
# encoding: utf-8
import os
import numpy as np
import matplotlib.pyplot as plt




bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"
# TO change
pathPlot=f'{bin_dir}/plot/'
pathOutFile=pathPlot+'./min/'


pathFile=f'{bin_dir}/min/min_cost.txt'


 #   z_b_first_guess  z_b_target
ite,cost,grad_cost, grad_cost2=np.loadtxt(pathFile,unpack=True)

#Creation of outputfile 
if (os.path.isdir(pathPlot)==False):
   os.mkdir(pathPlot)

if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)



plt.plot(ite,cost     ,'b',linewidth=1,label=r'$J_{hy}$')
plt.plot(ite,grad_cost,'r',linewidth=1,label=r'$||\nabla_{Q_1} J_{hy}||$')
#plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')

plt.xlabel('Iterations')
plt.legend()
#plt.xscale('log')
plt.yscale('log')

#titlefig = r'Evolution of $J$ and $||\nabla J||$ in function of iterations.'
titlefig = r'$J_{hy}$ and $||\nabla_Q J_{hy}||$ vs iterations.' 
plt.title(titlefig)

pathfile=pathOutFile+'min_cost.png'
plt.savefig(pathfile)
os.system('eog '+str(pathfile) + '&' ) 




#####
# J, GradJ pour SMASH 
ite =  [0,1,2,3,4,5,6,7,8,9,10,11]
J=    [1,  0.212184 ,0.208835 , 0.016089, 0.006804, 0.000117, 0.000013, 0.000012,0.000012, 0.000008,0.000004,0    ]
gradJ=[1,   0.000189       , 0.000189, 0.000084, 0.000048, 0.000006, 0.000002, 0.000001,0.000001,0.000001, 0.000001,0 ]

ite = np.asarray(ite)
J = np.asarray(J)
gradJ = np.asarray(gradJ)



plt.plot(ite,J     ,'b',linewidth=1,label=r'$J_{rr}$')
plt.plot(ite,gradJ,'r',linewidth=1,label=r'$||\nabla J_{rr}||$')
#plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')

plt.xlabel('Iterations')
plt.legend()
#plt.xscale('log')
plt.yscale('log')

#titlefig = r'Evolution of $J$ and $||\nabla J||$ in function of iterations.'
titlefig = r'$J_{rr}$ and $||\nabla J_{rr}||$ vs iterations.' 
plt.title(titlefig)

pathfile=pathOutFile+'min_cost_rr.png'
plt.savefig(pathfile)
os.system('eog '+str(pathfile) + '&' ) 


###
# avoir la norme 2 du vecteur de controle  

# a) avoir le vecteur de controle true : qtrue smash
# b) avoir les hydrogrames inférés dassflow /min/hydrogram.txt
# c) calculer  la norme 2 des 2 vecteurs

q_true_smash = np.asarray([3.38090814e-07, 6.67026557e-04, 4.65896772e-03, 2.36780345e-02,
        7.84027725e-02, 2.10827976e-01, 4.86194402e-01, 9.97156918e-01,
        1.86506212e+00, 3.24599671e+00, 4.17270851e+00, 5.18576813e+00,
        6.87647581e+00, 9.74914932e+00, 1.31693268e+01, 1.53667479e+01,
        1.50729065e+01, 1.24433346e+01, 8.46863270e+00, 6.72360468e+00,
        5.53490877e+00, 4.67793989e+00, 4.03365040e+00, 3.53335285e+00])

print(dassflow_dir)
bin_dir = f"{dassflow_dir}/code/bin_A"
min_dir = f"{bin_dir}/min"

files = os.listdir(min_dir)
q1_file = [x for x in files if  "hydrograph_001" in x]
q2_file = [x for x in files if  "hydrograph_002" in x]

id_q1 = [int(x[15:18]) for x in q1_file]
id_q2 = [int(x[15:18]) for x in q2_file]


res = dict()

res["q1"] =  dict()
res["q2"] =  dict()


for i in range(len(id_q1)):
    q = np.loadtxt(f"{min_dir}/{q1_file[i]}")
    q[np.where(q[:,1] <0),1]=0
    res["q1"][f"{id_q1[i]}"] = q
    
for i in range(len(id_q2)):
    q = np.loadtxt(f"{min_dir}/{q2_file[i]}")
    res["q2"][f"{id_q2[i]}"] = q
    
norm =[]
for i in range(len(res["q1"])):
    norm.append(np.linalg.norm(q_true_smash - res["q1"][f"{id_q1[i]}"][:,1]))

plt.plot(res["q1"]["9"][:,1], label = "qinf")
plt.plot(q_true_smash, label = "qtrue")
plt.legend()
plt.show()


plt.plot(res["q2"]["9"][:,1], label = "qinf")
plt.plot(q_true_smash, label = "qtrue")
plt.legend()
plt.show()