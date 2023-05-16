import numpy as np
import matplotlib.pylab as plt


def hydro_reel(t,tend) :
   h=10.0
   Ttotal=4*3600

   wi=[1./Ttotal,5./Ttotal,10./Ttotal,20./Ttotal]
   a=[8,2,0.5,0.1]
   phi=[1.5*np.pi,0,0,0]
   for i in range(len(wi)):
      h=h+a[i]*np.sin(t*wi[i]*np.pi*2.0+phi[i])
   return abs(h)

def hydro_aca(t,tend) :
   Qn=100
   # Original gaussian 1100
   #h=Qn*(1+10*np.exp(-np.power(t-tend/2,2)/np.power(50000,2)))
   h=Qn*(1+2*np.exp(-np.power(t-tend/2,2)/np.power(50000,2)))
   return abs(h)


def hydro(t,tend) :
   return hydro_aca(t,tend)


def hydro1stGuess(t) :

   return 100.0


factorTime=1./(3600)





tstart=0.0
tend  =24*4*3600.0
step  =241
stepTrue=241
stepFirstGuess=241
nameHydro='hydrograph.txt'
Hydro=[]
linspace=np.linspace(float(tstart),float(tend),step)#
nameHydro         ='hydrograph.txt'
nameHydroTarget   ='hydrograph_target.txt'
nameHydro1stGuess ='hydrograph_first_guess.txt'

linspace=np.linspace(float(tstart),float(tend),stepTrue)
linspace1stGuess=np.linspace(float(tstart),float(tend),stepFirstGuess)
#print linspace
#linspace1stGuess=[0.0,4.320191187752307E+04,  8.640210869660787E+04,  1.296022682942512E+05,  1.728024278915321E+05,2.160025874888131E+05,2.592027470860940E+05,tend-60.0]
#linspace1stGuess=[0.0,45697.400000000001,47796.699999999997, 74209.399999999994,76308.699999999997, 102721.39999999999,104820.7, 131233.39999999999,133332.70000000001, 159745.39999999999,161844.70000000001, 188257.39999999999,190356.70000000001, 216769.39999999999,218868.70000000001, 245281.39999999999, 247380.70000000001,273793.40000000002, 275892.70000000001,3*24*3600]

#linspace1stGuess=[0,17185.400000000001, 19284.700000000001, 45697.4, 47796.7, 74209.4, 76308.7, 102721.4, 104820.7, 131233.4, 133332.7, 159745.4, 161844.7, 188257.4, 190356.7, 216769.4, 218868.7, 245281.4, 247380.7,3*24*3600-10]

Hydro=[]
file=open(nameHydro,'w+')
file.write('!===============================================!\n')
file.write('!  Number of hydrograph                         !\n')
file.write('!===============================================!\n')
file.write('1                                                \n')
file.write('!===============================================!\n')
file.write('!  Hydrograph                                   !\n')
file.write('!===============================================!\n')
file.write('      '+str(len(linspace))+'\n'                     )
for i in linspace:
   file.write(str(i) + ' '+str(hydro(i,tend))+'\n')
   Hydro.append(hydro(i,tend))
   
plt.plot(linspace*factorTime,Hydro)
plt.ylabel('$Q$ ($m^3.s^{-1}$)')
plt.xlabel('$Times$ ($hours$)')
ymin=min(Hydro)
ymax=max(Hydro)
delta=ymax-ymin
plt.ylim(ymin-0.1*delta,ymax+0.1*delta)
plt.xlim([linspace[0]*factorTime,linspace[-1]*factorTime])
plt.savefig('hydro.png')


file=open(nameHydroTarget,'w+')
file.write('!===============================================!\n')
file.write('!  Number of hydrograph                         !\n')
file.write('!===============================================!\n')
file.write('1                                                \n')
file.write('!===============================================!\n')
file.write('!  Hydrograph                                   !\n')
file.write('!===============================================!\n')
file.write('      '+str(len(linspace))+'\n'                     )
for i in linspace:
   file.write(str(i) + ' '+str(hydro(i,tend))+'\n')
   Hydro.append(hydro(i,tend))




Hydro1stGuess=[]
file=open(nameHydro1stGuess,'w+')
file.write('!===============================================!\n')
file.write('!  Number of hydrograph                         !\n')
file.write('!===============================================!\n')
file.write('1                                                \n')
file.write('!===============================================!\n')
file.write('!  Hydrograph                                   !\n')
file.write('!===============================================!\n')
file.write('      '+str(len(linspace1stGuess))+'\n'                     )
for i in linspace1stGuess:
   file.write(str(i) + ' '+str(hydro1stGuess(i))+'\n')
   Hydro1stGuess.append(hydro1stGuess(i))
   

plt.plot(linspace1stGuess*factorTime,Hydro1stGuess)
plt.savefig('hydro_true_1stguess.png')
plt.clf()

plt.plot(linspace1stGuess*factorTime,Hydro1stGuess)
plt.savefig('hydro_1stguess.png')
plt.clf()




