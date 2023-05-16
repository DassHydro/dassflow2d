import numpy as np
import matplotlib.pylab as plt


data=np.loadtxt('data')
Day=data[:,0]
hours=data[:,1]
minu=data[:,2]
Q=data[:,3]

ToConserve=[]
NDay=[Day[0]]
prec=Day[0]
for i in range(1,len(Day)):
   if (Day[i]==prec):
      NDay.append(NDay[-1])
   elif (Day[i]==prec+1):
      NDay.append(NDay[-1]+1)
      prec=prec+1
   else :
      NDay.append(NDay[-1]+1)
      prec=Day[i]


NDay=np.array(NDay)

Day=data[:,0]-data[0,0]
Day=NDay[:]-NDay[0]

Time=3600*24*Day+hours*3600+minu*60
Time=Time-Time[0]
Time=Time/(3600*24)

plt.plot(Time*(3600*24),Q)
plt.title('Evolution of discharge at Toulouse\'s station')
plt.xlabel(r'Time ($s$)')
plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
plt.savefig('hydro_seconds.png')

plt.clf()
plt.plot(Time,Q,label='Q real')
plt.title('Evolution of discharge at Toulouse\'s station')
plt.xlabel(r'Time ($days$)')
plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
plt.savefig('hydro_days.png')


#plt.savefig('hydro.png')



nameHydro         ='hydrograph.txt'
nameHydro1stGuess ='hydrograph_first_guess'
hydro1stGuess     =100.0
linspace1stGuess=[0,17185.400000000001, 19284.700000000001, 45697.4, 47796.7, 74209.4, 76308.7, 102721.4, 104820.7, 131233.4, 133332.7, 159745.4, 161844.7,171780.0-10]
linspace=len(Time)

file=open(nameHydro,'w+')
file.write('!===============================================!\n')
file.write('!  Number of hydrograph                         !\n')
file.write('!===============================================!\n')
file.write('1                                                \n')
file.write('!===============================================!\n')
file.write('!  Hydrograph                                   !\n')
file.write('!===============================================!\n')
file.write('      '+str((linspace))+'\n'                     )
for i in range(linspace):
   file.write(str(Time[i]*(3600*24)) + ' '+str(Q[i])+'\n')


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
   file.write(str(i) + ' '+str(hydro1stGuess)+'\n')


y=Q
period=Time[-1]
T=Time[-1]
Nh=Time[0:70]
Nh=range(1,20)

dt=Time[1]
def F(t):
   return np.interp(t, Time, Q)

x=np.linspace(Time[0],Time[-1],1000)
#y=F(x)


def cn(n):
   c=0
   for i in range(len(Time)):
      c=c + y[i]*np.exp(-1j*2*n*np.pi*Time[i]/T)*dt
   return c/T
   #c = y*np.exp(-1j*2*n*np.pi*Time/T)
   #return c.sum()/c.size

def f(x):
   f=0
   for i in Nh:
      f= f+ 2*cn(i)*np.exp(1j*2*i*np.pi*x/T)
   return f

   #f = np.array([2*cn(i)*np.exp(1j*2*i*np.pi*x/T) for i in Nh])
   #return f.sum()

y2 = np.array([f(t).real for t in Time])+cn(0).real

plt.clf()
plt.plot(Time,Q, label= 'Q_ori')
plt.plot(Time,y2,label = 'Q fourier')

#plt.plot(Time,y, label='Q')
#plt.plot(Time,y2,label = 'Q fourier')
plt.savefig('plot.png')
#plt.show()















