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
period=50
T=Time[-1]
#Nh=[0.1,0.25,1,5,10,100]
Nt=Time 
Nh=range(1,20)
#Nh=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,100]
#Nh=np.linspace(0,T,280)
def f(t):
   return np.interp(t, Time, Q)



def an(n):
   a=0
   for i in range(len(Nt)-1):
      dt=Nt[i+1]-Nt[i]
      ts=Nt[i]
      tend=Nt[i+1]
      tm=ts+0.5*(tend-ts)
      a=a+ ((tend-ts)/6)*(f(ts)*np.cos((2.*n*ts*np.pi)/T) + 4*f(tm)*np.cos((2.*n*tm*np.pi)/T) +f(tend)*np.cos((2.*n*tend*np.pi)/T))
 #f(i)*np.cos((2.*n*i*np.pi)/T)*dt
   return (2./T)*a


def bn(n):
   b=0
   for i in range(len(Nt)-1):
      dt=Nt[i+1]-Nt[i]
      ts=Nt[i]
      tend=Nt[i+1]
      tm=ts+0.5*(tend-ts)
      b=b+ ((tend-ts)/6)*(f(ts)*np.sin((2.*n*ts*np.pi)/T) + 4*f(tm)*np.sin((2.*n*tm*np.pi)/T) +f(tend)*np.sin((2.*n*tend*np.pi)/T))
 #f(i)*np.cos((2.*n*i*np.pi)/T)*dt
   return (2./T)*b


   #for i in Nt:
   #   b=b+f(i)*np.sin((2.*n*i*np.pi)/T)*dt
   #return (2./T)*b


def f_f(x):
   f_f = 0.5*an(0)
   for i in Nh:
      f_f=f_f+(an(i)*np.cos((2*i*x*np.pi)/T)+bn(i)*np.sin((2*i*x*np.pi)/T))
   return f_f



#y2 = np.array([f_f(t).real for t in Time])
#erreur=np.power(y2-Q,2)
#erreur=sum(erreur)/sum(np.power(Q,2))
#print '\t- ErrL2='+str(erreur)
#plt.clf()


#plt.plot(Time,Q, label='Q')
#plt.plot(Time, y2,label = 'Q fourier')
#plt.savefig('plot2.png')
#print A

Erreur,ErreurM=[],[]
erreur,erreurM=[],[]
fichier=open('resultSF.txt','w+')
if True:
   for i in range(2,150):
      Nh=range(1,i)
      y2 = np.array([f_f(t).real for t in Time])
      erreur=np.power(y2-Q,2)
      erreur=np.sqrt(sum(erreur)/sum(np.power(Q,2)))
      erreurM=max(y2-Q)/max(Q)
      print 'i :' +str(i)
      print '\t- ErrL2='+str(erreur)
      print '\t- Errmax='+str(erreurM)
      print ''
      plt.clf()
      plt.plot(Time,Q, label='$Q(t)$')
      plt.plot(Time, y2,label = '$Q(t)$ Fourier series')
      plt.xlabel(r'Time ($days$)')
      plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
      plt.legend()
      #plt.title('ErrL2='+str(erreur))
      plt.savefig('plotSF/hydro_days_'+str(i)+'.png')
      Erreur.append(erreur*100)
      ErreurM.append(erreurM*100)
      fichier.write(str(erreur) + '\t\t'+str(erreurM) +'\n')
   


plt.clf()
plt.plot(Erreur,label='$L_2(e)$')
plt.plot(ErreurM,label='$L_\infty(e)$')
plt.xlabel('Number of N used')
plt.ylabel('Error ($\%$)')
plt.title('Error on $Q$ in function of number of Fourier coeff used')
plt.legend()
plt.savefig('plotSF/CVresult_number.png')



fichier.close()
#time=np.linspace(0,2,10)




#plt.show()













