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


def f(t):
   return np.interp(t, Time, Q)

Time_n=np.linspace(0,Time[-1],10)
Q_n=f(Time_n)

plt.clf()
plt.plot(Time,Q)
plt.plot(Time_n,Q_n)
#plt.show()


A=np.fft.fft(Q)  # Coeff
sortCoeff=np.argsort(abs(A)) #Trie des coeffs en fonction de abs(coeff)
sortCoeff=list(sortCoeff)    # list
sortCoeff.reverse()          # Inverse ( plus gros au plus petit)
A_bis=abs(A[sortCoeff])
np.savetxt('temp.txt',np.c_[np.array(sortCoeff),A_bis])


plt.clf()
plt.plot(abs(A),label='$|Y^i|$')
plt.xlabel('$i$')
plt.ylabel('$|Y^i|$')
plt.title('$|Y^i|$ in function of i')
plt.legend()
plt.savefig('plot.png')

#nbrCoeffToConserve=len(A)-1
#for i in range(nbrCoeffToConserve,len(A)):
#   A[sortCoeff[i]]=0
#Y=np.fft.ifft(A)

#erreur=np.power(Y.real-Q,2)
#erreur=sum(erreur)/sum(np.power(Q,2))

#plt.clf()
#plt.plot(Time,Q,label='Q')
#plt.plot(Time,Y.real,label='fft')
#plt.legend()
#plt.show()


Erreur,ErreurM=[],[]
PercentageCoeffConserved=[]
fichier=open('resultTF.txt','w+')
if False:
   Aref=np.fft.fft(Q)  # Coeff
   sortCoeff=np.argsort(abs(Aref)) #Trie des coeffs en fonction de abs(coeff)
   sortCoeff=list(sortCoeff)    # list
   sortCoeff.reverse()          # Inverse ( plus gros au plus petit)
   nbrCoeff=len(Aref)
   print Aref
   for j in range(1,nbrCoeff+1):
      nbrCoeffToConserve=j
      A=np.copy(Aref)
      #print A
      for i in range(nbrCoeffToConserve,nbrCoeff):
         A[sortCoeff[i]]=0
      Y=np.fft.ifft(A)
      y2=Y.real
      #print y2
      erreur=np.power(y2-Q,2)
      erreur=np.sqrt(sum(erreur)/sum(np.power(Q,2)))
      erreurM=max(y2-Q)/max(Q)
      percentageCoeffConserved=float(nbrCoeffToConserve)/nbrCoeff
      print 'i :' +str(j)
      print '\t- ErrL2='+str(erreur)
      print '\t- Errmax='+str(erreurM)
      print ''
      plt.clf()
      plt.plot(Time,Q, label='$Q(t)$')
      plt.plot(Time, y2,label = '$Q(t)$ Fourier transform')
      plt.xlabel(r'Time ($days$)')
      plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
      plt.legend()
      plt.savefig('plotTF/hydro_days_'+str(j)+'.png')
      Erreur.append(erreur*100)
      ErreurM.append(erreurM*100)
      PercentageCoeffConserved.append(percentageCoeffConserved*100)
      fichier.write(str(percentageCoeffConserved) + '\t\t'+str(erreur) + '\t\t'+str(erreurM) +'\n')
   
   plt.clf()
   print len(PercentageCoeffConserved)
   print len(Erreur)

   plt.plot(PercentageCoeffConserved,Erreur,label='$L_2(e)$')
   plt.plot(PercentageCoeffConserved,ErreurM,label='$L_\infty(e)$')
   plt.xlabel('Number of conserved coeff ($\%$)')
   plt.ylabel('Error ($\%$)')
   plt.title('Error on $Q$ in function of percent of Fourier coeff conserved')
   plt.legend()
   plt.savefig('plotTF/CVresult_percent.png')

   plt.clf()
   plt.plot(Erreur,label='$L_2(e)$')
   plt.plot(ErreurM,label='$L_\infty(e)$')
   plt.xlabel('Number of coeff Fourier conserved')
   plt.ylabel('Error ($\%$)')
   plt.title('Error on $Q$ in function of number of Fourier coeff conserved')
   plt.legend()
   plt.savefig('plotTF/CVresult_nbr.png')
   fichier.close()














