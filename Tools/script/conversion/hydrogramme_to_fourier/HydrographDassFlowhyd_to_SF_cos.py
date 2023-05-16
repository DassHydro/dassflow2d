import numpy as np
import matplotlib.pylab as plt


#Load data
Time,Q=np.loadtxt('hydrograph/hydrograph_serie_0.6j.txt',skiprows=4,unpack=True)


outputDir='plotSF_2'



y=Q
T=Time[-1]
Nt=Time 
Nh=range(1,20)


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
fichier=open(outputDir+'/resultSF.txt','w+')
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
      plt.savefig(outputDir+'/hydro_days_'+str(i)+'.png')
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
plt.savefig(outputDir+'/CVresult_number.png')



fichier.close()
#time=np.linspace(0,2,10)




#plt.show()













