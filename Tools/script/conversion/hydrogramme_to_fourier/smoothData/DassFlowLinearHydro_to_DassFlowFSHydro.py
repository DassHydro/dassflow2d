import numpy as np
import matplotlib.pylab as plt
import os

numberN=40
numberStepSmooth=1000
pathHydro='./hydrograph.txt'
TimeOri,QOri=np.loadtxt(pathHydro,skiprows=4,unpack=True)

Time=np.linspace(0,TimeOri[-1],numberStepSmooth)



def f(t,Time,Q):
   return np.interp(t, Time, Q)

def an(n):
   a=0
   for i in range(len(Nt)):
      a=a+ Q[i]*np.cos((2.*n*Time[i]*np.pi)/T)*dt
   return (2./T)*a


def an_tra(n):
   a=0
   for i in range(len(Nt)-1):
      dt=Nt[i+1]-Nt[i]
      ts=Nt[i]
      tend=Nt[i+1]
      tm=ts+0.5*(tend-ts)
      a=a+ ((tend-ts)/6)*(f(ts)*np.cos((2.*n*ts*np.pi)/T) + 4*f(tm)*np.cos((2.*n*tm*np.pi)/T) +f(tend)*np.cos((2.*n*tend*np.pi)/T))
   return (2./T)*a


def bn(n):
   b=0
   for i in range(len(Nt)):
      b=b+Q[i]*np.sin((2.*n*Time[i]*np.pi)/T)*dt
   return (2./T)*b


def bn_tra(n):
   b=0
   for i in range(len(Nt)-1):
      dt=Nt[i+1]-Nt[i]
      ts=Nt[i]
      tend=Nt[i+1]
      tm=ts+0.5*(tend-ts)
      b=b+ ((tend-ts)/6)*(f(ts)*np.sin((2.*n*ts*np.pi)/T) + 4*f(tm)*np.sin((2.*n*tm*np.pi)/T) +f(tend)*np.sin((2.*n*tend*np.pi)/T))
   return (2./T)*b





def f_f(x):
   f_f = 0.5*an(0)
   for i in Nh:
      f_f=f_f+(an_tra(i)*np.cos((2*i*x*np.pi)/T)+bn_tra(i)*np.sin((2*i*x*np.pi)/T))
      #f_f=f_f+(an(i)*np.cos((2*i*x*np.pi)/T)+bn(i)*np.sin((2*i*x*np.pi)/T))
   return f_f


def f_f_n(x,n,An,Bn):
   f_f_n=(An*np.cos((2*n*x*np.pi)/T)+Bn*np.sin((2*n*x*np.pi)/T))
   return f_f_n



Q=f(Time,TimeOri,QOri)


T=Time[-1]
Nt=Time 
dt=Time[1]


Erreur,ErreurM=[],[]
erreur,erreurM=[],[]
A,B=[],[]
outputDir='./plotSFHydro/'
if (os.path.isdir(outputDir)==False):
   os.mkdir(outputDir)
fichier=open(str(outputDir)+'resultSF.txt','w+')

fichierInputDassFlowFS=open(str(outputDir)+'hydrograph_FS.txt','w+')
fichierInputDassFlowFS.write('!===============================================!\n')
fichierInputDassFlowFS.write('!  Hydrograph SF                                 \n')
fichierInputDassFlowFS.write('!===============================================!\n')
fichierInputDassFlowFS.write('\t'+str(int(numberN)) +'\n')

fichierInputDassFlowFSFG=open(str(outputDir)+'hydrograph_FS_first_guess.txt','w+')
fichierInputDassFlowFSFG.write('!===============================================!\n')
fichierInputDassFlowFSFG.write('!  Hydrograph SF                                 \n')
fichierInputDassFlowFSFG.write('!===============================================!\n')
fichierInputDassFlowFSFG.write('\t'+str(int(numberN)) +'\n')


fichierInputDassFlow=open(str(outputDir)+'hydrograph.txt','w+')
fichierInputDassFlow.write('!===============================================!\n')
fichierInputDassFlow.write('!  Hydrograph                                    \n')
fichierInputDassFlow.write('!===============================================!\n')

fichierInputDassFlow.write('\t'+str(int(numberStepSmooth)) +'\n')
for i in range(numberStepSmooth):
   fichierInputDassFlow.write('\t'+str("{:14.7E}".format(Time[i]  )) +'\t' +str("{:13.7E}".format(Q[i]  ))+'\n')




if True:
   y2=0.5*an(0)
   fichierInputDassFlowFS  .write('\t'+str("{:13.7E}".format(y2*2 )) +'\n')
   fichierInputDassFlowFSFG.write('\t'+str("{:13.7E}".format(150.0 )) +'\n')
   for i in range(1,numberN):
      #Nh=range(1,i)
      #y2 = np.array([f_f(t).real for t in Time])
      An=an(i)
      Bn=bn(i)
      y2=y2+ np.array([f_f_n(t,i,An,Bn).real for t in Time])
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
      plt.savefig(str(outputDir)+'hydro_days_'+str(i)+'.png')
      Erreur.append(erreur*100)
      ErreurM.append(erreurM*100)
      fichier.write(str(erreur) + '\t\t'+str(erreurM) +'\n')
      A.append(An)
      B.append(Bn)
      fichierInputDassFlowFS.write('\t'+str(i)+'\t'+str("{:14.7E}".format(A[i-1]  )) +'\t' +str("{:13.7E}".format(B[i-1]  ))+'\n')
      fichierInputDassFlowFSFG.write('\t'+str(i)+'\t'+str("{:14.7E}".format(0  )) +'\t' +str("{:13.7E}".format(0  ))+'\n')


      plt.clf()
      plt.plot(Erreur,label='$e_2(Q)$')
      plt.plot(ErreurM,label='$e_\infty(Q)$')
      plt.xlabel('Number of N used')
      plt.ylabel('Error ($\%$)')
      plt.title('Error on $Q$ in function of number of Fourier coeff used')
      plt.legend()
      plt.savefig(str(outputDir)+'CVresult_number.png')

    

fichier.close()
#time=np.linspace(0,2,10)




#plt.show()













