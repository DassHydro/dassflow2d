import numpy as np
import matplotlib.pylab as plt

data=np.loadtxt('data')
Day=data[:,0]
hours=data[:,1]
minu=data[:,2]
Qori=data[:,3]

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


data=np.loadtxt('data_1')
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


#
#plt.plot(Time,Q)
#plt.show()


T=np.linspace(0,Time[-1],1000)

def f(t):
   return np.interp(t, Time, Q)


def smooth(tp,l_ts,te_smooth,tp_smooth):
    te_smooth,tp_smooth =[],[]
    for i in range(1,len(tp)):
        im = np.mod(i,l_ts)
        if im==0:
            te_smooth.append(int(i-l_ts/2))
            tp_smooth.append( np.mean( tp[i-l_ts:i] ) )
    # convert to float array     
    te_smooth = np.array( map(float,te_smooth) )
    tp_smooth = np.array( map(float,tp_smooth) )
    
    return [te_smooth,tp_smooth]

def smooth(T,Q,dt):
   Q_smooth=[]
   dtd2=float(dt)/2
   for i in range(0,len(T)):
      bsup=(T[i]+dtd2)
      binf=(T[i]-dtd2)

      meanQ=np.mean(Q[(T >= binf) & (T <= bsup)])
      Q_smooth.append(meanQ)

      

   return Q_smooth





dt=T[1]

Q_dt_cst=f(T)
dt_smooth_h=4.0
dt_smooth=dt_smooth_h/(24.0)
Q_dt_cst_smooth4=smooth(T,Q_dt_cst,dt_smooth)


dt_smooth_h=8.0
dt_smooth=dt_smooth_h/(24.0)
Q_dt_cst_smooth8=smooth(T,Q_dt_cst,dt_smooth)


dt_smooth_h=16.0
dt_smooth=dt_smooth_h/(24.0)
Q_dt_cst_smooth16=smooth(T,Q_dt_cst,dt_smooth)

plt.clf()


plt.plot(Time,Qori, label  ='$Q(t)$: original data')

plt.plot(T,Q_dt_cst_smooth4,label ='$Q(t)$ mean $4$ $h$')
plt.plot(T,Q_dt_cst_smooth8,label ='$Q(t)$ mean $8$ $h$')
plt.plot(T,Q_dt_cst_smooth16,label ='$Q(t)$ mean $16$ $h$')


plt.xlabel(r'Time ($days$)')
plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
plt.legend()
plt.savefig('plotTraitement.png')


np.savetxt('data_1_post',np.c_[T,Q_dt_cst_smooth8])
#plt.show()





