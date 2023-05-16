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
plt.title('Evolution of discharge ')#at Toulouse\'s station')
plt.xlabel(r'Time ($s$)')
plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
plt.savefig('hydro_seconds.png')

plt.clf()
plt.plot(Time,Q,label='Q')
plt.title('Evolution of discharge ')#at Toulouse\'s station')
plt.xlabel(r'Time ($days$)')
plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
plt.savefig('hydro_days.png')


plt.clf()
plt.plot(Time/(365.25),Q,label='Q')
plt.title('Evolution of discharge ')#at Toulouse\'s station')
plt.xlabel(r'Time ($years$)')
plt.ylabel(r'Discharge ($m^3.s^{-1}$)')
plt.xlim(Time[0]/(365.25),Time[-1]/(365.25))
plt.savefig('hydro_years.png')

print 'Annee : '+str(Time[-1]/365.25)

a=Time[1:]-Time[:-1]
print 'Moyenne entre 2 pts : ' +str(np.mean(a))
print 'Nombre pts : ' + str(len(a))
print 'Min :  ' + str(np.min(Q))
print 'Max : '+ str(np.max(Q))
print 'Mean : '+ str(np.mean(Q))



#plt.legend()
#plt.savefig('hydro_days.png')

#plt.show()













