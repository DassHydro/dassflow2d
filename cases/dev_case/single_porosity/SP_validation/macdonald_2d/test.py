import numpy as np
import matplotlib.pyplot as plt

g = 9.81

def h_ex(x, nx) :
    return 1+0.5*np.exp(-16*(x/nx - 0.5)**2)

nx = 100
x = np.linspace(0, nx)
h = h_ex(x, nx)

qref = 2
cref = 0.5*(qref/h[-1])**2 + g*h[-1]

zb = cref / g - qref**2/(2*g*h**2) - h
zb = np.zeros(len(x))
for i in range(len(x)) : 
    zb[i] = h[-1] - h[i] + qref**2/(2*g)*(1/h[-1]**2 - 1/h[i]**2)
print(h)
print(x)

plt.plot(x, h)
plt.plot(x, zb)
plt.show()

a = np.linspace(start=0, stop = 2000, num = 2)

print(a)