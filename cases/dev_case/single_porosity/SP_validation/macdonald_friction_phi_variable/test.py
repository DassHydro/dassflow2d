import numpy as np
import matplotlib.pyplot as plt

###########################
## PARAMETERS DEFINITION ##
###########################

def int_friction(x, lx):
    N = 10000
    t = np.linspace(0, x, N)
    integrand = h_ex(t, lx)**(-10/3)
    result = np.trapz(integrand, t)
    return result

def h_ex(x, lx) :
    return 1+0.5*np.exp(-16*(x/lx - 0.5)**2)

def bathy(x, lx, n, cref, zref, qref, xref) :
    if x == xref:
        return zref
    else :
        h = h_ex(x, lx)
        return cref - qref**2/(2*g*h**2) - h - (qref**2)*(n**2) * int_friction(x, lx)

g = 9.81
lx = 100
n = 0.05
xref = 0
zref = 0
qref = 2
href = h_ex(xref, lx)
cref = qref**2 / (2*g*href**2) + href + zref
c = 11.868197654947547
print(cref)
print(c/g)

###########################
## TEST SIMULATIONS ##
###########################
x = np.linspace(0, 100, 100)   
h = h_ex(x, lx)
z = np.zeros_like(x)
for i in range(len(x)) :
    z[i] = bathy(x[i], lx, n, cref, zref, qref, xref)
plt.plot(x, z, label='Bathy modifiée')
plt.plot(x, h+z, label='Cote eau modifiée')
plt.xlabel('x (m)')
plt.ylabel('Elevation (m)')
plt.title('Bathymétrie et cote eau modifiées selon MacDonald')
plt.legend()
plt.show()

plt.plot(x, h, label='Hauteur eau modifiée')
plt.xlabel('x (m)')
plt.ylabel('Hauteur eau (m)')
plt.title('Hauteur d\'eau modifiée selon MacDonald')
plt.legend()
plt.show()

###########################
## CALCUL Z INLET OUTLET ##
###########################
x = -0.5
z = c/9.81 - 2**2/(2*g*h_ex(x, 100)**2) - h_ex(x, 100) - n**2*qref**2 * int_friction(x, lx)
z2 = bathy(x, lx, n, cref, zref, qref, xref)
print("z inlet = ", z)

x = 100.5
z = c/9.81 - 2**2/(2*g*h_ex(x, 100)**2) - h_ex(x, 100) - n**2*qref**2 * int_friction(x, lx)
z2 = bathy(x, lx, n, cref, zref, qref, xref)
print("z outlet = ", z)

###########################
## HAUTEUR EAU EN SORTIE ##
###########################
print("h_fin = ", h_ex(100, 100))