import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg


def read_dassflow_result(filename, y_min, y_max):
    bathy_list = []
    u_list = []
    h_list = []
    y_list = []
    x_list = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("#") or line.startswith("Gnuplot"):
                continue
            parts = line.strip().split()
            if len(parts) >= 8 :
                try :
                    y = float(parts[2])
                    if y <= y_max and y >= y_min :
                        x = float(parts[1])
                        bathy = float(parts[3])
                        h = float(parts[4])
                        u = float(parts[7])
                        x_list.append(x)
                        y_list.append(y)
                        bathy_list.append(bathy)
                        u_list.append(u)
                        h_list.append(h)
                except ValueError :
                    continue
    return np.array(x_list), np.array(bathy_list), np.array(u_list), np.array(h_list)


def get_time_from_filename(filename):
    if "initial" in filename: return 0.0
    if "final" in filename: return float('inf')
    try:
        time_str = filename.replace('result_', '').replace('.dat', '')
        return float(time_str)
    except (ValueError, IndexError): return -1.0


img_paths = []
input_dir = os.getcwd()
res_dir = os.path.join(input_dir, "bin_A//res")

# Récupérer les fichiers
files = [f for f in os.listdir(res_dir) if f.startswith("result_") and f.endswith(".dat")]


filepath = os.path.join(res_dir, files[-2])
g = 9.81
x_m, bathy_m, u_m, h_m = read_dassflow_result(filepath, 1.7, 8)
x_g, bathy_g, u_g, h_g = read_dassflow_result(filepath, 0, 1.7)
x_d, bathy_d, u_d, h_d = read_dassflow_result(filepath, 8, 10)
H_m = h_m + bathy_m
H_g = h_g + bathy_g
H_d = h_d + bathy_d
#h_true  =  (4._8/g)**0.33333333333 * ( 1 + 0.5 * exp( - 16._8 * ( x / lx - 0.5 )**2 ) )

lx=100

def h_exacte(x):
    return (4/g)**(1/3)*(1+0.5*np.exp(-16*x/lx - 1))**2
    #return (4/g)**(1/3)*(1+0.5*np.exp(-16*(x/lx-0.5)**2))

h_supp = h_exacte(x_m)
print("hauteur d'eau à imposer à la fin : ", h_exacte(lx))

q0 = 2
u0 = q0 / h_supp[0]
zb0 = bathy_m[0]
c0 = 0.5*u0**2 + g*(h_supp[0] + zb0)
bathy_sup = c0/g - q0**2 / (2*g*h_supp**2) - h_supp
q_m = h_m * u_m
q_g = h_g * u_g
q_d = h_d * u_d


u0 = q0 / h_supp[0]
z0 = bathy_m[0]
new_alt = z0
c0 = 0.5*u0**2 + g*(h_supp[0] + z0)
bathy_theo =  c0/g - q0**2/(2*g*h_supp**2) - h_supp

print("barycentres des cellules : ", x_m)
print("hauteur d'eau théorique : ", h_supp)
print("hauteur d'eau calculée : ", h_m)
print("débit calculé : ", q_m)



plt.plot(x_m, h_supp + bathy_m, label = "Surface libre théorique", color='red')
plt.plot(x_m, H_m, label = "Surface libre au temps final", linestyle='--', color='orange')
plt.plot(x_m, bathy_m, label = "Bathymétrie du modèle", linestyle='--', color='black')
#plt.plot(x_m, bathy_theo)
plt.xlabel('Distance (en m)')
plt.ylabel('Hauteur d\'eau (m)')
plt.legend()
plt.show()

plt.plot(x_m, q_m, label="débit au milieu", color='blue')
#plt.plot(x_g, q_g, label = "débit à gauche", color='lightblue')
#plt.plot(x_d, q_d, label = "débit à droite", color='lightpink')
plt.title("Débit")
plt.xlabel('Distance (en m)')
plt.ylabel('Débit (en m²/s)')
plt.legend()
plt.show()