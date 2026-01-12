import os
import matplotlib.pyplot as plt
import numpy as np


#################
### FONCTIONS ###
#################

def read_dassflow_result(filename, y_min, y_max):
    bathy_list = []
    u_list = []
    h_list = []
    y_list = []
    x_list = []
    phi_list = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("#") or line.startswith("Gnuplot"):
                continue
            parts = line.strip().split()
            if len(parts) >= 9 :
                try :
                    y = float(parts[2])
                    if y <= y_max and y >= y_min :
                        x = float(parts[1])
                        bathy = float(parts[3])
                        h = float(parts[4])
                        u = float(parts[7])
                        phi = float(parts[9])
                        x_list.append(x)
                        y_list.append(y)
                        bathy_list.append(bathy)
                        u_list.append(u)
                        h_list.append(h)
                        phi_list.append(phi)
                except ValueError :
                    continue
    return np.array(x_list), np.array(bathy_list), np.array(u_list), np.array(h_list), np.array(phi_list)

def get_time_from_filename(filename):
    if "initial" in filename: return 0.0
    if "final" in filename: return float('inf')
    try:
        time_str = filename.replace('result_', '').replace('.dat', '')
        return float(time_str)
    except (ValueError, IndexError): return -1.0

def h_exacte(x, lx):
    return (1+0.5*np.exp(-16*(x/lx-0.5)**2))


#################
### VARIABLES ###
#################
img_paths = []
input_dir = os.getcwd()
res_dir = os.path.join(input_dir, "res")

lx = 100
ly = 10
nx = 401
ny = 7
g = 9.81
q0 = 2
zb0 = 0


###################
### RECUP INFOS ###
###################

# Récupérer les fichiers
files = [f for f in os.listdir(res_dir) if f.startswith("result_") and f.endswith(".dat")]

filepath = os.path.join(res_dir, files[-2])

# Récupérer les valeurs pour la ligne milieu
x_m, bathy_m, u_m, h_m, phi_m = read_dassflow_result(filepath, 4.57, 4.59)
H_m = h_m + bathy_m
q_m = h_m * u_m * phi_m

# Hauteur d'eau supposée
x = lx + (lx / (nx - 1))/2
h_supp = h_exacte(x_m, lx)

e = np.linalg.norm(h_supp - h_m, ord=2) / np.linalg.norm(h_supp, ord=2)

########################
### PRINT POUR DEBUG ###
########################

print("barycentres des cellules : ", x_m)
print("hauteur d'eau théorique : ", h_supp)
print("hauteur d'eau calculée : ", h_m)
print("débit calculé : ", q_m)
print("bathy réelle : ", bathy_m)
print("Porosité : ", phi_m)


###############
### FIGURES ###
###############

# Première figure : hauteur d'eau et bathymétried'un côté, vitesse de l'eau de l'autre
#plt.plot(x_m, h_supp)
#plt.show()

fig, ax1 = plt.subplots()

ax1.plot(x_m, h_supp + bathy_m, label="Exact free surface", linestyle='--', color='lightblue')
ax1.plot(x_m, H_m, label="Computed free surface", color='lightblue')
ax1.plot(x_m, bathy_m, label="Bathymetry", color='black')
ax1.set_xlabel('Distance (m)')
ax1.set_ylabel("Altitude (m)")
ax1.legend(loc="upper right", fontsize='x-small')

ax2 = ax1.twinx()
ax2.plot(x_m, h_supp, label="Exact water height", linestyle='--', color='orange')
ax2.plot(x_m, h_m, label="Computed water height", color='orange')
ax2.set_ylabel("Water height (m)", color='orange')
ax2.tick_params(axis='y', labelcolor='orange')
ax2.legend(loc="upper right", fontsize = 'x-small', bbox_to_anchor=(1, 0.5)  )

plt.title("Evolution of water height and velocity along the channel at steady state")
plt.savefig('water_height.png')
plt.close()

plt.plot(x_m, abs(h_supp - h_m)/h_supp)
plt.title(f"Relative error on the water height along the channel (e={e:.4e})")
plt.xlabel("Distance (m)")
plt.ylabel("Relative error on the water height")
plt.savefig('error_water_height.png')
plt.close()

fig, ax1 = plt.subplots()
ax1.plot(x_m, q_m, label="Discharge*phi", color='blue')
ax1.set_ylim(0.2,1)
ax1.set_xlabel('Distance (m)')
ax1.set_ylabel('Lineic discharge*phi (m2/s)')
ax1.legend(loc="upper left")

ax2 = ax1.twinx()
ax2.plot(x_m, u_m, label='Velocity', linestyle='-.', color='blue')
ax2.set_ylabel('Velocity (m/s)', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
ax2.legend(loc="upper right")
plt.title("Evolution of the discharge and velocity along the channel")
plt.savefig('discharge_velocity.png')
plt.close()