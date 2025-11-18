import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg


def read_dassflow_result(filename):
    bathy_list = []
    speed_list = []
    H_list = []
    poro_list = []
    yn_list = []
    x_list = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("#") or line.startswith("Gnuplot"):
                continue
            parts = line.strip().split()
            if len(parts) >= 8 :
                try :
                    x = float(parts[1])
                    bathy = float(parts[3])
                    H = float(parts[5])
                    u = float(parts[7])
                    v = float(parts[8])
                    x_list.append(x)
                    speed = np.sqrt(u**2 + v**2)
                    bathy_list.append(bathy)
                    speed_list.append(speed)
                    H_list.append(H)
                except ValueError :
                    continue
            if len(parts) >= 10 :
                try : 
                    poro = float(parts[9])
                    yn = float(parts[10])
                    poro_list.append(poro)
                    yn_list.append(yn)
                except ValueError :
                    continue
    return np.array(x_list), np.array(bathy_list), np.array(speed_list), np.array(H_list), np.array(poro_list), np.array(yn_list)


def read_dassflow_parabola(filename):
    a_list = []
    gamma_list = []
    w_list = []
    hbanks_list = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("#") or line.startswith("Gnuplot"):
                continue
            parts = line.strip().split()
            if len(parts) >= 4 :
                try:
                    a = float(parts[1])
                    gamma = float(parts[2])
                    w = float(parts[3])
                    hbanks = float(parts[4])
                    a_list.append(a)
                    gamma_list.append(gamma)
                    w_list.append(w)
                    hbanks_list.append(hbanks)
                except ValueError:
                    continue
    return np.array(a_list), np.array(gamma_list), np.array(w_list), np.array(hbanks_list)


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

'''for f in range(0, len(files)) :
    filepath = os.path.join(res_dir, files[f])
    x, bathy, speed, H, poro, yn = read_dassflow_result(filepath)
    h = H - bathy
    h_theo = (4/10)**(1/3) * (1+0.5*np.exp(-16*x/1000 - 1)**2)
    time_sec = get_time_from_filename(files[f])
    if time_sec == float('inf') : 
        time_label = "final"
    else : 
        time_label = f"Heure {time_sec / 3600:.1f}"    

    plt.plot(H, label = "temps : " + time_label)'''


filepath = os.path.join(res_dir, files[-2])
g = 9.81
x, bathy, speed, H, poro, yn = read_dassflow_result(filepath)
h = H - bathy
time_sec = get_time_from_filename(files[-2])
if time_sec == float('inf') : 
    time_label = "final"
else : 
    time_label = f"Heure {time_sec / 3600:.1f}"   



h_supp = (4/g)**(1/3)*(1+0.5*np.exp(-16*x/1000 - 1))**2
print(x)
q0 = 2
u0 = q0 / h_supp[0]
zb0 = bathy[0]
c0 = 0.5*u0**2 + g*(h_supp[0] + zb0)
bathy_sup = c0/g - q0**2 / (2*g*h_supp**2) - h_supp
print("u0 =", u0, "c0 =", c0)

plt.plot(h_supp + bathy, label = "Surface libre théorique")
plt.plot(H, label = "Surface libre au temps " + time_label)
plt.plot(bathy_sup, label = "Bathymétrie supposée", linestyle='--', color='red')
plt.plot(bathy, label = "Bathymétrie du modèle", linestyle='--', color='black')
plt.legend()
plt.xlabel("Distance (m)")  
plt.ylabel("Hauteur d'eau (m)")
plt.title("Evolution de la hauteur d'eau le long du canal")
plt.show() 