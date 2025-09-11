import dassflow2d as df2d
import matplotlib.pyplot as plt
import shutil, os, sys
import numpy as np
import csv
import random
from mpi4py import MPI

path_to_SP = os.path.abspath(os.path.join(os.path.dirname(__file__),'../..'))

path_to_Outils = os.path.join(path_to_SP,'Outils')

sys.path.append(path_to_SP)

from Outils.out import *
from Outils.mesh_generator import *

# TODO : adapt this into ? model.meshing object + mesh_extent_box from a real dassflow mesh/setup (with wrapped object)
class mesh_box:
    def __init__(self, xmin, ymax, ncol, nrow, dx, dy, x_cell, y_cell):
        self.xmin = xmin
        self.ymax = ymax
        self.ncol = ncol
        self.nrow = nrow
        self.dx = dx
        self.dy = dy
        self.x_cell = x_cell
        self.y_cell = y_cell

##############################################################################################################
# --- CONFIGURATION DES CHEMINS ---
initial_code_dir = os.getcwd() # Sauvegarde le répertoire de travail initial
bin_dir = os.path.join(initial_code_dir, "bin_A")
res_dir = os.path.join(bin_dir, 'res') # Définir res_dir ici, basé sur le bin_dir correct

##############################################################################################################

##########
# Initialise bin
##########

os.chdir(initial_code_dir) # Revenir au répertoire initial avant les "make clean"

os.system("make cleanres")        # removes all in bin_dir/res directory
os.system("make cleanmsh")        # removes all in bin_dir/msh directory
os.system("make cleanmin")        # removes all in bin_dir/min directory


os.chdir(bin_dir)

##########
# MPI
##########

df2d.wrapping.m_mpi.init_mpi()
rank = df2d.wrapping.m_mpi.get_proc() # get the rank and number of processors
nproc = df2d.wrapping.m_mpi.get_np()
mpi = [rank, nproc]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

##########
# Mesh
##########

L = 100
dx = 1
type = 'channel'

mesh_name = "automaticaly_generated_mesh.txt "

# Remark : If you choose type = 'box' or 'flat_square' you must comment the call to the function plot_dat and add '_box' to the name function plot_vtk


##########
# Model
##########

ts = 20000
dtw= 2000 

# REMARK : In this test-case, you must set nland equal to the number of cells in your mesh

use_porosity = 1
df2d.wrapping.read_input("input.txt")
input_params={ "mesh_name": mesh_name ,
               "ts": ts ,
			   "use_obs":'0',
			   "use_UVobs":'0',
			   "use_Zobs":'0',
               "use_porosity": use_porosity,
            
			   "w_obs":'0',
               "w_vtk":'1',
               "w_gnuplot":'1',
               "w_tecplot":'0',
               "adapt_dt":'1',

               "dt":'0.1',

               "dtw": dtw ,
               "dta":"100",
               
               "bc_infil":"0",
			   "bc_rain":"0",
                "mesh_type": "'dassflow'" }



df2d.wrapping.read_input(os.path.join(bin_dir,"input.txt"))
my_model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = os.path.join(bin_dir,"res","simu.hdf5") , run_type = "direct", clean = True, custom_config=input_params)

#Re initialize mesh from new updated geo file
my_model.init_mesh()

# Input parameters for the generation of observations

Config = df2d.core.config.Config()
Config.set(custom_config = input_params)
nc = my_model.kernel.mesh.nc    # number of cells 
nland = nc
##########
# Friction
##########

nc = my_model.kernel.mesh.nc

nland = nc

#Create Python class by calling wrapped initialise routines
my_model.kernel.my_friction =  df2d.wrapping.m_model.friction_data(my_model.kernel.mesh)
#Allocate and get initial values from Fortran
my_model.kernel.my_friction.nland = nland
df2d.wrapping.call_model.init_friction(my_model.kernel)

#Provide values, on top of initial ones from Fortran initialization routine, in Python structure

my_model.kernel.my_friction.manning[:] = 0.0
my_model.kernel.my_friction.manning_beta[:] = 0
my_model.kernel.my_friction.land[:] = 1

##########
# Porosity
##########

# On initialise la structure porosité côté Python
my_model.kernel.my_porosity = df2d.wrapping.m_model.porosity_data(my_model.kernel.mesh)
my_model.kernel.my_porosity.nland = nc
df2d.wrapping.call_model.init_porosity(my_model.kernel)

# Allocation des tableaux (Python → Fortran)
my_model.kernel.my_porosity.a = np.zeros(nc, dtype=float)
my_model.kernel.my_porosity.gamma = np.zeros(nc, dtype=float)
my_model.kernel.my_porosity.width = np.zeros(nc, dtype=float)
my_model.kernel.my_porosity.hbanks = np.zeros(nc, dtype=float)
my_model.kernel.my_porosity.land[:] = np.arange(1, nc + 1)
my_model.kernel.my_porosity.phi = np.zeros(nc, dtype=float)

my_model.kernel.my_porosity.gamma = 2.0 
my_model.kernel.my_porosity.hbanks = 25.0 

##########
# Hydraulic states
##########

my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
my_model.kernel.dof0 = my_model.kernel.dof
my_model.kernel.dof0.h[:] = 10.0 # Altitude de surface libre Zs constante
my_model.kernel.dof0.u[:] = 0.0
my_model.kernel.dof0.v[:] = 0.0

##############################################################################################################

##########################################
# Run
##########################################

#my_model.meshing.plot_dev("cell")

df2d.wrapping.call_model.init_fortran(my_model.kernel)

df2d.wrapping.call_model.run(my_model.kernel, arg = "direct")

if (os.path.isdir("./obs")):
    shutil.rmtree('./obs')
#shutil.copytree("./res/obs", "./obs")



##############################################################################################################
def read_dassflow_result(filename):
    """Lit un fichier de résultat temporel (result_xxx.dat)."""
    speed_list, h_list = [], []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith("#"): continue
            parts = line.strip().split()
            if len(parts) >= 11:
                try:
                    h = float(parts[4]); 
                    u = float(parts[7]); 
                    v = float(parts[8]); 
                    speed = np.sqrt(u**2 + v**2)
                    h_list.append(h); speed_list.append(speed); 
                except (ValueError, IndexError): continue
    return np.array(speed_list), np.array(h_list)

def get_time_from_filename(filename):
    """Extrait la valeur numérique du temps depuis un nom de fichier."""
    if "initial" in filename: return 0.0
    if "final" in filename: return float('inf')
    try:
        time_str = filename.replace('result_', '').replace('.dat', '')
        return float(time_str)
    except (ValueError, IndexError): return -1.0

# --- GRAPHIQUE DE L'ÉVOLUTION VERS L'ÉTAT STATIONNAIRE ---

dat_files = sorted(
    [f for f in os.listdir(res_dir) if f.startswith("result_") and f.endswith(".dat") and "final" not in f],
    key=get_time_from_filename)

if len(dat_files) < 2:
    print("Pas assez de fichiers pour calculer des erreurs relatives (il faut au moins 2 fichiers).")
else:
    times_for_plot = []
    h_rel_errors, speed_rel_errors = [], []
    speed = []

    # Lecture du premier état (t=0)
    speed_curr, h_curr = read_dassflow_result(os.path.join(res_dir, dat_files[0]))
    
    # Boucle sur les états suivants (t=1, t=2, ...)
    for f_idx, f in enumerate(dat_files[1:]): # Utilise f_idx pour suivre l'itération
        time_curr = get_time_from_filename(f)
        times_for_plot.append(time_curr)

        speed_prev, h_prev = speed_curr, h_curr

        speed_curr, h_curr = read_dassflow_result(os.path.join(res_dir, f))
        speed.append(max(speed_curr))

        # print("speed_prev = ", speed_prev, "speed_curr = ", speed_curr)

        # --- Calcul de l'erreur relative maximale pour h ---
        h_diff = np.abs(h_curr - h_prev)
        errors_h_per_cell = h_diff / h_prev
        h_rel = np.max(errors_h_per_cell[errors_h_per_cell != -1.0])
        h_rel_errors.append(h_rel)

        speed_diff = np.abs(speed_curr - speed_prev)
        errors_speed_per_cell = speed_diff / speed_prev
        speed_rel = np.max(errors_speed_per_cell[errors_speed_per_cell != -1.0])
        speed_rel_errors.append(speed_rel)

    #print(h_rel_errors)
    #print(speed_rel_errors)

times_for_plot = np.array(times_for_plot)
h_rel_errors = np.array(h_rel_errors)
speed_rel_errors = np.array(speed_rel_errors)
speed = np.array(speed)

# print("times_for_plot = ", times_for_plot)

plt.figure(figsize=(12, 7))
plt.yscale("log")
plt.plot(times_for_plot, h_rel_errors, "o-", label="Erreur relative h")
plt.axhline(y=1e-7, color='r', linestyle='--', label="Seuil de convergence = 1e-7")
plt.xlabel("Temps (s)")
plt.ylabel("Erreur relative max (échelle log)")
plt.title("Convergence vers l'état stationnaire (erreur relative sur la hauteur d'eau)")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
save_path = os.path.join(res_dir, "evolution_erreur_relative_h.png")
plt.savefig(save_path)
plt.show()

plt.figure(figsize=(12, 7))
#plt.yscale("log")
plt.plot(times_for_plot, speed_rel_errors, "o-", label="Erreur relative vitesse")
#plt.axhline(y=1e-7, color='r', linestyle='--', label="Seuil de convergence = 1e-7")
plt.xlabel("Temps (s)")
plt.ylabel("Erreur relative max") 
plt.title("Convergence vers l'état stationnaire (erreur relative sur la vitesse)")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
save_path = os.path.join(res_dir, "evolution_erreur_relative_v.png")
plt.savefig(save_path)
plt.show()

plt.figure(figsize=(12, 7))
plt.yscale("log")
plt.plot(times_for_plot, speed, "o-", label="Evolution de la vitesse")
plt.axhline(y=1e-7, color='r', linestyle='--', label="Seuil de convergence = 1e-7")
plt.xlabel("Temps (s)")
plt.ylabel("Vitesse (m/s)")
plt.title("Convergence vers l'état stationnaire (évolution de la vitesse)")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
save_path = os.path.join(res_dir, "evolution_vitesse.png")
plt.savefig(save_path)
plt.show()



'''
        u_plot_times = times_arr[u_errors_arr != -1.0]
        u_plot_errors = u_errors_arr[u_errors_arr != -1.0]
        if len(u_plot_errors) > 0:
            plt.plot(u_plot_times, u_plot_errors, "s-", label="Erreur relative u")
        else:
            print("Aucune erreur relative u valide à tracer.")

        v_plot_times = times_arr[v_errors_arr != -1.0]
        v_plot_errors = v_errors_arr[v_errors_arr != -1.0]
        if len(v_plot_errors) > 0:
            plt.plot(v_plot_times, v_plot_errors, "^-", label="Erreur relative v")
        else:
            print("Aucune erreur relative v valide à tracer.")

        epsilon = 1e-7 # Un seuil de convergence typique
        plt.axhline(y=epsilon, color='r', linestyle='--', label=f"Seuil de convergence = {epsilon:.0e}")
        
        

        # Gestion des limites Y
        all_tracable_errors = []
        if len(h_plot_errors) > 0: all_tracable_errors.extend(h_plot_errors)
        if len(u_plot_errors) > 0: all_tracable_errors.extend(u_plot_errors)
        if len(v_plot_errors) > 0: all_tracable_errors.extend(v_plot_errors)

        if all_tracable_errors:
            min_plot_val = np.min(all_tracable_errors)
            max_plot_val = np.max(all_tracable_errors)
            # S'assurer que les limites Y sont toujours > 0 pour l'échelle log
            plt.ylim(max(1e-15, min_plot_val / 5.0), max_plot_val * 5.0)
        else:
            plt.ylim(1e-15, 1.0) # Limites par défaut

        
        plt.close()
        print(f"Graphique sauvegardé : {save_path}")
        '''

########################
# Outputs from python
########################

# Meaning of graphe = [ 'h' , 'u' , 'v' , 'qx' , 'qy' ]

graphe = [1,0,0,0,0]

# Remark : reference not usable here

display_ref = 0

# Récupération des tranches initiales (pour h0, u0, v0)
h0 = my_model.kernel.dof0.h[:nc]
u0 = my_model.kernel.dof0.u[:nc]
v0 = my_model.kernel.dof0.v[:nc]

# Récupération de la tranche active de 'h' à l'état final (pour le calcul de la moyenne)
h_active_final = my_model.kernel.dof.h[:nc] # Nommé différemment pour plus de clarté
u_active_final = my_model.kernel.dof.u[:nc]
v_active_final = my_model.kernel.dof.v[:nc]


# Calcul de la moyenne et de la déviation
h_final_mean = np.mean(h_active_final)
h_deviation_from_mean = h_active_final - h_final_mean


# --- Partie cruciale pour le "hack" VTK ---

# 1. SAUVEGARDER la tranche originale de 'h' AVANT de la modifier
original_h_slice = np.copy(my_model.kernel.dof.h[:nc]) # Utilisation de np.copy() pour une copie réelle

# 2. MODIFIER la tranche active de 'h' du modèle avec les déviations
my_model.kernel.dof.h[:nc] = h_deviation_from_mean

#plot_dat(graphe,display_ref,mesh_name,ts,h,u,v,h0,u0,v0) # Si décommenté, utilisera le h modifié

# To save pictures : save = 1
save = 1

# Put : 'initial' or 'final'
time = 'final'

# 3. TRACER le champ 'h'. Il tracera maintenant les déviations.
plot_vtk(initial_code_dir ,'h',time,ts,save)

# 4. RESTAURER la tranche originale de 'h' APRÈS le tracé
my_model.kernel.dof.h[:nc] = original_h_slice


# Les autres tracés (u, v, zs) sont effectués normalement, avec les VRAIES valeurs
# (car h a été restauré, et u,v,zs n'ont pas été modifiés par ce hack).
plot_vtk(initial_code_dir ,'u',time,ts,save)
plot_vtk(initial_code_dir ,'v',time,ts,save)
plot_vtk(initial_code_dir ,'zs','initial',ts,save)
plot_vtk(initial_code_dir ,'zs',time,ts,save)

df2d.wrapping.call_model.clean_model(my_model.kernel)
