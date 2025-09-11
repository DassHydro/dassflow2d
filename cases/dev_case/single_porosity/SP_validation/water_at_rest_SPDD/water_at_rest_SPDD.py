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
code_dir =  os.getcwd() #os.path.join(dassflow_dir,"code")
bin_dir = os.path.join(code_dir,"bin_A")

##############################################################################################################

##########
# Initialise bin
##########

os.chdir(code_dir)

os.system("make cleanres")			   # removes all in bin_dir/res directory
os.system("make cleanmsh")			   # removes all in bin_dir/msh directory
os.system("make cleanmin")			   # removes all in bin_dir/min directory

#if os.path.isfile(f"rm {bin_dir}/restart.bin"):
#	os.system(f"rm {bin_dir}/restart.bin")   # removes all in bin_dir/msh directory

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

ts = 13000
dtw= 1300 

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

#plot_vtk(code_dir,'porosity',time,ts,save)
def read_dassflow_result(filename):
    """Lit un fichier de résultat temporel (result_xxx.dat)."""
    u_list, H_list, v_list = [], [], []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith("#"): continue
            parts = line.strip().split()
            if len(parts) >= 11:
                try:
                    H = float(parts[4]); 
                    u = float(parts[7]); 
                    v = float(parts[8]); 
                    u_list.append(u); H_list.append(H);v_list.append(v); 
                except (ValueError, IndexError): continue
    return np.array(u_list),np.array(v_list), np.array(H_list)

def get_time_from_filename(filename):
    """Extrait la valeur numérique du temps depuis un nom de fichier."""
    if "initial" in filename: return 0.0
    if "final" in filename: return float('inf')
    try:
        time_str = filename.replace('result_', '').replace('.dat', '')
        return float(time_str)
    except (ValueError, IndexError): return -1.0

# --- GRAPHIQUE DE L'ÉVOLUTION VERS L'ÉTAT STATIONNAIRE ---
if rank == 0:
    print("\n--- Génération du graphique d'évolution vers l'état stationnaire (erreurs relatives) ---")

    res_dir = os.path.join(bin_dir, 'res')
    dat_files = sorted(
        [f for f in os.listdir(res_dir) if f.startswith("result_") and f.endswith(".dat") and "final" not in f],
        key=get_time_from_filename
    )

    if len(dat_files) < 2:
        print("Pas assez de fichiers pour calculer des erreurs relatives (il faut au moins 2 fichiers).")
    else:
        times_for_plot = []
        h_rel_errors, u_rel_errors, v_rel_errors = [], [], []

        # Lecture du premier état (t=0)
        u_prev, v_prev , h_prev = read_dassflow_result(os.path.join(res_dir, dat_files[0]))

        # Boucle sur les états suivants (t=1, t=2, ...)
        for f_idx, f in enumerate(dat_files[1:]): # Utilise f_idx pour suivre l'itération
            time_curr = get_time_from_filename(f)
            u_curr, v_curr, h_curr  = read_dassflow_result(os.path.join(res_dir, f))

            # --- Calcul de l'erreur relative maximale pour h ---
            with np.errstate(divide='ignore', invalid='ignore'):
                h_diff = np.abs(h_curr - h_prev)
                h_denom = np.abs(h_prev)
                zero_h_mask = np.isclose(h_denom, 0.0, atol=1e-12)
                errors_h_per_cell = np.where(zero_h_mask, -1.0, h_diff / h_denom)
                
                if np.all(errors_h_per_cell == -1.0):
                    h_rel = -1.0
                else:
                    h_rel = np.max(errors_h_per_cell[errors_h_per_cell != -1.0])

            # --- Calcul de l'erreur relative maximale pour u ---
            with np.errstate(divide='ignore', invalid='ignore'):
                u_diff = np.abs(u_curr - u_prev)
                u_denom = np.abs(u_prev)
                zero_u_mask = np.isclose(u_denom, 0.0, atol=1e-12)
                errors_u_per_cell = np.where(zero_u_mask, -1.0, u_diff / u_denom)
                
                if np.all(errors_u_per_cell == -1.0):
                    u_rel = -1.0
                else:
                    u_rel = np.max(errors_u_per_cell[errors_u_per_cell != -1.0])

            # --- Calcul de l'erreur relative maximale pour v ---
            with np.errstate(divide='ignore', invalid='ignore'):
                v_diff = np.abs(v_curr - v_prev)
                v_denom = np.abs(v_prev)
                zero_v_mask = np.isclose(v_denom, 0.0, atol=1e-12)
                errors_v_per_cell = np.where(zero_v_mask, -1.0, v_diff / v_denom)
                
                if np.all(errors_v_per_cell == -1.0):
                    v_rel = -1.0
                else:
                    v_rel = np.max(errors_v_per_cell[errors_v_per_cell != -1.0])

            times_for_plot.append(time_curr)
            h_rel_errors.append(h_rel)
            u_rel_errors.append(u_rel)
            v_rel_errors.append(v_rel)
            
            # Message de débogage amélioré
            print(f"Time: {time_curr:.1f}s | h_prev[0]={h_prev[0]:.2e}, u_prev[0]={u_prev[0]:.2e}, v_prev[0]={v_prev[0]:.2e} | h_curr[0]={h_curr[0]:.2e}, u_curr[0]={u_curr[0]:.2e}, v_curr[0]={v_curr[0]:.2e} | h_rel={h_rel:.4e}, u_rel={u_rel:.4e}, v_rel={v_rel:.4e}")

            # Mise à jour pour l'itération suivante
            h_prev, u_prev, v_prev = h_curr, u_curr, v_curr

        # --- Tracé du graphique ---
        plt.figure(figsize=(12, 7))
        plt.yscale("log") # Déplacé ici, avant le plot, c'est la meilleure pratique
        
        times_arr = np.array(times_for_plot)
        h_errors_arr = np.array(h_rel_errors)
        u_errors_arr = np.array(u_rel_errors)
        v_errors_arr = np.array(v_rel_errors)
        
        # Filtrer les points avec -1.0 avant de les tracer sur l'échelle log
        # car log(-1) n'est pas défini et causerait des problèmes.
        
        h_plot_times = times_arr[h_errors_arr != -1.0]
        h_plot_errors = h_errors_arr[h_errors_arr != -1.0]
        if len(h_plot_errors) > 0:
            plt.plot(h_plot_times, h_plot_errors, "o-", label="Erreur relative h")
        else:
            print("Aucune erreur relative h valide à tracer.")

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
        
        plt.xlabel("Temps (s)")
        plt.ylabel("Erreur relative max (échelle log)") # Libellé mis à jour pour être clair
        plt.title("Convergence vers l'état stationnaire (erreurs relatives)")
        plt.grid(True, which="both")
        plt.legend()
        plt.tight_layout()

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

        save_path = os.path.join(res_dir, "evolution_erreur_relative.png")
        plt.savefig(save_path)
        plt.close()
        print(f"Graphique sauvegardé : {save_path}")

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
plot_vtk(code_dir,'h',time,ts,save)

# 4. RESTAURER la tranche originale de 'h' APRÈS le tracé
my_model.kernel.dof.h[:nc] = original_h_slice


# Les autres tracés (u, v, zs) sont effectués normalement, avec les VRAIES valeurs
# (car h a été restauré, et u,v,zs n'ont pas été modifiés par ce hack).
plot_vtk(code_dir,'u',time,ts,save)
plot_vtk(code_dir,'v',time,ts,save)
plot_vtk(code_dir,'zs','initial',ts,save)
plot_vtk(code_dir,'zs',time,ts,save)

df2d.wrapping.call_model.clean_model(my_model.kernel)