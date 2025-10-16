import dassflow2d as df2d
import matplotlib.pyplot as plt
import shutil, os, sys
import numpy as np
import csv
import random
from mpi4py import MPI

# Assurez-vous que ces imports sont nécessaires et disponibles
# from Outils.mesh_generator import * # Plus nécessaire si on utilise un maillage statique
# from Outils.out import * # Si non utilisé, peut être supprimé

path_to_SP = os.path.abspath(os.path.join(os.path.dirname(__file__),'../..'))
path_to_Outils = os.path.join(path_to_SP,'Outils')
sys.path.append(path_to_SP)

# La classe mesh_box n'est pas utilisée, vous pouvez la supprimer.
# class mesh_box:
#     def __init__(self, xmin, ymax, ncol, nrow, dx, dy, x_cell, y_cell):
#         self.xmin = xmin
#         self.ymax = ymax
#         self.ncol = ncol
#         self.nrow = nrow
#         self.dx = dx
#         self.dy = dy
#         self.x_cell = x_cell
#         self.y_cell = y_cell

##############################################################################################################

initial_code_dir = os.getcwd() # Sauvegarde le répertoire de travail initial
bin_dir = os.path.join(initial_code_dir, "bin_A")
# res_dir est défini plus bas par DassFlow via my_model, pas besoin ici pour l'initialisation

##############################################################################################################

##########
# Initialise bin (et assure les dossiers nécessaires)
##########

os.chdir(initial_code_dir)

# Il est normal que ces "make clean" échouent si les répertoires sont déjà vides.
# L'important est que si des fichiers existent, ils soient nettoyés.
print("Cleaning old results (if makefiles are configured)...")
os.system("make cleanres")        # removes all in bin_dir/res directory
os.system("make cleanmsh")        # removes all in bin_dir/msh directory
os.system("make cleanmin")        # removes all in bin_dir/min directory

# if os.path.isfile(f"rm {bin_dir}/restart.bin"):
#    os.system(f"rm {bin_dir}/restart.bin")    # removes all in bin_dir/msh directory

os.chdir(bin_dir)
print("ok, changed to bin_A directory.")

# CRÉER LES DOSSIERS 'msh' ET 'res' S'ILS N'EXISTENT PAS
# DassFlow en aura besoin pour écrire ses fichiers intermédiaires et ses résultats.
os.makedirs(os.path.join(bin_dir, 'msh'), exist_ok=True)
os.makedirs(os.path.join(bin_dir, 'res'), exist_ok=True)
print(f"Ensured '{os.path.join(bin_dir, 'msh')}' and '{os.path.join(bin_dir, 'res')}' directories exist within bin_A.")


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

# COMMENTEZ OU SUPPRIMEZ LES PARAMÈTRES POUR LA GÉNÉRATION DE MAILLAGE DYNAMIQUE
# L_mesh = 100.0 # Longueur du canal
# dx_mesh = 1.0  # Taille de cellule en x
# mesh_type_for_gen ='channel'
# generated_mesh_filename = '{}_dx={}_L={}.geo'.format(mesh_type_for_gen, dx_mesh, L_mesh)
# mesh_full_path_generated = os.path.join(bin_dir, generated_mesh_filename)


# DÉFINIR LE NOM DU MAILLAGE STATIQUE DIRECTEMENT
mesh_name_for_dassflow = 'channel.geo' # C'est le fichier que vous avez copié !

# Nous utilisons maintenant un maillage statique, donc pas de génération.
# Assurez-vous que le fichier existe dans bin_A.
mesh_full_path_for_dassflow_read = os.path.join(bin_dir, mesh_name_for_dassflow)
if not os.path.exists(mesh_full_path_for_dassflow_read):
    print(f"ERROR: Static mesh file '{mesh_name_for_dassflow}' not found in '{bin_dir}'. Please ensure it was copied correctly.")
    sys.exit(1) # Quitter le script si le fichier n'est pas là
else:
    print(f"Using static mesh file '{mesh_name_for_dassflow}' from '{bin_dir}'. Skipping generation.")

# mesh_name_for_dassflow est déjà défini ci-dessus avec 'channel.geo'
##########
# Model
##########

ts = 10000
dtw = ts*0.2

use_porosity = 1

input_params={ "mesh_name": mesh_name_for_dassflow ,
                "ts": ts ,
                "use_obs":'0',
                "use_UVobs":'0',
                "use_Zobs":'0',
                "use_porosity": use_porosity,
                
                "w_obs":'1',
                "w_vtk":'1', # Gardez à 1 pour que DassFlow écrive les fichiers VTK nécessaires
                "w_gnuplot":'1',
                "w_tecplot":'0',
                "adapt_dt":'1',

                "dt":'0.1',

                "dtw": dtw ,
                "dta":"100",
                
                "bc_infil":"0",
                "bc_rain":"0",}

# Lisez le fichier input.txt pour être sûr que toutes les configs sont là, puis surchargez avec custom_config
df2d.wrapping.read_input(os.path.join(bin_dir,"input.txt"))
my_model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = os.path.join(bin_dir,"res","simu.hdf5") , run_type = "direct", clean = True, custom_config=input_params)


# Reinitialize mesh from new updated geo file
my_model.init_mesh()





# my_model.meshing.plot() # Cette fonction peut ouvrir une fenêtre graphique PyVista, à commenter si non souhaité

# Input parameters for the generation of observations
Config = df2d.core.config.Config()
Config.set(custom_config = input_params)
nc = my_model.kernel.mesh.nc    # number of cells 
nland = nc # nland est souvent égal à nc pour ce type de maillage

##########
# Friction
##########

# nc est déjà défini

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

# nc est déjà défini

#Create Python class by calling wrapped initialise routines
my_model.kernel.my_porosity = df2d.wrapping.m_model.porosity_data(my_model.kernel.mesh)

#Allocate and get initial values from Fortran
my_model.kernel.my_porosity.nland = nland # Utilisez nland défini précédemment
df2d.wrapping.call_model.init_porosity(my_model.kernel)

#Provide values, on top of initial ones from Fortran initialization routine, in Python structure
for i in range(nland) :
    my_model.kernel.my_porosity.phi[i] = 0.5
    my_model.kernel.my_porosity.hbanks[i] = 18.0 # Assurez-vous que c'est un float
    my_model.kernel.my_porosity.gamma[i] = 2.0  # Assurez-vous que c'est un float
my_model.kernel.my_porosity.land[:] = np.arange(1, nland + 1) # Plus pythonic que range

##########
# Hydraulic states
##########

my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
my_model.kernel.dof0 = my_model.kernel.dof

for i in range(nc) :
    my_model.kernel.dof0.h[i] = 1.0 + random.randint(0,1)*0.1 # Utilisez 1.0 pour float

my_model.kernel.dof0.u[:] = 0.0
my_model.kernel.dof0.v[:] = 0.0

# Définition de la bathymétrie (zb).
# C'est ici que vous pouvez définir le profil du fond du canal.
# Exemple d'une pente simple :
# Définition de la bathymétrie (zb).
# C'est ici que vous pouvez définir le profil du fond du canal.
# Exemple d'une pente simple :

print("Setting flat bathymetry (zb = 0.0 m for all cells).")
zb = np.zeros(nc, dtype=float)
zb[:] = 0.0 # Fond plat à une altitude de 0.0 mètre

# Assigner la bathymétrie au maillage Fortran
my_model.kernel.mesh.z[:] = zb

##############################################################################################################

##########################################
# Run
##########################################

df2d.wrapping.call_model.init_fortran(my_model.kernel)
df2d.wrapping.call_model.run(my_model.kernel, arg = "direct")

# IMPORTANT : Après le run, utilisez my_model.save_all() pour que les outputs soient générés
my_model.save_all()

if (os.path.isdir("./obs")):
    shutil.rmtree('./obs')

##############################################################################################################

########################
# Outputs from python (Visualisation avec plot_field)
########################

# Instancier explicitement l'objet Outputs et le lier au modèle
try:
    # Tentative 1: df2d.Outputs
    my_model.outputs = df2d.Outputs(my_model)
    my_model.outputs.load_outputs(custom_config=my_model.config.get_config())
    print("Outputs object successfully created and loaded from df2d.Outputs.")
except AttributeError:
    # Tentative 2: df2d.postprocess.Outputs (la plus probable pour votre version)
    try:
        import dassflow2d.postprocess # Importez le module postprocess
        my_model.outputs = dassflow2d.postprocess.Outputs(my_model)
        my_model.outputs.load_outputs(custom_config=my_model.config.get_config())
        print("Outputs object successfully created and loaded from dassflow2d.postprocess.Outputs.")
    except Exception as e_postprocess:
        print(f"Error creating/loading Outputs object from dassflow2d.postprocess: {e_postprocess}")
        print("Please check your DassFlow2D installation and module structure.")
        sys.exit("Cannot proceed with plotting without Outputs object.")
except Exception as e:
    print(f"Error creating/loading Outputs object (initial attempt df2d.Outputs): {e}")
    print("Please check your DassFlow2D installation and module structure.")
    sys.exit("Cannot proceed with plotting without Outputs object.")

# Calcul de la moyenne et de la déviation (si nécessaire pour un tracé spécifique)
# Utilisez my_model.outputs.result.h pour obtenir les données du dernier pas de temps
h_final_data = my_model.outputs.result.h[:, -1] # Prend la dernière colonne (dernier temps)
h_final_mean = np.mean(h_final_data)
h_deviation_from_mean = h_final_data - h_final_mean


# --- Visualisation des champs avec my_model.outputs.result.plot_field ---

# 1. Tracé de la Bathymétrie (zb)
plotter = my_model.outputs.result.plot_field(
    my_mesh = my_model.meshing.mesh_pyvista,
    what = "bathy", # Ou "zb" selon comment c'est nommé dans vos fichiers VTK/HDF5
    title_plot = "Bathymetry elevation",
    title_scale_bar ="z [m]",
    notebook = False
)
plotter.show()
print("Bathymetry plot displayed.")


# 2. Tracé de la Porosité (phi)
# Le champ 'phi' n'est souvent pas un champ VTK standard écrit par DassFlow.
# Nous le passons donc via 'my_scalar'.
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             my_scalar = my_model.kernel.my_porosity.phi,
                                             title_plot = "Porosity (from Python kernel)",
                                             title_scale_bar ="phi [-]",
                                             notebook = False )
plotter.show()
print("Porosity plot displayed.")


# 3. Tracé des déviations de hauteur d'eau
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             my_scalar = h_deviation_from_mean,
                                             title_plot = "Water depth deviations from mean (final state)",
                                             title_scale_bar ="h_dev [m]",
                                             notebook = False )
plotter.show()
print("Water depth deviations plot displayed.")


# 4. Tracé de la Hauteur d'eau (h) à l'état final
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "h",
                                             when = -1, # -1 pour le dernier pas de temps
                                             title_scale_bar ="h [m]",
                                             title_plot = "Water depth (final state)",
                                             notebook = False )
plotter.show()
print("Water depth (final state) plot displayed.")


# 5. Tracé de la vitesse u à l'état final
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "u",
                                             when = -1,
                                             title_scale_bar ="u [m/s]",
                                             title_plot = "Velocity u (final state)",
                                             notebook = False )
plotter.show()
print("Velocity u (final state) plot displayed.")


# 6. Tracé de la vitesse v à l'état final
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "v",
                                             when = -1,
                                             title_scale_bar ="v [m/s]",
                                             title_plot = "Velocity v (final state)",
                                             notebook = False )
plotter.show()
print("Velocity v (final state) plot displayed.")


df2d.wrapping.call_model.clean_model(my_model.kernel)
print("Model cleaned.")

# COMMENTEZ OU SUPPRIMEZ LES LIGNES DE SUPPRESSION DE MAILLAGE.
# Puisque nous utilisons un maillage statique, nous ne voulons pas le supprimer après l'exécution.
# delete_mesh(generated_mesh_filename, bin_dir=bin_dir) # Ancienne ligne
# print(f"Mesh file {generated_mesh_filename} deleted.")
