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

initial_code_dir = os.getcwd()
bin_dir = os.path.join(initial_code_dir, "bin_A")

##########
# Initialise bin (et assure les dossiers nécessaires)
##########

os.chdir(initial_code_dir)

print("Cleaning old results (if makefiles are configured)...")
os.system("make cleanres")        # removes all in bin_dir/res directory
os.system("make cleanmsh")        # removes all in bin_dir/msh directory
os.system("make cleanmin")        # removes all in bin_dir/min directory

os.chdir(bin_dir)

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

mesh_name = 'channel.geo'
# Nous utilisons maintenant un maillage statique, donc pas de génération.
# Assurez-vous que le fichier existe dans bin_A.
mesh_full_path_for_dassflow_read = os.path.join(bin_dir, mesh_name)
if not os.path.exists(mesh_full_path_for_dassflow_read):
    print(f"ERROR: Static mesh file '{mesh_name}' not found in '{bin_dir}'. Please ensure it was copied correctly.")
    sys.exit(1) 
else:
    print(f"Using static mesh file '{mesh_name}' from '{bin_dir}'. Skipping generation.")

##########
# Model
##########

ts = 1400
dtw = ts*0.2

use_porosity = 1

input_params={ "mesh_name": mesh_name,
                "ts": ts ,
                "use_obs":'0',
                "use_UVobs":'0',
                "use_Zobs":'0',
                "use_porosity": use_porosity,
                
                "w_obs":'1',
                "w_vtk":'1', 
                "w_gnuplot":'1',
                "w_tecplot":'0',
                "adapt_dt":'1',

                "dt":'0.1',

                "dtw": dtw ,
                "dta":"100",
                
                "bc_infil":"0",
                "bc_rain":"0",}

df2d.wrapping.read_input(os.path.join(bin_dir,"input.txt"))
my_model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = os.path.join(bin_dir,"res","simu.hdf5") , run_type = "direct", clean = True, custom_config=input_params)
# Reinitialize mesh from new updated geo file
my_model.init_mesh()
# my_model.meshing.plot() 

# Input parameters for the generation of observations
Config = df2d.core.config.Config()
Config.set(custom_config = input_params)
nc = my_model.kernel.mesh.nc    # number of cells 
nland = nc

##########
# Friction
##########

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

#Create Python class by calling wrapped initialise routines
my_model.kernel.my_porosity = df2d.wrapping.m_model.porosity_data(my_model.kernel.mesh)

#Allocate and get initial values from Fortran
my_model.kernel.my_porosity.nland = nland
df2d.wrapping.call_model.init_porosity(my_model.kernel)

#Provide values, on top of initial ones from Fortran initialization routine, in Python structure
for i in range(nland) :
    my_model.kernel.my_porosity.phi[i] = 0.5
    my_model.kernel.my_porosity.hbanks[i] = 18.0 
    my_model.kernel.my_porosity.gamma[i] = 2.0 
my_model.kernel.my_porosity.land[:] = np.arange(1, nland + 1)

##########
# Hydraulic states
##########

my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
my_model.kernel.dof0 = my_model.kernel.dof

for i in range(nc) :
    my_model.kernel.dof0.h[i] = 1.0 + random.randint(0,1)*0.1 

my_model.kernel.dof0.u[:] = 0.0
my_model.kernel.dof0.v[:] = 0.0


##########################################
# Run
##########################################

df2d.wrapping.call_model.init_fortran(my_model.kernel)
df2d.wrapping.call_model.run(my_model.kernel, arg = "direct")

########################
# Outputs from python (Gestion et Visualisation)
########################

# 1. Instancier explicitement l'objet Outputs et le lier au modèle.
#    Ceci CRÉE l'objet 'outputs' et l'attache à 'my_model', mais ne CHARGE PAS encore de données.
try:
    # Tentative 1: df2d.Outputs
    my_model.outputs = df2d.Outputs(my_model)
    print("Outputs object successfully created and linked to model.")
except AttributeError:
    # Tentative 2: df2d.postprocess.Outputs (la plus probable pour votre version)
    try:
        import dassflow2d.postprocess # Importez le module postprocess
        my_model.outputs = dassflow2d.postprocess.Outputs(my_model)
        print("Outputs object successfully created and linked from dassflow2d.postprocess.Outputs.")
    except Exception as e_postprocess:
        print(f"Error creating Outputs object from dassflow2d.postprocess: {e_postprocess}")
        print("Please check your DassFlow2D installation and module structure.")
        sys.exit("Cannot proceed with plotting without Outputs object.")
except Exception as e:
    print(f"Error creating Outputs object (initial attempt df2d.Outputs): {e}")
    print("Please check your DassFlow2D installation and module structure.")
    sys.exit("Cannot proceed with plotting without Outputs object.")


# 2. IMPORTANT : Maintenant que l'objet Outputs est créé (my_model.outputs existe),
#    on peut appeler my_model.save_all(). C'est cette fonction qui va écrire les
#    résultats Fortran dans les fichiers de sortie (ex: simu.hdf5).
my_model.save_all()
print("Model results saved to HDF5/VTK files.")


# 3. Nettoyage des fichiers 'obs' (si présents et non désirés)
#    Ceci devrait se faire APRÈS la sauvegarde mais AVANT le chargement si ces
#    fichiers ne sont pas ceux que load_outputs doit lire.
if (os.path.isdir("./obs")): # Utilisez le chemin exact où 'obs' est créé
    shutil.rmtree('./obs')
    print("Removed temporary './obs' directory.")


# 4. Charger les résultats dans l'objet Outputs.
#    Ceci lit les fichiers de sortie qui viennent d'être écrits par my_model.save_all().
my_model.outputs.load_outputs(custom_config=my_model.config.get_config())
print("Outputs loaded into the Outputs object for plotting.")
##############################################################################################################

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
