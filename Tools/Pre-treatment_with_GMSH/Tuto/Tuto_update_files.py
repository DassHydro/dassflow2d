import os
from Update_files_tools import apply_bathymetry, update_boundary_files

bin_folder = "meshtuto_bin" # Doit être dans le dossier courant
geo_file = os.path.join(bin_folder, "mesh_tuto.geo")

# -------------------------------------------------------------------------- #
#           EXEMPLE D'UTILISATION DE LA FONCTION apply_bathymetry
# -------------------------------------------------------------------------- #

def f(x, y):
    return (-0.5 * x) + 10

print(f"--- Application de la nouvelle bathymétrie depuis la fonction f(x,y) ---")
apply_bathymetry(f, geo_file)
print(f"--- Bathymétrie appliquée et fichier {geo_file} mis à jour ---")

# ---------------------------------------------------------------------------------- #
#           EXEMPLE D'UTILISATION DE  LA FONCTION update_boundary_files
# ---------------------------------------------------------------------------------- #


times = [0, 3600, 7200, 10800, 14400]  # en secondes
flows = [0.0, 5.0, 150.0, 30.0, 0.0]      # en m3/s

print(f"--- Mise à jour de l'hydrographe avec les nouvelles valeurs ---")
bc_hydro_file = os.path.join(bin_folder, "hydrograph.txt")
update_boundary_files(bc_hydro_file, target_group_id=1, times=times, values=flows)
print(f"--- Hydrographe {bc_hydro_file} mis à jour ---")

hpresc_values = [0.0, 5.0, 30.0, 20.0, 15.0]      # en mètres
print(f"--- Mise à jour de hpresc avec les nouvelles valeurs ---")
bc_hpresc_file = os.path.join(bin_folder, "hpresc.txt")
update_boundary_files(bc_hpresc_file, target_group_id=1, times=times, values=hpresc_values)
print(f"--- hpresc {bc_hpresc_file} mis à jour ---")

# ---------------------------------------------------------------------------------- #


