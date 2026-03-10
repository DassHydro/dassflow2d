import gmsh
from GMSH_to_df2d_librairy import convert_gmsh_to_df2d

Input_mesh = "mesh_tuto.msh" # Doit être dans le dossier courant
folder_output="meshtuto_bin"

#---------------------------------------- #
# Visualisation du maillage avec Gmsh
gmsh.initialize()
gmsh.open(Input_mesh)
gmsh.fltk.run()
gmsh.finalize()
#---------------------------------------- #

print(f"=== CONVERSION Test for {Input_mesh} ===")
convert_gmsh_to_df2d(Input_MSH=Input_mesh, Sim_duration=14400, Default_Q=10.0, Default_H=1.0, Default_Z=1.0, folder_output=folder_output)
print(f"=== CONVERSION Terminee ===")
