import gmsh
import dassflow2d as df2d


gmsh.initialize()
# ---------------------------------------------------------- #
#   Génération d'une géométrie simple (carré tordu)
# ---------------------------------------------------------- #
# Génération de Points géométriques avec coordonnée Z
p1 = gmsh.model.occ.addPoint(0, 0, 15)
p2 = gmsh.model.occ.addPoint(1000, 0, 5)
p3 = gmsh.model.occ.addPoint(1000, 100, 5)
p4 = gmsh.model.occ.addPoint(0, 100, 15)

# Génération de Splines géométriques
l_bas = gmsh.model.occ.addSpline([p1, p2])   # Mur : par défault
l_droit = gmsh.model.occ.addSpline([p2, p3]) # Sortie
l_haut = gmsh.model.occ.addSpline([p3, p4])  # Mur : par défault
l_gauche = gmsh.model.occ.addSpline([p4, p1]) # Entrée

# Génération de Surface géométrique
cl = gmsh.model.occ.addCurveLoop([l_bas, l_droit, l_haut, l_gauche])
s = gmsh.model.occ.addSurfaceFilling(cl)
gmsh.model.occ.synchronize()

# Réglages d'affichage avant ouverture de l'UI
gmsh.option.setNumber("Geometry.Points", 0)
gmsh.option.setNumber("Geometry.Lines", 1)
gmsh.option.setNumber("Geometry.Surfaces", 1)
gmsh.option.setNumber("Geometry.SurfaceType", 2)
gmsh.fltk.run() # Affichage intermédiaire sur l'interface GMSH

# ---------------------------------------------------------- #
#   /!\ IMPORTANT /!\
#  Définition des groupes physiques selon les conventions
# ---------------------------------------------------------- #

# 1. Attribution d'un groupe physique "DOMAIN" pour l'ensemble des surfaces ([s1,s2,etc], ici : s est la seule surface)

gmsh.model.addPhysicalGroup(2, [s], tag=1, name="DOMAIN")

gmsh.model.addPhysicalGroup(1, [l_gauche], tag=2, name="INLET:discharg1")

gmsh.model.addPhysicalGroup(1, [l_droit], tag=4, name="OUTLET:ratcurve")


# Option 2 : maillage structuré
gmsh.model.mesh.setTransfiniteSurface(s)
gmsh.model.mesh.setRecombine(2, s) # Optionnel : recombinaison en quadrangles
gmsh.model.mesh.setTransfiniteCurve(1, 100) # Optionnel : raffinement sur les courbes imposé
gmsh.model.mesh.setTransfiniteCurve(3, 100) # --> Le nombre de noeuds sur les courbes opposées doit être identique
gmsh.model.mesh.setTransfiniteCurve(2, 4) # Optionnel : recombinaison en quadrangles
gmsh.model.mesh.setTransfiniteCurve(4, 4) # /
gmsh.model.mesh.generate(2)


gmsh.option.setNumber("Geometry.Points", 0)
gmsh.option.setNumber("Geometry.Lines", 0)
gmsh.option.setNumber("Geometry.Surfaces", 0)
gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
gmsh.fltk.run()


gmsh.write("mesh_simple_channel_lineaire/mesh_tuto.msh")
gmsh.finalize()

from GMSH_to_df2d_librairy import convert_gmsh_to_df2d

Input_mesh = "mesh_simple_channel_lineaire/mesh_tuto.msh" # Doit être dans le dossier courant
folder_output="mesh_simple_channel_lineaire"

print(f"=== CONVERSION Test for {Input_mesh} ===")
convert_gmsh_to_df2d(Input_MSH=Input_mesh, Sim_duration=14400, Default_Q=10.0, Default_H=1.0, Default_Z=1.0, folder_output=folder_output)
print(f"=== CONVERSION Terminee ===")



folder_output="mesh_simple_channel_lineaire"

from Update_files_tools import apply_bathymetry, update_boundary_files
import os
# -------------------------------------------------------------------------- #
#           EXEMPLE D'UTILISATION DE LA FONCTION apply_bathymetry
# -------------------------------------------------------------------------- #
geo_file = os.path.join(folder_output, "mesh_simple_lineaire.geo")

def f(x, y):
    return (-0.25 * x) + 40 

print(f"--- Application de la nouvelle bathymétrie depuis la fonction f(x,y) ---")
apply_bathymetry(f, geo_file)
print(f"--- Bathymétrie appliquée et fichier {geo_file} mis à jour ---")


