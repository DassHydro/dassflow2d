#
# Exemple de génération d'un maillage GMSH respectant les conventions pour être converti 
# en fichier DASSFLOW2D par la librairie GMSH_to_df2d_librairy.py
# 
# Conventions principales :
# - Groupe physique "DOMAIN" pour le domaine de calcul : voir ci-dessous
# - Groupes physiques pour les conditions limites nommés INLET_*, OUTLET_* : voir ci-dessous
# - Maillage 2D en triangles ou quadrangles ou mélange des deux
# - Coordonnée Z des entités géométriques utilisée pour la bathymétrie initiale imposée
#    --> Par défault, dans GMSH_to_df2d_librairy.py :
#       - La bathymétrie d'une cellule est calculée comme la moyenne des Z de ses noeuds (A vérifier)
#       - La bathymétrie d'une cellule ghost est égale à la bathymétrie de sa cellule réelle associée
#    --> Bathymétrie modifiable ultérieurement avec modif_bathy.py si besoin
#       - A partir d'une fonction bathy(x,y) s'appliquant aux centres des cellules réelles et ghost
# - Noms des CL à respecter :
#       discharg1, discharg2, hpresc, zpresc, ratcurve, etc.
import gmsh

gmsh.initialize()
# ---------------------------------------------------------- #
#   Génération d'une géométrie simple (carré tordu)
# ---------------------------------------------------------- #
# Génération de Points géométriques avec coordonnée Z
p1 = gmsh.model.occ.addPoint(0, 0, 10)
p2 = gmsh.model.occ.addPoint(10, 0, 10)
p3 = gmsh.model.occ.addPoint(10, 10, 5)
p4 = gmsh.model.occ.addPoint(0, 10, 8)

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

# 2. Attribution de groupes physiques sur des courbes géométriques pour l'application de CL
# --- DEFINITION DES CONDITIONS AUX LIMITES --- 
#
# Conventions de nommage des groupes physiques :
#        name : "CATEGORIE:TYPE" ou "CATEGORIE:TYPE_id"* (* uniquement pour discharg_ et ratcurve)
#  --> CATEGORIE  : INLET ou OUTLET
#  --> TYPE       : type de condition limite (discharg1, hpresc, zpresc, ratcurve, etc.)
#  --> identifiant unique (optionnel) : _1, _2, etc. (si plusieurs conditions de même type)
#       --> /!\ En maintenance : vérifier la génération du fichier bc.txt dans la conversion
#
# Les conditions limites de type 'WALL' n'ont pas besoin d'être spécifiées (par défaut sur les autres bords non définis)
#

# Exemples :
# Entrée avec débit imposé (discharg1)
gmsh.model.addPhysicalGroup(1, [l_gauche], tag=2, name="INLET:discharg1")

#gmsh.model.addPhysicalGroup(1, [l_gauche], tag=2, name="INLET:discharg1_1") # /!\ A VERIFIER SI CA FONCTIONNE
#gmsh.model.addPhysicalGroup(1, [l_bas], tag=3, name="INLET:discharg1_2")

# Sortie avec hauteur imposée (hpresc)
gmsh.model.addPhysicalGroup(1, [l_droit], tag=4, name="OUTLET:hpresc")

# ---------------------------------------------------------- #
#       Génération automatique du maillage 2D
# Par défaut : triangles ; possibilité de recombiner en quadrangles
# ---------------------------------------------------------- #

# Option 1 : maillage simple (direct)
#gmsh.model.mesh.generate(2)
#gmsh.model.mesh.recombine() # Optionnel : recombinaison en quadrangles lorsque possible

# Option 2 : maillage structuré
gmsh.model.mesh.setTransfiniteSurface(s)
gmsh.model.mesh.setRecombine(2, s) # Optionnel : recombinaison en quadrangles
gmsh.model.mesh.setTransfiniteCurve(1, 5) # Optionnel : raffinement sur les courbes imposé
gmsh.model.mesh.setTransfiniteCurve(3, 5) # --> Le nombre de noeuds sur les courbes opposées doit être identique
gmsh.model.mesh.setTransfiniteCurve(2, 10) # Optionnel : recombinaison en quadrangles
gmsh.model.mesh.setTransfiniteCurve(4, 10) # /
gmsh.model.mesh.generate(2)


gmsh.option.setNumber("Geometry.Points", 0)
gmsh.option.setNumber("Geometry.Lines", 0)
gmsh.option.setNumber("Geometry.Surfaces", 0)
gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
gmsh.fltk.run()


gmsh.write("mesh_tuto.msh")
gmsh.finalize()