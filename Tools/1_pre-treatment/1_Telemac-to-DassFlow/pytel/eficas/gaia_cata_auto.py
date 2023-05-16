
# -*- coding: latin-1 -*-

from Accas import *
class DateJJMMAAAA:
  def __init__(self):
    self.ntuple=3

  def __convert__(self,valeur):
    if type(valeur) == types.StringType: return None
    if len(valeur) != self.ntuple: return None
    return valeur

  def info(self):
    return "Date : jj/mm/aaaa "

  __repr__=info
  __str__=info

class grma(GEOM):
  pass

import types
class Tuple:
  def __init__(self,ntuple):
    self.ntuple=ntuple

  def __convert__(self,valeur):
    if type(valeur) == types.StringType:
      return None
    if len(valeur) != self.ntuple:
      return None
    return valeur

  def info(self):
    return "Tuple de %s elements" % self.ntuple



JdC = JDC_CATA (code = 'GAIA',
                execmodul = None,
                )
# =======================================================================
# Catalog entry for the MAP function : c_pre_interfaceBody_mesh
# =======================================================================

VERSION_CATALOGUE="V8P2"
# -----------------------------------------------------------------------
COMPUTATION_ENVIRONMENT = PROC(nom= "COMPUTATION_ENVIRONMENT",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    GLOBAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TITLE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Titre du cas etudie.
Ce titre sera inscrit dans les sorties.""",
            ang = """Title of the case being considered.""",
        ),
#       -----------------------------------
        PARALLEL_PROCESSORS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [0],
            fr = """NOMBRE DE PROCESSEURS EN CALCUL PARALLELE
0 : 1 machine, compilation sans bibliotheque de parallelisme
1 : 1 machine, compilation avec bibliotheque de parallelisme
2 : 2 processeurs ou machines en parallele
etc...""",
            ang = """NUMBER OF PROCESSORS FOR PARALLEL PROCESSING
0 : 1 machine, compiling without parallel library
1 : 1 machine, compiling with a parallel library
2 : 2 processors or machines in parallel
etc....""",
        ),
#       -----------------------------------
        CHECKING_THE_MESH = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Si oui on appelle le sous-programme checkmesh qui verifie
la coherence du maillage, points superposes, etc.""",
            ang = """if this key word is equal to yes, a call to subroutine
checkmesh will look for errors in the mesh, superimposed points, etc.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 30,
            fr = """nombre maximal de frontieres differentes dans le maillage.
Sert au dimensionnement de la memoire, a augmenter si necessaire""",
            ang = """maximal number of boundaries in the mesh.
Used for dimensioning arrays. Can be increased if needed""",
        ),
    ),
#   -----------------------------------
    INPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DATA = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            VALIDATION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Option utilisee principalement pour le dossier de validation.
Si ce mot-cle vaut OUI, les resultats du calcul vont alors etre
compares aux valeurs du fichier de reference.
La comparaison est effectuee par le sous-programme VALIDA qui peut
etre modifie pour realiser, par exemple, une comparaison avec
une solution exacte.""",
                ang = """This option is primarily used for the validation
documents. If this keyword is equal to YES, the REFERENCE FILE
is then considered as a reference which the computation is
going to be compared with.
The comparison is made by the subroutine VALIDA, which can be
modified so as to include,for example,a comparison with an
exact solution.""",
            ),
#           -----------------------------------
            b_VALIDATIONG = BLOC(condition="VALIDATION == True",
#           -----------------------------------
#               -----------------------------------
                REFERENCE_FILE_FORMAT = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'TXM',
                    into = ['SERAFIN','SERAFIND','MED'],
                    defaut = 'SERAFIN',
                    fr = """Format du fichier de resultats du calcul precedent.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                    ang = """Previous computation results file format.
Possible values are:
- SERAFIN : classical single precision format in Telemac;
- SERAFIND: classical double precision format in Telemac;
- MED     : MED format based on HDF5""",
                ),
#               -----------------------------------
                REFERENCE_FILE = SIMP(statut ='f',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)'), max='**',
                    defaut = '',
                    fr = """Nom du fichier servant a valider le calcul.
Si VALIDATION = OUI, les resultats du calcul vont etre
comparees aux valeurs contenues dans ce fichier.
La comparaison est effectuee par le sous-programme VALIDA.""",
                    ang = """Name of the file used to validate the computation.
If VALIDATION = YES, the results of the computation will be
compared with the values of this file. The comparison is
made by the subroutine VALIDA.""",
                ),
            ),
        ),
    ),
#   -----------------------------------
    OUTPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        RESULTS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            NUMBER_OF_PRIVATE_ARRAYS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = [1],
                fr = """Nombre de tableaux mis a disposition de l utilisateur""",
                ang = """Number of arrays for own user programming""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
INTERNAL = PROC(nom= "INTERNAL",op = None,
# -----------------------------------------------------------------------
    UIinfo = {"groupes": ("CACHE")},
#   -----------------------------------
    RELEASE = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        defaut = 'TRUNK',
        fr = """Numero de version des bibliotheques utilisees par GAIA.""",
        ang = """Release of the libraries used by GAIA.""",
    ),
#   -----------------------------------
    DICTIONARY = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = 'gaia.dico',
        fr = """Dictionnaire des mots cles.""",
        ang = """Key word dictionary.""",
    ),
)
# -----------------------------------------------------------------------
GENERAL_PARAMETERS = PROC(nom= "GENERAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    DEBUGGER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [0],
        fr = """Pour imprimer la sequence des appels, mettre 1""",
        ang = """If 1, calls of subroutines will be printed in the listing""",
    ),
#   -----------------------------------
    LOCATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ORIGIN_COORDINATES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I', min= 2, max= 2,
            defaut = [0,0],
            fr = """Valeur en metres, utilise pour eviter les trop grands nombres,
transmis dans le format Selafin mais pas d''autre traitement pour
l''instant""",
            ang = """Value in metres, used to avoid large real numbers,
added in Selafin format, but so far no other treatment""",
        ),
    ),
)
# -----------------------------------------------------------------------
HYDRODYNAMICS = PROC(nom= "HYDRODYNAMICS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    TIDAL_FLATS_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Utilise si ''BANCS DECOUVRANTS'' est vrai
   1 : EQUATIONS RESOLUES PARTOUT AVEC CORRECTION
       SUR LES BANCS DECOUVRANTS
   2 : GEL DES ELEMENTS DECOUVRANTS
Il est conseille de choisir l''option 1 car elle permet de
conserver la masse.""",
            ang = """Used if ''TIDAL FLATS'' is true
   1 : EQUATIONS SOLVED EVERYWHERE WITH CORRECTION ON TIDAL FLATS
   2 : DRY ELEMENTS FROZEN
It is recommended to choose 1 since it ensures mass conservation.""",
        ),
#       -----------------------------------
        MINIMAL_VALUE_OF_THE_WATER_HEIGHT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-3,
            fr = """Fixe la valeur minimale de la hauteur d''eau.
Est utilise lorsque le mot cle BANCS DECOUVRANTS est egal a oui.""",
            ang = """Sets the minimum value of the water depth.
Is used when the keyword TIDAL FLATS is equal to yes.""",
        ),
#       -----------------------------------
        TIDAL_FLATS = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """permet de supprimer les tests sur les bancs decouvrants, dans
les cas ou l''on est certain qu''il n''y en aura pas.
En cas de doute : oui""",
            ang = """When no, the specific treatments for tidal flats
are by-passed.
This spares time, but of course you must be sure that you
have no tidal flats""",
        ),
    ),
)
# -----------------------------------------------------------------------
MISCELLANEOUS = PROC(nom= "MISCELLANEOUS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NESTOR = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Pour le couplage avec NESTOR""",
        ang = """For coupling with NESTOR""",
    ),
#   -----------------------------------
    NESTOR_ACTION_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier de commandes de nestor""",
        ang = """Name of the Nestor steering file""",
    ),
#   -----------------------------------
    NESTOR_POLYGON_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier de polygons de Nestor""",
        ang = """Name of the Nestor polygon file""",
    ),
#   -----------------------------------
    NESTOR_RESTART_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier de phydef-cf.cfg.ds de Nestor""",
        ang = """Name of the Nestor file phydef-cf.cfg.ds""",
    ),
#   -----------------------------------
    NESTOR_SURFACE_REFERENCE_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier de reference surface de Nestor""",
        ang = """Name of the Nestor file which contains the reference
         water surface""",
    ),
)
# -----------------------------------------------------------------------
GENERAL = PROC(nom= "GENERAL",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    FLUXLINE_INPUT_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Nom du fichier de fluxline""",
        ang = """Name of the Fluxline file""",
    ),
#   -----------------------------------
    FLUXLINE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """FLUXLINE""",
        ang = """Use Fluxline to compute flux over lines""",
    ),
#   -----------------------------------
    CONTROL_SECTIONS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I', min= 3, max= 3,
        fr = """Couples de points (numeros globaux dans le maillage) entre
lesquels les debits instantanes et cumules seront donnes.""",
        ang = """Couples of points (global numbers in the mesh) defining sections
 where the instantaneous and cumulated discharges will be given""",
    ),
#   -----------------------------------
    SECONDARY_CURRENTS = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Pour prendre en compte les courants secondaires""",
        ang = """using the parametrisation for secondary currents""",
    ),
#   -----------------------------------
    SECONDARY_CURRENTS_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Pour prendre en compte les courants secondaires avec FICHIER""",
        ang = """The radii needed for the parametrisation of secondary currents
are read from SELAFIN file""",
    ),
#   -----------------------------------
    EFFECT_OF_WAVES = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Prend en compte l''effet de la houle sur le transport solide""",
        ang = """Takes into account the effect of waves""",
    ),
#   -----------------------------------
    TYPE_OF_WAVES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 2,
        fr = """Utilise pour calculer Uw
 Si Uw est calcule par Tomawac, choisir 2 (defaut)
 1= vagues regulieres (monochormatique)
 2= vagues irregulieres (spectre)""",
        ang = """is used to calculate Uw
 if Uw is calculated with Tomawac, choose 2 (default)
 1= regular (monochromatic) waves
 2= irregular (spectral) waves""",
    ),
)
# -----------------------------------------------------------------------
COHESIVE_AND_NON_COHESIVE = PROC(nom= "COHESIVE_AND_NON_COHESIVE",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    BED_MODEL = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """3 types de modeles de lit sont disponibles :
1 : multicouche (couche active automatique si plusieurs classes)
2 : multicouche avec consolidation
3 : consolidation basee sur la theorie de Gibson""",
        ang = """3 kinds of bed model are available:
1 : multilayers (automatic active layer if several classes)
2 : multilayer with consolidation
3 : consolidation model based on Gibson theory""",
    ),
#   -----------------------------------
    CHARRIAGE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        CLASSES_SHIELDS_PARAMETERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [-9.,-9.],
            fr = """Utilise pour determiner la valeur de la contrainte critique
d''entrainement (sediments non-cohesifs). En multiclasse, specifier la
valeur pour chaque classe.  Donner valeur negative pour les sediments
cohesifs.""",
            ang = """Used to determine the critical bed shear stress value
(non-cohesive sediments).
For multi grain size,
the shields parameter needs to be specified for each class.
It is necessary to give a negative value
in the parameter file for cohesive sediments.""",
        ),
    ),
#   -----------------------------------
    SUSPENSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        CLASSES_SETTLING_VELOCITIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [-9,-9.],
            fr = """Fixe la valeur de vitesse de chute pour chaque sediment. Donner une
valuer negative pour utiliser les formules de Stokes, Zanke ou Van
Rijn (appliquées selon la taille des grains).""",
            ang = """Sets the value of settling velocity for every sediment. Give a
negative value to use the Stokes, Zanke or Van Rijn formulae
(depending on the grain size)""",
        ),
#       -----------------------------------
        EQUILIBRIUM_INFLOW_CONCENTRATION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Impose la concentration en entree du domaine dans les cas 2D.
Pour les sediments non-cohesifs, la concentration d''equilibre est
calculee en accord avec la formule de transport pour tous les sable.""",
            ang = """Imposes the equilibrium concentration at the inlet boundaries in 2D
cases.
For non cohesive sediments, the equilibrium near bed concentration is
computed with respect to the suspension transport formula for all
sands.""",
        ),
    ),
#   -----------------------------------
    BED_MATERIAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        CLASSES_SEDIMENT_DIAMETERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [.01,.01],
            fr = """Sets value of diameter dm for particular size class.""",
            ang = """Sets value of diameter dm for particular size class.""",
        ),
#       -----------------------------------
        CLASSES_HIDING_FACTOR = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [1.,1.],
            fr = """Fixe la valeur du facteur de pavage par classe
granulometrique""",
            ang = """Sets value of hiding factor for particular size class.""",
        ),
#       -----------------------------------
        CLASSES_INITIAL_FRACTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [1.,0.],
            fr = """Fixe la fraction initiale de chaque classe sédimentologique
dans le mélange. Attention, la somme de AVA0 sur toutes les
classes doit valoir 1.""",
            ang = """Sets the value of the initial fraction of each sediment class.
Beware that the sum over all classes must be equal to 1.""",
        ),
#       -----------------------------------
        LAYERS_INITIAL_THICKNESS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Epaisseurs initiales des sediments (m).""",
            ang = """Sediment layers thickness (m) for initialisation.""",
        ),
#       -----------------------------------
        CLASSES_TYPE_OF_SEDIMENT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            fr = """Liste des types de sediment: cohesif (CO) ou non cohesif (NCO).""",
            ang = """ Liste of types of sediment: cohesive (CO) ou non cohesive (NCO).""",
        ),
    ),
)
# -----------------------------------------------------------------------
BED_STRUCTURE = PROC(nom= "BED_STRUCTURE",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    INITIALIZATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_LAYERS_FOR_INITIAL_STRATIFICATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Nombre de couche du lit initial, defaut NUMSTRAT=1""",
            ang = """Number of layers for initial stratification, default NUMSTRAT=1""",
        ),
    ),
)
# -----------------------------------------------------------------------
USELESS = PROC(nom= "USELESS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    VECTOR_LENGTH = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """longueur du vecteur pour les machines vectorielles.""",
        ang = """vector length on vector machines.""",
    ),
#   -----------------------------------
    STEERING_FILE = SIMP(statut ='o',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier contenant les parametres du calcul
a realiser. Il peut-etre ecrit par l''utilisateur avec EDAMOX.""",
        ang = """Name of the file containing the parameters
of the computation. Could be written by the user with EDAMOX.""",
    ),
)
# -----------------------------------------------------------------------
INPUT_OUTPUT__FILES = PROC(nom= "INPUT_OUTPUT__FILES",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    GEOMETRY_FILE_FORMAT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ['SERAFIN','SERAFIND','MED'],
        defaut = 'SERAFIN',
        fr = """Format du fichier de geometrie.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
        ang = """Geometry file format.
Possible values are:
- SERAFIN : classical single precision format in Telemac;
- SERAFIND: classical double precision format in Telemac;
- MED     : MED format based on HDF5""",
    ),
#   -----------------------------------
    NAMES = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        GEOMETRY_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'), max='**',
            fr = """Nom du fichier contenant le maillage du calcul a realiser.""",
            ang = """Name of the file containing the mesh. This file may also
contain the topography and the friction coefficients.""",
        ),
#       -----------------------------------
        NAMES_OF_PRIVATE_VARIABLES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min= 2, max= 2,
            fr = """Noms des variables privees en 32 caracteres, 16 pour le nom
         16 pour l''unite. Elles correspondent au bloc PRIVE
         et peuvent etre lues dans le fichier de geometrie si elles
         y sont presentes avec leur nom""",
            ang = """Name of private variables in 32 characters, 16 for the name,
         16 for the unit. They are stored in the block PRIVE and
         can be read in the geometry file if they are here with their
         name""",
        ),
#       -----------------------------------
        BEDLOAD_BOUNDARIES_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'), max='**',
            defaut = '',
            fr = """Fichier de variations en temps des conditions aux limites en charriage.
Les donnees de ce fichier sont sur le canal GAI\_FILES(GAILIQ)%LU.""",
            ang = """Variations in time of boundary conditions in bedload. Data of this file
are read on channel GAI\_FILES(GAILIQ)%LU.""",
        ),
    ),
)
# -----------------------------------------------------------------------
RESULTS = PROC(nom= "RESULTS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    RESULTS_FILE_FORMAT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ['SERAFIN','SERAFIND','MED'],
        defaut = 'SERAFIN',
        fr = """Format du fichier de resultats.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
        ang = """Results file format. Possible values are:
- SERAFIN : classical single precision format in Telemac;
- SERAFIND: classical double precision format in Telemac;
- MED     : MED format based on HDF5""",
    ),
#   -----------------------------------
    RESULTS_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
        fr = """Nom du fichier dans lequel seront ecrits les resultats avec
une periodicite donnee par le mot cle PERIODE DE SORTIE GRAPHIQUE
dans telemac2d ou telemac3d.""",
        ang = """Name of the file into wich the computation results shall be
written, the periodicity being given by the keyword
GRAPHIC PRINTOUT PERIOD in telemac2d or telemac3d.""",
    ),
#   -----------------------------------
    VARIABLES_TO_BE_PRINTED = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM', max='**',
        into = ["velocity along x axis (m/s)","velocity along y axis (m/s)","wawe celerity (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","Froude number","scalar flowrate of fluid (m2/s)","flowrate along x axis (m2/s)","flowrate along y axis (m2/s)","bed-load discharge (m2/s)","bed-load discharge along x axis (m2/s)","bed-load discharge along y axis (m2/s)","bottom evolution (m)","non erodable bottom","bed friction coefficient (m if Nikuradse)","mean bottom friction (N/m2)","wave angle with axis Oy (deg)","wave height","wave period","fraction of non cohesive sediment of class*,*n* layer","fraction of non cohesive sediment of class1,1n* layer","fraction of non cohesive sediment of class2,2n* layer","fraction of non cohesive sediment of class3,3n* layer","fraction of non cohesive sediment of class4,4n* layer","fraction of non cohesive sediment of class5,5n* layer","fraction of non cohesive sediment of class6,6n* layer","fraction of non cohesive sediment of class7,7n* layer","fraction of non cohesive sediment of class8,8n* layer","fraction of non cohesive sediment of class9,9n* layer","fraction of non cohesive sediment of class10,10n* layer","fraction of non cohesive sediment of class11,11n* layer","fraction of non cohesive sediment of class12,12n* layer","fraction of non cohesive sediment of class13,13n* layer","fraction of non cohesive sediment of class14,14n* layer","fraction of non cohesive sediment of class15,15n* layer","fraction of non cohesive sediment of class16,16n* layer","fraction of non cohesive sediment of class17,17n* layer","fraction of non cohesive sediment of class18,18n* layer","fraction of non cohesive sediment of class19,19n* layer","fraction of non cohesive sediment of class**,**n* layer","fraction of non cohesive sediment of class*,*n1 layer","fraction of non cohesive sediment of class1,1n1 layer","fraction of non cohesive sediment of class2,2n1 layer","fraction of non cohesive sediment of class3,3n1 layer","fraction of non cohesive sediment of class4,4n1 layer","fraction of non cohesive sediment of class5,5n1 layer","fraction of non cohesive sediment of class6,6n1 layer","fraction of non cohesive sediment of class7,7n1 layer","fraction of non cohesive sediment of class8,8n1 layer","fraction of non cohesive sediment of class9,9n1 layer","fraction of non cohesive sediment of class10,10n1 layer","fraction of non cohesive sediment of class11,11n1 layer","fraction of non cohesive sediment of class12,12n1 layer","fraction of non cohesive sediment of class13,13n1 layer","fraction of non cohesive sediment of class14,14n1 layer","fraction of non cohesive sediment of class15,15n1 layer","fraction of non cohesive sediment of class16,16n1 layer","fraction of non cohesive sediment of class17,17n1 layer","fraction of non cohesive sediment of class18,18n1 layer","fraction of non cohesive sediment of class19,19n1 layer","fraction of non cohesive sediment of class**,**n1 layer","fraction of non cohesive sediment of class*,*n2 layer","fraction of non cohesive sediment of class1,1n2 layer","fraction of non cohesive sediment of class2,2n2 layer","fraction of non cohesive sediment of class3,3n2 layer","fraction of non cohesive sediment of class4,4n2 layer","fraction of non cohesive sediment of class5,5n2 layer","fraction of non cohesive sediment of class6,6n2 layer","fraction of non cohesive sediment of class7,7n2 layer","fraction of non cohesive sediment of class8,8n2 layer","fraction of non cohesive sediment of class9,9n2 layer","fraction of non cohesive sediment of class10,10n2 layer","fraction of non cohesive sediment of class11,11n2 layer","fraction of non cohesive sediment of class12,12n2 layer","fraction of non cohesive sediment of class13,13n2 layer","fraction of non cohesive sediment of class14,14n2 layer","fraction of non cohesive sediment of class15,15n2 layer","fraction of non cohesive sediment of class16,16n2 layer","fraction of non cohesive sediment of class17,17n2 layer","fraction of non cohesive sediment of class18,18n2 layer","fraction of non cohesive sediment of class19,19n2 layer","fraction of non cohesive sediment of class**,**n2 layer","fraction of non cohesive sediment of class*,*n3 layer","fraction of non cohesive sediment of class1,1n3 layer","fraction of non cohesive sediment of class2,2n3 layer","fraction of non cohesive sediment of class3,3n3 layer","fraction of non cohesive sediment of class4,4n3 layer","fraction of non cohesive sediment of class5,5n3 layer","fraction of non cohesive sediment of class6,6n3 layer","fraction of non cohesive sediment of class7,7n3 layer","fraction of non cohesive sediment of class8,8n3 layer","fraction of non cohesive sediment of class9,9n3 layer","fraction of non cohesive sediment of class10,10n3 layer","fraction of non cohesive sediment of class11,11n3 layer","fraction of non cohesive sediment of class12,12n3 layer","fraction of non cohesive sediment of class13,13n3 layer","fraction of non cohesive sediment of class14,14n3 layer","fraction of non cohesive sediment of class15,15n3 layer","fraction of non cohesive sediment of class16,16n3 layer","fraction of non cohesive sediment of class17,17n3 layer","fraction of non cohesive sediment of class18,18n3 layer","fraction of non cohesive sediment of class19,19n3 layer","fraction of non cohesive sediment of class**,**n3 layer","fraction of non cohesive sediment of class*,*n4 layer","fraction of non cohesive sediment of class1,1n4 layer","fraction of non cohesive sediment of class2,2n4 layer","fraction of non cohesive sediment of class3,3n4 layer","fraction of non cohesive sediment of class4,4n4 layer","fraction of non cohesive sediment of class5,5n4 layer","fraction of non cohesive sediment of class6,6n4 layer","fraction of non cohesive sediment of class7,7n4 layer","fraction of non cohesive sediment of class8,8n4 layer","fraction of non cohesive sediment of class9,9n4 layer","fraction of non cohesive sediment of class10,10n4 layer","fraction of non cohesive sediment of class11,11n4 layer","fraction of non cohesive sediment of class12,12n4 layer","fraction of non cohesive sediment of class13,13n4 layer","fraction of non cohesive sediment of class14,14n4 layer","fraction of non cohesive sediment of class15,15n4 layer","fraction of non cohesive sediment of class16,16n4 layer","fraction of non cohesive sediment of class17,17n4 layer","fraction of non cohesive sediment of class18,18n4 layer","fraction of non cohesive sediment of class19,19n4 layer","fraction of non cohesive sediment of class**,**n4 layer","fraction of non cohesive sediment of class*,*n5 layer","fraction of non cohesive sediment of class1,1n5 layer","fraction of non cohesive sediment of class2,2n5 layer","fraction of non cohesive sediment of class3,3n5 layer","fraction of non cohesive sediment of class4,4n5 layer","fraction of non cohesive sediment of class5,5n5 layer","fraction of non cohesive sediment of class6,6n5 layer","fraction of non cohesive sediment of class7,7n5 layer","fraction of non cohesive sediment of class8,8n5 layer","fraction of non cohesive sediment of class9,9n5 layer","fraction of non cohesive sediment of class10,10n5 layer","fraction of non cohesive sediment of class11,11n5 layer","fraction of non cohesive sediment of class12,12n5 layer","fraction of non cohesive sediment of class13,13n5 layer","fraction of non cohesive sediment of class14,14n5 layer","fraction of non cohesive sediment of class15,15n5 layer","fraction of non cohesive sediment of class16,16n5 layer","fraction of non cohesive sediment of class17,17n5 layer","fraction of non cohesive sediment of class18,18n5 layer","fraction of non cohesive sediment of class19,19n5 layer","fraction of non cohesive sediment of class**,**n5 layer","fraction of non cohesive sediment of class*,*n6 layer","fraction of non cohesive sediment of class1,1n6 layer","fraction of non cohesive sediment of class2,2n6 layer","fraction of non cohesive sediment of class3,3n6 layer","fraction of non cohesive sediment of class4,4n6 layer","fraction of non cohesive sediment of class5,5n6 layer","fraction of non cohesive sediment of class6,6n6 layer","fraction of non cohesive sediment of class7,7n6 layer","fraction of non cohesive sediment of class8,8n6 layer","fraction of non cohesive sediment of class9,9n6 layer","fraction of non cohesive sediment of class10,10n6 layer","fraction of non cohesive sediment of class11,11n6 layer","fraction of non cohesive sediment of class12,12n6 layer","fraction of non cohesive sediment of class13,13n6 layer","fraction of non cohesive sediment of class14,14n6 layer","fraction of non cohesive sediment of class15,15n6 layer","fraction of non cohesive sediment of class16,16n6 layer","fraction of non cohesive sediment of class17,17n6 layer","fraction of non cohesive sediment of class18,18n6 layer","fraction of non cohesive sediment of class19,19n6 layer","fraction of non cohesive sediment of class**,**n6 layer","fraction of non cohesive sediment of class*,*n7 layer","fraction of non cohesive sediment of class1,1n7 layer","fraction of non cohesive sediment of class2,2n7 layer","fraction of non cohesive sediment of class3,3n7 layer","fraction of non cohesive sediment of class4,4n7 layer","fraction of non cohesive sediment of class5,5n7 layer","fraction of non cohesive sediment of class6,6n7 layer","fraction of non cohesive sediment of class7,7n7 layer","fraction of non cohesive sediment of class8,8n7 layer","fraction of non cohesive sediment of class9,9n7 layer","fraction of non cohesive sediment of class10,10n7 layer","fraction of non cohesive sediment of class11,11n7 layer","fraction of non cohesive sediment of class12,12n7 layer","fraction of non cohesive sediment of class13,13n7 layer","fraction of non cohesive sediment of class14,14n7 layer","fraction of non cohesive sediment of class15,15n7 layer","fraction of non cohesive sediment of class16,16n7 layer","fraction of non cohesive sediment of class17,17n7 layer","fraction of non cohesive sediment of class18,18n7 layer","fraction of non cohesive sediment of class19,19n7 layer","fraction of non cohesive sediment of class**,**n7 layer","fraction of non cohesive sediment of class*,*n8 layer","fraction of non cohesive sediment of class1,1n8 layer","fraction of non cohesive sediment of class2,2n8 layer","fraction of non cohesive sediment of class3,3n8 layer","fraction of non cohesive sediment of class4,4n8 layer","fraction of non cohesive sediment of class5,5n8 layer","fraction of non cohesive sediment of class6,6n8 layer","fraction of non cohesive sediment of class7,7n8 layer","fraction of non cohesive sediment of class8,8n8 layer","fraction of non cohesive sediment of class9,9n8 layer","fraction of non cohesive sediment of class10,10n8 layer","fraction of non cohesive sediment of class11,11n8 layer","fraction of non cohesive sediment of class12,12n8 layer","fraction of non cohesive sediment of class13,13n8 layer","fraction of non cohesive sediment of class14,14n8 layer","fraction of non cohesive sediment of class15,15n8 layer","fraction of non cohesive sediment of class16,16n8 layer","fraction of non cohesive sediment of class17,17n8 layer","fraction of non cohesive sediment of class18,18n8 layer","fraction of non cohesive sediment of class19,19n8 layer","fraction of non cohesive sediment of class**,**n8 layer","fraction of non cohesive sediment of class*,*n9 layer","fraction of non cohesive sediment of class1,1n9 layer","fraction of non cohesive sediment of class2,2n9 layer","fraction of non cohesive sediment of class3,3n9 layer","fraction of non cohesive sediment of class4,4n9 layer","fraction of non cohesive sediment of class5,5n9 layer","fraction of non cohesive sediment of class6,6n9 layer","fraction of non cohesive sediment of class7,7n9 layer","fraction of non cohesive sediment of class8,8n9 layer","fraction of non cohesive sediment of class9,9n9 layer","fraction of non cohesive sediment of class10,10n9 layer","fraction of non cohesive sediment of class11,11n9 layer","fraction of non cohesive sediment of class12,12n9 layer","fraction of non cohesive sediment of class13,13n9 layer","fraction of non cohesive sediment of class14,14n9 layer","fraction of non cohesive sediment of class15,15n9 layer","fraction of non cohesive sediment of class16,16n9 layer","fraction of non cohesive sediment of class17,17n9 layer","fraction of non cohesive sediment of class18,18n9 layer","fraction of non cohesive sediment of class19,19n9 layer","fraction of non cohesive sediment of class**,**n9 layer","fraction of non cohesive sediment of class*,*n10 layer","fraction of non cohesive sediment of class1,1n10 layer","fraction of non cohesive sediment of class2,2n10 layer","fraction of non cohesive sediment of class3,3n10 layer","fraction of non cohesive sediment of class4,4n10 layer","fraction of non cohesive sediment of class5,5n10 layer","fraction of non cohesive sediment of class6,6n10 layer","fraction of non cohesive sediment of class7,7n10 layer","fraction of non cohesive sediment of class8,8n10 layer","fraction of non cohesive sediment of class9,9n10 layer","fraction of non cohesive sediment of class10,10n10 layer","fraction of non cohesive sediment of class11,11n10 layer","fraction of non cohesive sediment of class12,12n10 layer","fraction of non cohesive sediment of class13,13n10 layer","fraction of non cohesive sediment of class14,14n10 layer","fraction of non cohesive sediment of class15,15n10 layer","fraction of non cohesive sediment of class16,16n10 layer","fraction of non cohesive sediment of class17,17n10 layer","fraction of non cohesive sediment of class18,18n10 layer","fraction of non cohesive sediment of class19,19n10 layer","fraction of non cohesive sediment of class**,**n10 layer","fraction of non cohesive sediment of class*,*n11 layer","fraction of non cohesive sediment of class1,1n11 layer","fraction of non cohesive sediment of class2,2n11 layer","fraction of non cohesive sediment of class3,3n11 layer","fraction of non cohesive sediment of class4,4n11 layer","fraction of non cohesive sediment of class5,5n11 layer","fraction of non cohesive sediment of class6,6n11 layer","fraction of non cohesive sediment of class7,7n11 layer","fraction of non cohesive sediment of class8,8n11 layer","fraction of non cohesive sediment of class9,9n11 layer","fraction of non cohesive sediment of class10,10n11 layer","fraction of non cohesive sediment of class11,11n11 layer","fraction of non cohesive sediment of class12,12n11 layer","fraction of non cohesive sediment of class13,13n11 layer","fraction of non cohesive sediment of class14,14n11 layer","fraction of non cohesive sediment of class15,15n11 layer","fraction of non cohesive sediment of class16,16n11 layer","fraction of non cohesive sediment of class17,17n11 layer","fraction of non cohesive sediment of class18,18n11 layer","fraction of non cohesive sediment of class19,19n11 layer","fraction of non cohesive sediment of class**,**n11 layer","fraction of non cohesive sediment of class*,*n12 layer","fraction of non cohesive sediment of class1,1n12 layer","fraction of non cohesive sediment of class2,2n12 layer","fraction of non cohesive sediment of class3,3n12 layer","fraction of non cohesive sediment of class4,4n12 layer","fraction of non cohesive sediment of class5,5n12 layer","fraction of non cohesive sediment of class6,6n12 layer","fraction of non cohesive sediment of class7,7n12 layer","fraction of non cohesive sediment of class8,8n12 layer","fraction of non cohesive sediment of class9,9n12 layer","fraction of non cohesive sediment of class10,10n12 layer","fraction of non cohesive sediment of class11,11n12 layer","fraction of non cohesive sediment of class12,12n12 layer","fraction of non cohesive sediment of class13,13n12 layer","fraction of non cohesive sediment of class14,14n12 layer","fraction of non cohesive sediment of class15,15n12 layer","fraction of non cohesive sediment of class16,16n12 layer","fraction of non cohesive sediment of class17,17n12 layer","fraction of non cohesive sediment of class18,18n12 layer","fraction of non cohesive sediment of class19,19n12 layer","fraction of non cohesive sediment of class**,**n12 layer","fraction of non cohesive sediment of class*,*n13 layer","fraction of non cohesive sediment of class1,1n13 layer","fraction of non cohesive sediment of class2,2n13 layer","fraction of non cohesive sediment of class3,3n13 layer","fraction of non cohesive sediment of class4,4n13 layer","fraction of non cohesive sediment of class5,5n13 layer","fraction of non cohesive sediment of class6,6n13 layer","fraction of non cohesive sediment of class7,7n13 layer","fraction of non cohesive sediment of class8,8n13 layer","fraction of non cohesive sediment of class9,9n13 layer","fraction of non cohesive sediment of class10,10n13 layer","fraction of non cohesive sediment of class11,11n13 layer","fraction of non cohesive sediment of class12,12n13 layer","fraction of non cohesive sediment of class13,13n13 layer","fraction of non cohesive sediment of class14,14n13 layer","fraction of non cohesive sediment of class15,15n13 layer","fraction of non cohesive sediment of class16,16n13 layer","fraction of non cohesive sediment of class17,17n13 layer","fraction of non cohesive sediment of class18,18n13 layer","fraction of non cohesive sediment of class19,19n13 layer","fraction of non cohesive sediment of class**,**n13 layer","fraction of non cohesive sediment of class*,*n14 layer","fraction of non cohesive sediment of class1,1n14 layer","fraction of non cohesive sediment of class2,2n14 layer","fraction of non cohesive sediment of class3,3n14 layer","fraction of non cohesive sediment of class4,4n14 layer","fraction of non cohesive sediment of class5,5n14 layer","fraction of non cohesive sediment of class6,6n14 layer","fraction of non cohesive sediment of class7,7n14 layer","fraction of non cohesive sediment of class8,8n14 layer","fraction of non cohesive sediment of class9,9n14 layer","fraction of non cohesive sediment of class10,10n14 layer","fraction of non cohesive sediment of class11,11n14 layer","fraction of non cohesive sediment of class12,12n14 layer","fraction of non cohesive sediment of class13,13n14 layer","fraction of non cohesive sediment of class14,14n14 layer","fraction of non cohesive sediment of class15,15n14 layer","fraction of non cohesive sediment of class16,16n14 layer","fraction of non cohesive sediment of class17,17n14 layer","fraction of non cohesive sediment of class18,18n14 layer","fraction of non cohesive sediment of class19,19n14 layer","fraction of non cohesive sediment of class**,**n14 layer","fraction of non cohesive sediment of class*,*n15 layer","fraction of non cohesive sediment of class1,1n15 layer","fraction of non cohesive sediment of class2,2n15 layer","fraction of non cohesive sediment of class3,3n15 layer","fraction of non cohesive sediment of class4,4n15 layer","fraction of non cohesive sediment of class5,5n15 layer","fraction of non cohesive sediment of class6,6n15 layer","fraction of non cohesive sediment of class7,7n15 layer","fraction of non cohesive sediment of class8,8n15 layer","fraction of non cohesive sediment of class9,9n15 layer","fraction of non cohesive sediment of class10,10n15 layer","fraction of non cohesive sediment of class11,11n15 layer","fraction of non cohesive sediment of class12,12n15 layer","fraction of non cohesive sediment of class13,13n15 layer","fraction of non cohesive sediment of class14,14n15 layer","fraction of non cohesive sediment of class15,15n15 layer","fraction of non cohesive sediment of class16,16n15 layer","fraction of non cohesive sediment of class17,17n15 layer","fraction of non cohesive sediment of class18,18n15 layer","fraction of non cohesive sediment of class19,19n15 layer","fraction of non cohesive sediment of class**,**n15 layer","fraction of non cohesive sediment of class*,*n16 layer","fraction of non cohesive sediment of class1,1n16 layer","fraction of non cohesive sediment of class2,2n16 layer","fraction of non cohesive sediment of class3,3n16 layer","fraction of non cohesive sediment of class4,4n16 layer","fraction of non cohesive sediment of class5,5n16 layer","fraction of non cohesive sediment of class6,6n16 layer","fraction of non cohesive sediment of class7,7n16 layer","fraction of non cohesive sediment of class8,8n16 layer","fraction of non cohesive sediment of class9,9n16 layer","fraction of non cohesive sediment of class10,10n16 layer","fraction of non cohesive sediment of class11,11n16 layer","fraction of non cohesive sediment of class12,12n16 layer","fraction of non cohesive sediment of class13,13n16 layer","fraction of non cohesive sediment of class14,14n16 layer","fraction of non cohesive sediment of class15,15n16 layer","fraction of non cohesive sediment of class16,16n16 layer","fraction of non cohesive sediment of class17,17n16 layer","fraction of non cohesive sediment of class18,18n16 layer","fraction of non cohesive sediment of class19,19n16 layer","fraction of non cohesive sediment of class**,**n16 layer","fraction of non cohesive sediment of class*,*n17 layer","fraction of non cohesive sediment of class1,1n17 layer","fraction of non cohesive sediment of class2,2n17 layer","fraction of non cohesive sediment of class3,3n17 layer","fraction of non cohesive sediment of class4,4n17 layer","fraction of non cohesive sediment of class5,5n17 layer","fraction of non cohesive sediment of class6,6n17 layer","fraction of non cohesive sediment of class7,7n17 layer","fraction of non cohesive sediment of class8,8n17 layer","fraction of non cohesive sediment of class9,9n17 layer","fraction of non cohesive sediment of class10,10n17 layer","fraction of non cohesive sediment of class11,11n17 layer","fraction of non cohesive sediment of class12,12n17 layer","fraction of non cohesive sediment of class13,13n17 layer","fraction of non cohesive sediment of class14,14n17 layer","fraction of non cohesive sediment of class15,15n17 layer","fraction of non cohesive sediment of class16,16n17 layer","fraction of non cohesive sediment of class17,17n17 layer","fraction of non cohesive sediment of class18,18n17 layer","fraction of non cohesive sediment of class19,19n17 layer","fraction of non cohesive sediment of class**,**n17 layer","fraction of non cohesive sediment of class*,*n18 layer","fraction of non cohesive sediment of class1,1n18 layer","fraction of non cohesive sediment of class2,2n18 layer","fraction of non cohesive sediment of class3,3n18 layer","fraction of non cohesive sediment of class4,4n18 layer","fraction of non cohesive sediment of class5,5n18 layer","fraction of non cohesive sediment of class6,6n18 layer","fraction of non cohesive sediment of class7,7n18 layer","fraction of non cohesive sediment of class8,8n18 layer","fraction of non cohesive sediment of class9,9n18 layer","fraction of non cohesive sediment of class10,10n18 layer","fraction of non cohesive sediment of class11,11n18 layer","fraction of non cohesive sediment of class12,12n18 layer","fraction of non cohesive sediment of class13,13n18 layer","fraction of non cohesive sediment of class14,14n18 layer","fraction of non cohesive sediment of class15,15n18 layer","fraction of non cohesive sediment of class16,16n18 layer","fraction of non cohesive sediment of class17,17n18 layer","fraction of non cohesive sediment of class18,18n18 layer","fraction of non cohesive sediment of class19,19n18 layer","fraction of non cohesive sediment of class**,**n18 layer","fraction of non cohesive sediment of class*,*n19 layer","fraction of non cohesive sediment of class1,1n19 layer","fraction of non cohesive sediment of class2,2n19 layer","fraction of non cohesive sediment of class3,3n19 layer","fraction of non cohesive sediment of class4,4n19 layer","fraction of non cohesive sediment of class5,5n19 layer","fraction of non cohesive sediment of class6,6n19 layer","fraction of non cohesive sediment of class7,7n19 layer","fraction of non cohesive sediment of class8,8n19 layer","fraction of non cohesive sediment of class9,9n19 layer","fraction of non cohesive sediment of class10,10n19 layer","fraction of non cohesive sediment of class11,11n19 layer","fraction of non cohesive sediment of class12,12n19 layer","fraction of non cohesive sediment of class13,13n19 layer","fraction of non cohesive sediment of class14,14n19 layer","fraction of non cohesive sediment of class15,15n19 layer","fraction of non cohesive sediment of class16,16n19 layer","fraction of non cohesive sediment of class17,17n19 layer","fraction of non cohesive sediment of class18,18n19 layer","fraction of non cohesive sediment of class19,19n19 layer","fraction of non cohesive sediment of class**,**n19 layer","fraction of non cohesive sediment of class*,*n** layer","fraction of non cohesive sediment of class1,1n** layer","fraction of non cohesive sediment of class2,2n** layer","fraction of non cohesive sediment of class3,3n** layer","fraction of non cohesive sediment of class4,4n** layer","fraction of non cohesive sediment of class5,5n** layer","fraction of non cohesive sediment of class6,6n** layer","fraction of non cohesive sediment of class7,7n** layer","fraction of non cohesive sediment of class8,8n** layer","fraction of non cohesive sediment of class9,9n** layer","fraction of non cohesive sediment of class10,10n** layer","fraction of non cohesive sediment of class11,11n** layer","fraction of non cohesive sediment of class12,12n** layer","fraction of non cohesive sediment of class13,13n** layer","fraction of non cohesive sediment of class14,14n** layer","fraction of non cohesive sediment of class15,15n** layer","fraction of non cohesive sediment of class16,16n** layer","fraction of non cohesive sediment of class17,17n** layer","fraction of non cohesive sediment of class18,18n** layer","fraction of non cohesive sediment of class19,19n** layer","fraction of non cohesive sediment of class**,**n** layer","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","thicknes of bed layer k","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration for class i","supplementary variable A","supplementary variable G","supplementary variable L","supplementary variable O"],
        defaut = '',
        fr = """Nom des variables que l''utilisateur desire ecrire sur
le lisring. Meme possibilites que pour les sorties graphiques.""",
        ang = """Names of variables the user wants to write on the listing.
Each variable is represented by a letter in the same manner as
it is done in the graphic results file.""",
    ),
#   -----------------------------------
    LISTING_PRINTOUT_PERIOD = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Determine la periode en nombre de pas de temps d''impression
des ''VARIABLES A IMPRIMER'' (voir ce mot-cle).
La sortie des resultats est effectuee sur le fichier listing
(fichier cas\_numerodeprocessus.sortie sur station de travail).""",
        ang = """Determines, in number of time steps, the printout period of
the ''VARIABLES TO BE PRINTED''.
The results are printed out on the listing file
(file cas\_numerodeprocessus.sortie on a workstation).""",
    ),
#   -----------------------------------
    MASS_BALANCE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Determine si oui ou non le bilan de masse est realise.""",
        ang = """Determines whether a check of the mass-balance over the domain
is made or not""",
    ),
#   -----------------------------------
    SECTIONS_OUTPUT_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)','Sauvegarde'),
        defaut = '',
        fr = """sections output file, written by the master""",
        ang = """sections output file, written by the master""",
    ),
#   -----------------------------------
    C_VSM_RESULTS_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
        defaut = '',
        fr = """Nom du fichier dans lequel seront ecrits les resultats C-VSM du
calcul avec la periodicite donnee par le mot cle
\telkey{C-VSM FULL PRINTOUT PERIOD}.""",
        ang = """Name of the file into which the C-VSM results of the computation
are written, the periodicity being given by the keyword:
\telkey{C-VSM FULL PRINTOUT PERIOD}.""",
    ),
#   -----------------------------------
    C_VSM_RESULTS_FILE_FORMAT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ['SERAFIN','SERAFIND','MED'],
        defaut = 'SERAFIN',
        fr = """Format du \telkey{FICHIER DES C-VSM}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
        ang = """Format of the \telkey{C-VSM RESULT FILE}. Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
    ),
)
# -----------------------------------------------------------------------
DATA_FILES = PROC(nom= "DATA_FILES",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    WAVE_FILE_FORMAT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ['SERAFIN','SERAFIND','MED'],
        defaut = 'SERAFIN',
        fr = """Format du fichier de houle.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
        ang = """Wave file format.
Possible values are:
- SERAFIN : classical single precision format in Telemac;
- SERAFIND: classical double precision format in Telemac;
- MED     : MED format based on HDF5""",
    ),
#   -----------------------------------
    FORTRAN_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = 'FichierOuRepertoire',
        defaut = '',
        fr = """Nom du fichier FORTRAN a soumettre.""",
        ang = """Name of FORTRAN file to be submitted.""",
    ),
#   -----------------------------------
    BOUNDARY_CONDITIONS_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        fr = """Nom du fichier contenant les types de conditions aux limites.
Ce fichier est rempli de facon automatique par le mailleur au moyen de
couleurs affectees aux noeuds des frontieres du domaine de calcul.""",
        ang = """Name of the file containing the types of boundary conditions.
This file is filled automatically by the mesh generator through
colours that are assigned to the computation domain boundary nodes.""",
    ),
#   -----------------------------------
    WAVE_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom d''un fichier contenant les resultats d''un calcul
precedent TOMAWAC realise sur le meme maillage. La houle sera donnee,
par le dernier pas de temps
de ce fichier . Seules les donnees de houle de ce fichier
seront conservees (hauteur de houle, periode de
houle, angle de la houle).
Les donnees hydrodynamiques(hauteur d''eau,vitesses) seront lues
dans le fichier hydrodynamique (Verifier la compatibilite !!!)
Les donnees de houle peuvent aussi etre imposees
par l''utilisateur dans le sous programme USER\_FORCING\_GAIA.
ou encore lues dans le fichier hydrodynamique.""",
        ang = """Name of a file containing the results a previous
TOMAWAC computation  made on the same mesh. The wave data (wave height,
 wave period, wave angle ) will
be given by the last record of the file.
The user has to verify that both informations (wave and current data)
are consistent.
Remark :The wave data can also be specified in
 the hydrodynamic file. the user has also
 the possibility to give the values
of the wave data in the subroutine USER\_FORCING\_GAIA.
This is recommended for non-steady flow simulation.""",
    ),
#   -----------------------------------
    BOTTOM_TOPOGRAPHY_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier facultatif contenant la bathymetrie
associee au maillage.""",
        ang = """Name of the possible file containing the bathymetric data.""",
    ),
#   -----------------------------------
    SECTIONS_INPUT_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """sections input file, partitioned""",
        ang = """sections input file, partitioned""",
    ),
)
# -----------------------------------------------------------------------
INITIAL_CONDITIONS = PROC(nom= "INITIAL_CONDITIONS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    PREVIOUS_SEDIMENTOLOGICAL_COMPUTATION_FILE_FORMAT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ['SERAFIN','SERAFIND','MED'],
        defaut = 'SERAFIN',
        fr = """Format du fichier de resultats du calcul precedent.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
        ang = """Previous computation results file format.
Possible values are:
- SERAFIN : classical single precision format in Telemac;
- SERAFIND: classical double precision format in Telemac;
- MED     : MED format based on HDF5""",
    ),
#   -----------------------------------
    PREVIOUS_SEDIMENTOLOGICAL_COMPUTATION_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom d''un fichier contenant les resultats d''un
calcul precedent sedimentologique realise sur le meme maillage
et dont le dernier pas de temps enregistre va fournir les
conditions initiales pour une suite de de calcul.""",
        ang = """Name of a file containing the results of an
earlier  sedimentological computation which was made
on the same mesh. The last recorded time step will provide
the initial conditions for the new computation.""",
    ),
#   -----------------------------------
    COMPUTATION_CONTINUED = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Determine si le calcul en cours est independant de tout autre
resultat ou est une reprise effectuee a partir du resultat d''un calcul
precedent.
NON : Il s''agit du premier passage pour ce calcul et il est necessaire
de definir un jeu complet de conditions initiales
OUI : Il s''agit d''une reprise de calcul :
les conditions initiales sont constituees par le dernier pas de temps du
FICHIER PRECEDENT SEDIMENTOLOGIQUE du fichier des parametres utilise
pour soumettre le calcul.
Par contre, l''ensemble des donnees du fichier des parametres
peuvent etre redefinies
De meme, il est necessaire de definir des conditions aux limites""",
        ang = """Determines whether the computation under way is an independent
result or is following an earlier result.
NO: It is the first run for this computation and a whole set of
initial conditions should be defined.
YES: It follows a former computation:
the initial conditions consist in the last time step of the
PREVIOUS COMPUTATION FILE
in the steering file used for submitting the computation.
All the data from the steering file may be defined once again, which
provides an opportunity to change, for example, the time step.
It is also possible to define new boundary conditions.""",
    ),
)
# -----------------------------------------------------------------------
INPUT_OUTPUT__GRAPHICS_AND_LISTING = PROC(nom= "INPUT_OUTPUT__GRAPHICS_AND_LISTING",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    VARIABLES_FOR_GRAPHIC_PRINTOUTS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM', min=0, max='**',
        into = ["velocity along x axis (m/s)","velocity along y axis (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","scalar flowrate of fluid (m2/s)","flowrate along x axis (m2/s)","flowrate along y axis (m2/s)","non erodable bottom","Bed Shear stress (Totalfriction) (N/m2)","wave height","wave period","wave angle with axis Oy (deg)","bed-load discharge (kg/(m*s))","bed-load discharge along x axis (kg/(m*s))","bed-load discharge along y axis (kg/(m*s))","bottom evolution (m)","total bed roughness (m)","Skin friction correction factor","Mean grain diameter","wave orbital velocity (m/s)","fraction of non cohesive sediment of class*,*n* layer","fraction of non cohesive sediment of class1,1n* layer","fraction of non cohesive sediment of class2,2n* layer","fraction of non cohesive sediment of class3,3n* layer","fraction of non cohesive sediment of class4,4n* layer","fraction of non cohesive sediment of class5,5n* layer","fraction of non cohesive sediment of class6,6n* layer","fraction of non cohesive sediment of class7,7n* layer","fraction of non cohesive sediment of class8,8n* layer","fraction of non cohesive sediment of class9,9n* layer","fraction of non cohesive sediment of class10,10n* layer","fraction of non cohesive sediment of class11,11n* layer","fraction of non cohesive sediment of class12,12n* layer","fraction of non cohesive sediment of class13,13n* layer","fraction of non cohesive sediment of class14,14n* layer","fraction of non cohesive sediment of class15,15n* layer","fraction of non cohesive sediment of class16,16n* layer","fraction of non cohesive sediment of class17,17n* layer","fraction of non cohesive sediment of class18,18n* layer","fraction of non cohesive sediment of class19,19n* layer","fraction of non cohesive sediment of class**,**n* layer","fraction of non cohesive sediment of class*,*n1 layer","fraction of non cohesive sediment of class1,1n1 layer","fraction of non cohesive sediment of class2,2n1 layer","fraction of non cohesive sediment of class3,3n1 layer","fraction of non cohesive sediment of class4,4n1 layer","fraction of non cohesive sediment of class5,5n1 layer","fraction of non cohesive sediment of class6,6n1 layer","fraction of non cohesive sediment of class7,7n1 layer","fraction of non cohesive sediment of class8,8n1 layer","fraction of non cohesive sediment of class9,9n1 layer","fraction of non cohesive sediment of class10,10n1 layer","fraction of non cohesive sediment of class11,11n1 layer","fraction of non cohesive sediment of class12,12n1 layer","fraction of non cohesive sediment of class13,13n1 layer","fraction of non cohesive sediment of class14,14n1 layer","fraction of non cohesive sediment of class15,15n1 layer","fraction of non cohesive sediment of class16,16n1 layer","fraction of non cohesive sediment of class17,17n1 layer","fraction of non cohesive sediment of class18,18n1 layer","fraction of non cohesive sediment of class19,19n1 layer","fraction of non cohesive sediment of class**,**n1 layer","fraction of non cohesive sediment of class*,*n2 layer","fraction of non cohesive sediment of class1,1n2 layer","fraction of non cohesive sediment of class2,2n2 layer","fraction of non cohesive sediment of class3,3n2 layer","fraction of non cohesive sediment of class4,4n2 layer","fraction of non cohesive sediment of class5,5n2 layer","fraction of non cohesive sediment of class6,6n2 layer","fraction of non cohesive sediment of class7,7n2 layer","fraction of non cohesive sediment of class8,8n2 layer","fraction of non cohesive sediment of class9,9n2 layer","fraction of non cohesive sediment of class10,10n2 layer","fraction of non cohesive sediment of class11,11n2 layer","fraction of non cohesive sediment of class12,12n2 layer","fraction of non cohesive sediment of class13,13n2 layer","fraction of non cohesive sediment of class14,14n2 layer","fraction of non cohesive sediment of class15,15n2 layer","fraction of non cohesive sediment of class16,16n2 layer","fraction of non cohesive sediment of class17,17n2 layer","fraction of non cohesive sediment of class18,18n2 layer","fraction of non cohesive sediment of class19,19n2 layer","fraction of non cohesive sediment of class**,**n2 layer","fraction of non cohesive sediment of class*,*n3 layer","fraction of non cohesive sediment of class1,1n3 layer","fraction of non cohesive sediment of class2,2n3 layer","fraction of non cohesive sediment of class3,3n3 layer","fraction of non cohesive sediment of class4,4n3 layer","fraction of non cohesive sediment of class5,5n3 layer","fraction of non cohesive sediment of class6,6n3 layer","fraction of non cohesive sediment of class7,7n3 layer","fraction of non cohesive sediment of class8,8n3 layer","fraction of non cohesive sediment of class9,9n3 layer","fraction of non cohesive sediment of class10,10n3 layer","fraction of non cohesive sediment of class11,11n3 layer","fraction of non cohesive sediment of class12,12n3 layer","fraction of non cohesive sediment of class13,13n3 layer","fraction of non cohesive sediment of class14,14n3 layer","fraction of non cohesive sediment of class15,15n3 layer","fraction of non cohesive sediment of class16,16n3 layer","fraction of non cohesive sediment of class17,17n3 layer","fraction of non cohesive sediment of class18,18n3 layer","fraction of non cohesive sediment of class19,19n3 layer","fraction of non cohesive sediment of class**,**n3 layer","fraction of non cohesive sediment of class*,*n4 layer","fraction of non cohesive sediment of class1,1n4 layer","fraction of non cohesive sediment of class2,2n4 layer","fraction of non cohesive sediment of class3,3n4 layer","fraction of non cohesive sediment of class4,4n4 layer","fraction of non cohesive sediment of class5,5n4 layer","fraction of non cohesive sediment of class6,6n4 layer","fraction of non cohesive sediment of class7,7n4 layer","fraction of non cohesive sediment of class8,8n4 layer","fraction of non cohesive sediment of class9,9n4 layer","fraction of non cohesive sediment of class10,10n4 layer","fraction of non cohesive sediment of class11,11n4 layer","fraction of non cohesive sediment of class12,12n4 layer","fraction of non cohesive sediment of class13,13n4 layer","fraction of non cohesive sediment of class14,14n4 layer","fraction of non cohesive sediment of class15,15n4 layer","fraction of non cohesive sediment of class16,16n4 layer","fraction of non cohesive sediment of class17,17n4 layer","fraction of non cohesive sediment of class18,18n4 layer","fraction of non cohesive sediment of class19,19n4 layer","fraction of non cohesive sediment of class**,**n4 layer","fraction of non cohesive sediment of class*,*n5 layer","fraction of non cohesive sediment of class1,1n5 layer","fraction of non cohesive sediment of class2,2n5 layer","fraction of non cohesive sediment of class3,3n5 layer","fraction of non cohesive sediment of class4,4n5 layer","fraction of non cohesive sediment of class5,5n5 layer","fraction of non cohesive sediment of class6,6n5 layer","fraction of non cohesive sediment of class7,7n5 layer","fraction of non cohesive sediment of class8,8n5 layer","fraction of non cohesive sediment of class9,9n5 layer","fraction of non cohesive sediment of class10,10n5 layer","fraction of non cohesive sediment of class11,11n5 layer","fraction of non cohesive sediment of class12,12n5 layer","fraction of non cohesive sediment of class13,13n5 layer","fraction of non cohesive sediment of class14,14n5 layer","fraction of non cohesive sediment of class15,15n5 layer","fraction of non cohesive sediment of class16,16n5 layer","fraction of non cohesive sediment of class17,17n5 layer","fraction of non cohesive sediment of class18,18n5 layer","fraction of non cohesive sediment of class19,19n5 layer","fraction of non cohesive sediment of class**,**n5 layer","fraction of non cohesive sediment of class*,*n6 layer","fraction of non cohesive sediment of class1,1n6 layer","fraction of non cohesive sediment of class2,2n6 layer","fraction of non cohesive sediment of class3,3n6 layer","fraction of non cohesive sediment of class4,4n6 layer","fraction of non cohesive sediment of class5,5n6 layer","fraction of non cohesive sediment of class6,6n6 layer","fraction of non cohesive sediment of class7,7n6 layer","fraction of non cohesive sediment of class8,8n6 layer","fraction of non cohesive sediment of class9,9n6 layer","fraction of non cohesive sediment of class10,10n6 layer","fraction of non cohesive sediment of class11,11n6 layer","fraction of non cohesive sediment of class12,12n6 layer","fraction of non cohesive sediment of class13,13n6 layer","fraction of non cohesive sediment of class14,14n6 layer","fraction of non cohesive sediment of class15,15n6 layer","fraction of non cohesive sediment of class16,16n6 layer","fraction of non cohesive sediment of class17,17n6 layer","fraction of non cohesive sediment of class18,18n6 layer","fraction of non cohesive sediment of class19,19n6 layer","fraction of non cohesive sediment of class**,**n6 layer","fraction of non cohesive sediment of class*,*n7 layer","fraction of non cohesive sediment of class1,1n7 layer","fraction of non cohesive sediment of class2,2n7 layer","fraction of non cohesive sediment of class3,3n7 layer","fraction of non cohesive sediment of class4,4n7 layer","fraction of non cohesive sediment of class5,5n7 layer","fraction of non cohesive sediment of class6,6n7 layer","fraction of non cohesive sediment of class7,7n7 layer","fraction of non cohesive sediment of class8,8n7 layer","fraction of non cohesive sediment of class9,9n7 layer","fraction of non cohesive sediment of class10,10n7 layer","fraction of non cohesive sediment of class11,11n7 layer","fraction of non cohesive sediment of class12,12n7 layer","fraction of non cohesive sediment of class13,13n7 layer","fraction of non cohesive sediment of class14,14n7 layer","fraction of non cohesive sediment of class15,15n7 layer","fraction of non cohesive sediment of class16,16n7 layer","fraction of non cohesive sediment of class17,17n7 layer","fraction of non cohesive sediment of class18,18n7 layer","fraction of non cohesive sediment of class19,19n7 layer","fraction of non cohesive sediment of class**,**n7 layer","fraction of non cohesive sediment of class*,*n8 layer","fraction of non cohesive sediment of class1,1n8 layer","fraction of non cohesive sediment of class2,2n8 layer","fraction of non cohesive sediment of class3,3n8 layer","fraction of non cohesive sediment of class4,4n8 layer","fraction of non cohesive sediment of class5,5n8 layer","fraction of non cohesive sediment of class6,6n8 layer","fraction of non cohesive sediment of class7,7n8 layer","fraction of non cohesive sediment of class8,8n8 layer","fraction of non cohesive sediment of class9,9n8 layer","fraction of non cohesive sediment of class10,10n8 layer","fraction of non cohesive sediment of class11,11n8 layer","fraction of non cohesive sediment of class12,12n8 layer","fraction of non cohesive sediment of class13,13n8 layer","fraction of non cohesive sediment of class14,14n8 layer","fraction of non cohesive sediment of class15,15n8 layer","fraction of non cohesive sediment of class16,16n8 layer","fraction of non cohesive sediment of class17,17n8 layer","fraction of non cohesive sediment of class18,18n8 layer","fraction of non cohesive sediment of class19,19n8 layer","fraction of non cohesive sediment of class**,**n8 layer","fraction of non cohesive sediment of class*,*n9 layer","fraction of non cohesive sediment of class1,1n9 layer","fraction of non cohesive sediment of class2,2n9 layer","fraction of non cohesive sediment of class3,3n9 layer","fraction of non cohesive sediment of class4,4n9 layer","fraction of non cohesive sediment of class5,5n9 layer","fraction of non cohesive sediment of class6,6n9 layer","fraction of non cohesive sediment of class7,7n9 layer","fraction of non cohesive sediment of class8,8n9 layer","fraction of non cohesive sediment of class9,9n9 layer","fraction of non cohesive sediment of class10,10n9 layer","fraction of non cohesive sediment of class11,11n9 layer","fraction of non cohesive sediment of class12,12n9 layer","fraction of non cohesive sediment of class13,13n9 layer","fraction of non cohesive sediment of class14,14n9 layer","fraction of non cohesive sediment of class15,15n9 layer","fraction of non cohesive sediment of class16,16n9 layer","fraction of non cohesive sediment of class17,17n9 layer","fraction of non cohesive sediment of class18,18n9 layer","fraction of non cohesive sediment of class19,19n9 layer","fraction of non cohesive sediment of class**,**n9 layer","fraction of non cohesive sediment of class*,*n10 layer","fraction of non cohesive sediment of class1,1n10 layer","fraction of non cohesive sediment of class2,2n10 layer","fraction of non cohesive sediment of class3,3n10 layer","fraction of non cohesive sediment of class4,4n10 layer","fraction of non cohesive sediment of class5,5n10 layer","fraction of non cohesive sediment of class6,6n10 layer","fraction of non cohesive sediment of class7,7n10 layer","fraction of non cohesive sediment of class8,8n10 layer","fraction of non cohesive sediment of class9,9n10 layer","fraction of non cohesive sediment of class10,10n10 layer","fraction of non cohesive sediment of class11,11n10 layer","fraction of non cohesive sediment of class12,12n10 layer","fraction of non cohesive sediment of class13,13n10 layer","fraction of non cohesive sediment of class14,14n10 layer","fraction of non cohesive sediment of class15,15n10 layer","fraction of non cohesive sediment of class16,16n10 layer","fraction of non cohesive sediment of class17,17n10 layer","fraction of non cohesive sediment of class18,18n10 layer","fraction of non cohesive sediment of class19,19n10 layer","fraction of non cohesive sediment of class**,**n10 layer","fraction of non cohesive sediment of class*,*n11 layer","fraction of non cohesive sediment of class1,1n11 layer","fraction of non cohesive sediment of class2,2n11 layer","fraction of non cohesive sediment of class3,3n11 layer","fraction of non cohesive sediment of class4,4n11 layer","fraction of non cohesive sediment of class5,5n11 layer","fraction of non cohesive sediment of class6,6n11 layer","fraction of non cohesive sediment of class7,7n11 layer","fraction of non cohesive sediment of class8,8n11 layer","fraction of non cohesive sediment of class9,9n11 layer","fraction of non cohesive sediment of class10,10n11 layer","fraction of non cohesive sediment of class11,11n11 layer","fraction of non cohesive sediment of class12,12n11 layer","fraction of non cohesive sediment of class13,13n11 layer","fraction of non cohesive sediment of class14,14n11 layer","fraction of non cohesive sediment of class15,15n11 layer","fraction of non cohesive sediment of class16,16n11 layer","fraction of non cohesive sediment of class17,17n11 layer","fraction of non cohesive sediment of class18,18n11 layer","fraction of non cohesive sediment of class19,19n11 layer","fraction of non cohesive sediment of class**,**n11 layer","fraction of non cohesive sediment of class*,*n12 layer","fraction of non cohesive sediment of class1,1n12 layer","fraction of non cohesive sediment of class2,2n12 layer","fraction of non cohesive sediment of class3,3n12 layer","fraction of non cohesive sediment of class4,4n12 layer","fraction of non cohesive sediment of class5,5n12 layer","fraction of non cohesive sediment of class6,6n12 layer","fraction of non cohesive sediment of class7,7n12 layer","fraction of non cohesive sediment of class8,8n12 layer","fraction of non cohesive sediment of class9,9n12 layer","fraction of non cohesive sediment of class10,10n12 layer","fraction of non cohesive sediment of class11,11n12 layer","fraction of non cohesive sediment of class12,12n12 layer","fraction of non cohesive sediment of class13,13n12 layer","fraction of non cohesive sediment of class14,14n12 layer","fraction of non cohesive sediment of class15,15n12 layer","fraction of non cohesive sediment of class16,16n12 layer","fraction of non cohesive sediment of class17,17n12 layer","fraction of non cohesive sediment of class18,18n12 layer","fraction of non cohesive sediment of class19,19n12 layer","fraction of non cohesive sediment of class**,**n12 layer","fraction of non cohesive sediment of class*,*n13 layer","fraction of non cohesive sediment of class1,1n13 layer","fraction of non cohesive sediment of class2,2n13 layer","fraction of non cohesive sediment of class3,3n13 layer","fraction of non cohesive sediment of class4,4n13 layer","fraction of non cohesive sediment of class5,5n13 layer","fraction of non cohesive sediment of class6,6n13 layer","fraction of non cohesive sediment of class7,7n13 layer","fraction of non cohesive sediment of class8,8n13 layer","fraction of non cohesive sediment of class9,9n13 layer","fraction of non cohesive sediment of class10,10n13 layer","fraction of non cohesive sediment of class11,11n13 layer","fraction of non cohesive sediment of class12,12n13 layer","fraction of non cohesive sediment of class13,13n13 layer","fraction of non cohesive sediment of class14,14n13 layer","fraction of non cohesive sediment of class15,15n13 layer","fraction of non cohesive sediment of class16,16n13 layer","fraction of non cohesive sediment of class17,17n13 layer","fraction of non cohesive sediment of class18,18n13 layer","fraction of non cohesive sediment of class19,19n13 layer","fraction of non cohesive sediment of class**,**n13 layer","fraction of non cohesive sediment of class*,*n14 layer","fraction of non cohesive sediment of class1,1n14 layer","fraction of non cohesive sediment of class2,2n14 layer","fraction of non cohesive sediment of class3,3n14 layer","fraction of non cohesive sediment of class4,4n14 layer","fraction of non cohesive sediment of class5,5n14 layer","fraction of non cohesive sediment of class6,6n14 layer","fraction of non cohesive sediment of class7,7n14 layer","fraction of non cohesive sediment of class8,8n14 layer","fraction of non cohesive sediment of class9,9n14 layer","fraction of non cohesive sediment of class10,10n14 layer","fraction of non cohesive sediment of class11,11n14 layer","fraction of non cohesive sediment of class12,12n14 layer","fraction of non cohesive sediment of class13,13n14 layer","fraction of non cohesive sediment of class14,14n14 layer","fraction of non cohesive sediment of class15,15n14 layer","fraction of non cohesive sediment of class16,16n14 layer","fraction of non cohesive sediment of class17,17n14 layer","fraction of non cohesive sediment of class18,18n14 layer","fraction of non cohesive sediment of class19,19n14 layer","fraction of non cohesive sediment of class**,**n14 layer","fraction of non cohesive sediment of class*,*n15 layer","fraction of non cohesive sediment of class1,1n15 layer","fraction of non cohesive sediment of class2,2n15 layer","fraction of non cohesive sediment of class3,3n15 layer","fraction of non cohesive sediment of class4,4n15 layer","fraction of non cohesive sediment of class5,5n15 layer","fraction of non cohesive sediment of class6,6n15 layer","fraction of non cohesive sediment of class7,7n15 layer","fraction of non cohesive sediment of class8,8n15 layer","fraction of non cohesive sediment of class9,9n15 layer","fraction of non cohesive sediment of class10,10n15 layer","fraction of non cohesive sediment of class11,11n15 layer","fraction of non cohesive sediment of class12,12n15 layer","fraction of non cohesive sediment of class13,13n15 layer","fraction of non cohesive sediment of class14,14n15 layer","fraction of non cohesive sediment of class15,15n15 layer","fraction of non cohesive sediment of class16,16n15 layer","fraction of non cohesive sediment of class17,17n15 layer","fraction of non cohesive sediment of class18,18n15 layer","fraction of non cohesive sediment of class19,19n15 layer","fraction of non cohesive sediment of class**,**n15 layer","fraction of non cohesive sediment of class*,*n16 layer","fraction of non cohesive sediment of class1,1n16 layer","fraction of non cohesive sediment of class2,2n16 layer","fraction of non cohesive sediment of class3,3n16 layer","fraction of non cohesive sediment of class4,4n16 layer","fraction of non cohesive sediment of class5,5n16 layer","fraction of non cohesive sediment of class6,6n16 layer","fraction of non cohesive sediment of class7,7n16 layer","fraction of non cohesive sediment of class8,8n16 layer","fraction of non cohesive sediment of class9,9n16 layer","fraction of non cohesive sediment of class10,10n16 layer","fraction of non cohesive sediment of class11,11n16 layer","fraction of non cohesive sediment of class12,12n16 layer","fraction of non cohesive sediment of class13,13n16 layer","fraction of non cohesive sediment of class14,14n16 layer","fraction of non cohesive sediment of class15,15n16 layer","fraction of non cohesive sediment of class16,16n16 layer","fraction of non cohesive sediment of class17,17n16 layer","fraction of non cohesive sediment of class18,18n16 layer","fraction of non cohesive sediment of class19,19n16 layer","fraction of non cohesive sediment of class**,**n16 layer","fraction of non cohesive sediment of class*,*n17 layer","fraction of non cohesive sediment of class1,1n17 layer","fraction of non cohesive sediment of class2,2n17 layer","fraction of non cohesive sediment of class3,3n17 layer","fraction of non cohesive sediment of class4,4n17 layer","fraction of non cohesive sediment of class5,5n17 layer","fraction of non cohesive sediment of class6,6n17 layer","fraction of non cohesive sediment of class7,7n17 layer","fraction of non cohesive sediment of class8,8n17 layer","fraction of non cohesive sediment of class9,9n17 layer","fraction of non cohesive sediment of class10,10n17 layer","fraction of non cohesive sediment of class11,11n17 layer","fraction of non cohesive sediment of class12,12n17 layer","fraction of non cohesive sediment of class13,13n17 layer","fraction of non cohesive sediment of class14,14n17 layer","fraction of non cohesive sediment of class15,15n17 layer","fraction of non cohesive sediment of class16,16n17 layer","fraction of non cohesive sediment of class17,17n17 layer","fraction of non cohesive sediment of class18,18n17 layer","fraction of non cohesive sediment of class19,19n17 layer","fraction of non cohesive sediment of class**,**n17 layer","fraction of non cohesive sediment of class*,*n18 layer","fraction of non cohesive sediment of class1,1n18 layer","fraction of non cohesive sediment of class2,2n18 layer","fraction of non cohesive sediment of class3,3n18 layer","fraction of non cohesive sediment of class4,4n18 layer","fraction of non cohesive sediment of class5,5n18 layer","fraction of non cohesive sediment of class6,6n18 layer","fraction of non cohesive sediment of class7,7n18 layer","fraction of non cohesive sediment of class8,8n18 layer","fraction of non cohesive sediment of class9,9n18 layer","fraction of non cohesive sediment of class10,10n18 layer","fraction of non cohesive sediment of class11,11n18 layer","fraction of non cohesive sediment of class12,12n18 layer","fraction of non cohesive sediment of class13,13n18 layer","fraction of non cohesive sediment of class14,14n18 layer","fraction of non cohesive sediment of class15,15n18 layer","fraction of non cohesive sediment of class16,16n18 layer","fraction of non cohesive sediment of class17,17n18 layer","fraction of non cohesive sediment of class18,18n18 layer","fraction of non cohesive sediment of class19,19n18 layer","fraction of non cohesive sediment of class**,**n18 layer","fraction of non cohesive sediment of class*,*n19 layer","fraction of non cohesive sediment of class1,1n19 layer","fraction of non cohesive sediment of class2,2n19 layer","fraction of non cohesive sediment of class3,3n19 layer","fraction of non cohesive sediment of class4,4n19 layer","fraction of non cohesive sediment of class5,5n19 layer","fraction of non cohesive sediment of class6,6n19 layer","fraction of non cohesive sediment of class7,7n19 layer","fraction of non cohesive sediment of class8,8n19 layer","fraction of non cohesive sediment of class9,9n19 layer","fraction of non cohesive sediment of class10,10n19 layer","fraction of non cohesive sediment of class11,11n19 layer","fraction of non cohesive sediment of class12,12n19 layer","fraction of non cohesive sediment of class13,13n19 layer","fraction of non cohesive sediment of class14,14n19 layer","fraction of non cohesive sediment of class15,15n19 layer","fraction of non cohesive sediment of class16,16n19 layer","fraction of non cohesive sediment of class17,17n19 layer","fraction of non cohesive sediment of class18,18n19 layer","fraction of non cohesive sediment of class19,19n19 layer","fraction of non cohesive sediment of class**,**n19 layer","fraction of non cohesive sediment of class*,*n** layer","fraction of non cohesive sediment of class1,1n** layer","fraction of non cohesive sediment of class2,2n** layer","fraction of non cohesive sediment of class3,3n** layer","fraction of non cohesive sediment of class4,4n** layer","fraction of non cohesive sediment of class5,5n** layer","fraction of non cohesive sediment of class6,6n** layer","fraction of non cohesive sediment of class7,7n** layer","fraction of non cohesive sediment of class8,8n** layer","fraction of non cohesive sediment of class9,9n** layer","fraction of non cohesive sediment of class10,10n** layer","fraction of non cohesive sediment of class11,11n** layer","fraction of non cohesive sediment of class12,12n** layer","fraction of non cohesive sediment of class13,13n** layer","fraction of non cohesive sediment of class14,14n** layer","fraction of non cohesive sediment of class15,15n** layer","fraction of non cohesive sediment of class16,16n** layer","fraction of non cohesive sediment of class17,17n** layer","fraction of non cohesive sediment of class18,18n** layer","fraction of non cohesive sediment of class19,19n** layer","fraction of non cohesive sediment of class**,**n** layer","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","solid transport load of class*","solid transport load of class1","solid transport load of class2","solid transport load of class3","solid transport load of class4","solid transport load of class5","solid transport load of class6","solid transport load of class7","solid transport load of class8","solid transport load of class9","solid transport load of class10","solid transport load of class11","solid transport load of class12","solid transport load of class13","solid transport load of class14","solid transport load of class15","solid transport load of class16","solid transport load of class17","solid transport load of class18","solid transport load of class19","solid transport load of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class*","mass concentration of class1","mass concentration of class2","mass concentration of class3","mass concentration of class4","mass concentration of class5","mass concentration of class6","mass concentration of class7","mass concentration of class8","mass concentration of class9","mass concentration of class10","mass concentration of class11","mass concentration of class12","mass concentration of class13","mass concentration of class14","mass concentration of class15","mass concentration of class16","mass concentration of class17","mass concentration of class18","mass concentration of class19","mass concentration of class**","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","mass concentration of class* for 2D graphic printouts","mass concentration of class1 for 2D graphic printouts","mass concentration of class2 for 2D graphic printouts","mass concentration of class3 for 2D graphic printouts","mass concentration of class4 for 2D graphic printouts","mass concentration of class5 for 2D graphic printouts","mass concentration of class6 for 2D graphic printouts","mass concentration of class7 for 2D graphic printouts","mass concentration of class8 for 2D graphic printouts","mass concentration of class9 for 2D graphic printouts","mass concentration of class10 for 2D graphic printouts","mass concentration of class11 for 2D graphic printouts","mass concentration of class12 for 2D graphic printouts","mass concentration of class13 for 2D graphic printouts","mass concentration of class14 for 2D graphic printouts","mass concentration of class15 for 2D graphic printouts","mass concentration of class16 for 2D graphic printouts","mass concentration of class17 for 2D graphic printouts","mass concentration of class18 for 2D graphic printouts","mass concentration of class19 for 2D graphic printouts","mass concentration of class** for 2D graphic printouts","sediment viscosity along x axis (m2/s) - only 3D","sediment viscosity along y axis (m2/s) - only 3D","sediment viscosity along z axis (m2/s) - only 3D","bed load transport rate (kg/(m*s))","bed load transport rate x axis","bed load transport rate y axis","bedload transport rate of class i","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","supplementary variable A","supplementary variable G","supplementary variable L","supplementary variable O","fraction of cohesive sediment of class*,*n* layer","fraction of cohesive sediment of class1,1n* layer","fraction of cohesive sediment of class2,2n* layer","fraction of cohesive sediment of class3,3n* layer","fraction of cohesive sediment of class4,4n* layer","fraction of cohesive sediment of class5,5n* layer","fraction of cohesive sediment of class6,6n* layer","fraction of cohesive sediment of class7,7n* layer","fraction of cohesive sediment of class8,8n* layer","fraction of cohesive sediment of class9,9n* layer","fraction of cohesive sediment of class10,10n* layer","fraction of cohesive sediment of class11,11n* layer","fraction of cohesive sediment of class12,12n* layer","fraction of cohesive sediment of class13,13n* layer","fraction of cohesive sediment of class14,14n* layer","fraction of cohesive sediment of class15,15n* layer","fraction of cohesive sediment of class16,16n* layer","fraction of cohesive sediment of class17,17n* layer","fraction of cohesive sediment of class18,18n* layer","fraction of cohesive sediment of class19,19n* layer","fraction of cohesive sediment of class**,**n* layer","fraction of cohesive sediment of class*,*n1 layer","fraction of cohesive sediment of class1,1n1 layer","fraction of cohesive sediment of class2,2n1 layer","fraction of cohesive sediment of class3,3n1 layer","fraction of cohesive sediment of class4,4n1 layer","fraction of cohesive sediment of class5,5n1 layer","fraction of cohesive sediment of class6,6n1 layer","fraction of cohesive sediment of class7,7n1 layer","fraction of cohesive sediment of class8,8n1 layer","fraction of cohesive sediment of class9,9n1 layer","fraction of cohesive sediment of class10,10n1 layer","fraction of cohesive sediment of class11,11n1 layer","fraction of cohesive sediment of class12,12n1 layer","fraction of cohesive sediment of class13,13n1 layer","fraction of cohesive sediment of class14,14n1 layer","fraction of cohesive sediment of class15,15n1 layer","fraction of cohesive sediment of class16,16n1 layer","fraction of cohesive sediment of class17,17n1 layer","fraction of cohesive sediment of class18,18n1 layer","fraction of cohesive sediment of class19,19n1 layer","fraction of cohesive sediment of class**,**n1 layer","fraction of cohesive sediment of class*,*n2 layer","fraction of cohesive sediment of class1,1n2 layer","fraction of cohesive sediment of class2,2n2 layer","fraction of cohesive sediment of class3,3n2 layer","fraction of cohesive sediment of class4,4n2 layer","fraction of cohesive sediment of class5,5n2 layer","fraction of cohesive sediment of class6,6n2 layer","fraction of cohesive sediment of class7,7n2 layer","fraction of cohesive sediment of class8,8n2 layer","fraction of cohesive sediment of class9,9n2 layer","fraction of cohesive sediment of class10,10n2 layer","fraction of cohesive sediment of class11,11n2 layer","fraction of cohesive sediment of class12,12n2 layer","fraction of cohesive sediment of class13,13n2 layer","fraction of cohesive sediment of class14,14n2 layer","fraction of cohesive sediment of class15,15n2 layer","fraction of cohesive sediment of class16,16n2 layer","fraction of cohesive sediment of class17,17n2 layer","fraction of cohesive sediment of class18,18n2 layer","fraction of cohesive sediment of class19,19n2 layer","fraction of cohesive sediment of class**,**n2 layer","fraction of cohesive sediment of class*,*n3 layer","fraction of cohesive sediment of class1,1n3 layer","fraction of cohesive sediment of class2,2n3 layer","fraction of cohesive sediment of class3,3n3 layer","fraction of cohesive sediment of class4,4n3 layer","fraction of cohesive sediment of class5,5n3 layer","fraction of cohesive sediment of class6,6n3 layer","fraction of cohesive sediment of class7,7n3 layer","fraction of cohesive sediment of class8,8n3 layer","fraction of cohesive sediment of class9,9n3 layer","fraction of cohesive sediment of class10,10n3 layer","fraction of cohesive sediment of class11,11n3 layer","fraction of cohesive sediment of class12,12n3 layer","fraction of cohesive sediment of class13,13n3 layer","fraction of cohesive sediment of class14,14n3 layer","fraction of cohesive sediment of class15,15n3 layer","fraction of cohesive sediment of class16,16n3 layer","fraction of cohesive sediment of class17,17n3 layer","fraction of cohesive sediment of class18,18n3 layer","fraction of cohesive sediment of class19,19n3 layer","fraction of cohesive sediment of class**,**n3 layer","fraction of cohesive sediment of class*,*n4 layer","fraction of cohesive sediment of class1,1n4 layer","fraction of cohesive sediment of class2,2n4 layer","fraction of cohesive sediment of class3,3n4 layer","fraction of cohesive sediment of class4,4n4 layer","fraction of cohesive sediment of class5,5n4 layer","fraction of cohesive sediment of class6,6n4 layer","fraction of cohesive sediment of class7,7n4 layer","fraction of cohesive sediment of class8,8n4 layer","fraction of cohesive sediment of class9,9n4 layer","fraction of cohesive sediment of class10,10n4 layer","fraction of cohesive sediment of class11,11n4 layer","fraction of cohesive sediment of class12,12n4 layer","fraction of cohesive sediment of class13,13n4 layer","fraction of cohesive sediment of class14,14n4 layer","fraction of cohesive sediment of class15,15n4 layer","fraction of cohesive sediment of class16,16n4 layer","fraction of cohesive sediment of class17,17n4 layer","fraction of cohesive sediment of class18,18n4 layer","fraction of cohesive sediment of class19,19n4 layer","fraction of cohesive sediment of class**,**n4 layer","fraction of cohesive sediment of class*,*n5 layer","fraction of cohesive sediment of class1,1n5 layer","fraction of cohesive sediment of class2,2n5 layer","fraction of cohesive sediment of class3,3n5 layer","fraction of cohesive sediment of class4,4n5 layer","fraction of cohesive sediment of class5,5n5 layer","fraction of cohesive sediment of class6,6n5 layer","fraction of cohesive sediment of class7,7n5 layer","fraction of cohesive sediment of class8,8n5 layer","fraction of cohesive sediment of class9,9n5 layer","fraction of cohesive sediment of class10,10n5 layer","fraction of cohesive sediment of class11,11n5 layer","fraction of cohesive sediment of class12,12n5 layer","fraction of cohesive sediment of class13,13n5 layer","fraction of cohesive sediment of class14,14n5 layer","fraction of cohesive sediment of class15,15n5 layer","fraction of cohesive sediment of class16,16n5 layer","fraction of cohesive sediment of class17,17n5 layer","fraction of cohesive sediment of class18,18n5 layer","fraction of cohesive sediment of class19,19n5 layer","fraction of cohesive sediment of class**,**n5 layer","fraction of cohesive sediment of class*,*n6 layer","fraction of cohesive sediment of class1,1n6 layer","fraction of cohesive sediment of class2,2n6 layer","fraction of cohesive sediment of class3,3n6 layer","fraction of cohesive sediment of class4,4n6 layer","fraction of cohesive sediment of class5,5n6 layer","fraction of cohesive sediment of class6,6n6 layer","fraction of cohesive sediment of class7,7n6 layer","fraction of cohesive sediment of class8,8n6 layer","fraction of cohesive sediment of class9,9n6 layer","fraction of cohesive sediment of class10,10n6 layer","fraction of cohesive sediment of class11,11n6 layer","fraction of cohesive sediment of class12,12n6 layer","fraction of cohesive sediment of class13,13n6 layer","fraction of cohesive sediment of class14,14n6 layer","fraction of cohesive sediment of class15,15n6 layer","fraction of cohesive sediment of class16,16n6 layer","fraction of cohesive sediment of class17,17n6 layer","fraction of cohesive sediment of class18,18n6 layer","fraction of cohesive sediment of class19,19n6 layer","fraction of cohesive sediment of class**,**n6 layer","fraction of cohesive sediment of class*,*n7 layer","fraction of cohesive sediment of class1,1n7 layer","fraction of cohesive sediment of class2,2n7 layer","fraction of cohesive sediment of class3,3n7 layer","fraction of cohesive sediment of class4,4n7 layer","fraction of cohesive sediment of class5,5n7 layer","fraction of cohesive sediment of class6,6n7 layer","fraction of cohesive sediment of class7,7n7 layer","fraction of cohesive sediment of class8,8n7 layer","fraction of cohesive sediment of class9,9n7 layer","fraction of cohesive sediment of class10,10n7 layer","fraction of cohesive sediment of class11,11n7 layer","fraction of cohesive sediment of class12,12n7 layer","fraction of cohesive sediment of class13,13n7 layer","fraction of cohesive sediment of class14,14n7 layer","fraction of cohesive sediment of class15,15n7 layer","fraction of cohesive sediment of class16,16n7 layer","fraction of cohesive sediment of class17,17n7 layer","fraction of cohesive sediment of class18,18n7 layer","fraction of cohesive sediment of class19,19n7 layer","fraction of cohesive sediment of class**,**n7 layer","fraction of cohesive sediment of class*,*n8 layer","fraction of cohesive sediment of class1,1n8 layer","fraction of cohesive sediment of class2,2n8 layer","fraction of cohesive sediment of class3,3n8 layer","fraction of cohesive sediment of class4,4n8 layer","fraction of cohesive sediment of class5,5n8 layer","fraction of cohesive sediment of class6,6n8 layer","fraction of cohesive sediment of class7,7n8 layer","fraction of cohesive sediment of class8,8n8 layer","fraction of cohesive sediment of class9,9n8 layer","fraction of cohesive sediment of class10,10n8 layer","fraction of cohesive sediment of class11,11n8 layer","fraction of cohesive sediment of class12,12n8 layer","fraction of cohesive sediment of class13,13n8 layer","fraction of cohesive sediment of class14,14n8 layer","fraction of cohesive sediment of class15,15n8 layer","fraction of cohesive sediment of class16,16n8 layer","fraction of cohesive sediment of class17,17n8 layer","fraction of cohesive sediment of class18,18n8 layer","fraction of cohesive sediment of class19,19n8 layer","fraction of cohesive sediment of class**,**n8 layer","fraction of cohesive sediment of class*,*n9 layer","fraction of cohesive sediment of class1,1n9 layer","fraction of cohesive sediment of class2,2n9 layer","fraction of cohesive sediment of class3,3n9 layer","fraction of cohesive sediment of class4,4n9 layer","fraction of cohesive sediment of class5,5n9 layer","fraction of cohesive sediment of class6,6n9 layer","fraction of cohesive sediment of class7,7n9 layer","fraction of cohesive sediment of class8,8n9 layer","fraction of cohesive sediment of class9,9n9 layer","fraction of cohesive sediment of class10,10n9 layer","fraction of cohesive sediment of class11,11n9 layer","fraction of cohesive sediment of class12,12n9 layer","fraction of cohesive sediment of class13,13n9 layer","fraction of cohesive sediment of class14,14n9 layer","fraction of cohesive sediment of class15,15n9 layer","fraction of cohesive sediment of class16,16n9 layer","fraction of cohesive sediment of class17,17n9 layer","fraction of cohesive sediment of class18,18n9 layer","fraction of cohesive sediment of class19,19n9 layer","fraction of cohesive sediment of class**,**n9 layer","fraction of cohesive sediment of class*,*n10 layer","fraction of cohesive sediment of class1,1n10 layer","fraction of cohesive sediment of class2,2n10 layer","fraction of cohesive sediment of class3,3n10 layer","fraction of cohesive sediment of class4,4n10 layer","fraction of cohesive sediment of class5,5n10 layer","fraction of cohesive sediment of class6,6n10 layer","fraction of cohesive sediment of class7,7n10 layer","fraction of cohesive sediment of class8,8n10 layer","fraction of cohesive sediment of class9,9n10 layer","fraction of cohesive sediment of class10,10n10 layer","fraction of cohesive sediment of class11,11n10 layer","fraction of cohesive sediment of class12,12n10 layer","fraction of cohesive sediment of class13,13n10 layer","fraction of cohesive sediment of class14,14n10 layer","fraction of cohesive sediment of class15,15n10 layer","fraction of cohesive sediment of class16,16n10 layer","fraction of cohesive sediment of class17,17n10 layer","fraction of cohesive sediment of class18,18n10 layer","fraction of cohesive sediment of class19,19n10 layer","fraction of cohesive sediment of class**,**n10 layer","fraction of cohesive sediment of class*,*n11 layer","fraction of cohesive sediment of class1,1n11 layer","fraction of cohesive sediment of class2,2n11 layer","fraction of cohesive sediment of class3,3n11 layer","fraction of cohesive sediment of class4,4n11 layer","fraction of cohesive sediment of class5,5n11 layer","fraction of cohesive sediment of class6,6n11 layer","fraction of cohesive sediment of class7,7n11 layer","fraction of cohesive sediment of class8,8n11 layer","fraction of cohesive sediment of class9,9n11 layer","fraction of cohesive sediment of class10,10n11 layer","fraction of cohesive sediment of class11,11n11 layer","fraction of cohesive sediment of class12,12n11 layer","fraction of cohesive sediment of class13,13n11 layer","fraction of cohesive sediment of class14,14n11 layer","fraction of cohesive sediment of class15,15n11 layer","fraction of cohesive sediment of class16,16n11 layer","fraction of cohesive sediment of class17,17n11 layer","fraction of cohesive sediment of class18,18n11 layer","fraction of cohesive sediment of class19,19n11 layer","fraction of cohesive sediment of class**,**n11 layer","fraction of cohesive sediment of class*,*n12 layer","fraction of cohesive sediment of class1,1n12 layer","fraction of cohesive sediment of class2,2n12 layer","fraction of cohesive sediment of class3,3n12 layer","fraction of cohesive sediment of class4,4n12 layer","fraction of cohesive sediment of class5,5n12 layer","fraction of cohesive sediment of class6,6n12 layer","fraction of cohesive sediment of class7,7n12 layer","fraction of cohesive sediment of class8,8n12 layer","fraction of cohesive sediment of class9,9n12 layer","fraction of cohesive sediment of class10,10n12 layer","fraction of cohesive sediment of class11,11n12 layer","fraction of cohesive sediment of class12,12n12 layer","fraction of cohesive sediment of class13,13n12 layer","fraction of cohesive sediment of class14,14n12 layer","fraction of cohesive sediment of class15,15n12 layer","fraction of cohesive sediment of class16,16n12 layer","fraction of cohesive sediment of class17,17n12 layer","fraction of cohesive sediment of class18,18n12 layer","fraction of cohesive sediment of class19,19n12 layer","fraction of cohesive sediment of class**,**n12 layer","fraction of cohesive sediment of class*,*n13 layer","fraction of cohesive sediment of class1,1n13 layer","fraction of cohesive sediment of class2,2n13 layer","fraction of cohesive sediment of class3,3n13 layer","fraction of cohesive sediment of class4,4n13 layer","fraction of cohesive sediment of class5,5n13 layer","fraction of cohesive sediment of class6,6n13 layer","fraction of cohesive sediment of class7,7n13 layer","fraction of cohesive sediment of class8,8n13 layer","fraction of cohesive sediment of class9,9n13 layer","fraction of cohesive sediment of class10,10n13 layer","fraction of cohesive sediment of class11,11n13 layer","fraction of cohesive sediment of class12,12n13 layer","fraction of cohesive sediment of class13,13n13 layer","fraction of cohesive sediment of class14,14n13 layer","fraction of cohesive sediment of class15,15n13 layer","fraction of cohesive sediment of class16,16n13 layer","fraction of cohesive sediment of class17,17n13 layer","fraction of cohesive sediment of class18,18n13 layer","fraction of cohesive sediment of class19,19n13 layer","fraction of cohesive sediment of class**,**n13 layer","fraction of cohesive sediment of class*,*n14 layer","fraction of cohesive sediment of class1,1n14 layer","fraction of cohesive sediment of class2,2n14 layer","fraction of cohesive sediment of class3,3n14 layer","fraction of cohesive sediment of class4,4n14 layer","fraction of cohesive sediment of class5,5n14 layer","fraction of cohesive sediment of class6,6n14 layer","fraction of cohesive sediment of class7,7n14 layer","fraction of cohesive sediment of class8,8n14 layer","fraction of cohesive sediment of class9,9n14 layer","fraction of cohesive sediment of class10,10n14 layer","fraction of cohesive sediment of class11,11n14 layer","fraction of cohesive sediment of class12,12n14 layer","fraction of cohesive sediment of class13,13n14 layer","fraction of cohesive sediment of class14,14n14 layer","fraction of cohesive sediment of class15,15n14 layer","fraction of cohesive sediment of class16,16n14 layer","fraction of cohesive sediment of class17,17n14 layer","fraction of cohesive sediment of class18,18n14 layer","fraction of cohesive sediment of class19,19n14 layer","fraction of cohesive sediment of class**,**n14 layer","fraction of cohesive sediment of class*,*n15 layer","fraction of cohesive sediment of class1,1n15 layer","fraction of cohesive sediment of class2,2n15 layer","fraction of cohesive sediment of class3,3n15 layer","fraction of cohesive sediment of class4,4n15 layer","fraction of cohesive sediment of class5,5n15 layer","fraction of cohesive sediment of class6,6n15 layer","fraction of cohesive sediment of class7,7n15 layer","fraction of cohesive sediment of class8,8n15 layer","fraction of cohesive sediment of class9,9n15 layer","fraction of cohesive sediment of class10,10n15 layer","fraction of cohesive sediment of class11,11n15 layer","fraction of cohesive sediment of class12,12n15 layer","fraction of cohesive sediment of class13,13n15 layer","fraction of cohesive sediment of class14,14n15 layer","fraction of cohesive sediment of class15,15n15 layer","fraction of cohesive sediment of class16,16n15 layer","fraction of cohesive sediment of class17,17n15 layer","fraction of cohesive sediment of class18,18n15 layer","fraction of cohesive sediment of class19,19n15 layer","fraction of cohesive sediment of class**,**n15 layer","fraction of cohesive sediment of class*,*n16 layer","fraction of cohesive sediment of class1,1n16 layer","fraction of cohesive sediment of class2,2n16 layer","fraction of cohesive sediment of class3,3n16 layer","fraction of cohesive sediment of class4,4n16 layer","fraction of cohesive sediment of class5,5n16 layer","fraction of cohesive sediment of class6,6n16 layer","fraction of cohesive sediment of class7,7n16 layer","fraction of cohesive sediment of class8,8n16 layer","fraction of cohesive sediment of class9,9n16 layer","fraction of cohesive sediment of class10,10n16 layer","fraction of cohesive sediment of class11,11n16 layer","fraction of cohesive sediment of class12,12n16 layer","fraction of cohesive sediment of class13,13n16 layer","fraction of cohesive sediment of class14,14n16 layer","fraction of cohesive sediment of class15,15n16 layer","fraction of cohesive sediment of class16,16n16 layer","fraction of cohesive sediment of class17,17n16 layer","fraction of cohesive sediment of class18,18n16 layer","fraction of cohesive sediment of class19,19n16 layer","fraction of cohesive sediment of class**,**n16 layer","fraction of cohesive sediment of class*,*n17 layer","fraction of cohesive sediment of class1,1n17 layer","fraction of cohesive sediment of class2,2n17 layer","fraction of cohesive sediment of class3,3n17 layer","fraction of cohesive sediment of class4,4n17 layer","fraction of cohesive sediment of class5,5n17 layer","fraction of cohesive sediment of class6,6n17 layer","fraction of cohesive sediment of class7,7n17 layer","fraction of cohesive sediment of class8,8n17 layer","fraction of cohesive sediment of class9,9n17 layer","fraction of cohesive sediment of class10,10n17 layer","fraction of cohesive sediment of class11,11n17 layer","fraction of cohesive sediment of class12,12n17 layer","fraction of cohesive sediment of class13,13n17 layer","fraction of cohesive sediment of class14,14n17 layer","fraction of cohesive sediment of class15,15n17 layer","fraction of cohesive sediment of class16,16n17 layer","fraction of cohesive sediment of class17,17n17 layer","fraction of cohesive sediment of class18,18n17 layer","fraction of cohesive sediment of class19,19n17 layer","fraction of cohesive sediment of class**,**n17 layer","fraction of cohesive sediment of class*,*n18 layer","fraction of cohesive sediment of class1,1n18 layer","fraction of cohesive sediment of class2,2n18 layer","fraction of cohesive sediment of class3,3n18 layer","fraction of cohesive sediment of class4,4n18 layer","fraction of cohesive sediment of class5,5n18 layer","fraction of cohesive sediment of class6,6n18 layer","fraction of cohesive sediment of class7,7n18 layer","fraction of cohesive sediment of class8,8n18 layer","fraction of cohesive sediment of class9,9n18 layer","fraction of cohesive sediment of class10,10n18 layer","fraction of cohesive sediment of class11,11n18 layer","fraction of cohesive sediment of class12,12n18 layer","fraction of cohesive sediment of class13,13n18 layer","fraction of cohesive sediment of class14,14n18 layer","fraction of cohesive sediment of class15,15n18 layer","fraction of cohesive sediment of class16,16n18 layer","fraction of cohesive sediment of class17,17n18 layer","fraction of cohesive sediment of class18,18n18 layer","fraction of cohesive sediment of class19,19n18 layer","fraction of cohesive sediment of class**,**n18 layer","fraction of cohesive sediment of class*,*n19 layer","fraction of cohesive sediment of class1,1n19 layer","fraction of cohesive sediment of class2,2n19 layer","fraction of cohesive sediment of class3,3n19 layer","fraction of cohesive sediment of class4,4n19 layer","fraction of cohesive sediment of class5,5n19 layer","fraction of cohesive sediment of class6,6n19 layer","fraction of cohesive sediment of class7,7n19 layer","fraction of cohesive sediment of class8,8n19 layer","fraction of cohesive sediment of class9,9n19 layer","fraction of cohesive sediment of class10,10n19 layer","fraction of cohesive sediment of class11,11n19 layer","fraction of cohesive sediment of class12,12n19 layer","fraction of cohesive sediment of class13,13n19 layer","fraction of cohesive sediment of class14,14n19 layer","fraction of cohesive sediment of class15,15n19 layer","fraction of cohesive sediment of class16,16n19 layer","fraction of cohesive sediment of class17,17n19 layer","fraction of cohesive sediment of class18,18n19 layer","fraction of cohesive sediment of class19,19n19 layer","fraction of cohesive sediment of class**,**n19 layer","fraction of cohesive sediment of class*,*n** layer","fraction of cohesive sediment of class1,1n** layer","fraction of cohesive sediment of class2,2n** layer","fraction of cohesive sediment of class3,3n** layer","fraction of cohesive sediment of class4,4n** layer","fraction of cohesive sediment of class5,5n** layer","fraction of cohesive sediment of class6,6n** layer","fraction of cohesive sediment of class7,7n** layer","fraction of cohesive sediment of class8,8n** layer","fraction of cohesive sediment of class9,9n** layer","fraction of cohesive sediment of class10,10n** layer","fraction of cohesive sediment of class11,11n** layer","fraction of cohesive sediment of class12,12n** layer","fraction of cohesive sediment of class13,13n** layer","fraction of cohesive sediment of class14,14n** layer","fraction of cohesive sediment of class15,15n** layer","fraction of cohesive sediment of class16,16n** layer","fraction of cohesive sediment of class17,17n** layer","fraction of cohesive sediment of class18,18n** layer","fraction of cohesive sediment of class19,19n** layer","fraction of cohesive sediment of class**,**n** layer","porosity of k layer","mass of non cohesive sediment of class*,*n* layer","mass of non cohesive sediment of class1,1n* layer","mass of non cohesive sediment of class2,2n* layer","mass of non cohesive sediment of class3,3n* layer","mass of non cohesive sediment of class4,4n* layer","mass of non cohesive sediment of class5,5n* layer","mass of non cohesive sediment of class6,6n* layer","mass of non cohesive sediment of class7,7n* layer","mass of non cohesive sediment of class8,8n* layer","mass of non cohesive sediment of class9,9n* layer","mass of non cohesive sediment of class10,10n* layer","mass of non cohesive sediment of class11,11n* layer","mass of non cohesive sediment of class12,12n* layer","mass of non cohesive sediment of class13,13n* layer","mass of non cohesive sediment of class14,14n* layer","mass of non cohesive sediment of class15,15n* layer","mass of non cohesive sediment of class16,16n* layer","mass of non cohesive sediment of class17,17n* layer","mass of non cohesive sediment of class18,18n* layer","mass of non cohesive sediment of class19,19n* layer","mass of non cohesive sediment of class**,**n* layer","mass of non cohesive sediment of class*,*n1 layer","mass of non cohesive sediment of class1,1n1 layer","mass of non cohesive sediment of class2,2n1 layer","mass of non cohesive sediment of class3,3n1 layer","mass of non cohesive sediment of class4,4n1 layer","mass of non cohesive sediment of class5,5n1 layer","mass of non cohesive sediment of class6,6n1 layer","mass of non cohesive sediment of class7,7n1 layer","mass of non cohesive sediment of class8,8n1 layer","mass of non cohesive sediment of class9,9n1 layer","mass of non cohesive sediment of class10,10n1 layer","mass of non cohesive sediment of class11,11n1 layer","mass of non cohesive sediment of class12,12n1 layer","mass of non cohesive sediment of class13,13n1 layer","mass of non cohesive sediment of class14,14n1 layer","mass of non cohesive sediment of class15,15n1 layer","mass of non cohesive sediment of class16,16n1 layer","mass of non cohesive sediment of class17,17n1 layer","mass of non cohesive sediment of class18,18n1 layer","mass of non cohesive sediment of class19,19n1 layer","mass of non cohesive sediment of class**,**n1 layer","mass of non cohesive sediment of class*,*n2 layer","mass of non cohesive sediment of class1,1n2 layer","mass of non cohesive sediment of class2,2n2 layer","mass of non cohesive sediment of class3,3n2 layer","mass of non cohesive sediment of class4,4n2 layer","mass of non cohesive sediment of class5,5n2 layer","mass of non cohesive sediment of class6,6n2 layer","mass of non cohesive sediment of class7,7n2 layer","mass of non cohesive sediment of class8,8n2 layer","mass of non cohesive sediment of class9,9n2 layer","mass of non cohesive sediment of class10,10n2 layer","mass of non cohesive sediment of class11,11n2 layer","mass of non cohesive sediment of class12,12n2 layer","mass of non cohesive sediment of class13,13n2 layer","mass of non cohesive sediment of class14,14n2 layer","mass of non cohesive sediment of class15,15n2 layer","mass of non cohesive sediment of class16,16n2 layer","mass of non cohesive sediment of class17,17n2 layer","mass of non cohesive sediment of class18,18n2 layer","mass of non cohesive sediment of class19,19n2 layer","mass of non cohesive sediment of class**,**n2 layer","mass of non cohesive sediment of class*,*n3 layer","mass of non cohesive sediment of class1,1n3 layer","mass of non cohesive sediment of class2,2n3 layer","mass of non cohesive sediment of class3,3n3 layer","mass of non cohesive sediment of class4,4n3 layer","mass of non cohesive sediment of class5,5n3 layer","mass of non cohesive sediment of class6,6n3 layer","mass of non cohesive sediment of class7,7n3 layer","mass of non cohesive sediment of class8,8n3 layer","mass of non cohesive sediment of class9,9n3 layer","mass of non cohesive sediment of class10,10n3 layer","mass of non cohesive sediment of class11,11n3 layer","mass of non cohesive sediment of class12,12n3 layer","mass of non cohesive sediment of class13,13n3 layer","mass of non cohesive sediment of class14,14n3 layer","mass of non cohesive sediment of class15,15n3 layer","mass of non cohesive sediment of class16,16n3 layer","mass of non cohesive sediment of class17,17n3 layer","mass of non cohesive sediment of class18,18n3 layer","mass of non cohesive sediment of class19,19n3 layer","mass of non cohesive sediment of class**,**n3 layer","mass of non cohesive sediment of class*,*n4 layer","mass of non cohesive sediment of class1,1n4 layer","mass of non cohesive sediment of class2,2n4 layer","mass of non cohesive sediment of class3,3n4 layer","mass of non cohesive sediment of class4,4n4 layer","mass of non cohesive sediment of class5,5n4 layer","mass of non cohesive sediment of class6,6n4 layer","mass of non cohesive sediment of class7,7n4 layer","mass of non cohesive sediment of class8,8n4 layer","mass of non cohesive sediment of class9,9n4 layer","mass of non cohesive sediment of class10,10n4 layer","mass of non cohesive sediment of class11,11n4 layer","mass of non cohesive sediment of class12,12n4 layer","mass of non cohesive sediment of class13,13n4 layer","mass of non cohesive sediment of class14,14n4 layer","mass of non cohesive sediment of class15,15n4 layer","mass of non cohesive sediment of class16,16n4 layer","mass of non cohesive sediment of class17,17n4 layer","mass of non cohesive sediment of class18,18n4 layer","mass of non cohesive sediment of class19,19n4 layer","mass of non cohesive sediment of class**,**n4 layer","mass of non cohesive sediment of class*,*n5 layer","mass of non cohesive sediment of class1,1n5 layer","mass of non cohesive sediment of class2,2n5 layer","mass of non cohesive sediment of class3,3n5 layer","mass of non cohesive sediment of class4,4n5 layer","mass of non cohesive sediment of class5,5n5 layer","mass of non cohesive sediment of class6,6n5 layer","mass of non cohesive sediment of class7,7n5 layer","mass of non cohesive sediment of class8,8n5 layer","mass of non cohesive sediment of class9,9n5 layer","mass of non cohesive sediment of class10,10n5 layer","mass of non cohesive sediment of class11,11n5 layer","mass of non cohesive sediment of class12,12n5 layer","mass of non cohesive sediment of class13,13n5 layer","mass of non cohesive sediment of class14,14n5 layer","mass of non cohesive sediment of class15,15n5 layer","mass of non cohesive sediment of class16,16n5 layer","mass of non cohesive sediment of class17,17n5 layer","mass of non cohesive sediment of class18,18n5 layer","mass of non cohesive sediment of class19,19n5 layer","mass of non cohesive sediment of class**,**n5 layer","mass of non cohesive sediment of class*,*n6 layer","mass of non cohesive sediment of class1,1n6 layer","mass of non cohesive sediment of class2,2n6 layer","mass of non cohesive sediment of class3,3n6 layer","mass of non cohesive sediment of class4,4n6 layer","mass of non cohesive sediment of class5,5n6 layer","mass of non cohesive sediment of class6,6n6 layer","mass of non cohesive sediment of class7,7n6 layer","mass of non cohesive sediment of class8,8n6 layer","mass of non cohesive sediment of class9,9n6 layer","mass of non cohesive sediment of class10,10n6 layer","mass of non cohesive sediment of class11,11n6 layer","mass of non cohesive sediment of class12,12n6 layer","mass of non cohesive sediment of class13,13n6 layer","mass of non cohesive sediment of class14,14n6 layer","mass of non cohesive sediment of class15,15n6 layer","mass of non cohesive sediment of class16,16n6 layer","mass of non cohesive sediment of class17,17n6 layer","mass of non cohesive sediment of class18,18n6 layer","mass of non cohesive sediment of class19,19n6 layer","mass of non cohesive sediment of class**,**n6 layer","mass of non cohesive sediment of class*,*n7 layer","mass of non cohesive sediment of class1,1n7 layer","mass of non cohesive sediment of class2,2n7 layer","mass of non cohesive sediment of class3,3n7 layer","mass of non cohesive sediment of class4,4n7 layer","mass of non cohesive sediment of class5,5n7 layer","mass of non cohesive sediment of class6,6n7 layer","mass of non cohesive sediment of class7,7n7 layer","mass of non cohesive sediment of class8,8n7 layer","mass of non cohesive sediment of class9,9n7 layer","mass of non cohesive sediment of class10,10n7 layer","mass of non cohesive sediment of class11,11n7 layer","mass of non cohesive sediment of class12,12n7 layer","mass of non cohesive sediment of class13,13n7 layer","mass of non cohesive sediment of class14,14n7 layer","mass of non cohesive sediment of class15,15n7 layer","mass of non cohesive sediment of class16,16n7 layer","mass of non cohesive sediment of class17,17n7 layer","mass of non cohesive sediment of class18,18n7 layer","mass of non cohesive sediment of class19,19n7 layer","mass of non cohesive sediment of class**,**n7 layer","mass of non cohesive sediment of class*,*n8 layer","mass of non cohesive sediment of class1,1n8 layer","mass of non cohesive sediment of class2,2n8 layer","mass of non cohesive sediment of class3,3n8 layer","mass of non cohesive sediment of class4,4n8 layer","mass of non cohesive sediment of class5,5n8 layer","mass of non cohesive sediment of class6,6n8 layer","mass of non cohesive sediment of class7,7n8 layer","mass of non cohesive sediment of class8,8n8 layer","mass of non cohesive sediment of class9,9n8 layer","mass of non cohesive sediment of class10,10n8 layer","mass of non cohesive sediment of class11,11n8 layer","mass of non cohesive sediment of class12,12n8 layer","mass of non cohesive sediment of class13,13n8 layer","mass of non cohesive sediment of class14,14n8 layer","mass of non cohesive sediment of class15,15n8 layer","mass of non cohesive sediment of class16,16n8 layer","mass of non cohesive sediment of class17,17n8 layer","mass of non cohesive sediment of class18,18n8 layer","mass of non cohesive sediment of class19,19n8 layer","mass of non cohesive sediment of class**,**n8 layer","mass of non cohesive sediment of class*,*n9 layer","mass of non cohesive sediment of class1,1n9 layer","mass of non cohesive sediment of class2,2n9 layer","mass of non cohesive sediment of class3,3n9 layer","mass of non cohesive sediment of class4,4n9 layer","mass of non cohesive sediment of class5,5n9 layer","mass of non cohesive sediment of class6,6n9 layer","mass of non cohesive sediment of class7,7n9 layer","mass of non cohesive sediment of class8,8n9 layer","mass of non cohesive sediment of class9,9n9 layer","mass of non cohesive sediment of class10,10n9 layer","mass of non cohesive sediment of class11,11n9 layer","mass of non cohesive sediment of class12,12n9 layer","mass of non cohesive sediment of class13,13n9 layer","mass of non cohesive sediment of class14,14n9 layer","mass of non cohesive sediment of class15,15n9 layer","mass of non cohesive sediment of class16,16n9 layer","mass of non cohesive sediment of class17,17n9 layer","mass of non cohesive sediment of class18,18n9 layer","mass of non cohesive sediment of class19,19n9 layer","mass of non cohesive sediment of class**,**n9 layer","mass of non cohesive sediment of class*,*n10 layer","mass of non cohesive sediment of class1,1n10 layer","mass of non cohesive sediment of class2,2n10 layer","mass of non cohesive sediment of class3,3n10 layer","mass of non cohesive sediment of class4,4n10 layer","mass of non cohesive sediment of class5,5n10 layer","mass of non cohesive sediment of class6,6n10 layer","mass of non cohesive sediment of class7,7n10 layer","mass of non cohesive sediment of class8,8n10 layer","mass of non cohesive sediment of class9,9n10 layer","mass of non cohesive sediment of class10,10n10 layer","mass of non cohesive sediment of class11,11n10 layer","mass of non cohesive sediment of class12,12n10 layer","mass of non cohesive sediment of class13,13n10 layer","mass of non cohesive sediment of class14,14n10 layer","mass of non cohesive sediment of class15,15n10 layer","mass of non cohesive sediment of class16,16n10 layer","mass of non cohesive sediment of class17,17n10 layer","mass of non cohesive sediment of class18,18n10 layer","mass of non cohesive sediment of class19,19n10 layer","mass of non cohesive sediment of class**,**n10 layer","mass of non cohesive sediment of class*,*n11 layer","mass of non cohesive sediment of class1,1n11 layer","mass of non cohesive sediment of class2,2n11 layer","mass of non cohesive sediment of class3,3n11 layer","mass of non cohesive sediment of class4,4n11 layer","mass of non cohesive sediment of class5,5n11 layer","mass of non cohesive sediment of class6,6n11 layer","mass of non cohesive sediment of class7,7n11 layer","mass of non cohesive sediment of class8,8n11 layer","mass of non cohesive sediment of class9,9n11 layer","mass of non cohesive sediment of class10,10n11 layer","mass of non cohesive sediment of class11,11n11 layer","mass of non cohesive sediment of class12,12n11 layer","mass of non cohesive sediment of class13,13n11 layer","mass of non cohesive sediment of class14,14n11 layer","mass of non cohesive sediment of class15,15n11 layer","mass of non cohesive sediment of class16,16n11 layer","mass of non cohesive sediment of class17,17n11 layer","mass of non cohesive sediment of class18,18n11 layer","mass of non cohesive sediment of class19,19n11 layer","mass of non cohesive sediment of class**,**n11 layer","mass of non cohesive sediment of class*,*n12 layer","mass of non cohesive sediment of class1,1n12 layer","mass of non cohesive sediment of class2,2n12 layer","mass of non cohesive sediment of class3,3n12 layer","mass of non cohesive sediment of class4,4n12 layer","mass of non cohesive sediment of class5,5n12 layer","mass of non cohesive sediment of class6,6n12 layer","mass of non cohesive sediment of class7,7n12 layer","mass of non cohesive sediment of class8,8n12 layer","mass of non cohesive sediment of class9,9n12 layer","mass of non cohesive sediment of class10,10n12 layer","mass of non cohesive sediment of class11,11n12 layer","mass of non cohesive sediment of class12,12n12 layer","mass of non cohesive sediment of class13,13n12 layer","mass of non cohesive sediment of class14,14n12 layer","mass of non cohesive sediment of class15,15n12 layer","mass of non cohesive sediment of class16,16n12 layer","mass of non cohesive sediment of class17,17n12 layer","mass of non cohesive sediment of class18,18n12 layer","mass of non cohesive sediment of class19,19n12 layer","mass of non cohesive sediment of class**,**n12 layer","mass of non cohesive sediment of class*,*n13 layer","mass of non cohesive sediment of class1,1n13 layer","mass of non cohesive sediment of class2,2n13 layer","mass of non cohesive sediment of class3,3n13 layer","mass of non cohesive sediment of class4,4n13 layer","mass of non cohesive sediment of class5,5n13 layer","mass of non cohesive sediment of class6,6n13 layer","mass of non cohesive sediment of class7,7n13 layer","mass of non cohesive sediment of class8,8n13 layer","mass of non cohesive sediment of class9,9n13 layer","mass of non cohesive sediment of class10,10n13 layer","mass of non cohesive sediment of class11,11n13 layer","mass of non cohesive sediment of class12,12n13 layer","mass of non cohesive sediment of class13,13n13 layer","mass of non cohesive sediment of class14,14n13 layer","mass of non cohesive sediment of class15,15n13 layer","mass of non cohesive sediment of class16,16n13 layer","mass of non cohesive sediment of class17,17n13 layer","mass of non cohesive sediment of class18,18n13 layer","mass of non cohesive sediment of class19,19n13 layer","mass of non cohesive sediment of class**,**n13 layer","mass of non cohesive sediment of class*,*n14 layer","mass of non cohesive sediment of class1,1n14 layer","mass of non cohesive sediment of class2,2n14 layer","mass of non cohesive sediment of class3,3n14 layer","mass of non cohesive sediment of class4,4n14 layer","mass of non cohesive sediment of class5,5n14 layer","mass of non cohesive sediment of class6,6n14 layer","mass of non cohesive sediment of class7,7n14 layer","mass of non cohesive sediment of class8,8n14 layer","mass of non cohesive sediment of class9,9n14 layer","mass of non cohesive sediment of class10,10n14 layer","mass of non cohesive sediment of class11,11n14 layer","mass of non cohesive sediment of class12,12n14 layer","mass of non cohesive sediment of class13,13n14 layer","mass of non cohesive sediment of class14,14n14 layer","mass of non cohesive sediment of class15,15n14 layer","mass of non cohesive sediment of class16,16n14 layer","mass of non cohesive sediment of class17,17n14 layer","mass of non cohesive sediment of class18,18n14 layer","mass of non cohesive sediment of class19,19n14 layer","mass of non cohesive sediment of class**,**n14 layer","mass of non cohesive sediment of class*,*n15 layer","mass of non cohesive sediment of class1,1n15 layer","mass of non cohesive sediment of class2,2n15 layer","mass of non cohesive sediment of class3,3n15 layer","mass of non cohesive sediment of class4,4n15 layer","mass of non cohesive sediment of class5,5n15 layer","mass of non cohesive sediment of class6,6n15 layer","mass of non cohesive sediment of class7,7n15 layer","mass of non cohesive sediment of class8,8n15 layer","mass of non cohesive sediment of class9,9n15 layer","mass of non cohesive sediment of class10,10n15 layer","mass of non cohesive sediment of class11,11n15 layer","mass of non cohesive sediment of class12,12n15 layer","mass of non cohesive sediment of class13,13n15 layer","mass of non cohesive sediment of class14,14n15 layer","mass of non cohesive sediment of class15,15n15 layer","mass of non cohesive sediment of class16,16n15 layer","mass of non cohesive sediment of class17,17n15 layer","mass of non cohesive sediment of class18,18n15 layer","mass of non cohesive sediment of class19,19n15 layer","mass of non cohesive sediment of class**,**n15 layer","mass of non cohesive sediment of class*,*n16 layer","mass of non cohesive sediment of class1,1n16 layer","mass of non cohesive sediment of class2,2n16 layer","mass of non cohesive sediment of class3,3n16 layer","mass of non cohesive sediment of class4,4n16 layer","mass of non cohesive sediment of class5,5n16 layer","mass of non cohesive sediment of class6,6n16 layer","mass of non cohesive sediment of class7,7n16 layer","mass of non cohesive sediment of class8,8n16 layer","mass of non cohesive sediment of class9,9n16 layer","mass of non cohesive sediment of class10,10n16 layer","mass of non cohesive sediment of class11,11n16 layer","mass of non cohesive sediment of class12,12n16 layer","mass of non cohesive sediment of class13,13n16 layer","mass of non cohesive sediment of class14,14n16 layer","mass of non cohesive sediment of class15,15n16 layer","mass of non cohesive sediment of class16,16n16 layer","mass of non cohesive sediment of class17,17n16 layer","mass of non cohesive sediment of class18,18n16 layer","mass of non cohesive sediment of class19,19n16 layer","mass of non cohesive sediment of class**,**n16 layer","mass of non cohesive sediment of class*,*n17 layer","mass of non cohesive sediment of class1,1n17 layer","mass of non cohesive sediment of class2,2n17 layer","mass of non cohesive sediment of class3,3n17 layer","mass of non cohesive sediment of class4,4n17 layer","mass of non cohesive sediment of class5,5n17 layer","mass of non cohesive sediment of class6,6n17 layer","mass of non cohesive sediment of class7,7n17 layer","mass of non cohesive sediment of class8,8n17 layer","mass of non cohesive sediment of class9,9n17 layer","mass of non cohesive sediment of class10,10n17 layer","mass of non cohesive sediment of class11,11n17 layer","mass of non cohesive sediment of class12,12n17 layer","mass of non cohesive sediment of class13,13n17 layer","mass of non cohesive sediment of class14,14n17 layer","mass of non cohesive sediment of class15,15n17 layer","mass of non cohesive sediment of class16,16n17 layer","mass of non cohesive sediment of class17,17n17 layer","mass of non cohesive sediment of class18,18n17 layer","mass of non cohesive sediment of class19,19n17 layer","mass of non cohesive sediment of class**,**n17 layer","mass of non cohesive sediment of class*,*n18 layer","mass of non cohesive sediment of class1,1n18 layer","mass of non cohesive sediment of class2,2n18 layer","mass of non cohesive sediment of class3,3n18 layer","mass of non cohesive sediment of class4,4n18 layer","mass of non cohesive sediment of class5,5n18 layer","mass of non cohesive sediment of class6,6n18 layer","mass of non cohesive sediment of class7,7n18 layer","mass of non cohesive sediment of class8,8n18 layer","mass of non cohesive sediment of class9,9n18 layer","mass of non cohesive sediment of class10,10n18 layer","mass of non cohesive sediment of class11,11n18 layer","mass of non cohesive sediment of class12,12n18 layer","mass of non cohesive sediment of class13,13n18 layer","mass of non cohesive sediment of class14,14n18 layer","mass of non cohesive sediment of class15,15n18 layer","mass of non cohesive sediment of class16,16n18 layer","mass of non cohesive sediment of class17,17n18 layer","mass of non cohesive sediment of class18,18n18 layer","mass of non cohesive sediment of class19,19n18 layer","mass of non cohesive sediment of class**,**n18 layer","mass of non cohesive sediment of class*,*n19 layer","mass of non cohesive sediment of class1,1n19 layer","mass of non cohesive sediment of class2,2n19 layer","mass of non cohesive sediment of class3,3n19 layer","mass of non cohesive sediment of class4,4n19 layer","mass of non cohesive sediment of class5,5n19 layer","mass of non cohesive sediment of class6,6n19 layer","mass of non cohesive sediment of class7,7n19 layer","mass of non cohesive sediment of class8,8n19 layer","mass of non cohesive sediment of class9,9n19 layer","mass of non cohesive sediment of class10,10n19 layer","mass of non cohesive sediment of class11,11n19 layer","mass of non cohesive sediment of class12,12n19 layer","mass of non cohesive sediment of class13,13n19 layer","mass of non cohesive sediment of class14,14n19 layer","mass of non cohesive sediment of class15,15n19 layer","mass of non cohesive sediment of class16,16n19 layer","mass of non cohesive sediment of class17,17n19 layer","mass of non cohesive sediment of class18,18n19 layer","mass of non cohesive sediment of class19,19n19 layer","mass of non cohesive sediment of class**,**n19 layer","mass of non cohesive sediment of class*,*n** layer","mass of non cohesive sediment of class1,1n** layer","mass of non cohesive sediment of class2,2n** layer","mass of non cohesive sediment of class3,3n** layer","mass of non cohesive sediment of class4,4n** layer","mass of non cohesive sediment of class5,5n** layer","mass of non cohesive sediment of class6,6n** layer","mass of non cohesive sediment of class7,7n** layer","mass of non cohesive sediment of class8,8n** layer","mass of non cohesive sediment of class9,9n** layer","mass of non cohesive sediment of class10,10n** layer","mass of non cohesive sediment of class11,11n** layer","mass of non cohesive sediment of class12,12n** layer","mass of non cohesive sediment of class13,13n** layer","mass of non cohesive sediment of class14,14n** layer","mass of non cohesive sediment of class15,15n** layer","mass of non cohesive sediment of class16,16n** layer","mass of non cohesive sediment of class17,17n** layer","mass of non cohesive sediment of class18,18n** layer","mass of non cohesive sediment of class19,19n** layer","mass of non cohesive sediment of class**,**n** layer","mass of cohesive sediment of class*,*n* layer","mass of cohesive sediment of class1,1n* layer","mass of cohesive sediment of class2,2n* layer","mass of cohesive sediment of class3,3n* layer","mass of cohesive sediment of class4,4n* layer","mass of cohesive sediment of class5,5n* layer","mass of cohesive sediment of class6,6n* layer","mass of cohesive sediment of class7,7n* layer","mass of cohesive sediment of class8,8n* layer","mass of cohesive sediment of class9,9n* layer","mass of cohesive sediment of class10,10n* layer","mass of cohesive sediment of class11,11n* layer","mass of cohesive sediment of class12,12n* layer","mass of cohesive sediment of class13,13n* layer","mass of cohesive sediment of class14,14n* layer","mass of cohesive sediment of class15,15n* layer","mass of cohesive sediment of class16,16n* layer","mass of cohesive sediment of class17,17n* layer","mass of cohesive sediment of class18,18n* layer","mass of cohesive sediment of class19,19n* layer","mass of cohesive sediment of class**,**n* layer","mass of cohesive sediment of class*,*n1 layer","mass of cohesive sediment of class1,1n1 layer","mass of cohesive sediment of class2,2n1 layer","mass of cohesive sediment of class3,3n1 layer","mass of cohesive sediment of class4,4n1 layer","mass of cohesive sediment of class5,5n1 layer","mass of cohesive sediment of class6,6n1 layer","mass of cohesive sediment of class7,7n1 layer","mass of cohesive sediment of class8,8n1 layer","mass of cohesive sediment of class9,9n1 layer","mass of cohesive sediment of class10,10n1 layer","mass of cohesive sediment of class11,11n1 layer","mass of cohesive sediment of class12,12n1 layer","mass of cohesive sediment of class13,13n1 layer","mass of cohesive sediment of class14,14n1 layer","mass of cohesive sediment of class15,15n1 layer","mass of cohesive sediment of class16,16n1 layer","mass of cohesive sediment of class17,17n1 layer","mass of cohesive sediment of class18,18n1 layer","mass of cohesive sediment of class19,19n1 layer","mass of cohesive sediment of class**,**n1 layer","mass of cohesive sediment of class*,*n2 layer","mass of cohesive sediment of class1,1n2 layer","mass of cohesive sediment of class2,2n2 layer","mass of cohesive sediment of class3,3n2 layer","mass of cohesive sediment of class4,4n2 layer","mass of cohesive sediment of class5,5n2 layer","mass of cohesive sediment of class6,6n2 layer","mass of cohesive sediment of class7,7n2 layer","mass of cohesive sediment of class8,8n2 layer","mass of cohesive sediment of class9,9n2 layer","mass of cohesive sediment of class10,10n2 layer","mass of cohesive sediment of class11,11n2 layer","mass of cohesive sediment of class12,12n2 layer","mass of cohesive sediment of class13,13n2 layer","mass of cohesive sediment of class14,14n2 layer","mass of cohesive sediment of class15,15n2 layer","mass of cohesive sediment of class16,16n2 layer","mass of cohesive sediment of class17,17n2 layer","mass of cohesive sediment of class18,18n2 layer","mass of cohesive sediment of class19,19n2 layer","mass of cohesive sediment of class**,**n2 layer","mass of cohesive sediment of class*,*n3 layer","mass of cohesive sediment of class1,1n3 layer","mass of cohesive sediment of class2,2n3 layer","mass of cohesive sediment of class3,3n3 layer","mass of cohesive sediment of class4,4n3 layer","mass of cohesive sediment of class5,5n3 layer","mass of cohesive sediment of class6,6n3 layer","mass of cohesive sediment of class7,7n3 layer","mass of cohesive sediment of class8,8n3 layer","mass of cohesive sediment of class9,9n3 layer","mass of cohesive sediment of class10,10n3 layer","mass of cohesive sediment of class11,11n3 layer","mass of cohesive sediment of class12,12n3 layer","mass of cohesive sediment of class13,13n3 layer","mass of cohesive sediment of class14,14n3 layer","mass of cohesive sediment of class15,15n3 layer","mass of cohesive sediment of class16,16n3 layer","mass of cohesive sediment of class17,17n3 layer","mass of cohesive sediment of class18,18n3 layer","mass of cohesive sediment of class19,19n3 layer","mass of cohesive sediment of class**,**n3 layer","mass of cohesive sediment of class*,*n4 layer","mass of cohesive sediment of class1,1n4 layer","mass of cohesive sediment of class2,2n4 layer","mass of cohesive sediment of class3,3n4 layer","mass of cohesive sediment of class4,4n4 layer","mass of cohesive sediment of class5,5n4 layer","mass of cohesive sediment of class6,6n4 layer","mass of cohesive sediment of class7,7n4 layer","mass of cohesive sediment of class8,8n4 layer","mass of cohesive sediment of class9,9n4 layer","mass of cohesive sediment of class10,10n4 layer","mass of cohesive sediment of class11,11n4 layer","mass of cohesive sediment of class12,12n4 layer","mass of cohesive sediment of class13,13n4 layer","mass of cohesive sediment of class14,14n4 layer","mass of cohesive sediment of class15,15n4 layer","mass of cohesive sediment of class16,16n4 layer","mass of cohesive sediment of class17,17n4 layer","mass of cohesive sediment of class18,18n4 layer","mass of cohesive sediment of class19,19n4 layer","mass of cohesive sediment of class**,**n4 layer","mass of cohesive sediment of class*,*n5 layer","mass of cohesive sediment of class1,1n5 layer","mass of cohesive sediment of class2,2n5 layer","mass of cohesive sediment of class3,3n5 layer","mass of cohesive sediment of class4,4n5 layer","mass of cohesive sediment of class5,5n5 layer","mass of cohesive sediment of class6,6n5 layer","mass of cohesive sediment of class7,7n5 layer","mass of cohesive sediment of class8,8n5 layer","mass of cohesive sediment of class9,9n5 layer","mass of cohesive sediment of class10,10n5 layer","mass of cohesive sediment of class11,11n5 layer","mass of cohesive sediment of class12,12n5 layer","mass of cohesive sediment of class13,13n5 layer","mass of cohesive sediment of class14,14n5 layer","mass of cohesive sediment of class15,15n5 layer","mass of cohesive sediment of class16,16n5 layer","mass of cohesive sediment of class17,17n5 layer","mass of cohesive sediment of class18,18n5 layer","mass of cohesive sediment of class19,19n5 layer","mass of cohesive sediment of class**,**n5 layer","mass of cohesive sediment of class*,*n6 layer","mass of cohesive sediment of class1,1n6 layer","mass of cohesive sediment of class2,2n6 layer","mass of cohesive sediment of class3,3n6 layer","mass of cohesive sediment of class4,4n6 layer","mass of cohesive sediment of class5,5n6 layer","mass of cohesive sediment of class6,6n6 layer","mass of cohesive sediment of class7,7n6 layer","mass of cohesive sediment of class8,8n6 layer","mass of cohesive sediment of class9,9n6 layer","mass of cohesive sediment of class10,10n6 layer","mass of cohesive sediment of class11,11n6 layer","mass of cohesive sediment of class12,12n6 layer","mass of cohesive sediment of class13,13n6 layer","mass of cohesive sediment of class14,14n6 layer","mass of cohesive sediment of class15,15n6 layer","mass of cohesive sediment of class16,16n6 layer","mass of cohesive sediment of class17,17n6 layer","mass of cohesive sediment of class18,18n6 layer","mass of cohesive sediment of class19,19n6 layer","mass of cohesive sediment of class**,**n6 layer","mass of cohesive sediment of class*,*n7 layer","mass of cohesive sediment of class1,1n7 layer","mass of cohesive sediment of class2,2n7 layer","mass of cohesive sediment of class3,3n7 layer","mass of cohesive sediment of class4,4n7 layer","mass of cohesive sediment of class5,5n7 layer","mass of cohesive sediment of class6,6n7 layer","mass of cohesive sediment of class7,7n7 layer","mass of cohesive sediment of class8,8n7 layer","mass of cohesive sediment of class9,9n7 layer","mass of cohesive sediment of class10,10n7 layer","mass of cohesive sediment of class11,11n7 layer","mass of cohesive sediment of class12,12n7 layer","mass of cohesive sediment of class13,13n7 layer","mass of cohesive sediment of class14,14n7 layer","mass of cohesive sediment of class15,15n7 layer","mass of cohesive sediment of class16,16n7 layer","mass of cohesive sediment of class17,17n7 layer","mass of cohesive sediment of class18,18n7 layer","mass of cohesive sediment of class19,19n7 layer","mass of cohesive sediment of class**,**n7 layer","mass of cohesive sediment of class*,*n8 layer","mass of cohesive sediment of class1,1n8 layer","mass of cohesive sediment of class2,2n8 layer","mass of cohesive sediment of class3,3n8 layer","mass of cohesive sediment of class4,4n8 layer","mass of cohesive sediment of class5,5n8 layer","mass of cohesive sediment of class6,6n8 layer","mass of cohesive sediment of class7,7n8 layer","mass of cohesive sediment of class8,8n8 layer","mass of cohesive sediment of class9,9n8 layer","mass of cohesive sediment of class10,10n8 layer","mass of cohesive sediment of class11,11n8 layer","mass of cohesive sediment of class12,12n8 layer","mass of cohesive sediment of class13,13n8 layer","mass of cohesive sediment of class14,14n8 layer","mass of cohesive sediment of class15,15n8 layer","mass of cohesive sediment of class16,16n8 layer","mass of cohesive sediment of class17,17n8 layer","mass of cohesive sediment of class18,18n8 layer","mass of cohesive sediment of class19,19n8 layer","mass of cohesive sediment of class**,**n8 layer","mass of cohesive sediment of class*,*n9 layer","mass of cohesive sediment of class1,1n9 layer","mass of cohesive sediment of class2,2n9 layer","mass of cohesive sediment of class3,3n9 layer","mass of cohesive sediment of class4,4n9 layer","mass of cohesive sediment of class5,5n9 layer","mass of cohesive sediment of class6,6n9 layer","mass of cohesive sediment of class7,7n9 layer","mass of cohesive sediment of class8,8n9 layer","mass of cohesive sediment of class9,9n9 layer","mass of cohesive sediment of class10,10n9 layer","mass of cohesive sediment of class11,11n9 layer","mass of cohesive sediment of class12,12n9 layer","mass of cohesive sediment of class13,13n9 layer","mass of cohesive sediment of class14,14n9 layer","mass of cohesive sediment of class15,15n9 layer","mass of cohesive sediment of class16,16n9 layer","mass of cohesive sediment of class17,17n9 layer","mass of cohesive sediment of class18,18n9 layer","mass of cohesive sediment of class19,19n9 layer","mass of cohesive sediment of class**,**n9 layer","mass of cohesive sediment of class*,*n10 layer","mass of cohesive sediment of class1,1n10 layer","mass of cohesive sediment of class2,2n10 layer","mass of cohesive sediment of class3,3n10 layer","mass of cohesive sediment of class4,4n10 layer","mass of cohesive sediment of class5,5n10 layer","mass of cohesive sediment of class6,6n10 layer","mass of cohesive sediment of class7,7n10 layer","mass of cohesive sediment of class8,8n10 layer","mass of cohesive sediment of class9,9n10 layer","mass of cohesive sediment of class10,10n10 layer","mass of cohesive sediment of class11,11n10 layer","mass of cohesive sediment of class12,12n10 layer","mass of cohesive sediment of class13,13n10 layer","mass of cohesive sediment of class14,14n10 layer","mass of cohesive sediment of class15,15n10 layer","mass of cohesive sediment of class16,16n10 layer","mass of cohesive sediment of class17,17n10 layer","mass of cohesive sediment of class18,18n10 layer","mass of cohesive sediment of class19,19n10 layer","mass of cohesive sediment of class**,**n10 layer","mass of cohesive sediment of class*,*n11 layer","mass of cohesive sediment of class1,1n11 layer","mass of cohesive sediment of class2,2n11 layer","mass of cohesive sediment of class3,3n11 layer","mass of cohesive sediment of class4,4n11 layer","mass of cohesive sediment of class5,5n11 layer","mass of cohesive sediment of class6,6n11 layer","mass of cohesive sediment of class7,7n11 layer","mass of cohesive sediment of class8,8n11 layer","mass of cohesive sediment of class9,9n11 layer","mass of cohesive sediment of class10,10n11 layer","mass of cohesive sediment of class11,11n11 layer","mass of cohesive sediment of class12,12n11 layer","mass of cohesive sediment of class13,13n11 layer","mass of cohesive sediment of class14,14n11 layer","mass of cohesive sediment of class15,15n11 layer","mass of cohesive sediment of class16,16n11 layer","mass of cohesive sediment of class17,17n11 layer","mass of cohesive sediment of class18,18n11 layer","mass of cohesive sediment of class19,19n11 layer","mass of cohesive sediment of class**,**n11 layer","mass of cohesive sediment of class*,*n12 layer","mass of cohesive sediment of class1,1n12 layer","mass of cohesive sediment of class2,2n12 layer","mass of cohesive sediment of class3,3n12 layer","mass of cohesive sediment of class4,4n12 layer","mass of cohesive sediment of class5,5n12 layer","mass of cohesive sediment of class6,6n12 layer","mass of cohesive sediment of class7,7n12 layer","mass of cohesive sediment of class8,8n12 layer","mass of cohesive sediment of class9,9n12 layer","mass of cohesive sediment of class10,10n12 layer","mass of cohesive sediment of class11,11n12 layer","mass of cohesive sediment of class12,12n12 layer","mass of cohesive sediment of class13,13n12 layer","mass of cohesive sediment of class14,14n12 layer","mass of cohesive sediment of class15,15n12 layer","mass of cohesive sediment of class16,16n12 layer","mass of cohesive sediment of class17,17n12 layer","mass of cohesive sediment of class18,18n12 layer","mass of cohesive sediment of class19,19n12 layer","mass of cohesive sediment of class**,**n12 layer","mass of cohesive sediment of class*,*n13 layer","mass of cohesive sediment of class1,1n13 layer","mass of cohesive sediment of class2,2n13 layer","mass of cohesive sediment of class3,3n13 layer","mass of cohesive sediment of class4,4n13 layer","mass of cohesive sediment of class5,5n13 layer","mass of cohesive sediment of class6,6n13 layer","mass of cohesive sediment of class7,7n13 layer","mass of cohesive sediment of class8,8n13 layer","mass of cohesive sediment of class9,9n13 layer","mass of cohesive sediment of class10,10n13 layer","mass of cohesive sediment of class11,11n13 layer","mass of cohesive sediment of class12,12n13 layer","mass of cohesive sediment of class13,13n13 layer","mass of cohesive sediment of class14,14n13 layer","mass of cohesive sediment of class15,15n13 layer","mass of cohesive sediment of class16,16n13 layer","mass of cohesive sediment of class17,17n13 layer","mass of cohesive sediment of class18,18n13 layer","mass of cohesive sediment of class19,19n13 layer","mass of cohesive sediment of class**,**n13 layer","mass of cohesive sediment of class*,*n14 layer","mass of cohesive sediment of class1,1n14 layer","mass of cohesive sediment of class2,2n14 layer","mass of cohesive sediment of class3,3n14 layer","mass of cohesive sediment of class4,4n14 layer","mass of cohesive sediment of class5,5n14 layer","mass of cohesive sediment of class6,6n14 layer","mass of cohesive sediment of class7,7n14 layer","mass of cohesive sediment of class8,8n14 layer","mass of cohesive sediment of class9,9n14 layer","mass of cohesive sediment of class10,10n14 layer","mass of cohesive sediment of class11,11n14 layer","mass of cohesive sediment of class12,12n14 layer","mass of cohesive sediment of class13,13n14 layer","mass of cohesive sediment of class14,14n14 layer","mass of cohesive sediment of class15,15n14 layer","mass of cohesive sediment of class16,16n14 layer","mass of cohesive sediment of class17,17n14 layer","mass of cohesive sediment of class18,18n14 layer","mass of cohesive sediment of class19,19n14 layer","mass of cohesive sediment of class**,**n14 layer","mass of cohesive sediment of class*,*n15 layer","mass of cohesive sediment of class1,1n15 layer","mass of cohesive sediment of class2,2n15 layer","mass of cohesive sediment of class3,3n15 layer","mass of cohesive sediment of class4,4n15 layer","mass of cohesive sediment of class5,5n15 layer","mass of cohesive sediment of class6,6n15 layer","mass of cohesive sediment of class7,7n15 layer","mass of cohesive sediment of class8,8n15 layer","mass of cohesive sediment of class9,9n15 layer","mass of cohesive sediment of class10,10n15 layer","mass of cohesive sediment of class11,11n15 layer","mass of cohesive sediment of class12,12n15 layer","mass of cohesive sediment of class13,13n15 layer","mass of cohesive sediment of class14,14n15 layer","mass of cohesive sediment of class15,15n15 layer","mass of cohesive sediment of class16,16n15 layer","mass of cohesive sediment of class17,17n15 layer","mass of cohesive sediment of class18,18n15 layer","mass of cohesive sediment of class19,19n15 layer","mass of cohesive sediment of class**,**n15 layer","mass of cohesive sediment of class*,*n16 layer","mass of cohesive sediment of class1,1n16 layer","mass of cohesive sediment of class2,2n16 layer","mass of cohesive sediment of class3,3n16 layer","mass of cohesive sediment of class4,4n16 layer","mass of cohesive sediment of class5,5n16 layer","mass of cohesive sediment of class6,6n16 layer","mass of cohesive sediment of class7,7n16 layer","mass of cohesive sediment of class8,8n16 layer","mass of cohesive sediment of class9,9n16 layer","mass of cohesive sediment of class10,10n16 layer","mass of cohesive sediment of class11,11n16 layer","mass of cohesive sediment of class12,12n16 layer","mass of cohesive sediment of class13,13n16 layer","mass of cohesive sediment of class14,14n16 layer","mass of cohesive sediment of class15,15n16 layer","mass of cohesive sediment of class16,16n16 layer","mass of cohesive sediment of class17,17n16 layer","mass of cohesive sediment of class18,18n16 layer","mass of cohesive sediment of class19,19n16 layer","mass of cohesive sediment of class**,**n16 layer","mass of cohesive sediment of class*,*n17 layer","mass of cohesive sediment of class1,1n17 layer","mass of cohesive sediment of class2,2n17 layer","mass of cohesive sediment of class3,3n17 layer","mass of cohesive sediment of class4,4n17 layer","mass of cohesive sediment of class5,5n17 layer","mass of cohesive sediment of class6,6n17 layer","mass of cohesive sediment of class7,7n17 layer","mass of cohesive sediment of class8,8n17 layer","mass of cohesive sediment of class9,9n17 layer","mass of cohesive sediment of class10,10n17 layer","mass of cohesive sediment of class11,11n17 layer","mass of cohesive sediment of class12,12n17 layer","mass of cohesive sediment of class13,13n17 layer","mass of cohesive sediment of class14,14n17 layer","mass of cohesive sediment of class15,15n17 layer","mass of cohesive sediment of class16,16n17 layer","mass of cohesive sediment of class17,17n17 layer","mass of cohesive sediment of class18,18n17 layer","mass of cohesive sediment of class19,19n17 layer","mass of cohesive sediment of class**,**n17 layer","mass of cohesive sediment of class*,*n18 layer","mass of cohesive sediment of class1,1n18 layer","mass of cohesive sediment of class2,2n18 layer","mass of cohesive sediment of class3,3n18 layer","mass of cohesive sediment of class4,4n18 layer","mass of cohesive sediment of class5,5n18 layer","mass of cohesive sediment of class6,6n18 layer","mass of cohesive sediment of class7,7n18 layer","mass of cohesive sediment of class8,8n18 layer","mass of cohesive sediment of class9,9n18 layer","mass of cohesive sediment of class10,10n18 layer","mass of cohesive sediment of class11,11n18 layer","mass of cohesive sediment of class12,12n18 layer","mass of cohesive sediment of class13,13n18 layer","mass of cohesive sediment of class14,14n18 layer","mass of cohesive sediment of class15,15n18 layer","mass of cohesive sediment of class16,16n18 layer","mass of cohesive sediment of class17,17n18 layer","mass of cohesive sediment of class18,18n18 layer","mass of cohesive sediment of class19,19n18 layer","mass of cohesive sediment of class**,**n18 layer","mass of cohesive sediment of class*,*n19 layer","mass of cohesive sediment of class1,1n19 layer","mass of cohesive sediment of class2,2n19 layer","mass of cohesive sediment of class3,3n19 layer","mass of cohesive sediment of class4,4n19 layer","mass of cohesive sediment of class5,5n19 layer","mass of cohesive sediment of class6,6n19 layer","mass of cohesive sediment of class7,7n19 layer","mass of cohesive sediment of class8,8n19 layer","mass of cohesive sediment of class9,9n19 layer","mass of cohesive sediment of class10,10n19 layer","mass of cohesive sediment of class11,11n19 layer","mass of cohesive sediment of class12,12n19 layer","mass of cohesive sediment of class13,13n19 layer","mass of cohesive sediment of class14,14n19 layer","mass of cohesive sediment of class15,15n19 layer","mass of cohesive sediment of class16,16n19 layer","mass of cohesive sediment of class17,17n19 layer","mass of cohesive sediment of class18,18n19 layer","mass of cohesive sediment of class19,19n19 layer","mass of cohesive sediment of class**,**n19 layer","mass of cohesive sediment of class*,*n** layer","mass of cohesive sediment of class1,1n** layer","mass of cohesive sediment of class2,2n** layer","mass of cohesive sediment of class3,3n** layer","mass of cohesive sediment of class4,4n** layer","mass of cohesive sediment of class5,5n** layer","mass of cohesive sediment of class6,6n** layer","mass of cohesive sediment of class7,7n** layer","mass of cohesive sediment of class8,8n** layer","mass of cohesive sediment of class9,9n** layer","mass of cohesive sediment of class10,10n** layer","mass of cohesive sediment of class11,11n** layer","mass of cohesive sediment of class12,12n** layer","mass of cohesive sediment of class13,13n** layer","mass of cohesive sediment of class14,14n** layer","mass of cohesive sediment of class15,15n** layer","mass of cohesive sediment of class16,16n** layer","mass of cohesive sediment of class17,17n** layer","mass of cohesive sediment of class18,18n** layer","mass of cohesive sediment of class19,19n** layer","mass of cohesive sediment of class**,**n** layer","reference level for Nestor"],
        defaut = ["velocity along x axis (m/s)","velocity along y axis (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","non erodable bottom","bottom evolution (m)"],
        fr = """Noms des variables que l''utilisateur veut ecrire dans
le fichier des resultats.
Chaque variable est representee par une lettre.
Le choix des separateurs est libre. Voir CHOIX ci-dessus.
 On peut utiliser *, *A* signifie : toutes les fractions""",
        ang = """Names of variables the user wants to write
into the graphic results file.
Each variable is represented by a letter. See CHOIX1 above.
 One can use *, *A* means all fractions""",
    ),
)
# -----------------------------------------------------------------------
TIME = PROC(nom= "TIME",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NUMBER_OF_SUB_ITERATIONS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """permet de realiser des sous iterations au sein de la
boucle en temps (En mode non permananet).
Peut etre utile en non permanent lorsque le
pas de temps qui est donne par le pas de sortie graphique du
FICHIER DE CALCUL PRECEDENT est trop grand.""",
        ang = """enable to realize sub-iteration inside a time step
(this key word is not used if the key word VARIABLE TIME-STEP
is set equal to yes). It could be useful for a non steady case
be useful for a non steady case when the time step which is fixed
by the graphic printout period of the HYDRODYNAMIC FILE
is too large.""",
    ),
#   -----------------------------------
    ORIGINAL_HOUR_OF_TIME = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I', min= 3, max= 3,
        defaut = [0,0,0],
        fr = """Permet de fixer l''heure d''origine des temps du modele lors
de la prise en compte de la force generatrice de la maree.""",
        ang = """Give the time of the time origin of the model when taking into
account of the tide generator force.""",
    ),
#   -----------------------------------
    ORIGINAL_DATE_OF_TIME = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I', min= 3, max= 3,
        defaut = [0,0,0],
        fr = """Permet de fixer la date d''origine des temps du modele lors
de la prise en compte de la force generatrice de la maree.""",
        ang = """Give the date of the time origin of the model when taking into
account the tide generating force.""",
    ),
)
# -----------------------------------------------------------------------
NON_COHESIVE = PROC(nom= "NON_COHESIVE",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    LAYERS_NON_COHESIVE_BED_POROSITY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        defaut = [0.4,0.4],
        fr = """la concentration volumique  du lit est definie par
CSF= (1-porosite)
Ce parametre est utilise pour les sediments non-cohesifs.""",
        ang = """The bed volume concentration CSF=(1-porosity) is used to
calculate the bed evolution of non-cohesive sand transport.""",
    ),
#   -----------------------------------
    CHARRIAGE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        BOUNDARY_CONDITIONS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            CLASSES_IMPOSED_SOLID_DISCHARGES_DISTRIBUTION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Donne la proportion du debit solide total impose par classe.
Donner un nombre par classe de non-cohesifs""",
                ang = """Gives the proportion of the imposed solid discharge for each class.
Give one numbre for each non-cohesive class""",
            ),
        ),
    ),
#   -----------------------------------
    BEDLOAD = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        MORPHOLOGICAL_FACTOR_ON_TIME_SCALE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """Coefficient d''amplification de l''echelle des temps""",
            ang = """amplification coefficient of time scale""",
        ),
#       -----------------------------------
        MORPHOLOGICAL_FACTOR_ON_BED_EVOLUTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """Coefficient d''amplification de l''evolution du lit""",
            ang = """amplification coefficient of bed evolution""",
        ),
#       -----------------------------------
        MINIMUM_DEPTH_FOR_BEDLOAD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-2,
            fr = """Suppression des flux de sediment de et vers les points secs""",
            ang = """To cancel sediment fluxes to and from dry points""",
        ),
#       -----------------------------------
        BED_LOAD_FOR_ALL_SANDS = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """calcul avec charriage""",
            ang = """""",
        ),
#       -----------------------------------
        BED_LOAD_TRANSPORT_FORMULA_FOR_ALL_SANDS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """10 formules de transport solide sont implementees dans GAIA.
Les formules Ne3, Ne30 et Ne9 ne doivent pas etre utilisees en cas de
couplage avec la suspension.
Les formules Ne4, Ne5, Ne8 and Ne9 calculent le taux de transport sous
l''action combinee de la houle et du courant :
                 1 : MEYER-PETER (charriage)
                 2 : EINSTEIN-BROWN (charriage)
                 3 : ENGELUND-HANSEN + CHOLLET ET CUNGE (VERSION 5.3)
                 30: ENGELUND-HANSEN (total)
                 4 : BIJKER (charriage + suspension)
                 5 : SOULSBY - VAN RIJN (charriage + suspension)
                 6 : HUNZIKER (uniquement granulometrie etendue)
                      DE MASQUAGE DE HUNZIKER APPLIQUEE
                     et mot-cle HIDING-FACTOR not used
                 7 : VAN RIJN (bed load)
                 8 : BAILARD (charriage + suspension)
                 9 : DIBAJNIA ET WATANABE (total)
            		10 : WILCOCK ET CROWE (granulometrie etendue)
L''utilisateur a aussi la possibilite d''imposer une autre formule de
transport solide (sous-programme bedload\_qb\_user.f) en fixant la
valeur du mot cle a 0 :
                 0 :  IMPOSEE PAR L''UTILISATEUR
Attention : dans ce cas, il n''est pas possible de choisir l''option
PAS DE TEMPS VARIABLE.""",
            ang = """10 bed-load or total load transport formulas are implemented in
GAIA.
The formula Ne3, Ne30 and Ne9 should not be used in the case of coupling
with the suspension.
The formula Ne4, Ne5, Ne8 and Ne9  model the transport under the
combined action of currents and waves :
                 1 : MEYER-PETER (bed load)
                 2 : EINSTEIN-BROWN (bed load)
                 3 : ENGELUND-HANSEN + CHOLLET AND CUNGE (VERSION 5.3)
                 30: ENGELUND-HANSEN (total)
                 4 : BIJKER (bed load + suspension)
                 5 : SOULSBY - VAN RIJN (bed load + suspension)
                 6 : HUNZIKER (only for sand grading)
                     IN THIS CASE HIDING FACTOR KEYWORD DISCARDED
                     And Hunziker formula used
                 7 : VAN RIJN (bed load)
                 8  : BAILARD (bed load + suspension)
                 9 : DIBAJNIA ET WATANABE (total)
            		10 : WILCOCK AND CROWE (graded sediment)
Users can also program other formulas (subroutine bedload\_qb\_user.f)
setting this key word to zero :
                 0 : FORMULA PROGRAMMED BY USER
Warning : it is not then possible to choose the option
VARIABLE TIME-STEP""",
        ),
#       -----------------------------------
        B_VALUE_FOR_THE_BIJKER_FORMULA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 2.E0,
            fr = """Coefficient b de la formule de Bijker""",
            ang = """b value for the Bijker formula""",
        ),
#       -----------------------------------
        MPM_COEFFICIENT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 8.0E-00,
            fr = """""",
            ang = """""",
        ),
#       -----------------------------------
        BOUNDARY_CONDITIONS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            PRESCRIBED_SOLID_DISCHARGES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min= 2, max= 2,
                fr = """ Valeurs des debits solides imposes aux frontieres
liquides entrantes (kg/s).
Une valeur par frontiere liquide""",
                ang = """Values of prescribed solid discharges
at the inflow boundaries (kg/s).
One value per liquid boundary""",
            ),
        ),
#       -----------------------------------
        SLOPE_INFLUENCE = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            FORMULA_FOR_SLOPE_EFFECT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """1 : formule de Koch et Flokstra, modifie le transport solide
            mot cle associe : BETA
        2 : formule de Soulsby, modifie la contrainte seuil, ne peut
            donc etre utilisee que avec une formule a seuil.
            mot cle associe : ANGLE DE REPOS DU SEDIMENT""",
                ang = """1 : formula of Koch et Flokstra, modification of bed load
             linked keyword : BETA
         2 : formula of Soulsby, modification critical shear stress,
             can only be used with a threshold fomula
             linked keyword : FRICTION ANGLE OF THE SEDIMENT""",
            ),
#           -----------------------------------
            FRICTION_ANGLE_OF_THE_SEDIMENT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 40.,
                fr = """Angle de repos du sediment, intervient pour la prise de compte
 de la pente sur la contrainte critique par la formule de Soulsby.
Utiliser si ...=2""",
                ang = """Angle of repose of the sediment. Used in the Soulsby formula to
take into account the influence of bed slope on critical shear stress.
Use if ...=2""",
            ),
#           -----------------------------------
            FORMULA_FOR_DEVIATION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """1 : Koch et Flokstra
	2 : formule de Talmon et al. 1995, JHR 33(4) formules (1) et
(17) mot cle associe : BETA2""",
                ang = """1: Koch and Flokstra
	 2: formula of Talmon et al. 1995, JHR 33(4) formulas (1) and
(17) linked keyword : BETA2""",
            ),
#           -----------------------------------
            PARAMETER_FOR_DEVIATION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.85,
                fr = """Parametre pour la deviation causee par effet de pente
 pour la formule de Talmon et al.
Une valeur elevee provoque une faible deviation""",
                ang = """Parameter pour la deviation pour la formule de Talmon et al.""",
            ),
#           -----------------------------------
            SLOPE_EFFECT = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [True ],
                fr = """Prise en compte de l''effet de pente :
deviation et modification du seuil critique.
NON supprime les mots-cles
 POUR EFFET DE PENTE et  POUR LA DEVIATION""",
                ang = """If yes, slope effect taken into account:
deviation + modification of critical shear stress.
NO will cancel the key-words
FORMULA FOR SLOPE EFFECT and FORMULA FOR DEVIATION""",
            ),
#           -----------------------------------
            BETA = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 1.3,
                fr = """Determine la valeur du coefficient beta qui intervient dans la
formulation de l''effet de pente de Koch et Flokstra.""",
                ang = """Specifies the value of the beta coefficient used in the Koch
and Flokstra slope effect formulation.""",
            ),
        ),
    ),
#   -----------------------------------
    FRICTION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        RATIO_BETWEEN_SKIN_FRICTION_AND_MEAN_DIAMETER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 3.0,
            fr = """Ratio pour le calcul du frottement de peau.
rugosite de peau = ratio * diametre moyen.
( pour la granulometrie etendue, le diametre moyen utilise
est une valeur moyenne par noeud calculee a partir de la fraction
et du diametre moyen de chaque sediment en chaque noeud du maillage)""",
            ang = """ Ratio for the computation of skin friction.
skin roughness = ratio * mean diameter
(for the mixture of sand, the mean diameter used is a value per node
which is computed thanks to the fraction and the mean diameter of each
sediment for each node of the mesh)
if KSPRATIO =0 : use skin friction prediction from Van Rijn (2007)
for currents and the Wiberg and Harris method for waves""",
        ),
#       -----------------------------------
        SKIN_FRICTION_CORRECTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Correction du frottement pour le transport solide (voir aussi KSPRATIO)
pour faire en sorte que le frottement de peau soit utilisé dans
les formules au lieu du frottement issu de l''hydrodynamique,
qui comprend souvent d''autres effets (surtout en 2D)
       0 : Pas de correction (TAUP= TOB), valide si la rugosite
           donnee a TELEMAC est physique - proche de la rugosite de peau
           estimee a 3d50
       1 : Correction pour fond plat (KSP= KSPRATIO * D50)
       2 : Prise en compte des rides""",
            ang = """Formula to modify the shear stress in sediment flow rate formulae
so they use the skin bed roughness (see also KSPRATIO)
       0 : No correction (TAUP= TOB), valid if the roughness provided
           to TELEMAC is physical - close to the skin roughness, usually
           estimated to 3d50
       1 : Correction for a flat bed (KSP= KSPRATIO * D50)
       2 : Ripple correction factor""",
        ),
#       -----------------------------------
        ADVANCED = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            BED_ROUGHNESS_PREDICTOR_OPTION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """1 : Fond plat ks=KSPRATIO D50,
2: Fond ride (methode de Wiberg et Harris),
3: Dunes et megarides (methode de Van Rijn)""",
                ang = """1: Flat bed, 2: Rippled bed,
3: Dunes and mega ripples (Method of Van Rijn)""",
            ),
#           -----------------------------------
            COMPUTE_BED_ROUGHNESS_AT_SEDIMENT_SCALE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Calcul de la rugosite de Nikuradse a l''echelle du sediment
- voir OPTION DU PREDICTEUR DE RUGOSITE -. Elle peut etre
differente de la rugosite utilisee pour l''hydrodynamique
dans les cas avec vagues ou avec une granulometrie variable en
temps et en espace. Le frottement
modifie pourra eventuellement etre envoye a Telemac mais ce n''est
pas encore supporte.""",
                ang = """Compute a bed roughness at the sediment scale - see
BED ROUGHNESS PREDICTOR OPTION -. It can be different from
the hydrodynamics roughness due to the action of waves or to
a space-time varying grain size distribution. This roughness
could be sent to Telemac but it is not supported yet.""",
            ),
        ),
    ),
#   -----------------------------------
    SUSPENSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SETTLING_LAG = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Uses the velocity profile based on the Miles approximation""",
            ang = """Uses the velocity profile based on the Miles approximation""",
        ),
#       -----------------------------------
        SUSPENSION_FOR_ALL_SANDS = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Active la suspension pour tous les sables. Il n''est pas possible
d''avoir un comportement different entre classe de sable. La vase
est toujours consideree en suspension.""",
            ang = """Activate suspension for all the sands in the simulation. It is not
possible to have a different behaviour between sand classes. Mud
is always considered in suspension.""",
        ),
#       -----------------------------------
        SUSPENSION_TRANSPORT_FORMULA_FOR_ALL_SANDS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Different formules sont proposees pour le calcul de la concentration
d''equilibre:
1 : formule de Zyserman et Fredsoe
2 : methode de Bijker. La concentration au fond
   est reliee au taux de transport par charriage
3 : formule de Van Rijn (1987)
4 : formule de Soulsy\_van Rijn""",
            ang = """Different choice to compute the equilibrium near-bed concentration
1 : Zysderman and Fredsoe, equilibrium formula
2: Bijker method. The near bed concentration
is related to the bedload . This option cannot be used
without bedload transport
3: Van Rijn formula
4: Soulsby\_van Rijn formula""",
        ),
    ),
#   -----------------------------------
    BED_MATERIAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ACTIVE_LAYER_THICKNESS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [10000.],
            fr = """Epaisseur de reference pour la stratification du lit. La
composition de la premiere couche sert a calculer le transport
solide. Utiliser une tres grande valeur pour ne pas avoir de
stratification.""",
            ang = """Thickness for bed stratification. Composition of first
layer is used to compute bed load transport rate. If you do not want
a stratification, use a large value""",
        ),
#       -----------------------------------
        HIDING_FACTOR_FORMULA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [0],
            fr = """4 formules pour le hiding factor sont programmees dans GAIA
     0: const => il faut donner le HIDING FACTOR PAR CLASSE
     1: Egiazaroff
     2: Ashida \& Michiue
      :
     4: Karim, Holly \& Yang""",
            ang = """4 hiding factor formulas are implemented in GAIA
     0: const => need to give CLASSES HIDING FACTOR
     1: Egiazaroff
     2: Ashida \& Michiue
      :
     4: Karim, Holly \& Yang""",
        ),
#       -----------------------------------
        CONSTANT_ACTIVE_LAYER_THICKNESS = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """epaisseur de couche active constante ou non""",
            ang = """constant active layer thickness or not""",
        ),
#       -----------------------------------
        D90_SAND_DIAMETER_FOR_ONLY_ONE_CLASS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = .01,
            fr = """Fixe la valeur de D90 pour des simulations avec une
 seule classe de sable. Pour des simulations avec plusieurs
 classes de sable, D90 est calcule par GAIA.""",
            ang = """Sets the value of diameter d90 for simulations with one sand class.
 With multiple sand classes, D90 is computed by GAIA.""",
        ),
    ),
#   -----------------------------------
    BED_STRUCTURE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ADVANCED = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SEDIMENT_SLIDE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Prise en compte de la pente d''equilibre du sable donnee par le
mot-cle ANGLE DE FROTTEMENT DU SEDIMENT""",
                ang = """If yes, the key-word FRICTION ANGLE OF THE SEDIMENT is taken
into account for slope stability""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
INITIAL_CONDITION = PROC(nom= "INITIAL_CONDITION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    SETTING = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUSPENSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            INITIAL_SUSPENDED_SEDIMENTS_CONCENTRATION_VALUES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                defaut = [0.,0.],
                fr = """Fixe les valeurs initiales de concentration des sediments en
suspension.""",
                ang = """Sets the initial values of suspended sediment concentration.""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
BOUNDARY_CONDITIONS = PROC(nom= "BOUNDARY_CONDITIONS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    SETTING = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUSPENSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            PRESCRIBED_SUSPENDED_SEDIMENTS_CONCENTRATION_VALUES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', max='**',
                fr = """Valeurs de concentration des sediment en suspension imposees aux
frontieres liquides entrantes. Les valeurs sont imposees a la premiere
frontiere, puis la deuxieme etc. suivant la meme logique que pour les
traceurs.""",
                ang = """Suspended sediment concentration values prescribed at the inflow
boundaries. Determines the imposed value of sediments at the first
boundary, then at the second and so on, with the same logic as tracers.""",
            ),
#           -----------------------------------
            VERTICAL_PROFILES_OF_SUSPENDED_SEDIMENTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["User defined","Constant","Constant or Rouse if sediment", "Normalized Rouse profile and imposed conc", "Modified Rouse profile accounting for molecular viscosity"],
                fr = """Permet de specifier le type de profil de concentration des
sediments sur la verticale (uniquement pour des simulations 3D).
Les choix possibles sont :
\begin{itemize}
\item 0 : Programmation utilisateur ;
\item 1 : Constant ;
\item 2 : Rouse equilibrium concentration ;
\item 3 : Rouse (normalise) et concentration imposee.
\item 4 : Rouse modifie avec viscosite moleculaire.
\end{itemize}""",
                ang = """Specifies the type of profiles of sediment concentration on the
vertical (only for 3D simulations). Possible choices are:
\begin{itemize}
\item 0: user defined,
\item 1: constant,
\item 2: Rouse equilibrium, constant (diluted sediment)
or Rouse (sediment),
\item 3: Rouse (normalized) and imposed concentration.
\item 4: Rouse modified with molecular viscosity.
\end{itemize}""",
            ),
#           -----------------------------------
            SUSPENDED_SEDIMENTS_CONCENTRATION_VALUES_AT_THE_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Valeurs des sediments en suspension a chacune des sources.
toutes les sources pour le premier sediment
puis toutes les sources du deuxieme sediment, etc.
Par exemple, s''il y a 3 sediments en suspension (SED1, SED2 et SED3)
et 2 sources (S1 et S2), la syntaxe suivante est utilisee :\\
S1\_SED1;S1\_SED2;S1\_SED3;S2\_SED1;S2\_SED2;S2\_SED3\\
10.0; 10.0; 0.0;  0.0; 10.0; 10.0""",
                ang = """Values of the suspended sediments at the sources.
All sources for the first suspended sediment, then
all sources for the second suspended sediment, etc.
For example, if there are 3 suspended sediments (SED1, SED2 and SED3)
and 2 sources (S1 and S2), the following syntax is used:\\
S1\_SED1;S1\_SED2;S1\_SED3;S2\_SED1;S2\_SED2;S2\_SED3\\
10.0; 10.0; 0.0;  0.0; 10.0; 10.0""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
NUMERICAL_PARAMETERS = PROC(nom= "NUMERICAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    MAXIMUM_NUMBER_OF_ITERATIONS_FOR_POSITIVE_THICKNESS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [50],
        fr = """Seulement pour le charriage elements finis.
Parametres de positive\_depths.
Ancienne valeur par defaut = 10 jusqu en version 8.1.""",
        ang = """Only for bedload in finite elements.
Parameter for positive\_depths.
Old default value = 10 until release 8.1.""",
    ),
#   -----------------------------------
    BEDLOAD = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        UPWINDING_FOR_BEDLOAD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.5E0,
            fr = """Parametre pour la résolution VF de equation d Exner,
0.5 = Centre (precis), 1 = Decentrage amont (stable car diffusif)""",
            ang = """Parameter for FV solving the Exner equation,
0.5 = Centered (precise), 1 = Upwind (stable because diffusive)""",
        ),
    ),
#   -----------------------------------
    ADVECTION_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUSPENSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_SUSPENDED_SEDIMENTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["NO ADVECTION","CHARACTERISTICS","EXPLICIT + SUPG","EXPLICIT LEO POSTMA","EXPLICIT + MURD SCHEME N","EXPLICIT + MURD SCHEME PSI","N-SCHEME FOR TIDAL FLATS LP","N-SCHEME FOR TIDAL FLATS","ERIA SCHEME - ONLY IN 2D"],
                defaut = ["EXPLICIT + MURD SCHEME PSI"],
                fr = """Choix du schema de convection pour les sediments en suspension,
ERIA fonctionne uniquement en 2D. A donner en suivant l''ordre
des sediments en suspension""",
                ang = """Choice of the advection scheme for the suspended sediments,
ERIA works only in 3D. The order of the chosen scheme must follow
the order of the suspended sediments.""",
            ),
#           -----------------------------------
            SCHEME_OPTION_FOR_ADVECTION_OF_SUSPENDED_SEDIMENTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I', min=0, max='**',
                defaut = [4,4],
                fr = """Si schema PSI ou N : 1=explicite 2=predicteur-correcteur
3=predicteur-correcteur deuxieme ordre en temps
4=implicite""",
                ang = """If N or PSI SCHEME: 1=explicit 2=predictor-corrector
3= predictor-corrector second-order in time 4= implicit""",
            ),
        ),
    ),
#   -----------------------------------
    SOLVER = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUSPENSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_DIFFUSION_OF_SUSPENDED_SEDIMENTS_IN_3D = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO DIFFUSION","IMPLICIT","VERTICAL DIFFUSION ONLY"],
                defaut = "IMPLICIT",
                fr = """Permet de specifier si l''on utilise ou non la diffusion
des sediments pour des simulations 3d.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de diffusion,
\item 1 : implicite;
\item 2 : diffusion verticale seulement.
\end{itemize}""",
                ang = """Monitors the choice of the diffusion scheme
for sediments in 3D simulations.
Possible choices are:
\begin{itemize}
\item 0: no diffusion,
\item 1: implicit,
\item 2: vertical diffusion only.
\end{itemize}""",
            ),
#           -----------------------------------
            SOLVER_FOR_DIFFUSION_OF_SUSPENSION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres","direct solver"],
                defaut = ["conjugate gradient"],
                fr = """Permet de choisir le solveur utilise pour la resolution de
la suspension.
Les choix possibles sont :
\begin{itemize}
\item 1 : gradient conjugue ;
\item 2 : residu conjugue ;
\item 3 : gradient conjugue sur equation normale ;
\item 4 : erreur minimale ;
\item 5 : gradient conjugue carre ;
\item 6 : CGSTAB ;
\item 7 : GMRES ;
\item 8 : solveur direct.
\end{itemize}""",
                ang = """Choice of the solver for suspension resolution.
Possible choices are:
\begin{itemize}
\item 1: conjugate gradient,
\item 2: conjugate residual,
\item 3: conjugate gradient on a normal equation,
\item 4: minimum error,
\item 5: squared conjugate gradient,
\item 6: CGSTAB,
\item 7: GMRES,
\item 8: direct solver.
\end{itemize}""",
            ),
#           -----------------------------------
            SOLVER_OPTION_FOR_DIFFUSION_OF_SUSPENSION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 5,
                fr = """Dimension de l''espace de Krylov pour la methode GMRES (7).""",
                ang = """Dimension of Krylov space for the GMRES method (7).""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_ITERATIONS_FOR_SOLVER_FOR_SUSPENSION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 60,
                fr = """Limite le nombre d''iterations du solveur a chaque pas
de temps pour le calcul de la diffusion de la suspension.""",
                ang = """Limits the number of solver iterations for the diffusion of
sediments.""",
            ),
#           -----------------------------------
            ACCURACY_FOR_DIFFUSION_OF_SUSPENSION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [1.E-8],
                fr = """Fixe la precision demandee pour le calcul de la diffusion
de la suspension. Une seule valeur est acceptee pour tous les
solveurs.""",
                ang = """Sets the accuracy needed for the computation of
the diffusion of suspension. It is not possible to set different
values for different solvers, only one is accepted.""",
            ),
#           -----------------------------------
            PRECONDITIONING_FOR_DIFFUSION_OF_SUSPENSION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["no preconditioning","diagonal","diagonal condensed","diagonal with absolute values","Crout","Gauss-Seidel EBE","Matrix defined by the user","diagonal and Crout","direct solver on the vertical","diagonal condensed and Crout","diagonal and direct solver on the vertical"],
                defaut = "diagonal",
                fr = """Permet de preconditionner le systeme relatif
a la diffusion de la suspension. Les choix possibles sont :
\begin{itemize}
\item 0 : aucun ;
\item 2 : diagonal ;
\item 3 : diagonal avec matrice condensee en 3D;
\item 5 : diagonal avec valeurs absolues en 3D ;
\item 7 : Crout ;
\item 11 : Gauss-Seidel EBE en 3D;
\item 13 : matrice fournie par l''utilisateur en 3D ;
\item 14 : diagonal et Crout ;
\item 17 : solveur direct sur la verticale en 3D;
\item 21 : diagonal condensee et Crout en 3D;
\item 34 : diagonal et solveur direct sur la verticale en 3D.
\end{itemize}""",
                ang = """Choice of preconditioning for the diffusion of sediments.
Possible choices are:
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal,
\item 3: diagonal with the condensed matrix in 3D,
\item 5: diagonal with absolute values in 3D,
\item 7: Crout,
\item 11: Gauss-Seidel EBE in 3D,
\item 13: matrix defined by the user in 3D,
\item 14: diagonal and Crout,
\item 17: direct solver on the vertical in 3D,
\item 21: diagonal condensed and Crout in 3D,
\item 34: diagonal and direct solver on the vertical in 3D.
\end{itemize}""",
            ),
        ),
    ),
#   -----------------------------------
    AUTOMATIC_DIFFERENTIATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        AD_NUMBER_OF_DERIVATIVES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [0],
            fr = """Definit le nombre de derivees utilisateurs, dans le cadre
de la differentiation algorithmique.""",
            ang = """Defines the number of user derivatives, within the framework
of the algorithmic differentiation.""",
        ),
#       -----------------------------------
        AD_NAMES_OF_DERIVATIVES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min= 2, max= 2,
            fr = """Noms des derivees utilisateurs en 32 caracteres,
         16 pour le nom, 16 pour l''unite""",
            ang = """Name of user derivatives in 32 characters,
         16 for the name, 16 for the unit.""",
        ),
#       -----------------------------------
        AD_NUMBER_OF_DIRECTIONS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """Definit le nombre de directions de differentiateurs.""",
            ang = """Defines the number of directions for the differentiators""",
        ),
#       -----------------------------------
        AD_SYMBOLIC_LINEAR_SOLVER = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Permet le solveur lineaire symbolique pour l AD.""",
            ang = """Enables the symbolic linear solver for AD.""",
        ),
#       -----------------------------------
        AD_LINEAR_SOLVER_RESET_DERIVATIVES = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Remet a zero les derivees pour l AD.""",
            ang = """Resets the derivatives for AD.""",
        ),
#       -----------------------------------
        AD_LINEAR_SOLVER_DERIVATIVE_CONVERGENCE = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Solveur lineaire iteratif : test de convergence des derivees
pour l AD.""",
            ang = """Iterative linear solvers: derivative convergence test for AD.""",
        ),
    ),
)
# -----------------------------------------------------------------------
SEDIMENT_INFO = PROC(nom= "SEDIMENT_INFO",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    SETTLING_VELOCITY = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ADVECTION_DIFFUSION_SCHEME_WITH_SETTLING_VELOCITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Choix de schema vertical pour la diffusion et le depot du
sediment (uniquement pour le 3D):
\begin{itemize}
\item 0 : Diffusion implicite ;
\item 1 : Schema implicite de convection-diffusion
(tridiagonal matrix solver) ;
\item 2 : Convection faible  \telfile{sed\_fall.f}
\end{itemize}""",
            ang = """Choice of the vertical scheme for diffusion and settling of
sediment (only in 3D):
\begin{itemize}
\item 0: Implicit-diffusion scheme,
\item 1: Implicit-convection scheme (Tridiagonal matrix solver),
\item 2: \telfile{set\_fall.f}
\end{itemize}""",
        ),
    ),
)
# -----------------------------------------------------------------------
COHESIVE = PROC(nom= "COHESIVE",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    INITIALIZATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUSPENSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            LAYERS_CRITICAL_EROSION_SHEAR_STRESS_OF_THE_MUD = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                defaut = [0.5,1.],
                fr = """Taux critique d erosion de la vase needs to be defined
        for each layer (N par m2)""",
                ang = """Critical erosion shear stress of the mud per layer (N per m2)""",
            ),
        ),
    ),
#   -----------------------------------
    SETTLING_VELOCITY = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        HINDERED_SETTLING = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide si la formulation entravee doit etre utilisee
pour calculer la vitesse de chute de la vase.""",
            ang = """Decides if hindered formulation is to be used to
compute settling velocity for mud.""",
        ),
#       -----------------------------------
        HINDERED_SETTLING_FORMULA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Type de vitesse de chute entravee :
\begin{itemize}
\item 1 : Whitehouse et al. (2000) - fonctionne ;
\item 2 : Winterwerp (1999) - ne fonctionne pas actuellement.
\end{itemize}""",
            ang = """Type of hindered settling:
\begin{itemize}
\item 1: Whitehouse et al. (2000) - working,
\item 2: Winterwerp (1999) - not currently working.
\end{itemize}""",
        ),
#       -----------------------------------
        WEAK_SOIL_CONCENTRATION_FOR_MUD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.0,
            fr = """Concentration limite en kg/m$^3$ au-dela de laquelle
la couche de vase fluide devient solide.
Cette valeur est demandee lorsque
\telkey{VITESSE DE CHUTE ENTRAVEE} = OUI.""",
            ang = """The sediment concentration at which sediment
forms a weak soil in kg/m$^3$. These values are needed when
\telkey{HINDERED SETTLING} = YES.""",
        ),
#       -----------------------------------
        THRESHOLD_CONCENTRATION_FOR_HINDERED_SETTLING = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.0,
            fr = """Concentration de sediment a laquelle la chute entravee est
initiee. Ces valeurs sont necessaires lorsque
\telkey{VITESSE DE CHUTE ENTRAVEE} = OUI.""",
            ang = """The sediment concentration at which hindered settling is
initiated. These values are needed when
\telkey{HINDERED SETTLING} = YES.""",
        ),
#       -----------------------------------
        FLOCCULATION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide si la formulation entravee doit etre utilisee
pour calculer la vitesse de chute pour la vase.""",
            ang = """Decides if hindered formulation is to be used to
compute settling velocity for mud.""",
        ),
#       -----------------------------------
        FLOCCULATION_FORMULA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Formule pour floculation :
\begin{itemize}
\item 1: Van Leussen ;
\item 2: Soulsby et  al. (2013).
\end{itemize}""",
            ang = """Type of flocculation formula:
\begin{itemize}
\item 1: Van Leussen,
\item 2: Soulsby et  al. (2013).
\end{itemize}""",
        ),
#       -----------------------------------
        FLOCCULATION_COEFFICIENT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.3,
            fr = """Coefficient intervenant dans la modelisation de l''influence de
la turbulence sur la floculation, il intervient plus precisement dans
le terme de formation des flocs par les contraintes turbulentes
(coefficient $a$ de la formule de Van Leussen).
Valeur a imposer si
\telkey{INFLUENCE DE LA TURBULENCE SUR LA VITESSE DE CHUTE}
= OUI.""",
            ang = """When the influence of turbulence on the settling velocity
is modelled, this coefficient traduces the formation of flocs by
turbulence (coefficient $a$ of Van Leussen formula).
Value to be imposed if
\telkey{INFLUENCE OF TURBULENCE ON SETTLING VELOCITY} = YES.""",
        ),
#       -----------------------------------
        COEFFICIENT_RELATIVE_TO_FLOC_DESTRUCTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.09,
            fr = """Coefficient intervenant dans la modelisation de l''influence de
la turbulence sur la floculation, il intervient plus precisement dans
le terme de destruction des flocs par les contraintes turbulentes
(coefficient $b$ de la formulede Van Leussen).
Valeur a imposer si
\telkey{INFLUENCE DE LA TURBULENCE SUR LA VITESSE DE CHUTE}
= OUI.""",
            ang = """When the influence of turbulence on the settling velocity
is modelled, this coefficient traduces the breaking of flocs by
turbulence (coefficient $b$ of Van Leussen formula).
Value to be imposed if
\telkey{INFLUENCE OF TURBULENCE ON SETTLING VELOCITY} = YES.""",
        ),
    ),
#   -----------------------------------
    SUSPENSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        CLASSES_CRITICAL_SHEAR_STRESS_FOR_MUD_DEPOSITION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [1000.,1000.],
            fr = """Contrainte critique de depot (Pa)""",
            ang = """Critical shear stress for deposition (Pa)""",
        ),
#       -----------------------------------
        LAYERS_PARTHENIADES_CONSTANT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [1.E-3,1.E-3],
            fr = """constante de la loi d''erosion de la vase (Kg/m2/s)""",
            ang = """constant of the Krone and Partheniades erosion law (Kg/m2/s)""",
        ),
    ),
#   -----------------------------------
    BED_STRUCTURE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        INITIALIZATION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            LAYERS_MUD_CONCENTRATION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Concentration du lit de vase en g/ l - defini par couches""",
                ang = """Concentrations of the mud-bed in g per l (per layer)""",
            ),
        ),
    ),
#   -----------------------------------
    CONSOLIDATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_LAYERS_OF_THE_CONSOLIDATION_MODEL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """Structure verticale du lit cohesif- le nombre de couche doit
       etre inferieur a 10""",
            ang = """Vertical bed structure - The number of layers should be less
       than 10""",
        ),
#       -----------------------------------
        LAYERS_MASS_TRANSFER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Coefficients de transfert de masse du modele
de tassement multicouche en s-1""",
            ang = """Mass transfert coefficients of
the multilayer consolidation model in s-1""",
        ),
    ),
)
# -----------------------------------------------------------------------
NUMERICAL = PROC(nom= "NUMERICAL",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    ZERO = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-10,
        fr = """Fixe le zero de GAIA utilise pour les clippings.""",
        ang = """Sets the zero of GAIA used for clipping values.""",
    ),
#   -----------------------------------
    FINITE_VOLUMES = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Methodes volumes finis ou pas""",
        ang = """Set finite volumes method or not""",
    ),
#   -----------------------------------
    MATRIX_VECTOR_PRODUCT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """""",
        ang = """""",
    ),
#   -----------------------------------
    MATRIX_STORAGE = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """""",
        ang = """""",
    ),
)
# -----------------------------------------------------------------------
PHYSICS = PROC(nom= "PHYSICS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    CLASSES_SEDIMENT_DENSITY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        defaut = [2650.,2650.],
        fr = """Fixe la valeur de la masse volumique du sediment par classe en Kg/m3""",
        ang = """Sets the value of the sediment density for each class en Kg/m3""",
    ),
#   -----------------------------------
    WATER_VISCOSITY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-6,
        fr = """Definit la viscosite cinematique de l''eau.
         M/S2""",
        ang = """Specifies the water kinematic viscosity.
         M/S2""",
    ),
)
# -----------------------------------------------------------------------
SUSPENSION = PROC(nom= "SUSPENSION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    CORRECTION_ON_CONVECTION_VELOCITY = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Modification du champ convecteur 2D pour prise en compte du
gradient vertical de vitesse et concentration""",
        ang = """Modification of 2D convection velocities  to account for
velocity and concentration profiles""",
    ),
#   -----------------------------------
    THETA_IMPLICITATION_FOR_SUSPENSION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [1.],
        fr = """facteur implicitation du terme de depot et de la diffusion.
 si teta =0, tout le terme de depot est traite de maniere explicite.
 Valide seulement pour le modèle 2D.""",
        ang = """ implicitation factor for the deposition flux and the diffusion.
 for teta =0, the deposition flux is only explicit.
 Only valid for the 2D model.""",
    ),
#   -----------------------------------
    TURBULENCE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        COEFFICIENT_FOR_DIFFUSION_OF_SUSPENDED_SEDIMENTS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [1.E-6],
            fr = """Fixe la valeur du coefficient de diffusion des sediments en
suspension en 2D.
L''influence de ce parametre sur l''evolution des sediments dans
le temps est importante. C'' est un scalaire (une seule valeur pour
tous les sediments.""",
            ang = """Sets the value of the suspended sediments diffusivity in 2D.
This value may have a significant effect on the evolution of
sediments in time. It is a scalar (one value for all sediments).""",
        ),
#       -----------------------------------
        COEFFICIENT_FOR_HORIZONTAL_DIFFUSION_OF_SUSPENDED_SEDIMENTS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Fixe les valeurs de coefficients de diffusion horizontal des
sediments, en 3D.  L''influence de ce parametre sur l''evolution des
sediments dans le temps est importante.
C'' est un tableau avec une valeur par sediment en suspension,
 separation par un point virgule.""",
            ang = """Sets the values of the horizontal diffusion of sediments in 3D.
These values may have a significant effect on the evolution of
sediments in time.
It is an array, with one value per suspended sediment, separated by
 semicolons.""",
        ),
    ),
)
# -----------------------------------------------------------------------
SUSPENSIONS = PROC(nom= "SUSPENSIONS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    TURBULENCE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        COEFFICIENT_FOR_VERTICAL_DIFFUSION_OF_SUSPENDED_SEDIMENTS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Fixe les valeurs de coefficients de diffusion vertical des
sediments, en 3D.  L''influence de ce parametre sur l''evolution des
 sediments dans le temps est importante.
C'' est un tableau avec une valeur par sediemnts en suspension,
separation par un point virgule.""",
            ang = """Sets the values of the vertical diffusion of sediments in 3D.
These values may have a significant effect on the evolution of
sediments in time.
It is an array, with one value per suspended sediment, separated by
 semicolons.""",
        ),
    ),
)
# -----------------------------------------------------------------------
SEDIMENTOLOGY = PROC(nom= "SEDIMENTOLOGY",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    GENERAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SECONDARY_CURRENTS_ALPHA_COEFFICIENT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.0E-00,
            fr = """(-)""",
            ang = """Alpha coefficient of secondary current(-),
Should be chosen between 0.75 (rough bottom) and 1 (smooth bottom)""",
        ),
    ),
)
# -----------------------------------------------------------------------
BED_MATERIAL = PROC(nom= "BED_MATERIAL",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    C_VSM = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        VERTICAL_GRAIN_SORTING_MODEL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [0],
            fr = """(-)""",
            ang = """Defines the model of the vertical grain sorting:
        0 = HR-VSM = Layer Model (Classic Hirano / Ribberink approach)
        1 = C-VSM (Continous Vertical Grain Sorting Model)""",
        ),
#       -----------------------------------
        C_VSM_MAXIMUM_SECTIONS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [200],
            fr = """(-)""",
            ang = """Defines the maximum discretisation of the
         Continous Vertical Sorting Model:
         Should be bigger than 8xNumber of Fractions.
         The bigger the higher the RAM requirements,
         but the faster and accurater the
         bookkeeping of the sediments.""",
        ),
#       -----------------------------------
        C_VSM_FULL_PRINTOUT_PERIOD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [0],
            fr = """(-)""",
            ang = """Number of Timesteps to next printout
         of the full C-VSM. These printouts are highly
         time and disc consuming.
         0 = Coupled to GRAPHIC PRINTOUT PERIOD
         >0 = Own printout period for the C-VSM""",
        ),
#       -----------------------------------
        C_VSM_PRINTOUT_SELECTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            defaut = '0;0;0;0;0;0;0;0;0;0;0;0;0',
            fr = """(-)""",
            ang = """Printout the C-VSM for the whole model as 3D
         \telkey{C-VSM RESULTS FILE}
         or / and for some nodes as
         \telkey{C-VSM NODE FILE}
         Give Up to 100 INTEGER numbers separated by ";"
         0 = Full model .-> VSPRES
         N = 1,2...NPOINT; 2D-ID of a SELFIN MESH POINT ->*\_VSP.CSV""",
        ),
#       -----------------------------------
        ACTIVE_LAYER_THICKNESS_FORMULA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [0],
            fr = """(-)""",
            ang = """MODEL FOR ACTIVE LAYER THICKNESS
         0 = ELAY0 (Keyword: ACTIVE LAYER THICKNESS)
         1 = Hunziker \& G$\ddot{u}$nther
         2 = Fredsoe \& Deigaard (1992)
         3 = van RIJN (1993)
         4 = Wong (2006)
         5 = Malcherek (2003)
         6 = $3*d50$ within last time steps ALT""",
        ),
    ),
)
TEXTE_NEW_JDC = "\
COMPUTATION_ENVIRONMENT();\
GENERAL_PARAMETERS();\
HYDRODYNAMICS();\
MISCELLANEOUS();\
GENERAL();\
COHESIVE_AND_NON_COHESIVE();\
BED_STRUCTURE();\
USELESS();\
INPUT_OUTPUT__FILES();\
RESULTS();\
DATA_FILES();\
INITIAL_CONDITIONS();\
INPUT_OUTPUT__GRAPHICS_AND_LISTING();\
TIME();\
NON_COHESIVE();\
INITIAL_CONDITION();\
BOUNDARY_CONDITIONS();\
NUMERICAL_PARAMETERS();\
SEDIMENT_INFO();\
COHESIVE();\
NUMERICAL();\
PHYSICS();\
SUSPENSION();\
SUSPENSIONS();\
SEDIMENTOLOGY();\
BED_MATERIAL();\
"
Ordre_Des_Commandes = (
'COMPUTATION_ENVIRONMENT',
'INTERNAL',
'GENERAL_PARAMETERS',
'HYDRODYNAMICS',
'MISCELLANEOUS',
'GENERAL',
'COHESIVE_AND_NON_COHESIVE',
'BED_STRUCTURE',
'USELESS',
'INPUT_OUTPUT__FILES',
'RESULTS',
'DATA_FILES',
'INITIAL_CONDITIONS',
'INPUT_OUTPUT__GRAPHICS_AND_LISTING',
'TIME',
'NON_COHESIVE',
'INITIAL_CONDITION',
'BOUNDARY_CONDITIONS',
'NUMERICAL_PARAMETERS',
'SEDIMENT_INFO',
'COHESIVE',
'NUMERICAL',
'PHYSICS',
'SUSPENSION',
'SUSPENSIONS',
'SEDIMENTOLOGY',
'BED_MATERIAL')
try:
    import TelApy
    source = "eficas"
except Exception as excpt:
    source = "Telemac"
enum = source+'.gaia_enum_auto'
dicoCasEn = source+'.gaia_dicoCasEnToCata'
dicoCasFr = source+'.gaia_dicoCasFrToCata'
