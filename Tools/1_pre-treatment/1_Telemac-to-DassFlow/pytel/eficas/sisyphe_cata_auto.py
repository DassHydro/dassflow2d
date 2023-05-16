
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



JdC = JDC_CATA (code = 'SISYPHE',
                execmodul = None,
                )
# =======================================================================
# Catalog entry for the MAP function : c_pre_interfaceBody_mesh
# =======================================================================

VERSION_CATALOGUE="V8P2"
# -----------------------------------------------------------------------
INPUT_OUTPUT__INFORMATION = PROC(nom= "INPUT_OUTPUT__INFORMATION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    COMPUTATIONAL_INFORMATION = FACT(statut='o',
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
    ),
)
# -----------------------------------------------------------------------
MISCELLANEOUS = PROC(nom= "MISCELLANEOUS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    PARALLEL_PROCESSORS = SIMP(statut ='f',
#   -----------------------------------
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
#   -----------------------------------
    VALIDATION = SIMP(statut ='f',
#   -----------------------------------
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
#   -----------------------------------
    b_VALIDATIONG = BLOC(condition="VALIDATION == True",
#   -----------------------------------
#       -----------------------------------
        REFERENCE_FILE_FORMAT = SIMP(statut ='f',
#       -----------------------------------
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
#       -----------------------------------
        REFERENCE_FILE = SIMP(statut ='f',
#       -----------------------------------
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
#   -----------------------------------
    NUMBER_OF_PRIVATE_ARRAYS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """Nombre de tableaux mis a disposition de l utilisateur""",
        ang = """Number of arrays for own user programming""",
    ),
#   -----------------------------------
    ORIGIN_COORDINATES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I', min= 2, max= 2,
        defaut = [0,0],
        fr = """Valeur en metres, utilise pour eviter les trop grands nombres,
transmis dans le format Selafin mais pas d''autre traitement pour
l''instant""",
        ang = """Value in metres, used to avoid large real numbers,
added in Selafin format, but so far no other treatment""",
    ),
#   -----------------------------------
    DEBUGGER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [0],
        fr = """Pour imprimer la sequence des appels, mettre 1""",
        ang = """If 1, calls of subroutines will be printed in the listing""",
    ),
#   -----------------------------------
    OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATS = SIMP(statut ='f',
#   -----------------------------------
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
    CHECKING_THE_MESH = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Si oui on appelle le sous-programme checkmesh qui verifie
la coherence du maillage, points superposes, etc.""",
        ang = """if this key word is equal to yes, a call to subroutine
checkmesh will look for errors in the mesh, superimposed points, etc.""",
    ),
#   -----------------------------------
    MAXIMUM_NUMBER_OF_BOUNDARIES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 30,
        fr = """nombre maximal de frontieres differentes dans le maillage.
Sert au dimensionnement de la memoire, a augmenter si necessaire""",
        ang = """maximal number of boundaries in the mesh.
Used for dimensioning arrays. Can be increased if needed""",
    ),
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
    STATIONARY_MODE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """(-)""",
        ang = """(-)""",
    ),
#   -----------------------------------
    SHIELDS_PARAMETERS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        fr = """Utilise pour determiner la valeur de la contrainte critique
d''entrainement (sediments non-cohesifs). En multiclasse, specifier la
valeur pour chaque classe,sinon on prend la meme valeur.  Par defaut (si
aucune valeur donnee) le code calcule lui meme le parametre de Shields
en fonction du diametre.""",
        ang = """Used to determine the critical bed shear stress value
(non-cohesive sediments).
For multi grain size,
the shields parameter needs to be specified for each class. If only one
value is specified, the shields parameter will be considered constant.
The default option (no shields given in parameter file)
is to calculate the shields parameter as a function of
sand grain diameter (see logical CALAC).""",
    ),
#   -----------------------------------
    OPTION_FOR_THE_TREATMENT_OF_NON_ERODABLE_BEDS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [0],
        fr = """Ce parametre permet de fixer la methode retenue pour traiter
les fonds non erodables :
        0 = FONDS ERODABLES PARTOUT
        1 = MINIMISATION DU DEBIT SOLIDE POUR LES FONDS NON ERODABLES
        2 = DEBIT SOLIDE NUL POUR LES FONDS NON ERODABLES
        3 = MINIMISATION DU DEBIT SOLIDE EN ELEMENTS FINIS/MASS-LUMPING
        4 = MINIMISATION DU DEBIT SOLIDE EN VOLUMES FINIS
Quand le fond rigide peut etre atteint lors de la simulation, il est
conseille d utiliser la metode 3 ou bien 4""",
        ang = """This parameters determines the method used to treat the
non erodable bottoms :
        0 = ERODABLE BOTTOMS EVERYWHERE
        1 = MINIMISATION OF THE SOLID DISCHARGE
        2 = NUL SOLID DISCHARGE
        3 = MINIMISATION OF THE SOLID DISCHARGE IN FE / MASS-LUMPING
        4 = MINIMISATION OF THE SOLID DISCHARGE IN FINITE VOLUMES
When the rigid bed can be reached during the computation, it is advised
to use the method 3 or the method 4""",
    ),
#   -----------------------------------
    MINIMAL_VALUE_OF_THE_WATER_HEIGHT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-3,
        fr = """Fixe la valeur minimale de la hauteur d''eau.
Est utilise lorsque le mot cle BANCS DECOUVRANTS est egal a oui.""",
        ang = """Sets the minimum value of the water depth.
Is used when the keyword TIDAL FLATS is equal to yes.""",
    ),
#   -----------------------------------
    TIDAL_FLATS = SIMP(statut ='f',
#   -----------------------------------
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
    GRAIN_FEEDING = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Now suppressed""",
        ang = """Now suppressed""",
    ),
#   -----------------------------------
    STEADY_CASE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Definit si l''hydrodynamique a prendre en compte est
permanente ou non. Si la valeur OUI est donnee a ce mot-cle,
le dernier enregistrement du fichier du calcul precedent
constituera le champ (h,u,v et eventuellement hauteur et periode
de houle) a prendre en compte.""",
        ang = """Specifies steady or unsteady case.If this keyword
is equal to YES, the last record of the previous computation
file will give the values of h,u,v and eventually wave height
and period to be considered.""",
    ),
#   -----------------------------------
    CONSTANT_FLOW_DISCHARGE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """constant flow discharge or not""",
        ang = """constant flow discharge or not""",
    ),
#   -----------------------------------
    NUMBER_OF_ITERATIONS_FOR_TELEMAC = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 500,
        fr = """Nombre d''iteration a effecuter avec telemac pour obtenir
un nouvel ecoulement quasi stationnaire.
A utiliser avec l''option constant flow discharge""",
        ang = """Number of iteration to do wtih telemac in order to obtain a
new quasi-stationary flow.
To use with the option constant flow discharge""",
    ),
#   -----------------------------------
    CRITERION_TO_UPDATE_THE_FLOW = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.1,
        fr = """Critere (Hauteur du fond>CRIT\_CFD*Hauteur d''eau) pour mettre
a jour l''ecoulement.
A utiliser avec l''option constant flow discharge""",
        ang = """Criterion (Bottom height>CRIT\_CFD*Water depth) in order to
update the flow.
To use with the option constant flow discharge""",
    ),
#   -----------------------------------
    CRITICAL_EVOLUTION_RATIO = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.,
        fr = """Definit le moment ou l''extrapolation du champ de courant
par SISYPHE n''est plus valable. Ce nombre represente le rapport
maximum entre les evolutions et la hauteur d''eau admissible.
Generalement, on admet qu''une evolution inferieure a 0,1 fois
la hauteur d''eau ne modifie pas sensiblement la repartition
du champ de courant.""",
        ang = """Specifies the moment when the SISYPHE extrapolation current
filed is no more valid. This value set the maximum ratio between
evolutions and the water depth.
Generally, it is considered that an evolution lower than 0,1 time the
water depth does not perceptibly modify the current field distribution.""",
    ),
#   -----------------------------------
    NUMBER_OF_BED_LOAD_MODEL_LAYERS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 2,
        fr = """Comme son nom l''indique..., defaut NOMBLAY=2""",
        ang = """This is the given allocation limit, secure default NOMLAY=2""",
    ),
#   -----------------------------------
    MASS_CONCENTRATION = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """par defaut, le code calcule les concentrations volumiques.
Ce mot cle permet dimprimer et de relire des concentrations en g/l.
Les concentrations imposees aux limites
(fichier condim, concentration de condim\_susp)
sont alors aussi donnees par l utilisateur en g/l et
reconverties en  concentration volumique par le code.""",
        ang = """Determines if concentrations (input and output)
are  mass concentrations in g/l or adimensionnal volume concentrations
(default option).""",
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
    MIXED_SEDIMENT = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Melange sable vase: 2 classes seulement""",
        ang = """Mixture of cohesive and non cohesive sediment : 2 class only""",
    ),
)
# -----------------------------------------------------------------------
USELESS = PROC(nom= "USELESS",op = None,
# -----------------------------------------------------------------------
    UIinfo = {"groupes": ("CACHE")},
#   -----------------------------------
    VECTOR_LENGTH = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """longueur du vecteur pour les machines vectorielles.""",
        ang = """vector length on vector machines.""",
    ),
#   -----------------------------------
    STEERING_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier contenant les parametres du calcul
a realiser. Il peut-etre ecrit par l''utilisateur avec EDAMOX.""",
        ang = """Name of the file containing the parameters
of the computation. Could be written by the user with EDAMOX.""",
    ),
#   -----------------------------------
    MEAN_DIAMETER_OF_THE_SEDIMENT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        defaut = [.01,.01,.01,.01,.01,.01,.01,.01,.01,.01],
        fr = """Sets value of diameter dm for particular size class.""",
        ang = """Sets value of diameter dm for particular size class.""",
    ),
#   -----------------------------------
    GEOMETRY_FILE_BINARY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        defaut = 'STD',
        fr = """Type de binaire pour l''ecriture du fichier de geometrie.
Ce type depend de la machine ou ce fichier est genere.
Les differentes valeurs sont :
   IBM, pour un fichier sur IBM (depuis un CRAY)
   I3E, pour un fichier sur HP (depuis un CRAY)
   STD, pour le type de binaire de la machine ou l''utilisateur
        travaille. Les commandes standard READ et WRITE sont
        alors utilisees.""",
        ang = """Binary file type used for writing the geometry file.
This type depends on the machine on which the file was generated.
The possible values are as follows :
   IBM, for a file on an IBM (from a CRAY)
   I3E, for a file on an HP (from a CRAY)
   STD, binary type of the machine on which the user is working.
        The normal READ and WRITE commands are then used.""",
    ),
#   -----------------------------------
    HYDRODYNAMIC_FILE_BINARY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        defaut = 'STD',
        fr = """obsolete""",
        ang = """obsolete""",
    ),
#   -----------------------------------
    BINARY_OF_THE_PREVIOUS_SEDIMENTOLOGICAL_COMPUTATION_FILE = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        defaut = 'STD',
        fr = """Type de binaire pour l''ecriture du fichier precedent
sedimentologique. Ce type depend de la machine ou le fichier
precedent sedimentologique est genere. Les differentes valeurs
sont identiques a celles du fichier de geometrie.""",
        ang = """Binary file type used for writing the previous
sedimentological computation results file.
This type depends on the machine on which the file
was generated. The possible values are the same as for
the geometry file.""",
    ),
#   -----------------------------------
    RESULTS_FILE_BINARY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        defaut = 'STD',
        fr = """Type de binaire pour l''ecriture du fichier des resultats.
Ce type depend de la machine ou le fichier des resultats est genere.
Les differentes valeurs sont identiques a celles du fichier de
geometrie.""",
        ang = """Binary file type used for writing the results file.
This type depends on the machine on which the file was generated.
The possible values are the same as for the geometry file.""",
    ),
#   -----------------------------------
    REFERENCE_FILE_BINARY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        defaut = 'STD',
        fr = """Type de binaire pour l''ecriture du fichier de reference.
Ce type depend de la machine ou le fichier de reference est genere.
Les differentes valeurs sont identiques a celles du fichier
de geometrie.""",
        ang = """Binary file type used for writing the reference file.
This type depends on the machine on which the file was generated.
The possible values are the same as for the geometry file.""",
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
        HYDRODYNAMIC_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'), max='**',
            defaut = '',
            fr = """Nom d''un fichier contenant les resultats d''un calcul
precedent realise sur le meme maillage. L''hydrodynamique
du calcul sera donnee, soit par le dernier pas de temps
de ce fichier si le calcul est permanent,
soit par les pas de temps representant la maree
ou la crue consideree si le calcul est instationnaire
Remarque : Si l''on modelise le transport sous l''action
combinee du courant et de la houle (mot-cle  DE
TRANSPORT SOLIDE =4) ce fichier doit contenir non
seulement les donnees hydrodynamiques(hauteur d''eau,vitesses)
mais aussi les donnees de houle (hauteur de houle, periode de
houle).Les donnees de houle peuvent toutefois etre imposees
par l''utilisateur dans le sous programme CONDIM.""",
            ang = """Name of a file containing the results a previous
computation  made on the same mesh. The hydrodynamic will
be given by the last record of the file if the case is
steady or, if the case is unsteady, by the time steps
describing the tide or flood.
Remark :If the bed-load transport under the combined
action of currents and wave is modelled (keyword
BED-LOAD TRANSPORT FORMULA set equal to 4), this file
must contain not only the hydrodynamic data
(water height, velocities) but also the wave data
(wave height, wave period).However, the user has also
 the possibility to give the values
of the wave data in the subroutine CONDIM.""",
        ),
#       -----------------------------------
        NAMES_OF_PRIVATE_VARIABLES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
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
        LIQUID_BOUNDARIES_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'), max='**',
            defaut = '',
            fr = """Fichier de variations en temps des conditions aux limites.
Les donnees de ce fichier sont sur le canal SIS\_FILES(SISLIQ)%LU.""",
            ang = """Variations in time of boundary conditions. Data of this file
are read on channel SIS\_FILES(SISLIQ)%LU.""",
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
une periodicite donnee par le mot cle PERIODE DE SORTIE GRAPHIQUE.""",
        ang = """Name of the file into wich the computation results shall be
written, the periodicity being given by the keyword
GRAPHIC PRINTOUT PERIOD.""",
    ),
#   -----------------------------------
    VARIABLES_TO_BE_PRINTED = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM', max='**',
        into = ["velocity along x axis (m/s)","velocity along y axis (m/s)","wawe celerity (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","Froude number","scalar flowrate of fluid (m2/s)","flowrate along x axis (m2/s)","flowrate along y axis (m2/s)","bed-load discharge (m2/s)","bed-load discharge along x axis (m2/s)","bed-load discharge along y axis (m2/s)","bottom evolution (m)","non erodable bottom","bed friction coefficient (m if Nikuradse)","mean bottom friction (N/m2)","wave angle with axis Oy (deg)","wave height","wave period","fraction of sediment of class**n layer*","fraction of sediment of class11n layer*","fraction of sediment of class22n layer*","fraction of sediment of class33n layer*","fraction of sediment of class44n layer*","fraction of sediment of class55n layer*","fraction of sediment of class66n layer*","fraction of sediment of class77n layer*","fraction of sediment of class88n layer*","fraction of sediment of class99n layer*","fraction of sediment of class1010n layer*","fraction of sediment of class1111n layer*","fraction of sediment of class1212n layer*","fraction of sediment of class1313n layer*","fraction of sediment of class1414n layer*","fraction of sediment of class1515n layer*","fraction of sediment of class1616n layer*","fraction of sediment of class1717n layer*","fraction of sediment of class1818n layer*","fraction of sediment of class1919n layer*","fraction of sediment of class****n layer*","fraction of sediment of class**n layer1","fraction of sediment of class11n layer1","fraction of sediment of class22n layer1","fraction of sediment of class33n layer1","fraction of sediment of class44n layer1","fraction of sediment of class55n layer1","fraction of sediment of class66n layer1","fraction of sediment of class77n layer1","fraction of sediment of class88n layer1","fraction of sediment of class99n layer1","fraction of sediment of class1010n layer1","fraction of sediment of class1111n layer1","fraction of sediment of class1212n layer1","fraction of sediment of class1313n layer1","fraction of sediment of class1414n layer1","fraction of sediment of class1515n layer1","fraction of sediment of class1616n layer1","fraction of sediment of class1717n layer1","fraction of sediment of class1818n layer1","fraction of sediment of class1919n layer1","fraction of sediment of class****n layer1","fraction of sediment of class**n layer2","fraction of sediment of class11n layer2","fraction of sediment of class22n layer2","fraction of sediment of class33n layer2","fraction of sediment of class44n layer2","fraction of sediment of class55n layer2","fraction of sediment of class66n layer2","fraction of sediment of class77n layer2","fraction of sediment of class88n layer2","fraction of sediment of class99n layer2","fraction of sediment of class1010n layer2","fraction of sediment of class1111n layer2","fraction of sediment of class1212n layer2","fraction of sediment of class1313n layer2","fraction of sediment of class1414n layer2","fraction of sediment of class1515n layer2","fraction of sediment of class1616n layer2","fraction of sediment of class1717n layer2","fraction of sediment of class1818n layer2","fraction of sediment of class1919n layer2","fraction of sediment of class****n layer2","fraction of sediment of class**n layer3","fraction of sediment of class11n layer3","fraction of sediment of class22n layer3","fraction of sediment of class33n layer3","fraction of sediment of class44n layer3","fraction of sediment of class55n layer3","fraction of sediment of class66n layer3","fraction of sediment of class77n layer3","fraction of sediment of class88n layer3","fraction of sediment of class99n layer3","fraction of sediment of class1010n layer3","fraction of sediment of class1111n layer3","fraction of sediment of class1212n layer3","fraction of sediment of class1313n layer3","fraction of sediment of class1414n layer3","fraction of sediment of class1515n layer3","fraction of sediment of class1616n layer3","fraction of sediment of class1717n layer3","fraction of sediment of class1818n layer3","fraction of sediment of class1919n layer3","fraction of sediment of class****n layer3","fraction of sediment of class**n layer4","fraction of sediment of class11n layer4","fraction of sediment of class22n layer4","fraction of sediment of class33n layer4","fraction of sediment of class44n layer4","fraction of sediment of class55n layer4","fraction of sediment of class66n layer4","fraction of sediment of class77n layer4","fraction of sediment of class88n layer4","fraction of sediment of class99n layer4","fraction of sediment of class1010n layer4","fraction of sediment of class1111n layer4","fraction of sediment of class1212n layer4","fraction of sediment of class1313n layer4","fraction of sediment of class1414n layer4","fraction of sediment of class1515n layer4","fraction of sediment of class1616n layer4","fraction of sediment of class1717n layer4","fraction of sediment of class1818n layer4","fraction of sediment of class1919n layer4","fraction of sediment of class****n layer4","fraction of sediment of class**n layer5","fraction of sediment of class11n layer5","fraction of sediment of class22n layer5","fraction of sediment of class33n layer5","fraction of sediment of class44n layer5","fraction of sediment of class55n layer5","fraction of sediment of class66n layer5","fraction of sediment of class77n layer5","fraction of sediment of class88n layer5","fraction of sediment of class99n layer5","fraction of sediment of class1010n layer5","fraction of sediment of class1111n layer5","fraction of sediment of class1212n layer5","fraction of sediment of class1313n layer5","fraction of sediment of class1414n layer5","fraction of sediment of class1515n layer5","fraction of sediment of class1616n layer5","fraction of sediment of class1717n layer5","fraction of sediment of class1818n layer5","fraction of sediment of class1919n layer5","fraction of sediment of class****n layer5","fraction of sediment of class**n layer6","fraction of sediment of class11n layer6","fraction of sediment of class22n layer6","fraction of sediment of class33n layer6","fraction of sediment of class44n layer6","fraction of sediment of class55n layer6","fraction of sediment of class66n layer6","fraction of sediment of class77n layer6","fraction of sediment of class88n layer6","fraction of sediment of class99n layer6","fraction of sediment of class1010n layer6","fraction of sediment of class1111n layer6","fraction of sediment of class1212n layer6","fraction of sediment of class1313n layer6","fraction of sediment of class1414n layer6","fraction of sediment of class1515n layer6","fraction of sediment of class1616n layer6","fraction of sediment of class1717n layer6","fraction of sediment of class1818n layer6","fraction of sediment of class1919n layer6","fraction of sediment of class****n layer6","fraction of sediment of class**n layer7","fraction of sediment of class11n layer7","fraction of sediment of class22n layer7","fraction of sediment of class33n layer7","fraction of sediment of class44n layer7","fraction of sediment of class55n layer7","fraction of sediment of class66n layer7","fraction of sediment of class77n layer7","fraction of sediment of class88n layer7","fraction of sediment of class99n layer7","fraction of sediment of class1010n layer7","fraction of sediment of class1111n layer7","fraction of sediment of class1212n layer7","fraction of sediment of class1313n layer7","fraction of sediment of class1414n layer7","fraction of sediment of class1515n layer7","fraction of sediment of class1616n layer7","fraction of sediment of class1717n layer7","fraction of sediment of class1818n layer7","fraction of sediment of class1919n layer7","fraction of sediment of class****n layer7","fraction of sediment of class**n layer8","fraction of sediment of class11n layer8","fraction of sediment of class22n layer8","fraction of sediment of class33n layer8","fraction of sediment of class44n layer8","fraction of sediment of class55n layer8","fraction of sediment of class66n layer8","fraction of sediment of class77n layer8","fraction of sediment of class88n layer8","fraction of sediment of class99n layer8","fraction of sediment of class1010n layer8","fraction of sediment of class1111n layer8","fraction of sediment of class1212n layer8","fraction of sediment of class1313n layer8","fraction of sediment of class1414n layer8","fraction of sediment of class1515n layer8","fraction of sediment of class1616n layer8","fraction of sediment of class1717n layer8","fraction of sediment of class1818n layer8","fraction of sediment of class1919n layer8","fraction of sediment of class****n layer8","fraction of sediment of class**n layer9","fraction of sediment of class11n layer9","fraction of sediment of class22n layer9","fraction of sediment of class33n layer9","fraction of sediment of class44n layer9","fraction of sediment of class55n layer9","fraction of sediment of class66n layer9","fraction of sediment of class77n layer9","fraction of sediment of class88n layer9","fraction of sediment of class99n layer9","fraction of sediment of class1010n layer9","fraction of sediment of class1111n layer9","fraction of sediment of class1212n layer9","fraction of sediment of class1313n layer9","fraction of sediment of class1414n layer9","fraction of sediment of class1515n layer9","fraction of sediment of class1616n layer9","fraction of sediment of class1717n layer9","fraction of sediment of class1818n layer9","fraction of sediment of class1919n layer9","fraction of sediment of class****n layer9","fraction of sediment of class**n layer10","fraction of sediment of class11n layer10","fraction of sediment of class22n layer10","fraction of sediment of class33n layer10","fraction of sediment of class44n layer10","fraction of sediment of class55n layer10","fraction of sediment of class66n layer10","fraction of sediment of class77n layer10","fraction of sediment of class88n layer10","fraction of sediment of class99n layer10","fraction of sediment of class1010n layer10","fraction of sediment of class1111n layer10","fraction of sediment of class1212n layer10","fraction of sediment of class1313n layer10","fraction of sediment of class1414n layer10","fraction of sediment of class1515n layer10","fraction of sediment of class1616n layer10","fraction of sediment of class1717n layer10","fraction of sediment of class1818n layer10","fraction of sediment of class1919n layer10","fraction of sediment of class****n layer10","fraction of sediment of class**n layer11","fraction of sediment of class11n layer11","fraction of sediment of class22n layer11","fraction of sediment of class33n layer11","fraction of sediment of class44n layer11","fraction of sediment of class55n layer11","fraction of sediment of class66n layer11","fraction of sediment of class77n layer11","fraction of sediment of class88n layer11","fraction of sediment of class99n layer11","fraction of sediment of class1010n layer11","fraction of sediment of class1111n layer11","fraction of sediment of class1212n layer11","fraction of sediment of class1313n layer11","fraction of sediment of class1414n layer11","fraction of sediment of class1515n layer11","fraction of sediment of class1616n layer11","fraction of sediment of class1717n layer11","fraction of sediment of class1818n layer11","fraction of sediment of class1919n layer11","fraction of sediment of class****n layer11","fraction of sediment of class**n layer12","fraction of sediment of class11n layer12","fraction of sediment of class22n layer12","fraction of sediment of class33n layer12","fraction of sediment of class44n layer12","fraction of sediment of class55n layer12","fraction of sediment of class66n layer12","fraction of sediment of class77n layer12","fraction of sediment of class88n layer12","fraction of sediment of class99n layer12","fraction of sediment of class1010n layer12","fraction of sediment of class1111n layer12","fraction of sediment of class1212n layer12","fraction of sediment of class1313n layer12","fraction of sediment of class1414n layer12","fraction of sediment of class1515n layer12","fraction of sediment of class1616n layer12","fraction of sediment of class1717n layer12","fraction of sediment of class1818n layer12","fraction of sediment of class1919n layer12","fraction of sediment of class****n layer12","fraction of sediment of class**n layer13","fraction of sediment of class11n layer13","fraction of sediment of class22n layer13","fraction of sediment of class33n layer13","fraction of sediment of class44n layer13","fraction of sediment of class55n layer13","fraction of sediment of class66n layer13","fraction of sediment of class77n layer13","fraction of sediment of class88n layer13","fraction of sediment of class99n layer13","fraction of sediment of class1010n layer13","fraction of sediment of class1111n layer13","fraction of sediment of class1212n layer13","fraction of sediment of class1313n layer13","fraction of sediment of class1414n layer13","fraction of sediment of class1515n layer13","fraction of sediment of class1616n layer13","fraction of sediment of class1717n layer13","fraction of sediment of class1818n layer13","fraction of sediment of class1919n layer13","fraction of sediment of class****n layer13","fraction of sediment of class**n layer14","fraction of sediment of class11n layer14","fraction of sediment of class22n layer14","fraction of sediment of class33n layer14","fraction of sediment of class44n layer14","fraction of sediment of class55n layer14","fraction of sediment of class66n layer14","fraction of sediment of class77n layer14","fraction of sediment of class88n layer14","fraction of sediment of class99n layer14","fraction of sediment of class1010n layer14","fraction of sediment of class1111n layer14","fraction of sediment of class1212n layer14","fraction of sediment of class1313n layer14","fraction of sediment of class1414n layer14","fraction of sediment of class1515n layer14","fraction of sediment of class1616n layer14","fraction of sediment of class1717n layer14","fraction of sediment of class1818n layer14","fraction of sediment of class1919n layer14","fraction of sediment of class****n layer14","fraction of sediment of class**n layer15","fraction of sediment of class11n layer15","fraction of sediment of class22n layer15","fraction of sediment of class33n layer15","fraction of sediment of class44n layer15","fraction of sediment of class55n layer15","fraction of sediment of class66n layer15","fraction of sediment of class77n layer15","fraction of sediment of class88n layer15","fraction of sediment of class99n layer15","fraction of sediment of class1010n layer15","fraction of sediment of class1111n layer15","fraction of sediment of class1212n layer15","fraction of sediment of class1313n layer15","fraction of sediment of class1414n layer15","fraction of sediment of class1515n layer15","fraction of sediment of class1616n layer15","fraction of sediment of class1717n layer15","fraction of sediment of class1818n layer15","fraction of sediment of class1919n layer15","fraction of sediment of class****n layer15","fraction of sediment of class**n layer16","fraction of sediment of class11n layer16","fraction of sediment of class22n layer16","fraction of sediment of class33n layer16","fraction of sediment of class44n layer16","fraction of sediment of class55n layer16","fraction of sediment of class66n layer16","fraction of sediment of class77n layer16","fraction of sediment of class88n layer16","fraction of sediment of class99n layer16","fraction of sediment of class1010n layer16","fraction of sediment of class1111n layer16","fraction of sediment of class1212n layer16","fraction of sediment of class1313n layer16","fraction of sediment of class1414n layer16","fraction of sediment of class1515n layer16","fraction of sediment of class1616n layer16","fraction of sediment of class1717n layer16","fraction of sediment of class1818n layer16","fraction of sediment of class1919n layer16","fraction of sediment of class****n layer16","fraction of sediment of class**n layer17","fraction of sediment of class11n layer17","fraction of sediment of class22n layer17","fraction of sediment of class33n layer17","fraction of sediment of class44n layer17","fraction of sediment of class55n layer17","fraction of sediment of class66n layer17","fraction of sediment of class77n layer17","fraction of sediment of class88n layer17","fraction of sediment of class99n layer17","fraction of sediment of class1010n layer17","fraction of sediment of class1111n layer17","fraction of sediment of class1212n layer17","fraction of sediment of class1313n layer17","fraction of sediment of class1414n layer17","fraction of sediment of class1515n layer17","fraction of sediment of class1616n layer17","fraction of sediment of class1717n layer17","fraction of sediment of class1818n layer17","fraction of sediment of class1919n layer17","fraction of sediment of class****n layer17","fraction of sediment of class**n layer18","fraction of sediment of class11n layer18","fraction of sediment of class22n layer18","fraction of sediment of class33n layer18","fraction of sediment of class44n layer18","fraction of sediment of class55n layer18","fraction of sediment of class66n layer18","fraction of sediment of class77n layer18","fraction of sediment of class88n layer18","fraction of sediment of class99n layer18","fraction of sediment of class1010n layer18","fraction of sediment of class1111n layer18","fraction of sediment of class1212n layer18","fraction of sediment of class1313n layer18","fraction of sediment of class1414n layer18","fraction of sediment of class1515n layer18","fraction of sediment of class1616n layer18","fraction of sediment of class1717n layer18","fraction of sediment of class1818n layer18","fraction of sediment of class1919n layer18","fraction of sediment of class****n layer18","fraction of sediment of class**n layer19","fraction of sediment of class11n layer19","fraction of sediment of class22n layer19","fraction of sediment of class33n layer19","fraction of sediment of class44n layer19","fraction of sediment of class55n layer19","fraction of sediment of class66n layer19","fraction of sediment of class77n layer19","fraction of sediment of class88n layer19","fraction of sediment of class99n layer19","fraction of sediment of class1010n layer19","fraction of sediment of class1111n layer19","fraction of sediment of class1212n layer19","fraction of sediment of class1313n layer19","fraction of sediment of class1414n layer19","fraction of sediment of class1515n layer19","fraction of sediment of class1616n layer19","fraction of sediment of class1717n layer19","fraction of sediment of class1818n layer19","fraction of sediment of class1919n layer19","fraction of sediment of class****n layer19","fraction of sediment of class**n layer**","fraction of sediment of class11n layer**","fraction of sediment of class22n layer**","fraction of sediment of class33n layer**","fraction of sediment of class44n layer**","fraction of sediment of class55n layer**","fraction of sediment of class66n layer**","fraction of sediment of class77n layer**","fraction of sediment of class88n layer**","fraction of sediment of class99n layer**","fraction of sediment of class1010n layer**","fraction of sediment of class1111n layer**","fraction of sediment of class1212n layer**","fraction of sediment of class1313n layer**","fraction of sediment of class1414n layer**","fraction of sediment of class1515n layer**","fraction of sediment of class1616n layer**","fraction of sediment of class1717n layer**","fraction of sediment of class1818n layer**","fraction of sediment of class1919n layer**","fraction of sediment of class****n layer**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**", "thicknes of bed layer k", "concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","saturated concentration (kg/m3)","supplementary variable A","supplementary variable G","supplementary variable L","supplementary variable O","reference level for Nestor "],
        defaut = '',
        fr = """Nom des variables que l''utilisateur desire ecrire sur
le lisring. Meme possibilites que pour les sorties graphiques.""",
        ang = """Names of variables the user wants to write on the listing.
Each variable is represented by a letter in the same manner as
it is done in the graphic results file.""",
    ),
#   -----------------------------------
    GRAPHIC_PRINTOUT_PERIOD = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Determine la periode en nombre de pas de temps d''impression
des ''VARIABLES POUR LES SORTIES GRAPHIQUES'' (voir ce mot-cle)
dans le ''FICHIER DES RESULTATS''.""",
        ang = """Determines, in number of time steps, the printout period for
the ''VARIABLES FOR GRAPHIC PRINTOUTS'' in the ''RESULTS FILE''.""",
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
    HYDRODYNAMIC_FILE_FORMAT = SIMP(statut ='f',
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
        typ = ('Fichier','All Files (*)'), max='**',
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
par l''utilisateur dans le sous programme CONDIM\_SISYPHE.
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
of the wave data in the subroutine CONDIM.
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
        into = ["velocity along x axis (m/s)","velocity along y axis (m/s)","wawe celerity (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","Froude number","scalar flowrate of fluid (m2/s)","flowrate along x axis (m2/s)","flowrate along y axis (m2/s)","bed-load discharge (m2/s)","bed-load discharge along x axis (m2/s)","bed-load discharge along y axis (m2/s)","bottom evolution (m)","non erodable bottom","total bed roughness (m)","Bed Shear stress (Totalfriction) (N/m2)","Skin friction correction factor","Mean grain diameter","wave angle with axis Oy (deg)","suspended load transport rate (m2/s)","bed load transport rate (m2/s)","wave height","wave period","wave orbital velocity (m/s)","fraction of sediment of class**n layer*","fraction of sediment of class11n layer*","fraction of sediment of class22n layer*","fraction of sediment of class33n layer*","fraction of sediment of class44n layer*","fraction of sediment of class55n layer*","fraction of sediment of class66n layer*","fraction of sediment of class77n layer*","fraction of sediment of class88n layer*","fraction of sediment of class99n layer*","fraction of sediment of class1010n layer*","fraction of sediment of class1111n layer*","fraction of sediment of class1212n layer*","fraction of sediment of class1313n layer*","fraction of sediment of class1414n layer*","fraction of sediment of class1515n layer*","fraction of sediment of class1616n layer*","fraction of sediment of class1717n layer*","fraction of sediment of class1818n layer*","fraction of sediment of class1919n layer*","fraction of sediment of class****n layer*","fraction of sediment of class**n layer1","fraction of sediment of class11n layer1","fraction of sediment of class22n layer1","fraction of sediment of class33n layer1","fraction of sediment of class44n layer1","fraction of sediment of class55n layer1","fraction of sediment of class66n layer1","fraction of sediment of class77n layer1","fraction of sediment of class88n layer1","fraction of sediment of class99n layer1","fraction of sediment of class1010n layer1","fraction of sediment of class1111n layer1","fraction of sediment of class1212n layer1","fraction of sediment of class1313n layer1","fraction of sediment of class1414n layer1","fraction of sediment of class1515n layer1","fraction of sediment of class1616n layer1","fraction of sediment of class1717n layer1","fraction of sediment of class1818n layer1","fraction of sediment of class1919n layer1","fraction of sediment of class****n layer1","fraction of sediment of class**n layer2","fraction of sediment of class11n layer2","fraction of sediment of class22n layer2","fraction of sediment of class33n layer2","fraction of sediment of class44n layer2","fraction of sediment of class55n layer2","fraction of sediment of class66n layer2","fraction of sediment of class77n layer2","fraction of sediment of class88n layer2","fraction of sediment of class99n layer2","fraction of sediment of class1010n layer2","fraction of sediment of class1111n layer2","fraction of sediment of class1212n layer2","fraction of sediment of class1313n layer2","fraction of sediment of class1414n layer2","fraction of sediment of class1515n layer2","fraction of sediment of class1616n layer2","fraction of sediment of class1717n layer2","fraction of sediment of class1818n layer2","fraction of sediment of class1919n layer2","fraction of sediment of class****n layer2","fraction of sediment of class**n layer3","fraction of sediment of class11n layer3","fraction of sediment of class22n layer3","fraction of sediment of class33n layer3","fraction of sediment of class44n layer3","fraction of sediment of class55n layer3","fraction of sediment of class66n layer3","fraction of sediment of class77n layer3","fraction of sediment of class88n layer3","fraction of sediment of class99n layer3","fraction of sediment of class1010n layer3","fraction of sediment of class1111n layer3","fraction of sediment of class1212n layer3","fraction of sediment of class1313n layer3","fraction of sediment of class1414n layer3","fraction of sediment of class1515n layer3","fraction of sediment of class1616n layer3","fraction of sediment of class1717n layer3","fraction of sediment of class1818n layer3","fraction of sediment of class1919n layer3","fraction of sediment of class****n layer3","fraction of sediment of class**n layer4","fraction of sediment of class11n layer4","fraction of sediment of class22n layer4","fraction of sediment of class33n layer4","fraction of sediment of class44n layer4","fraction of sediment of class55n layer4","fraction of sediment of class66n layer4","fraction of sediment of class77n layer4","fraction of sediment of class88n layer4","fraction of sediment of class99n layer4","fraction of sediment of class1010n layer4","fraction of sediment of class1111n layer4","fraction of sediment of class1212n layer4","fraction of sediment of class1313n layer4","fraction of sediment of class1414n layer4","fraction of sediment of class1515n layer4","fraction of sediment of class1616n layer4","fraction of sediment of class1717n layer4","fraction of sediment of class1818n layer4","fraction of sediment of class1919n layer4","fraction of sediment of class****n layer4","fraction of sediment of class**n layer5","fraction of sediment of class11n layer5","fraction of sediment of class22n layer5","fraction of sediment of class33n layer5","fraction of sediment of class44n layer5","fraction of sediment of class55n layer5","fraction of sediment of class66n layer5","fraction of sediment of class77n layer5","fraction of sediment of class88n layer5","fraction of sediment of class99n layer5","fraction of sediment of class1010n layer5","fraction of sediment of class1111n layer5","fraction of sediment of class1212n layer5","fraction of sediment of class1313n layer5","fraction of sediment of class1414n layer5","fraction of sediment of class1515n layer5","fraction of sediment of class1616n layer5","fraction of sediment of class1717n layer5","fraction of sediment of class1818n layer5","fraction of sediment of class1919n layer5","fraction of sediment of class****n layer5","fraction of sediment of class**n layer6","fraction of sediment of class11n layer6","fraction of sediment of class22n layer6","fraction of sediment of class33n layer6","fraction of sediment of class44n layer6","fraction of sediment of class55n layer6","fraction of sediment of class66n layer6","fraction of sediment of class77n layer6","fraction of sediment of class88n layer6","fraction of sediment of class99n layer6","fraction of sediment of class1010n layer6","fraction of sediment of class1111n layer6","fraction of sediment of class1212n layer6","fraction of sediment of class1313n layer6","fraction of sediment of class1414n layer6","fraction of sediment of class1515n layer6","fraction of sediment of class1616n layer6","fraction of sediment of class1717n layer6","fraction of sediment of class1818n layer6","fraction of sediment of class1919n layer6","fraction of sediment of class****n layer6","fraction of sediment of class**n layer7","fraction of sediment of class11n layer7","fraction of sediment of class22n layer7","fraction of sediment of class33n layer7","fraction of sediment of class44n layer7","fraction of sediment of class55n layer7","fraction of sediment of class66n layer7","fraction of sediment of class77n layer7","fraction of sediment of class88n layer7","fraction of sediment of class99n layer7","fraction of sediment of class1010n layer7","fraction of sediment of class1111n layer7","fraction of sediment of class1212n layer7","fraction of sediment of class1313n layer7","fraction of sediment of class1414n layer7","fraction of sediment of class1515n layer7","fraction of sediment of class1616n layer7","fraction of sediment of class1717n layer7","fraction of sediment of class1818n layer7","fraction of sediment of class1919n layer7","fraction of sediment of class****n layer7","fraction of sediment of class**n layer8","fraction of sediment of class11n layer8","fraction of sediment of class22n layer8","fraction of sediment of class33n layer8","fraction of sediment of class44n layer8","fraction of sediment of class55n layer8","fraction of sediment of class66n layer8","fraction of sediment of class77n layer8","fraction of sediment of class88n layer8","fraction of sediment of class99n layer8","fraction of sediment of class1010n layer8","fraction of sediment of class1111n layer8","fraction of sediment of class1212n layer8","fraction of sediment of class1313n layer8","fraction of sediment of class1414n layer8","fraction of sediment of class1515n layer8","fraction of sediment of class1616n layer8","fraction of sediment of class1717n layer8","fraction of sediment of class1818n layer8","fraction of sediment of class1919n layer8","fraction of sediment of class****n layer8","fraction of sediment of class**n layer9","fraction of sediment of class11n layer9","fraction of sediment of class22n layer9","fraction of sediment of class33n layer9","fraction of sediment of class44n layer9","fraction of sediment of class55n layer9","fraction of sediment of class66n layer9","fraction of sediment of class77n layer9","fraction of sediment of class88n layer9","fraction of sediment of class99n layer9","fraction of sediment of class1010n layer9","fraction of sediment of class1111n layer9","fraction of sediment of class1212n layer9","fraction of sediment of class1313n layer9","fraction of sediment of class1414n layer9","fraction of sediment of class1515n layer9","fraction of sediment of class1616n layer9","fraction of sediment of class1717n layer9","fraction of sediment of class1818n layer9","fraction of sediment of class1919n layer9","fraction of sediment of class****n layer9","fraction of sediment of class**n layer10","fraction of sediment of class11n layer10","fraction of sediment of class22n layer10","fraction of sediment of class33n layer10","fraction of sediment of class44n layer10","fraction of sediment of class55n layer10","fraction of sediment of class66n layer10","fraction of sediment of class77n layer10","fraction of sediment of class88n layer10","fraction of sediment of class99n layer10","fraction of sediment of class1010n layer10","fraction of sediment of class1111n layer10","fraction of sediment of class1212n layer10","fraction of sediment of class1313n layer10","fraction of sediment of class1414n layer10","fraction of sediment of class1515n layer10","fraction of sediment of class1616n layer10","fraction of sediment of class1717n layer10","fraction of sediment of class1818n layer10","fraction of sediment of class1919n layer10","fraction of sediment of class****n layer10","fraction of sediment of class**n layer11","fraction of sediment of class11n layer11","fraction of sediment of class22n layer11","fraction of sediment of class33n layer11","fraction of sediment of class44n layer11","fraction of sediment of class55n layer11","fraction of sediment of class66n layer11","fraction of sediment of class77n layer11","fraction of sediment of class88n layer11","fraction of sediment of class99n layer11","fraction of sediment of class1010n layer11","fraction of sediment of class1111n layer11","fraction of sediment of class1212n layer11","fraction of sediment of class1313n layer11","fraction of sediment of class1414n layer11","fraction of sediment of class1515n layer11","fraction of sediment of class1616n layer11","fraction of sediment of class1717n layer11","fraction of sediment of class1818n layer11","fraction of sediment of class1919n layer11","fraction of sediment of class****n layer11","fraction of sediment of class**n layer12","fraction of sediment of class11n layer12","fraction of sediment of class22n layer12","fraction of sediment of class33n layer12","fraction of sediment of class44n layer12","fraction of sediment of class55n layer12","fraction of sediment of class66n layer12","fraction of sediment of class77n layer12","fraction of sediment of class88n layer12","fraction of sediment of class99n layer12","fraction of sediment of class1010n layer12","fraction of sediment of class1111n layer12","fraction of sediment of class1212n layer12","fraction of sediment of class1313n layer12","fraction of sediment of class1414n layer12","fraction of sediment of class1515n layer12","fraction of sediment of class1616n layer12","fraction of sediment of class1717n layer12","fraction of sediment of class1818n layer12","fraction of sediment of class1919n layer12","fraction of sediment of class****n layer12","fraction of sediment of class**n layer13","fraction of sediment of class11n layer13","fraction of sediment of class22n layer13","fraction of sediment of class33n layer13","fraction of sediment of class44n layer13","fraction of sediment of class55n layer13","fraction of sediment of class66n layer13","fraction of sediment of class77n layer13","fraction of sediment of class88n layer13","fraction of sediment of class99n layer13","fraction of sediment of class1010n layer13","fraction of sediment of class1111n layer13","fraction of sediment of class1212n layer13","fraction of sediment of class1313n layer13","fraction of sediment of class1414n layer13","fraction of sediment of class1515n layer13","fraction of sediment of class1616n layer13","fraction of sediment of class1717n layer13","fraction of sediment of class1818n layer13","fraction of sediment of class1919n layer13","fraction of sediment of class****n layer13","fraction of sediment of class**n layer14","fraction of sediment of class11n layer14","fraction of sediment of class22n layer14","fraction of sediment of class33n layer14","fraction of sediment of class44n layer14","fraction of sediment of class55n layer14","fraction of sediment of class66n layer14","fraction of sediment of class77n layer14","fraction of sediment of class88n layer14","fraction of sediment of class99n layer14","fraction of sediment of class1010n layer14","fraction of sediment of class1111n layer14","fraction of sediment of class1212n layer14","fraction of sediment of class1313n layer14","fraction of sediment of class1414n layer14","fraction of sediment of class1515n layer14","fraction of sediment of class1616n layer14","fraction of sediment of class1717n layer14","fraction of sediment of class1818n layer14","fraction of sediment of class1919n layer14","fraction of sediment of class****n layer14","fraction of sediment of class**n layer15","fraction of sediment of class11n layer15","fraction of sediment of class22n layer15","fraction of sediment of class33n layer15","fraction of sediment of class44n layer15","fraction of sediment of class55n layer15","fraction of sediment of class66n layer15","fraction of sediment of class77n layer15","fraction of sediment of class88n layer15","fraction of sediment of class99n layer15","fraction of sediment of class1010n layer15","fraction of sediment of class1111n layer15","fraction of sediment of class1212n layer15","fraction of sediment of class1313n layer15","fraction of sediment of class1414n layer15","fraction of sediment of class1515n layer15","fraction of sediment of class1616n layer15","fraction of sediment of class1717n layer15","fraction of sediment of class1818n layer15","fraction of sediment of class1919n layer15","fraction of sediment of class****n layer15","fraction of sediment of class**n layer16","fraction of sediment of class11n layer16","fraction of sediment of class22n layer16","fraction of sediment of class33n layer16","fraction of sediment of class44n layer16","fraction of sediment of class55n layer16","fraction of sediment of class66n layer16","fraction of sediment of class77n layer16","fraction of sediment of class88n layer16","fraction of sediment of class99n layer16","fraction of sediment of class1010n layer16","fraction of sediment of class1111n layer16","fraction of sediment of class1212n layer16","fraction of sediment of class1313n layer16","fraction of sediment of class1414n layer16","fraction of sediment of class1515n layer16","fraction of sediment of class1616n layer16","fraction of sediment of class1717n layer16","fraction of sediment of class1818n layer16","fraction of sediment of class1919n layer16","fraction of sediment of class****n layer16","fraction of sediment of class**n layer17","fraction of sediment of class11n layer17","fraction of sediment of class22n layer17","fraction of sediment of class33n layer17","fraction of sediment of class44n layer17","fraction of sediment of class55n layer17","fraction of sediment of class66n layer17","fraction of sediment of class77n layer17","fraction of sediment of class88n layer17","fraction of sediment of class99n layer17","fraction of sediment of class1010n layer17","fraction of sediment of class1111n layer17","fraction of sediment of class1212n layer17","fraction of sediment of class1313n layer17","fraction of sediment of class1414n layer17","fraction of sediment of class1515n layer17","fraction of sediment of class1616n layer17","fraction of sediment of class1717n layer17","fraction of sediment of class1818n layer17","fraction of sediment of class1919n layer17","fraction of sediment of class****n layer17","fraction of sediment of class**n layer18","fraction of sediment of class11n layer18","fraction of sediment of class22n layer18","fraction of sediment of class33n layer18","fraction of sediment of class44n layer18","fraction of sediment of class55n layer18","fraction of sediment of class66n layer18","fraction of sediment of class77n layer18","fraction of sediment of class88n layer18","fraction of sediment of class99n layer18","fraction of sediment of class1010n layer18","fraction of sediment of class1111n layer18","fraction of sediment of class1212n layer18","fraction of sediment of class1313n layer18","fraction of sediment of class1414n layer18","fraction of sediment of class1515n layer18","fraction of sediment of class1616n layer18","fraction of sediment of class1717n layer18","fraction of sediment of class1818n layer18","fraction of sediment of class1919n layer18","fraction of sediment of class****n layer18","fraction of sediment of class**n layer19","fraction of sediment of class11n layer19","fraction of sediment of class22n layer19","fraction of sediment of class33n layer19","fraction of sediment of class44n layer19","fraction of sediment of class55n layer19","fraction of sediment of class66n layer19","fraction of sediment of class77n layer19","fraction of sediment of class88n layer19","fraction of sediment of class99n layer19","fraction of sediment of class1010n layer19","fraction of sediment of class1111n layer19","fraction of sediment of class1212n layer19","fraction of sediment of class1313n layer19","fraction of sediment of class1414n layer19","fraction of sediment of class1515n layer19","fraction of sediment of class1616n layer19","fraction of sediment of class1717n layer19","fraction of sediment of class1818n layer19","fraction of sediment of class1919n layer19","fraction of sediment of class****n layer19","fraction of sediment of class**n layer**","fraction of sediment of class11n layer**","fraction of sediment of class22n layer**","fraction of sediment of class33n layer**","fraction of sediment of class44n layer**","fraction of sediment of class55n layer**","fraction of sediment of class66n layer**","fraction of sediment of class77n layer**","fraction of sediment of class88n layer**","fraction of sediment of class99n layer**","fraction of sediment of class1010n layer**","fraction of sediment of class1111n layer**","fraction of sediment of class1212n layer**","fraction of sediment of class1313n layer**","fraction of sediment of class1414n layer**","fraction of sediment of class1515n layer**","fraction of sediment of class1616n layer**","fraction of sediment of class1717n layer**","fraction of sediment of class1818n layer**","fraction of sediment of class1919n layer**","fraction of sediment of class****n layer**","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the* layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the1 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the2 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the3 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the4 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the5 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the6 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the7 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the8 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the9 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the10 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the11 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the12 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the13 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the14 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the15 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the16 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the17 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the18 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the19 layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","thickness of the** layer","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer*","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer1","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer2","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer3","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer4","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer5","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer6","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer7","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer8","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer9","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer10","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer11","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer12","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer13","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer14","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer15","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer16","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer17","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer18","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer19","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","concentration of bed layer**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","bed load transport rate of sediment of class*","bed load transport rate of sediment of class1","bed load transport rate of sediment of class2","bed load transport rate of sediment of class3","bed load transport rate of sediment of class4","bed load transport rate of sediment of class5","bed load transport rate of sediment of class6","bed load transport rate of sediment of class7","bed load transport rate of sediment of class8","bed load transport rate of sediment of class9","bed load transport rate of sediment of class10","bed load transport rate of sediment of class11","bed load transport rate of sediment of class12","bed load transport rate of sediment of class13","bed load transport rate of sediment of class14","bed load transport rate of sediment of class15","bed load transport rate of sediment of class16","bed load transport rate of sediment of class17","bed load transport rate of sediment of class18","bed load transport rate of sediment of class19","bed load transport rate of sediment of class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","concentration volumic or mass concentration for class*","concentration volumic or mass concentration for class1","concentration volumic or mass concentration for class2","concentration volumic or mass concentration for class3","concentration volumic or mass concentration for class4","concentration volumic or mass concentration for class5","concentration volumic or mass concentration for class6","concentration volumic or mass concentration for class7","concentration volumic or mass concentration for class8","concentration volumic or mass concentration for class9","concentration volumic or mass concentration for class10","concentration volumic or mass concentration for class11","concentration volumic or mass concentration for class12","concentration volumic or mass concentration for class13","concentration volumic or mass concentration for class14","concentration volumic or mass concentration for class15","concentration volumic or mass concentration for class16","concentration volumic or mass concentration for class17","concentration volumic or mass concentration for class18","concentration volumic or mass concentration for class19","concentration volumic or mass concentration for class**","saturated concentration (kg/m3)","supplementary variable A","supplementary variable G","supplementary variable L","supplementary variable O","reference level for Nestor "],
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
    STARTING_TIME_OF_THE_HYDROGRAM = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = -1000.,
        fr = """Ce mot cle fixe le temps auquel le calcul SISYPHE debute
sauf en cas de suite de calcul car le temps initial est
alors lu sur le FICHIER PRECEDENT SEDIMENTOLOGIQUE. En non
permanent, ce mot cle designe de plus le temps correspondant
au premier enregistrement a lire dans le "fichier du calcul
precedent" contenant les donnees hydrodynamiques.""",
        ang = """this key word specifies the time when SISYPHE computation
begins except when a computation is continued (the initial time
is then read on the "previous sendimentological file".
For an unsteady case, it moreover specifies the time which
corresponds to the 1st record to be read in the "previous
computation file"(the file which contains the hydrodynamic
 data).""",
    ),
#   -----------------------------------
    NUMBER_OF_TIME_STEPS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Definit, en permanent, le nombre de pas de temps effectues
lors de l''execution du code. En non permanent, ce mot-cle
n''est pas utilise.""",
        ang = """Specifies, for a steady case, the number of time steps
performed when running the code. For an unsteady case, this
keyword is not used.""",
    ),
#   -----------------------------------
    TIME_STEP = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.,
        fr = """Definit le pas de temps du calcul en secondes en
permanent. En non permanent, ce pas de temps est donne par le
pas de sortie graphique du fichier precedent. Mais si aucun
nom n''est donne pour le ''FICHIER DU CALCUL PRECEDENT'' (les
variables hydrodynamiques pouvant etre imposees dans CONDIM),
la valeur du pas de temps donnee dans le fichier des
parametres est la aussi consideree.
Remarque : si le mot-cle "PAS DE TEMPS VARIABLE EST EGAL A OUI"
le pas de temps necessaire a une bonne resolution est calcule
dans le code et des sous iterations sont realisees.""",
        ang = """Specifies the time step in seconds in steady case.
For an unsteady case, this time step is fixed by the graphic
printout period of the previous computation file, except if
no name is given for the ''HYDRODYNAMIC FILE'' in
the steering file.
Remark : If the keyword ''VARIABLE TIME STEP'' is set equal
to yes, the time step required for a correct resolution is
computed in the code and sub-iterations are performed""",
    ),
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
    NUMBER_OF_TIDES_OR_FLOODS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """En non permanent, definit le nombre de marees ou de crues
effectuees lors de l''execution du code.En permanent, ce mot-cle
n''est pas utilise.""",
        ang = """For an unsteady case, specifies the number of tides or
floods performed when running the code. For a steady case, this
keyword is not used.""",
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
#   -----------------------------------
    TIDE_PERIOD = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 44640.,
        fr = """Fixe la valeur de la periode de l''evenement
(maree ou crue) en non permanent.""",
        ang = """Sets the period of the event (tide or flood)
for an unsteady case.""",
    ),
)
# -----------------------------------------------------------------------
EQUATIONS__BOUNDARY_CONDITIONS = PROC(nom= "EQUATIONS__BOUNDARY_CONDITIONS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    PRESCRIBED_SOLID_DISCHARGES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        fr = """ Valeurs des debits solides imposes aux frontieres
liquides entrantes (m3/s sans les vides).
Une valeur par frontiere liquide""",
        ang = """Values of prescribed solid discharges
at the inflow boundaries (m3/s without voids).
One value per liquid boundary""",
    ),
)
# -----------------------------------------------------------------------
SOLVER = PROC(nom= "SOLVER",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    SOLVER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["conjugate gradient on a normal equation","conjugate gradient","conjugate residual","minimum error","cgstab","gmres","direct"],
        defaut = ["conjugate gradient on a normal equation"],
        fr = """Permet de choisir le solveur utilise pour la resolution de
l''equation de continuite sur le fond (ce parametre n''est utilise
que si le mot cle MASS LUMPING est egal a faux). Toutes les methodes
proposees actuellement s''apparentent au Gradient Conjugue. Ce sont :
  1 : gradient conjugue
  2 : residu conjugue
  3 : gradient conjugue sur equation normale
  4 : erreur minimale
  5 : gradient conjugue carre (non programme)
  6 : gradient conjugue carre stabilise (cgstab)
  7 : gmres (voir aussi option du solveur)""",
        ang = """Makes it possible to select the solver used for solving the
bottom evolution equation (Used only if the key-word MASS LUMPING
is equal to false). All the currently available methods are
variations of the Conjugate Gradient method. They are as follows:
1: conjugate gradient
2: conjugate residual
3: conjugate gradient on a normal equation
4: minimum error
5: conjugate gradient squared (not implemented)
6: conjugate gradient squared stabilised (cgstab)
7: gmres (see option for solver)""",
    ),
#   -----------------------------------
    SOLVER_OPTION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [2],
        fr = """si le solveur est GMRES (7) le mot cle est la dimension de
l''espace de KRILOV (valeurs conseillees entre 2 et 15). Ce parametre
n''est utilise que si le mot cle MASS LUMPING est egal a faux.""",
        ang = """WHEN GMRES (7) IS CHOSEN, DIMENSION OF THE KRYLOV SPACE
TRY VALUES BETWEEN 2 AND 15. Used only if the key-word MASS LUMPING
is equal to false""",
    ),
#   -----------------------------------
    PRECONDITIONING = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["diagonal","no preconditioning","diagonal condensee","crout","diagonal and crout","diagonal condensed  and crout"],
        defaut = ["diagonal"],
        fr = """Permet de preconditionner le systeme lineaire afin d''accelerer
la convergence lors de sa resolution (ce parametre n''est utilise que
si le mot cle MASS LUMPING est egal a faux).
 - 0 : pas de preconditionnement;
 - 2 : preconditionnement diagonal.
 - 3 : preconditionnement diagonal-bloc
 - 7 : preconditionnement de Crout par element.
Certains preconditionnements sont cumulables
(les diagonaux 2 ou 3 avec les autres)
Pour cette raison on ne retient que les nombres premiers pour
designer les preconditionnements. Si l''on souhaite en cumuler
plusieurs on formera le produit des options correspondantes.""",
        ang = """Choice of the preconditioning in the resolution of the linear
system that the convergence is speeded up when it is being solved (Used
only if the key-word MASS LUMPING  is equal to false).
 0: no preconditioning
 2: diagonal preconditioning
 3: diagonal preconditioning with the condensed matrix
 7: Crout''s preconditioning per element (not implemented).
Some operations (either 2 or 3 diagonal preconditioning) can be
performed concurrently with the others.
Only prime numbers are therefore kept to denote the preconditioning
operations. When several of them are to be performed concurrently,
the product of relevant options shall be made.""",
    ),
#   -----------------------------------
    MAXIMUM_NUMBER_OF_ITERATIONS_FOR_SOLVER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [60],
        fr = """Les algorithmes utilises pour la resolution de l''etape de
propagation etant iteratifs; il est necessaire de limiter le nombre
d''iterations autorisees.
Remarque : Ce parametre n''est utilise que si le mot cle MASS LUMPING
est egal a faux.""",
        ang = """Since the algorithms used for solving the propagation step are
iterative, the allowed number of iterations should be limited.
NOTE:Used only if the key-word MASS LUMPING is equal to false .""",
    ),
#   -----------------------------------
    SOLVER_OPTION_FOR_SUSPENSION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [2],
        fr = """Parametre supplementaire disponible pour le solveur.
Dans le cas du solveur gmres, il s''agit de la dimension
de l''espace de Krylov.""",
        ang = """""",
    ),
#   -----------------------------------
    MAXIMUM_NUMBER_OF_ITERATIONS_FOR_SOLVER_FOR_SUSPENSION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [50],
        fr = """""",
        ang = """""",
    ),
#   -----------------------------------
    SOLVER_ACCURACY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [1.E-7],
        fr = """Precision demandee pour la resolution du systeme
(ce parametre n''est utilise que si le mot cle MASS LUMPING
 est egal a faux).""",
        ang = """Required accuracy for solving the linear system
(used only if the key word MASS LUMPING is equal to false).""",
    ),
#   -----------------------------------
    SOLVER_ACCURACY_FOR_SUSPENSION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [1.E-8],
        fr = """""",
        ang = """""",
    ),
#   -----------------------------------
    PRECONDITIONING_FOR_SUSPENSION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["diagonal","no preconditioning","diagonal condensee","crout"],
        defaut = ["diagonal"],
        fr = """Permet de preconditionner le systeme de l''etape de propagation
afin d''accelerer la convergence lors de sa resolution.
   - 0 : pas de preconditionnement,
   - 2 : preconditionnement diagonal.
   - 3 : preconditionnement diagonal avec la matrice condensee.
   - 7 : preconditionnement de Crout par element.
 Certains preconditionnements sont cumulables
 (les diagonaux 2 ou 3 avec les autres)
 Pour cette raison on ne retient que les nombres premiers pour
 designer les preconditionnements. Si l''on souhaite en cumuler
 plusieurs on formera le produit des options correspondantes.""",
        ang = """""",
    ),
)
# -----------------------------------------------------------------------
NUMERICAL = PROC(nom= "NUMERICAL",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    MASS_LUMPING = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Si oui, on effectue du mass-lumping sur la cote du fond.
Ceci revient a ramener toute la matrice de masse sur sa
diagonale lors de la resolution du systeme. Cette technique permet
d''accelerer le code et de le rendre egalement plus stable. Cependant,
les solutions obtenues se trouvent lissees.""",
        ang = """If this key word is equal to yes, the mass matrix is then
condensed on its diagonal.This technique is used to accelerate the
computation and also to make it more stable.However, the solutions
obtained are smoothed.""",
    ),
#   -----------------------------------
    TETA = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Definit le coefficient d''implicitation du schema numerique.""",
        ang = """Specifies the implicitation coefficient of the numerical
scheme.""",
    ),
#   -----------------------------------
    ZERO = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-10,
        fr = """Fixe le zero du code.""",
        ang = """Sets the zero of the code.""",
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
    TYPE_OF_ADVECTION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["CHARACTERISTICS","SUPG","CONSERVATIVE N-SCHEME LP","CONSERVATIVE N-SCHEME","CONSERVATIVE PSI-SCHEME","NON CONSERVATIVE PSI SCHEME","IMPLICIT NON CONSERVATIVE N SCHEME","EDGE-BASED N-SCHEME LP","EDGE-BASED N-SCHEME","ERIA SCHEME"],
        defaut = ["CHARACTERISTICS"],
        fr = """Determine le schema utilise pour la convection :
1: caracteristiques
2: schema semi-implicite + supg
3 et 4: schema N
5: schema psi
6: schema psi non conservatif
7: schema N non conservatif
13 et 14: schema N par segment (recommande pour les bancs decouvrants)""",
        ang = """Scheme used for advection of suspended sediment :
1: characteristics
2: semi-implicit SUPG
3 et 4: N scheme
5: psi scheme
6: non conservative psi scheme
7: non conservative N scheme
13 et 14: Edge-based N scheme (recommended for tidal flats)""",
    ),
#   -----------------------------------
    SUPG_OPTION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [2],
        fr = """""",
        ang = """""",
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
#   -----------------------------------
    OPTION_FOR_THE_DIFFUSION_OF_TRACER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["div( nu grad(T) )","1/h div ( h nu grad(T)"],
        defaut = ["div( nu grad(T) )"],
        fr = """1: Diffusion de la forme div( nu grad(T) )
        2: Diffusion de la forme 1/h div ( h nu grad(T) )""",
        ang = """1: Diffusion in the form div( nu grad(T) )
        2: Diffusion in the form 1/h div ( h nu grad(T) )""",
    ),
)
# -----------------------------------------------------------------------
NUMERICAL_PARAMETERS = PROC(nom= "NUMERICAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    MAXIMUM_NUMBER_OF_ITERATIONS_FOR_ADVECTION_SCHEMES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [50],
        fr = """Seulement pour schemes 13 et 14.
Ancienne valeur par defaut = 10 jusqu en version 8.1.""",
        ang = """Only for schemes 13 and 14.
Old default value = 10 until release 8.1.""",
    ),
#   -----------------------------------
    GENERAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        PARTITIONING_TOOL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['METIS','SCOTCH','PARMETIS','PTSCOTCH'],
            defaut = 'METIS',
            fr = """CHOIX DU PARTITIONNEUR
1 : METIS
2 : SCOTCH
3 : PARMETIS
4 : PTSCOTCH
etc...""",
            ang = """PARTITIONING TOOL SELECTION
1 : METIS
2 : SCOTCH
3 : PARMETIS
4 : PTSCOTCH
etc...""",
        ),
#       -----------------------------------
        NUMBER_OF_CORRECTIONS_OF_DISTRIBUTIVE_SCHEMES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """Pour les options avec predicteur-correcteur""",
            ang = """For predictor-corrector options""",
        ),
#       -----------------------------------
        NUMBER_OF_SUB_STEPS_OF_DISTRIBUTIVE_SCHEMES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """Pour les options predicteur-correcteur
avec schema localement implicite""",
            ang = """Only for implicit scheme with predictor-corrector""",
        ),
#       -----------------------------------
        TREATMENT_OF_FLUXES_AT_THE_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["Priority to prescribed values","Priority to fluxes"],
            defaut = ["Priority to fluxes"],
            fr = """Utilise pour les schemas PSI et N, avec option 2, on ne retrouve
pas exactement les valeurs imposees des traceurs,
mais le flux est correct""",
            ang = """Used so far only with the PSI and N schemes.
With option 2, Dirichlet prescribed values are not obeyed,
but the fluxes are correct""",
        ),
#       -----------------------------------
        BED_ROUGHNESS_PREDICTOR_OPTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """1 : Fond plat ks=KSPRATIO D50,
2: Fond ride (methode de Wiberg et Harris),
3: Dunes et megarides (methode de Van Rijn)""",
            ang = """1: Flat bed, 2: Rippled bed,
3: Dunes and mega ripples (Method of Van Rijn)""",
        ),
    ),
#   -----------------------------------
    SOLVER = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUSPENSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SOLVER_FOR_SUSPENSION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","gmres (see option for the solver for tracer diffusion)","direct"],
                defaut = ["conjugate gradient on a normal equation"],
                fr = """Permet de choisir le solveur utilise pour la resolution
de l''etape de propagation. Toutes les methodes proposees
actuellement s''apparentent au Gradient Conjugue. Ce sont :
 1 : gradient conjugue
 2 : residu conjugue
 3 : gradient conjugue sur equation normale
 4 : erreur minimale
 7 : gmres""",
                ang = """""",
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
PHYSICS = PROC(nom= "PHYSICS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    WATER_DENSITY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1000.,
        fr = """Fixe la valeur de la masse volumique de l''eau.""",
        ang = """sets the value of water density.""",
    ),
#   -----------------------------------
    SEDIMENT_DENSITY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 2650.,
        fr = """Fixe la valeur de la masse volumique du sediment
           en Kg/m3""",
        ang = """sets the value of the sediment density
           en Kg/m3""",
    ),
#   -----------------------------------
    NON_COHESIVE_BED_POROSITY = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.4,
        fr = """la concentration volumique  du lit est definie par
CSF= (1-porosite)
Ce parametre est utilise pour les sediments non-cohesifs.""",
        ang = """The bed volume concentration CSF=(1-porosity) is used to
calculate the bed evolution of non-cohesive sand transport.""",
    ),
#   -----------------------------------
    GRAVITY_ACCELERATION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 9.81,
        fr = """Fixe la valeur de l''acceleration de la pesanteur.
         M/S2""",
        ang = """Sets the value of the acceleration due to gravity.
         M/S2""",
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
#   -----------------------------------
    SETTLING_LAG = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Uses the velocity profile based on the Miles approximation""",
        ang = """Uses the velocity profile based on the Miles approximation""",
    ),
)
# -----------------------------------------------------------------------
SUSPENSION = PROC(nom= "SUSPENSION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    SETTLING_VELOCITIES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=10, max=10,
        fr = """Pas de valeur par defaut
Si non donne par l utilisateur, on utilise la subroutine vitchu-sisyphe
:
formules de Stokes, Zanke ou Van Rijn, selon la taille des grains""",
        ang = """The default value is not given. If the user does not
give a value, the subroutine vitchu-sisyphe is used:
Stockes, Zanke or Van Rijn formulae depending on the grain size""",
    ),
#   -----------------------------------
    SUSPENSION = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """calcul avec suspension""",
        ang = """""",
    ),
#   -----------------------------------
    DISPERSION_ALONG_THE_FLOW = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [1.E-2],
        fr = """""",
        ang = """""",
    ),
#   -----------------------------------
    DISPERSION_ACROSS_THE_FLOW = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [1.E-2],
        fr = """""",
        ang = """""",
    ),
#   -----------------------------------
    EQUILIBRIUM_INFLOW_CONCENTRATION = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """impose la concentration en entree du domaine et les
concentrations initiales en utilisant la formule de Fredsoe pour les
sediments non-cohesifs""",
        ang = """impose the equilibrium concentration for the inflow and at t=0
in the whole domain thanks to the formula of Fredsoe for non cohesive
sediments""",
    ),
#   -----------------------------------
    REFERENCE_CONCENTRATION_FORMULA = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """1 : formule de Zyserman et Fredsoe
        2 : methode de Bijker. La concentration au fond
        est reliee au taux de transport par charriage
        3 : formule de Van Rijn (1987)
        4 : formule de Soulsy\_van Rijn""",
        ang = """1 : Zysderman and Fredsoe, equilibrium formula
         2: Bijker method. The near bed concentration
         is related to the bedload . This option cannot be used
         without bedload transport
         3: Van Rijn formula
         4: Soulsby\_van Rijn formula""",
    ),
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
    INITIAL_SUSPENSION_CONCENTRATIONS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min= 2, max= 2,
        fr = """Pour la suspension, sert a initialiser la valeur
 de la concentration volumique pour chaque classe. ne sera pas pris
 en compte si CONCENTRATION EN ENTREE IMPOSEE=OUI""",
        ang = """In case of suspension, will be used to initialize the value
of volume concentration for each class. Will not be used if
 EQUILIBRIUM INFLOW CONCENTRATION=YES""",
    ),
#   -----------------------------------
    CONCENTRATION_PER_CLASS_AT_BOUNDARIES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min= 2, max= 2,
        fr = """Pour la suspension, sert a initialiser la valeur
 de la concentration volumique de chaque classe pour chaque frontiere
 ordre : frontiere 1 (classe 1, classe 2, etc.) , puis frontiere 2, etc""",
        ang = """In case of suspension, will be used to initialize the value
of volume concentration for each class and each boundary
order: boundary 1 (class 1, class2, etc., then boundary 2, etc.""",
    ),
#   -----------------------------------
    DIFFUSION = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [True ],
        fr = """Prise en compte de la diffusion de la concentration de sediment
en suspension""",
        ang = """If yes, diffusion of the concentration of suspended
sediment is done""",
    ),
#   -----------------------------------
    OPTION_FOR_THE_DISPERSION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """ 1 les mots cles dispersion longitudinale
et dispersion transversale permettent d affecter une valeur constante,
2 K1=alphal u*h et K2=alphat u*h affectent les valeurs alphal et alphat
(par default alphal=6 et alphat=0.6, 3 dipersion fournie par
 telemac2d""",
        ang = """ 1 les mots cles dispersion longitudinale
et dispersion transversale permettent d affecter une valeur constante,
2 K1=alphal u*h et K2=alphat u*h affectent les valeurs alphal et alphat
(par default alphal=6 et alphat=0.6, 3 dipersion fournie par
 telemac2d""",
    ),
#   -----------------------------------
    TETA_SUSPENSION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [1.],
        fr = """facteur implicitation du terme de depot et de la diffusion.
 si teta =0, tout le terme de depot est traite de maniere explicite.""",
        ang = """ implicitation factor for the deposition flux and the diffusion.
 for teta =0, the deposition flux is only explicit.""",
    ),
#   -----------------------------------
    CRITICAL_SHEAR_VELOCITY_FOR_MUD_DEPOSITION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1000.,
        fr = """Vitesse critique de depot de la vase (m/s)""",
        ang = """Critical shear velocity for deposition (m/s)""",
    ),
#   -----------------------------------
    PARTHENIADES_CONSTANT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-03,
        fr = """constante de la loi d''erosion de la vase (Kg/m2/s)""",
        ang = """constant of the Krone and Partheniades erosion law (Kg/m2/s)""",
    ),
)
# -----------------------------------------------------------------------
BED_MATERIAL = PROC(nom= "BED_MATERIAL",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    D90 = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        defaut = [.01,.01,.01,.01,.01,.01,.01,.01,.01,.01],
        fr = """Sets value of diameter d90 for particular size class.
Si le mot cle n est pas entre, la valeur par defaut est
celle du diametre moyen des grains""",
        ang = """Sets value of diameter d90 for particular size class.
If the keyword is not in the sterring file, the default value
is the value of the mean diameter of the sediment.""",
    ),
#   -----------------------------------
    SEDIMENT_DIAMETERS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        defaut = [.01,.01,.01,.01,.01,.01,.01,.01,.01,.01],
        fr = """Sets value of diameter dm for particular size class.""",
        ang = """Sets value of diameter dm for particular size class.""",
    ),
#   -----------------------------------
    HIDING_FACTOR_FOR_PARTICULAR_SIZE_CLASS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=10, max=10,
        defaut = [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],
        fr = """Fixe la valeur du facteur de pavage par classe
granulometrique""",
        ang = """Sets value of hiding factor for particular size class.""",
    ),
#   -----------------------------------
    NUMBER_OF_SIZE_CLASSES_OF_BED_MATERIAL = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """Fixe la valeur du nombre de classes granulometriques
de grains considerees dans le calcul""",
        ang = """Sets value of number of size classes of bed materials.""",
    ),
#   -----------------------------------
    INITIAL_FRACTION_FOR_PARTICULAR_SIZE_CLASS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=0, max='**',
        defaut = [1.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
        fr = """Sets value of initial fraction for particular size class.""",
        ang = """Sets value of initial fraction for particular size class.""",
    ),
#   -----------------------------------
    ACTIVE_LAYER_THICKNESS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [10000],
        fr = """Epaisseur de reference pour la stratification du lit. La
composition de la premiere couche sert a calculer le transport
solide. Utiliser une tres grande valeur pour ne pas avoir de
stratification.""",
        ang = """Thickness for bed stratification. Composition of first
layer is used to compute bed load transport rate. If you do not want
a stratification, use a large value""",
    ),
#   -----------------------------------
    HIDING_FACTOR_FORMULA = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [0],
        fr = """4 formules pour le hiding factor sont programmees dans SISYPHE
     0: const => il faut donner le HIDING FACTOR PAR CLASSE GRANULO
     1: Egiazaroff
     2: Ashida \& Michiue
      :
     4: Karim, Holly \& Yang""",
        ang = """4 hiding factor formulas are implemented in SISYPHE
     0: const => need to give HIDING FACTOR FOR PARTICULAR SIZE CLASS
     1: Egiazaroff
     2: Ashida \& Michiue
      :
     4: Karim, Holly \& Yang""",
    ),
#   -----------------------------------
    CONSTANT_ACTIVE_LAYER_THICKNESS = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """epaisseur de couche active constante ou non""",
        ang = """constant active layer thickness or not""",
    ),
#   -----------------------------------
    COHESIVE_SEDIMENTS = SIMP(statut ='f',
#   -----------------------------------
        typ = bool, min=0, max='**',
        defaut = [False,False,False,False,False,False,False,False,False,False],
        fr = """""",
        ang = """""",
    ),
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
            typ = 'TXM', min=13, max=13,
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
        C_VSM_DYNAMIC_ALT_MODEL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [5],
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
# -----------------------------------------------------------------------
FRICTION = PROC(nom= "FRICTION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    RATIO_BETWEEN_SKIN_FRICTION_AND_MEAN_DIAMETER = SIMP(statut ='f',
#   -----------------------------------
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
#   -----------------------------------
    SKIN_FRICTION_CORRECTION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """prise en compte du frottement de peau (voir aussi KSPRATIO)
       0 : pas de correction (TAUP= TOB)
       1 : fond plat (KSP= KSPRATIO * D50)
       2 : prise en compte des rides""",
        ang = """formula to predict the skin bed roughness (see also KSPRATIO)
       0 : NO correction (TAUP= TOB)
       1 : Flat bed (KSP= KSPRATIO * D50)
       2 : Ripple correction factor""",
    ),
#   -----------------------------------
    FRICTION_COEFFICIENT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 50.,
        fr = """Fixe la valeur du coefficient de Strickler
intervenant dans le calcul de la contrainte de
frottement au fond. Sa signification depend de
la LOI DE FROTTEMENT SUR LE FOND.""",
        ang = """Sets the value of the friction coefficient
to calculate the bed shear stress.
Depends on the LAW OF BOTTOM FRICTION.""",
    ),
#   -----------------------------------
    LAW_OF_BOTTOM_FRICTION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["FLAT BOTTOM","EQUILIBRIUM SAND RIPPLES (WAVES ONLY) KS=(MAX 3D50,ETA)","CHEZY","STRICKLER","MANNING","NIKURADSE"],
        defaut = ["STRICKLER"],
        fr = """Selectionne le type de formulation utilisee pour le calcul
du frottement sur le fond, voir COEFFICIENT DE FROTTEMENT.
Pour les lois possibles, voir CHOIX ci-dessus (cf. Note de principe).
En cas de couplage avec Telemac, le choix du frottement est impose par
Telemac, sauf si le mot-cle :PREDICTION DE LA RUGOSITE est mis a OUI""",
        ang = """Selects the type of formulation used for the bottom friction.
To know the possible laws see CHOIX1 above. See FRICTION COEFFICIENT.
Beware: in the case of internal coupling with Telemac, the friction
coefficient is selected in the Telemac steering file, except when BED
ROUGHNESS PREDICTION is set to YES""",
    ),
)
# -----------------------------------------------------------------------
SLOPE_EFFECT = PROC(nom= "SLOPE_EFFECT",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    FORMULA_FOR_SLOPE_EFFECT = SIMP(statut ='f',
#   -----------------------------------
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
#   -----------------------------------
    FRICTION_ANGLE_OF_THE_SEDIMENT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 40.,
        fr = """Angle de repos du sediment, intervient pour la prise de compte
 de la pente sur la contrainte critique par la formule de Soulsby.
Utiliser si ...=2""",
        ang = """Angle of repose of the sediment. Used in the Soulsby formula to
take into account the influence of bed slope on critical shear stress.
Use if ...=2""",
    ),
#   -----------------------------------
    FORMULA_FOR_DEVIATION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """1 : Koch et Flokstra
	2 : formule de Talmon et al. 1995, JHR 33(4) formules (1) et
(17) mot cle associe : BETA2""",
        ang = """1: Koch and Flokstra
	 2: formula of Talmon et al. 1995, JHR 33(4) formulas (1) and
(17) linked keyword : BETA2""",
    ),
#   -----------------------------------
    PARAMETER_FOR_DEVIATION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.85,
        fr = """Parametre pour la deviation causee par effet de pente
 pour la formule de Talmon et al.
Une valeur elevee provoque une faible deviation""",
        ang = """Parameter pour la deviation pour la formule de Talmon et al.""",
    ),
#   -----------------------------------
    SEDIMENT_SLIDE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Prise en compte de la pente d''equilibre du sable donnee par le
mot-cle ANGLE DE FROTTEMENT DU SEDIMENT""",
        ang = """If yes, the key-word FRICTION ANGLE OF THE SEDIMENT is taken
into account for slope stability""",
    ),
#   -----------------------------------
    SLOPE_EFFECT = SIMP(statut ='f',
#   -----------------------------------
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
#   -----------------------------------
    BETA = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.3,
        fr = """Determine la valeur du coefficient beta qui intervient dans la
formulation de l''effet de pente de Koch et Flokstra.""",
        ang = """Specifies the value of the beta coefficient used in the Koch
and Flokstra slope effect formulation.""",
    ),
)
# -----------------------------------------------------------------------
SEDIMENT_TRANSPORT = PROC(nom= "SEDIMENT_TRANSPORT",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NONEQUILIBRIUM_BED_LOAD = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        BED_ROUGHNESS_PREDICTION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Calcul de la rugosite de Nikuradse
- voir OPTION DU PREDICTEUR DE RUGOSITE -
la loi de frottement est forcee a 5 et le coefficient de frottement
ne sont pas utilises.
En cas de couplage, le frottement est envoye a Telemac""",
            ang = """The bed roughness is predicted according to the selected
BED ROUGHNESS PREDICTOR OPTION. In case of coupling with Telemac2d,
the calculated bed roughness is sent to Telemac.
The FRICTION COEFFICIENT and FRICTION LAW are no longer
used (KFROT is set to 5)""",
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
COMPUTATIONAL_INFORMATION = PROC(nom= "COMPUTATIONAL_INFORMATION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    GENERAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        MORPHOLOGICAL_FACTOR = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """Coefficient de l''echelle des temps""",
            ang = """Amplification for the morphological time scale""",
        ),
#       -----------------------------------
        MINIMUM_DEPTH_FOR_BEDLOAD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-2,
            fr = """Suppression des flux de sediment de et vers les points secs""",
            ang = """To cancel sediment fluxes to and from dry points""",
        ),
    ),
)
# -----------------------------------------------------------------------
BED_LOAD = PROC(nom= "BED_LOAD",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    BED_LOAD = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """calcul avec charriage""",
        ang = """""",
    ),
#   -----------------------------------
    BED_LOAD_TRANSPORT_FORMULA = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """11 formules de transport solide sont implementees dans SISYPHE.
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
transport solide (sous-programme QSFORM.F) en fixant la valeur
du mot cle a 0 :
                 0 :  IMPOSEE PAR L''UTILISATEUR
Attention : dans ce cas, il n''est pas possible de choisir l''option
PAS DE TEMPS VARIABLE.""",
        ang = """11 bed-load or total load transport formulas are implemented in
SISYPHE.
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
Users can also program other formulas (subroutine QSFORM.f) setting
this key word to zero :
                 0 : FORMULA PROGRAMMED BY USER
Warning : it is not then possible to choose the option
VARIABLE TIME-STEP""",
    ),
#   -----------------------------------
    B_VALUE_FOR_THE_BIJKER_FORMULA = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 2.E0,
        fr = """Coefficient b de la formule de Bijker""",
        ang = """b value for the Bijker formula""",
    ),
#   -----------------------------------
    MPM_COEFFICIENT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 8.0E-00,
        fr = """(-)""",
        ang = """Meyer-Peter Mueller Coefficient(-)""",
    ),
)
# -----------------------------------------------------------------------
EQUATIONS__ADVECTION = PROC(nom= "EQUATIONS__ADVECTION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    GENERAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SCHEME_OPTION_FOR_ADVECTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Si present remplace et a priorite sur :
OPTION POUR LES CARACTERISTIQUES (pas programme)
OPTION DE SUPG
Avec schema PSI : 1=explicit 2=predicteur-correcteur
pour les traceurs""",
            ang = """If present replaces and has priority over:
OPTION FOR CHARACTERISTICS (not yet implemented)
SUPG OPTION
IF PSI SCHEME: 1=explicit 2=predictor-corrector
for tracers""",
        ),
    ),
)
# -----------------------------------------------------------------------
CONSOLIDATION = PROC(nom= "CONSOLIDATION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    CONSOLIDATION_MODEL = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """1:Modele multicouche de Walther,
2: Modele de Thiebot (base sur la theorie de Gibson)""",
        ang = """1: Multilayer model of Walther,
2: Model of Thiebot (based on Gibson theory)""",
    ),
#   -----------------------------------
    GEL_CONCENTRATION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 310.E0,
        fr = """Concentration de transition pour modele de Thiebot (Kg/m3)""",
        ang = """Gel Concentration (Kg/m3)""",
    ),
#   -----------------------------------
    MAXIMUM_CONCENTRATION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 364.E0,
        fr = """Concentration maximale du modele de Thiebot (Kg/m3)""",
        ang = """Maximum Concentration for Thiebot consolidation model(Kg/m3)""",
    ),
#   -----------------------------------
    PERMEABILITY_COEFFICIENT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 8.E0,
        fr = """Coefficient de permeabilite pour le modele de consolidation""",
        ang = """Coefficient of permeability for consolidation model""",
    ),
#   -----------------------------------
    MUD_CONSOLIDATION = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Prise en compte du tassement par un modele multi-couche""",
        ang = """consolidation of the mud or sand mud-mixture sediment bed
 accounted for""",
    ),
#   -----------------------------------
    NUMBER_OF_LAYERS_OF_THE_CONSOLIDATION_MODEL = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """Structure verticale du lit cohesif- le nombre de couche doit
       etre inferieur a 10""",
        ang = """Vertical bed structure - The number of layers should be less
       than 10""",
    ),
#   -----------------------------------
    MASS_TRANSFER_PER_LAYER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min= 2, max= 2,
        fr = """Coefficients de transfert de masse du modele
de tassement multicouche en s-1""",
        ang = """Mass transfert coefficients of
the multilayer consolidation model in s-1""",
    ),
)
# -----------------------------------------------------------------------
COHESIVE_SEDIMENT = PROC(nom= "COHESIVE_SEDIMENT",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    MUD_CONCENTRATION_PER_LAYER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min= 2, max= 2,
        fr = """Concentration du lit de vase en g/ l - defini par couches""",
        ang = """Concentrations of the mud-bed in g per l (per layer)""",
    ),
#   -----------------------------------
    CRITICAL_EROSION_SHEAR_STRESS_OF_THE_MUD = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min=10, max=10,
        defaut = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,1.],
        fr = """Taux critique d erosion de la vase needs to be defined
        for each layer (N par m2)""",
        ang = """Critical erosion shear stress of the mud per layer (N per m2)""",
    ),
)
# -----------------------------------------------------------------------
INTERNAL = PROC(nom= "INTERNAL",op = None,
# -----------------------------------------------------------------------
    UIinfo = {"groupes": ("CACHE")},
#   -----------------------------------
    CONCATENATE_PARTEL_OUTPUT = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Avec cette option partel ne genereras non plus un fichier (GEO/CLI/PAR)
par processeur mais une concaténation de ceux-ci, ainsi qu''un fichier
d''index associé. Ainsi plutot que d''avoir 3P fichiers, il n''y en a
plus que 6.""",
        ang = """With this option partel no more generates a file (GEO/CLI/PAR) per
process but a single concatenate file of them, associated to an index
file. Then instead of having partel generating 3P files, it only
generates 6 files.""",
    ),
#   -----------------------------------
    DICTIONARY = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = 'sisyphev6p2.dico',
        fr = """Dictionnaire des mots cles.""",
        ang = """Key word dictionary.""",
    ),
)
TEXTE_NEW_JDC = "\
INPUT_OUTPUT__INFORMATION();\
MISCELLANEOUS();\
GENERAL();\
INPUT_OUTPUT__FILES();\
RESULTS();\
DATA_FILES();\
INITIAL_CONDITIONS();\
INPUT_OUTPUT__GRAPHICS_AND_LISTING();\
TIME();\
EQUATIONS__BOUNDARY_CONDITIONS();\
SOLVER();\
NUMERICAL();\
NUMERICAL_PARAMETERS();\
PHYSICS();\
SUSPENSION();\
BED_MATERIAL();\
FRICTION();\
SLOPE_EFFECT();\
SEDIMENT_TRANSPORT();\
SEDIMENTOLOGY();\
COMPUTATIONAL_INFORMATION();\
BED_LOAD();\
EQUATIONS__ADVECTION();\
CONSOLIDATION();\
COHESIVE_SEDIMENT();\
"
Ordre_Des_Commandes = (
'INPUT_OUTPUT__INFORMATION',
'MISCELLANEOUS',
'GENERAL',
'USELESS',
'INPUT_OUTPUT__FILES',
'RESULTS',
'DATA_FILES',
'INITIAL_CONDITIONS',
'INPUT_OUTPUT__GRAPHICS_AND_LISTING',
'TIME',
'EQUATIONS__BOUNDARY_CONDITIONS',
'SOLVER',
'NUMERICAL',
'NUMERICAL_PARAMETERS',
'PHYSICS',
'SUSPENSION',
'BED_MATERIAL',
'FRICTION',
'SLOPE_EFFECT',
'SEDIMENT_TRANSPORT',
'SEDIMENTOLOGY',
'COMPUTATIONAL_INFORMATION',
'BED_LOAD',
'EQUATIONS__ADVECTION',
'CONSOLIDATION',
'COHESIVE_SEDIMENT',
'INTERNAL')
try:
    import TelApy
    source = "eficas"
except Exception as excpt:
    source = "Telemac"
enum = source+'.sisyphe_enum_auto'
dicoCasEn = source+'.sisyphe_dicoCasEnToCata'
dicoCasFr = source+'.sisyphe_dicoCasFrToCata'
