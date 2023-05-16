
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



JdC = JDC_CATA (code = 'ARTEMIS',
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
    INPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DATA = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            GEOMETRY_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                fr = """Nom du fichier contenant le maillage du calcul a realiser.""",
                ang = """Name of the file which contains the computational mesh.""",
            ),
#           -----------------------------------
            BOTTOM_TOPOGRAPHY_SMOOTHING = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Nombre de lissages effectues sur la topographie.
Chaque lissage, effectue a l aide d une matrice de masse,
est conservatif.
Utilise lorsque les donnees de bathymetrie donnent des resultats
trop irreguliers apres interpolation.""",
                ang = """Number of smoothings done on the topography.
Each smoothing, using a mass matrix, is conservative.
It is used when bathymetric data provide too irregular results
after interpolation.""",
            ),
#           -----------------------------------
            FORTRAN_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'FichierOuRepertoire',
                defaut = '',
                fr = """Nom du fichier FORTRAN a soumettre, contenant les
sous-programmes specifiques au modele.""",
                ang = """Name of the FORTRAN file to be submitted, including specific
subroutines of the model.""",
            ),
#           -----------------------------------
            BOUNDARY_CONDITIONS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                fr = """Nom du fichier contenant les types de conditions aux limites.
Ce fichier est construit de facon automatique par le mailleur et \stbtel
au moyen de couleurs affectees aux noeuds des frontieres du domaine
de calcul.""",
                ang = """Name of the boundary conditions file. It is automatically built
by \stbtel or by the mesh generator MATISSE.""",
            ),
#           -----------------------------------
            BOTTOM_TOPOGRAPHY_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Nom du fichier eventuel contenant la bathymetrie associee au
maillage.
Si ce mot-cle est utilise, c est cette bathymetrie qui sera
utilisee pour le calcul.""",
                ang = """Name of a potential bathymetry file. If this keyword is specified,
the bathymetry which it is defining is accounted for.""",
            ),
#           -----------------------------------
            BINARY_DATA_FILE_1 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Fichier de donnees, code en binaire, mis a la disposition de
l utilisateur.""",
                ang = """Data file, written in binary mode, at the disposal of the user.""",
            ),
#           -----------------------------------
            BINARY_DATA_FILE_2 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Fichier de donnees, code en binaire, mis a la disposition de
l utilisateur.""",
                ang = """Data file, written in binary mode, at the disposal of the user.""",
            ),
#           -----------------------------------
            FORMATTED_DATA_FILE_1 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Fichier de donnees formate mis a la disposition de l utilisateur.""",
                ang = """Data file, written in ASCII mode, at the disposal of the user.""",
            ),
#           -----------------------------------
            FORMATTED_DATA_FILE_2 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Fichier de donnees formate mis a la disposition de l utilisateur.""",
                ang = """Data file, written in ASCII mode, at the disposal of the user.""",
            ),
#           -----------------------------------
            VALIDATION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Option utilisee principalement pour le dossier de validation.
Le \telkey{FICHIER DE REFERENCE} est alors considere comme une
reference a laquelle on va comparer le calcul. La comparaison est
effectuee par le sous-programme \telfile{BIEF\_VALIDA} qui peut etre une
comparaison avec une solution exacte par exemple.""",
                ang = """This option is primarily used for the validation documents.
The \telkey{REFERENCE FILE} is then considered as a
reference which the computation is going to be compared with.
The comparison is made by the subroutine \telfile{BIEF\_VALIDA},
which can be modified as to
so as to include, for example, a comparison with an exact solution.""",
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
                    fr = """Format du \telkey{FICHIER DE REFERENCE}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                    ang = """\telkey{REFERENCE FILE} format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
                ),
#               -----------------------------------
                REFERENCE_FILE = SIMP(statut ='f',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)'), max='**',
                    defaut = '',
                    fr = """Fichier de resultats de reference pour la validation.""",
                    ang = """Binary-coded result file for validation.""",
                ),
            ),
#           -----------------------------------
            GEOMETRY_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE GEOMETRIE}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{GEOMETRY FILE}.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            BINARY_DATA_FILE_1_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE DONNEES BINAIRE 1}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """\telkey{BINARY DATA FILE 1} format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            BINARY_DATA_FILE_2_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE DONNEES BINAIRE 2}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """\telkey{BINARY DATA FILE 2} format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
        ),
    ),
#   -----------------------------------
    GLOBAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TITLE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Titre du cas etudie.""",
            ang = """Title of the studied case.""",
        ),
#       -----------------------------------
        VECTOR_LENGTH = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Longueur du vecteur pour les machines vectorielles.""",
            ang = """Vector length on vector machines.""",
        ),
#       -----------------------------------
        PARALLEL_PROCESSORS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Nombre de processeurs pour la decomposition en parallele :
\begin{itemize}
\item 0 : 1 machine, compilation sans bibliotheque de parallelisme ;
\item 1 : 1 machine, compilation avec bibliotheque de parallelisme ;
\item 2 : 2 processeurs ou machines en parallele.
etc...
\end{itemize}""",
            ang = """Number of processors for domain partition:
\begin{itemize}
\item 0: 1 machine, compiling without parallel library,
\item 1: 1 machine, compiling with a parallel library,
\item 2: 2 processors or machines in parallel.
etc...
\end{itemize}""",
        ),
#       -----------------------------------
        CHECKING_THE_MESH = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si OUI on appelle le sous-programme \telfile{CHECKMESH}
qui verifie la coherence du maillage, points superposes, etc.""",
            ang = """If this keyword is equal to YES, a call to subroutine
\telfile{CHECKMESH} will look for errors in the mesh,
superimposed points, etc.""",
        ),
    ),
#   -----------------------------------
    OUTPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        RESULTS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            GRAPHIC_PRINTOUT_PERIOD = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode, en nombre de periodes de houle,
d impression des \telkey{VARIABLES POUR LES SORTIES GRAPHIQUES}
(voir ce mot-cle) dans le \telkey{FICHIER DES RESULTATS}.""",
                ang = """Fixes the period, in number of wave periods, for the writing
of the \telkey{VARIABLES FOR GRAPHIC PRINTOUTS} (see this keyword)
in the \telkey{RESULTS FILE}.""",
            ),
#           -----------------------------------
            VARIABLES_FOR_GRAPHIC_PRINTOUTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["wave height","wave phase","velocity u (free surface)(t=0)","velocity v (free surface)(t=0)","free surface elevation (t=0)","bottom elevation","still water height","phase velocity","group velocity","wave number","real potential","imaginal potential","prive(1,1)","prive(1,2)","prive(1,3)","prive(1,4)","first mean spectral period","second mean spectral period","third mean spectral period","force along X","force along Y","wave incidence radian","breaking rate","SXX stress","SXY stress","SYY stress"],
                defaut = [],
                fr = """Noms des variables que l utilisateur veut ecrire dans
le \telkey{FICHIER DES RESULTATS}.
Le choix des separateurs est libre.
Les possibilites offertes sont les suivantes :
\begin{itemize}
\item HS hauteur de la houle ;
\item PHAS phase de la houle ;
\item U0 vitesse u en surface      (a $t$ = 0) ;
\item V0 vitesse v en surface      (a $t$ = 0) ;
\item ZS cote de la surface libre  (a $t$ = 0) ;
\item ZF fond ;
\item HW hauteur d eau au repos ;
\item C vitesse de phase ;
\item CG vitesse de groupe ;
\item K nombre d onde ;
\item PHIR potentiel reel ;
\item PHII potentiel imaginaire ;
\item D prive(1,1) (variable 13) ;
\item E prive(1,2) (variable 14) ;
\item F prive(1,3) (variable 15) ;
\item G prive(1,4) (variable 16) ;
\item T01 premiere periode moyenne spectrale ;
\item T02 deuxieme periode moyenne spectrale ;
\item TM troisieme periode moyenne spectrale ;
\item FX force en X ;
\item FY force en Y ;
\item INC incidence de la houle ;
\item QB taux de deferlement ;
\item SXX contrainte SXX ;
\item SXY contrainte SXY ;
\item SYY contrainte SYY.
\end{itemize}
L utilisateur dispose de 4 champs libres, qu il peut
utiliser pour ecrire dans le fichier des resultats des variables
qu il cree lui-meme. Ces variables propres a l utlisateur doivent
etre calculees dans le sous-programme \telfile{CALRES} et le nom que
l on desire leur donner doit etre ecrit dans le sous-programme
\telfile{NOMVAR}.
Ces 4 champs sont :
\telfile{D, E, F, G} qui correspondent aux tableaux
\telfile{PRIVE(1,1), PRIVE(1,2), PRIVE(1,3), PRIVE(1,4)}.
A la difference des variables
     precedentes, celles-ci sont conservees dans tout le programme, et
     peuvent donc etre reutilisees.
     Dans ce dernier cas ne pas oublier de donner une taille
     suffisante au tableau PRIVE,
en precisant le parametre \telfile{NPRIV}
(dans le programme principal).""",
                ang = """Names of the variables that the user wants to write in the
\telkey{RESULTS FILE}.
Separators between variable names can be choosen free.
 The allowable values are:
\begin{itemize}
\item HS=wave height,
\item PHAS=wave phase,
\item U0=velocity u (free surface at $t$ = 0),
\item V0=velocity v (free surface at $t$ = 0),
\item ZS=free surface elevation (at $t$ = 0),
\item ZF=bottom elevation,
\item HW=still water height,
\item C=phase velocity,
\item CG=group velocity,
\item K=wave number,
\item PHIR=real potential,
\item PHII=imaginal potential,
\item D=prive(1,1),
\item E=prive(1,2),
\item F=prive(1,3),
\item G=prive(1,4),
\item T01=first mean spectral period,
\item T02=second mean spectral period,
\item TM=third mean spectral period,
\item FX=force along X,
\item FY=force along Y,
\item INC=wave incidence radian,
\item QB=breaking rate,
\item SXX=SXX stress,
\item SXY=SXY stress,
\item SYY=SYY stress.
\end{itemize}
The user has 4 free variables at his/her disposal to create other
variables by him/herself. These variables have to be computed
in the \telfile{CALRES} subroutine, and the name that we want to
attribute has to be precibed in the \telfile{NOMVAR} subroutine.
The 4 free variable fields are:
\telfile{D, E, F, G} which corresponds to the private arrays
\telfile{PRIVE(1,1), PRIVE(1,2), PRIVE(1,3) and PRIVE (1,4)}.
Contrary to the previous
     variables, these are conserved all through the computation, and can
     be used again.
     Do not forget to specify the number of private arrays you want to
use in the principal programme (variable \telfile{NPRIV}).""",
            ),
#           -----------------------------------
            RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les resultats du calcul,
avec la periodicite donnee par le mot cle
\telkey{PERIODE DE SORTIE GRAPHIQUE}.
Sur IBM, ce fichier est alloue automatiquement s il n existe pas,
avec les caracteristiques suivantes :
  Format d enregistrement   : VBS ;
  Longueur d enregistrement : X ;
  Taille de bloc            : 6204 ;
  Nombre de pistes          : 50 en primaire, 10 en secondaire.
La place memoire ainsi reservee est suffisante pour la plupart des
calculs de dimension moyenne.""",
                ang = """Name of the results file corresponding to the computations and
which contains the variables specified by the keyword
\telkey{VARIABLES FOR GRAPHIC PRINTOUTS}.""",
            ),
#           -----------------------------------
            BINARY_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier des resultats, code en binaire, mis a la disposition de
l utilisateur.""",
                ang = """Results file, written in binary mode, at the disposal of the user.""",
            ),
#           -----------------------------------
            FORMATTED_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier des resultats formate mis a la disposition de l utilisateur.""",
                ang = """Results file, written in ASCII mode, at the disposal of the user.""",
            ),
#           -----------------------------------
            NUMBER_OF_PRIVATE_VARIABLES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Permet de fixer le nombre de variables privees.""",
                ang = """Give the number of private variables.""",
            ),
#           -----------------------------------
            RESULTS_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE RESULTATS}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{RESULTS FILE}. Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
        ),
#       -----------------------------------
        LISTING = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            LISTING_PRINTOUT_PERIOD = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode, en nombre de periodes de houle,
d impression des \telkey{VARIABLES A IMPRIMER} (voir ce mot-cle). Pour
la mise au point, il faut savoir que la sortie des resultats est
effectuee systematiquement sur le fichier de retour d execution du code
(actuellement accessible par le menu 3.e de SPF sur IBM).""",
                ang = """Fixes the period, in number of wave periods, for the writing
of the \telkey{VARIABLES TO BE PRINTED} (see this keyword).""",
            ),
#           -----------------------------------
            LISTING_PRINTOUT = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Sortie des resultats sur support papier.
Si l on met NON le listing ne contient que l en-tete et la mention
FIN NORMALE DU PROGRAMME.
Commande a eviter.""",
                ang = """If NO is specified for this keyword, the printout listing just
contains the head and the sentence END OF PROGRAM.
It is advised not to use this way.""",
            ),
#           -----------------------------------
            INFORMATIONS_ABOUT_SOLVER = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Donne le nombre d iterations necessaires a la convergence du solveur.""",
                ang = """Gives the iterations number which was necessary for the solver
to converge.""",
            ),
#           -----------------------------------
            VARIABLES_TO_BE_PRINTED = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["wave height","wave phase","velocity u (free surface)(t=0)","velocity v (free surface)(t=0)","free surface elevation (t=0)","bottom elevation","still water height","phase velocity","group velocity","wave number","real potential","imaginal potential","prive(1,1)","prive(1,2)","prive(1,3)","prive(1,4)","first mean spectral period","second mean spectral period","third mean spectral period","force along X","force along Y","wave incidence radian","breaking rate","SXX stress","SXY stress","SYY stress"],
                defaut = '',
                fr = """Nom des variables que l utilisateur desire ecrire a l ecran.
Memes possibilites que pour les sorties graphiques.""",
                ang = """Name of variables taht the user whishes to write on the screen.
Possibilities are the same as for graphic outputs.""",
            ),
        ),
#       -----------------------------------
        FREE_SURFACE = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            FREE_SURFACE_ANIMATION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Option utilisee si l on veut produire le fichier de phases
et amplitudes qui servira a calculer la position de la surface libre
en tout point du maillage.
A besoin des noms de \telkey{FICHIER DES PHASES ET AMPLITUDES} et
\telkey{FICHIER DE SURFACE LIBRE}.""",
                ang = """This option informs that the phase and amplitude file is to
be generated. This file will be used to determine the free surface
elevation throughout the model area.
Remember to give \telkey{AMPLITUDE AND PHASE FILE} and
\telkey{FREE SURFACE FILE}.""",
            ),
#           -----------------------------------
            AMPLITUDE_AND_PHASE_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les phases et amplitudes
pour toutes les periodes et directions simulees dans le fichier des
parametres.""",
                ang = """Name of the results file corresponding to the computations and
which contains the phase and amplitude at every point for all the
periods and directions specified in the steering file.""",
            ),
#           -----------------------------------
            AMPLITUDE_AND_PHASE_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DES PHASES ET AMPLITUDES}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{AMPLITUDE AND PHASE FILE}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            FREE_SURFACE_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', max='**',
                defaut = '',
                fr = """Nom du fichier dans lequel sera ecrite la surface libre fonction du
temps pour permettre une animation du fichier.""",
                ang = """Name of the results file corresponding to the computations and
which contains the free surface at every point as a function of time.
This file can be animated to check free surface variations with time.""",
            ),
#           -----------------------------------
            FREE_SURFACE_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE SURFACE LIBRE}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{FREE SURFACE FILE}.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            FIRST_TIME_IN_THE_FREE_SURFACE_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [10000.],
                fr = """Determine le temps a partir duquel debute l ecriture
des resultats dans le \telkey{FICHIER DE SURFACE LIBRE}.""",
                ang = """Determines the time from which the results are written
in the \telkey{FREE SURFACE FILE}.""",
            ),
#           -----------------------------------
            TIME_STEP = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.25,
                fr = """Definit le pas de temps en secondes pour la sortie graphique
\telkey{FICHIER DE SURFACE LIBRE}.""",
                ang = """Specifies the time step in seconds in
 \telkey{FREE SURFACE FILE}.""",
            ),
#           -----------------------------------
            NUMBER_OF_TIME_STEPS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Definit le nombre de pas de temps pour la sortie graphique
\telkey{FICHIER DE SURFACE LIBRE}.""",
                ang = """Specifies the number of time steps in
 \telkey{FREE SURFACE FILE}.""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
GENERAL_PARAMETERS = PROC(nom= "GENERAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    DEBUGGER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 0,
        fr = """Pour imprimer la sequence des appels, mettre 1.""",
        ang = """If 1, calls of subroutines will be printed in the listing.""",
    ),
#   -----------------------------------
    TIME = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ORIGINAL_DATE_OF_TIME = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I', min= 3, max= 3,
            defaut = [0,0,0],
            fr = """Permet de fixer la date d origine des temps du modele lors
de la prise en compte de la force generatrice de la maree.""",
            ang = """Give the date of the time origin of the model when taking into
account the tide generating force.""",
        ),
#       -----------------------------------
        ORIGINAL_HOUR_OF_TIME = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I', min= 3, max= 3,
            defaut = [0,0,0],
            fr = """Permet de fixer l heure d origine des temps du modele lors
de la prise en compte de la force generatrice de la maree.""",
            ang = """Give the time of the time origin of the model when taking into
account the tide generating force.""",
        ),
    ),
#   -----------------------------------
    LOCATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ORIGIN_COORDINATES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I', min= 2, max= 2,
            defaut = [0,0],
            fr = """Valeur en metres, utilise pour eviter les trop grands
nombres, transmis
dans le format SERAFIN mais pas d autre traitement pour l instant.""",
            ang = """Value in metres, used to avoid large real numbers,
added in SERAFIN format, but so far no other treatment.""",
        ),
    ),
)
# -----------------------------------------------------------------------
NUMERICAL_PARAMETERS = PROC(nom= "NUMERICAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    MATRIX_STORAGE = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["Classical EBE","Assembled EBE","Edge-based storage"],
        defaut = "Edge-based storage",
        fr = """Permet de definir la methode de stockage des matrices.
Les choix possibles sont :
\begin{itemize}
\item 1 : EBE classique ;
\item 2 : EBE assemble ;
\item 3 : stockage par segments.
\end{itemize}
Attention, avec 2, il faut une numerotation speciale des points.""",
        ang = """Defines the method to store matrices. The possible choices are:
\begin{itemize}
\item 1: classical EBE,
\item 2: assembled EBE,
\item 3: edge-based storage.
\end{itemize}
Beware, with option 2, a special numbering of points is required.""",
    ),
#   -----------------------------------
    MATRIX_VECTOR_PRODUCT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Les choix possibles sont :
\begin{itemize}
\item 1 : Ancien Produit ;
\item 2 : Nouveau Produit Frontal.
\end{itemize}""",
        ang = """Possible choices are:
\begin{itemize}
\item 1: Classical Product,
\item 2: New Frontal Product.
\end{itemize}""",
    ),
#   -----------------------------------
    ZERO = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [1.E-12],
        fr = """Non active pour l instant.""",
        ang = """Non active at the moment.""",
    ),
#   -----------------------------------
    SOLVER_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        MAXIMUM_NUMBER_OF_ITERATIONS_FOR_SOLVER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 60000,
            fr = """Les algorithmes utilises pour la resolution du systeme
matriciel etant iteratifs, il est necessaire de limiter le nombre
d iterations autorisees.""",
            ang = """Algorithms used for solving the matrix system are iterative.
It is then necessary to limit the maximum number of iterations.""",
        ),
#       -----------------------------------
        PRECONDITIONING = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["no preconditioning","diagonal preconditioning","block-diagonal preconditioning","absolute value diagonal preconditioning","Crout preconditioning"],
            defaut = "diagonal preconditioning",
            fr = """Permet de preconditionner le systeme de l etape de propagation afin
d accelerer la convergence lors de sa resolution.
\begin{itemize}
\item 0 : pas de preconditionnement ;
\item 2 : preconditionnement diagonal ;
\item 3 : preconditionnement bloc-diagonal ;
\item 5 : preconditionnement diagonal en valeur absolue ;
\item 7 : preconditionnement de Crout par element.
\end{itemize}
Certains preconditionnements sont cumulables
(les diagonaux 2 ou 3 avec les autres).
Pour cette raison on ne retient que les nombres premiers pour
designer les preconditionnements. Si l on souhaite en cumuler
plusieurs on formera le produit des options correspondantes.""",
            ang = """Enables to apply preconditionning the matrix system to accelerate
the convergence of the solver.
\begin{itemize}
\item 0: no preconditionning,
\item 2: diagonal preconditionning,
\item 3: block-diagonal preconditionning,
\item 5: diagonal preconditionning in absolute value,
\item 7: Element Crout preconditionning.
\end{itemize}
Few of them can be combined
(numbers 2 or 3 with the other).
To combine some preconditionnings, impose the product of the previous
numbers: example 6 means preconditioning 2 and 3 applied.""",
        ),
#       -----------------------------------
        SOLVER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","CGSTAB","GMRES","direct","parallel direct (MUMPS)"],
            defaut = "direct",
            fr = """Permet de choisir le solveur utilise pour la resolution de l etape de
propagation. Toutes les methodes proposees actuellement s apparentent
au Gradient Conjugue. Ce sont :
\begin{itemize}
\item 1 : gradient conjugue ;
\item 2 : residu conjugue ;
\item 3 : gradient conjugue sur equation normale ;
\item 4 : erreur minimale ;
\item 5 : gradient conjugue carre (non programme) ;
\item 6 : gradient conjugue de type CGSTAB ;
\item 7 : GMRES ;
\item 8 : solveur direct.
\end{itemize}""",
            ang = """Enables to choose the solver used for solving the matrix system.
They are:
\begin{itemize}
\item 1: conjugate gradient,
\item 2: conjugate residual,
\item 3: conjugate gradient on the normal equation,
\item 4: minimum error,
\item 5: squared conjugate gradient (not programmed),
\item 6: CGSTAB conjugate gradient,
\item 7: GMRES,
\item 8: direct solver.
\end{itemize}""",
        ),
#       -----------------------------------
        SOLVER_OPTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 3,
            fr = """Parametre definissant la dimension de l espace de Krylov
pour le solveur 7 (GMRES).""",
            ang = """Defines the dimension of the Krylov space when using
the solver 7 (GMRES).""",
        ),
#       -----------------------------------
        SOLVER_ACCURACY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-4,
            fr = """Precision demandee pour la resolution de l equation de Berkhoff.""",
            ang = """Accuracy requested for the linear system solver.""",
        ),
    ),
#   -----------------------------------
    DISSIPATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUB_ITERATIONS_ACCURACY_FOR_DISSIPATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1.E-2],
            fr = """Donne la precision requise pour les sous-iterations du calcul
du coefficient de dissipation.""",
            ang = """Fixes the accuracy requested for sub-iterations necessary to
determine the dissipation coefficients.""",
        ),
#       -----------------------------------
        MAXIMUM_OF_SUB_ITERATIONS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [15],
            fr = """Donne le nombre maximum admis de sous-iterations pour le calcul
du coefficient de dissipation.""",
            ang = """Fixes the maximum number of sub-iterations for the computation
of dissipation.""",
        ),
#       -----------------------------------
        DISSIPATION_RELAXATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [0.5],
            fr = """Donne le coefficient de relaxation entre deux sous-iterations
pour le calcul du coefficient de dissipation.""",
            ang = """Fixes the relaxation coefficient used between two sub-iterations
for the computation of the dissipation term.""",
        ),
    ),
)
# -----------------------------------------------------------------------
BOUNDARY_CONDITIONS = PROC(nom= "BOUNDARY_CONDITIONS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    PHASE_REFERENCE_COORDINATES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min= 2, max= 2,
        defaut = [0.,0.],
        fr = """Coordonnees pour l origine des phases. Ne change rien
aux hauteurs de vagues calculees.""",
        ang = """Coordinates of reference point for phase. Will
not change the wave height computed.""",
    ),
#   -----------------------------------
    RANDOM_WAVE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_PERIODS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [5],
            fr = """Valeur utilisee avec l option
\telkey{HOULE ALEATOIRE MONODIRECTIONNELLE} = OUI
  ou avec l option
\telkey{HOULE ALEATOIRE MULTIDIRECTIONNELLE} = OUI.
Pour un calcul en houle aleatoire monodirectionnelle ou
multidirectionnelle, nombre de bandes d egale energie servant a
discretiser le spectre d energie en frequence.""",
            ang = """Used with option
\telkey{MONODIRECTIONAL RANDOM WAVE} = YES
   or
\telkey{MULTIDIRECTIONAL RANDOM WAVE} = YES.
It fixes the number of iso-energy frequency bands which discretize
the energy spectrum.""",
        ),
#       -----------------------------------
        PEAK_PERIOD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [10.0],
            fr = """Valeur utilisee avec l option
\telkey{HOULE ALEATOIRE MONODIRECTIONNELLE} = OUI
  ou avec l option
\telkey{HOULE ALEATOIRE MULTIDIRECTIONNELLE} = OUI.
Fixe la periode de pic (en secondes) du spectre d energie.""",
            ang = """Used with option
\telkey{MONODIRECTIONAL RANDOM WAVE} = YES
   or
\telkey{MULTIDIRECTIONAL RANDOM WAVE} = YES.
Fixes the peak period (in seconds) of the energy spectrum.""",
        ),
#       -----------------------------------
        GAMMA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [3.3],
            fr = """Valeur utilisee avec l option
\telkey{HOULE ALEATOIRE MONODIRECTIONNELLE} = OUI
  ou avec l option
\telkey{HOULE ALEATOIRE MULTIDIRECTIONNELLE} = OUI.
Indique la valeur de $\gamma$ pour le spectre d energie :
\begin{itemize}
\item GAMMA = 1   : spectre de Pierson-Moskowitz ;
\item GAMMA = 3.3 : spectre de JONSWAP moyen (valeur par defaut).
\end{itemize}""",
            ang = """Used with option
\telkey{MONODIRECTIONAL RANDOM WAVE} = YES
   or
\telkey{MULTIDIRECTIONAL RANDOM WAVE} = YES.
Fixes the $\gamma$ value tor the JONSWAP wave energy spectrum:
\begin{itemize}
\item GAMMA = 1: Pierson-Moskowitz,
\item GAMMA = 3.3: mean JONSWAP spectrum (default value).
\end{itemize}""",
        ),
#       -----------------------------------
        MINIMUM_SPECTRAL_PERIOD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [0.02],
            fr = """Valeur de la periode minimum voulue en secondes
si on veut tronquer le spectre pour le calcul
des periodes en houle aleatoire (voir sous-programme \telfile{PERALE}).""",
            ang = """Minimum period value requested in seconds
if it is necessary to alter the energy spectrum
for the computation of the periods in the case
of random waves (see \telfile{PERALE} subroutine).""",
        ),
#       -----------------------------------
        MAXIMUM_SPECTRAL_PERIOD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [200.],
            fr = """Valeur de la periode maximum voulue en secondes
si on veut tronquer le spectre pour le calcul
des periodes en houle aleatoire (voir sous-programme \telfile{PERALE}).""",
            ang = """Maximum period value requested in seconds
if it is necessary to alter the energy spectrum
for the computation of the periods in the case
of random waves (see \telfile{PERALE} subroutine).""",
        ),
#       -----------------------------------
        MONODIRECTIONAL = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            MONODIRECTIONAL_RANDOM_WAVE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Oui, si l on veut effectuer un calcul en houle aleatoire
monodirectionnelle (voir reels index 12, 13 et entier index 10).""",
                ang = """Yes, if one wants to run computation in random monodirectional waves
(see real keywords of index 12, 13 and integer of index 10).""",
            ),
        ),
#       -----------------------------------
        MULTIDIRECTIONAL = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            MULTIDIRECTIONAL_RANDOM_WAVE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Oui, si l on veut effectuer un calcul en houle aleatoire
multidirectionnelle (voir les reels index 12, 13, 14, 15 et 16 et
les entiers index 10 et 11.""",
                ang = """Yes, if one wants to run computation in random multidirectional waves
(see real keywords of index 12, 13 and integer of index 10).""",
            ),
#           -----------------------------------
            NUMBER_OF_DIRECTIONS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = [5],
                fr = """Valeur utilisee avec l option
\telkey{HOULE ALEATOIRE MULTIDIRECTIONNELLE} = OUI.
Pour un calcul en houle aleatoire multidirectionnelle,
nombre de bandes d egale energie servant a discretiser le spectre
directionnel d energie.""",
                ang = """Used with the option
\telkey{MULTIDIRECTIONAL RANDOM WAVE} = YES.
It fixes the number of iso-energy bands which discretizes the wave
directional spectrum.""",
            ),
#           -----------------------------------
            MINIMUM_ANGLE_OF_PROPAGATION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [-180.],
                fr = """Valeur utilisee avec l option
\telkey{HOULE ALEATOIRE MULTIDIRECTIONNELLE} = OUI.
Indique la borne inferieure de l intervalle des directions de
  propagation dans le cas d une houle aleatoire multidirectionnelle
  L angle est donne en degres et est compte positivement dans le sens
  direct a partir de l axe $x$.""",
                ang = """Used with the option
\telkey{MULTIDIRECTIONAL RANDOM WAVE} = YES.
Fixes the minimum value (in degrees) of the directions range. It is
counted positively in the trigonometric sense relatively to the $x$
axis.""",
            ),
#           -----------------------------------
            MAXIMUM_ANGLE_OF_PROPAGATION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [180.],
                fr = """Valeur utilisee avec l option
\telkey{HOULE ALEATOIRE MULTIDIRECTIONNELLE} = OUI.
Indique la borne superieure de l intervalle des directions de
  propagation dans le cas d une houle aleatoire multidirectionnelle.
  L angle est donne en degres et est compte positivement dans le sens
  direct a partir de l axe $x$.""",
                ang = """Used with the option
\telkey{MULTIDIRECTIONAL RANDOM WAVE} = YES.
Fixes the maximum value (in degrees) of the directions range. It is
counted positively in the trigonometric sense relatively to the $x$
axis.""",
            ),
#           -----------------------------------
            S_EXPONENT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [20.],
                fr = """Valeur utilisee avec l option
\telkey{HOULE ALEATOIRE MULTIDIRECTIONNELLE} = OUI.
Indique la valeur maximale de l exposant $s$ dans l expression donnant
  la repartition directionnelle de la houle.
  Cette expression est celle donnee par Goda dans Random Seas and
  Design of Maritime Structures - University of Tokyo Press:
  G(f,teta) = G0 * (cos(teta/2))**2s. f est la frequence et teta est
  la direction de propagation de la houle.""",
                ang = """Used with the option
\telkey{MULTIDIRECTIONAL RANDOM WAVE} = YES.
Fixes the maximum value of exponent $S$ in the Goda formula used to
express the directional wave energy spreading.
See GODA Y., Random Seas and Design of Maritime Structures - Univ.
of Tokyo Press, 1987.""",
            ),
        ),
    ),
#   -----------------------------------
    PERIODS_SCANNING = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        PERIOD_SCANNING = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Oui, si l on veut effectuer plusieurs calculs en balayant un
intervalle de periodes (voir reels index 8,9 et 10).""",
            ang = """Yes, if one wants to run computations by scanning a period range
(resonance computations, see also reals of index 8, 9, and 10).""",
        ),
#       -----------------------------------
        DATA = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            BEGINNING_PERIOD_FOR_PERIOD_SCANNING = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Valeur utilisee avec l option
\telkey{BALAYAGE EN PERIODE} = OUI.
Indique la borne gauche de l intervalle de periodes a parcourir
  (pour par exemple rechercher les periodes de resonances).""",
                ang = """Used with the option
\telkey{PERIOD SCANNING} = YES.
Fixes the minimum value (in seconds) of the period range to be used for
the period scanning.""",
            ),
#           -----------------------------------
            ENDING_PERIOD_FOR_PERIOD_SCANNING = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Valeur utilisee avec l option
\telkey{BALAYAGE EN PERIODE} = OUI.
Indique la borne droite de l intervalle de periodes a parcourir
  (pour par exemple rechercher les periodes de resonances).""",
                ang = """Used with the option
\telkey{PERIOD SCANNING} = YES.
Fixes the maximum value (in seconds) of the period range to be used for
the period scanning.""",
            ),
#           -----------------------------------
            STEP_FOR_PERIOD_SCANNING = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Valeur utilisee avec l option
\telkey{BALAYAGE EN PERIODE} = OUI.
Indique le pas a prendre pour effectuer le balayage en periodes
  (pour par exemple rechercher les periodes de resonances).""",
                ang = """Used with the option
\telkey{PERIOD SCANNING} = YES.
Fixes the value of the period step (in seconds) to be used for
the period scanning.""",
            ),
        ),
    ),
#   -----------------------------------
    MONOCHROMATIC_WAVE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WAVE_PERIOD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [10.],
            fr = """Definit la periode de la houle en mode monochromatique.""",
            ang = """Defines the wave period for monochromatic mode.""",
        ),
#       -----------------------------------
        DIRECTION_OF_WAVE_PROPAGATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [0.0],
            fr = """Donne la direction du vecteur d onde de la houle incidente. L angle
est donne en degres et est compte positivement dans le sens direct
a partir de l axe des $x$.
Il s agit de la direction principale de propagation.
Cette direction est la meme a toutes les frontieres maritimes.
Si l utilisateur veut specifier des directions differentes sur
differentes frontieres, il doit le faire dans son FORTRAN dans le
sous-programme \telfile{BORH} en specifiant la variable
\telfile{TETAB}.""",
            ang = """Fixes the direction towards the incident waves at boundaries go to.
It is counted in degrees and positively in the trigonometric sense
relatively to the $x$ axis.
This value is prescribed as a constant value along all the wave
incident type boundaries. If one wants to specify a non uniform value,
the user has to specify the value \telfile{TETAB} in the subroutine
 \telfile{BORH}.""",
        ),
#       -----------------------------------
        WAVE_HEIGHTS_SMOOTHING = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """OUI si on souhaite lisser les hauteurs de houle
pour ameliorer le calcul des contraintes de radiation
(actif uniquement en houle reguliere).
Valeur par defaut = NON.""",
            ang = """YES when one wants to smooth the wave heights
to improve the radiation stresses computation
(only used in regular wave mode).
Default value = NO.""",
        ),
    ),
#   -----------------------------------
    AUTOMATIC_TETAP_ANGLE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUB_ITERATIONS_ACCURACY_FOR_TETAP = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1.E-2],
            fr = """Donne la precision requise pour les sous-iterations du calcul
automatique de cos($\theta_p$).""",
            ang = """Fixes the accuracy requested for sub-iterations necessary to
determine value of $\theta_p$ (criterion on cos($\theta_p$)).""",
        ),
#       -----------------------------------
        AUTOMATIC_TETAP_CALCULATION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Si = OUI : calcul automatique des angles $\theta_p$
(basee sur la direction de la vitesse).""",
            ang = """If = TRUE: automatic calculation of $\theta_p$
(based on velocity direction).""",
        ),
#       -----------------------------------
        MAXIMUM_OF_SUB_ITERATIONS_FOR_TETAP = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [15],
            fr = """Donne le nombre maximum admis de sous-iterations pour le calcul
automatique de $\theta_p$.""",
            ang = """Fixes the maximum number of sub-iterations for the automatic
computation of $\theta_p$.""",
        ),
#       -----------------------------------
        RELAXATION_ON_TETAP = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1.],
            fr = """Donne le coefficient de relaxation entre deux sous-iterations
pour le calcul de l angle d incidence automatique.""",
            ang = """Fixes the relaxation coefficient used between two sub-iterations
for the computation of automatic $\theta_p$.""",
        ),
    ),
#   -----------------------------------
    AUTOMATIC_PHASE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        AUTOMATIC_CALCULATION_OF_PHASE = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """OUI : calcul automatique de la phase
(basee sur une profondeur de reference).""",
            ang = """TRUE: automatic calculation of incident phase
(based on reference water depth).""",
        ),
#       -----------------------------------
        REFERENCE_WATER_DEPTH_FOR_AUTOMATIC_PHASE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [-1.0],
            fr = """Profondeur de reference pour le calcul de la phase.
Essayez de mettre la frontiere incidente sur une zone
de bathymetrie homogene.
La profondeur a renseigner doit etre representative de la profondeur
d eau sur la frontiere.""",
            ang = """Water depth for automatic incident phase calculation.
Try to put the incident wave boundary on a regular topography zone.
The reference water depth should be representative of the water depth
on the boundary.""",
        ),
    ),
)
# -----------------------------------------------------------------------
PHYSICAL_PARAMETERS = PROC(nom= "PHYSICAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    GRAVITY_ACCELERATION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 9.81,
        fr = """Fixe la valeur de l acceleration de la pesanteur en m/s$^2$.""",
        ang = """Sets the value of the acceleration due to gravity in m/s$^2$.""",
    ),
#   -----------------------------------
    INITIALIZATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        INITIAL_WATER_LEVEL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Valeur utilisee avec l option
\telkey{CONDITIONS INITIALES} : ''COTE CONSTANTE''.""",
            ang = """Value to be used with the option
\telkey{INITIAL CONDITIONS}: ''CONSTANT ELEVATION''.""",
        ),
#       -----------------------------------
        INITIAL_DEPTH = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [0.],
            fr = """Valeur utilisee avec l option
\telkey{CONDITIONS INITIALES} : ''HAUTEUR CONSTANTE''.""",
            ang = """Value to be used along with the option
\telkey{INITIAL CONDITIONS}: ''CONSTANT DEPTH''.""",
        ),
#       -----------------------------------
        INITIAL_CONDITIONS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['ZERO ELEVATION','CONSTANT ELEVATION','ZERO DEPTH','CONSTANT DEPTH','SPECIAL'],
            defaut = 'ZERO ELEVATION',
            fr = """Permet de definir les conditions initiales sur les hauteurs d eau.
Les valeurs possibles sont :
\begin{itemize}
\item COTE NULLE : Initialise la cote de surface libre a 0.
           Les hauteurs d eau initiales sont alors retrouvees en
           faisant la difference entre les cotes de surface libre
           et du fond ;
\item COTE CONSTANTE : Initialise la cote de surface libre a la
valeur donnee par le mot-cle \telkey{COTE INITIALE}. Les hauteurs
d eau initiales sont calculees comme precedemment.
\item HAUTEUR NULLE : Initialise les hauteurs d eau a 0.
\item HAUTEUR CONSTANTE : Initialise les hauteurs d eau a la valeur
donnee par le mot-cle \telkey{HAUTEUR INITIALE}.
\item PARTICULIERES : Les conditions initiales sur la hauteur d eau
doivent etre precisees dans le sous-programme \telfile{CONDIH}.
\end{itemize}""",
            ang = """Enables to define the initial conditions on water depths.
The possible values are as follows:
\begin{itemize}
\item ZERO ELEVATION: fixes the free surface level to 0.
Water depths are then equal to the difference between
free surface level and bottom level,
\item CONSTANT ELEVATION: fixes the free surface level to the value
specified by the keyword \telkey{INITIAL WATER LEVEL}. Water
level are then computed as before,
\item ZERO DEPTH: initializes the water depths to 0,
\item CONSTANT DEPTH: initializes the water depths to the value
specified by the keyword \telkey{INITIAL DEPTH},
\item SPECIAL: initial conditions on water depths are to be
precised in the subroutine \telfile{CONDIH}.
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    DISSIPATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        BATHYMETRIC_BREAKING = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            BREAKING = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Oui, si l on souhaite integrer le processus de deferlement
bathymetrique (voir reels index 18, 19, 20, 21, 22, 23
et entiers index 12, 13).""",
                ang = """Yes, if one wants to account for breaking process (see also
reals of index 18, 19, 20, 21, 22, 23, and integer of index
12, 13).""",
            ),
#           -----------------------------------
            KDALLY = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.1],
                fr = """Donne le coefficient $K$ dans la formulation de la dissipation
par deferlement d apres Dally et al., 1984.""",
                ang = """Fixes the coefficient $K$ used in the formulation of the dissipation
coefficient proposed by Dally et al. 1984.""",
            ),
#           -----------------------------------
            BREAKING_LAW = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["BATTJES \& JANSSEN","DALLY"],
                defaut = ["BATTJES \& JANSSEN"],
                fr = """Specifie la formulation choisie pour le coefficient de dissipation
par deferlement. N est effectif qu en houle reguliere.
\begin{itemize}
\item 1 : Formulation de Battjes \& Jansen, 1978 ;
\item 2 : Formulation de Dally et al., 1984.
\end{itemize}
En houle aleatoire, la seule formulation utilisee est celle de
Battjes \& Janssen, 1978.""",
                ang = """Specifies the formulation chosen for calculating the dissipation
coefficient through breaking. Only effective for monochromatic wave
mode.
\begin{itemize}
\item 1: Formulation of Battjes \& Janssen, 1978,
\item 2: Formulation of Dally et al., 1984.
\end{itemize}
In random wave mode, the formulation of B \& J, 1978 is the only one
to be used.""",
            ),
#           -----------------------------------
            ALPHA = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [1.0],
                fr = """Donne le coefficient $\alpha$ dans la formulation de la dissipation
par deferlement en houle aleatoire d apres Battjes \& Janssen.""",
                ang = """Fixes the coefficient $\alpha$ used in the formulation of the
dissipation coefficient through breaking proposed by Battjes \& Janssen,
1978 for random waves.""",
            ),
#           -----------------------------------
            GAMMAS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.88],
                fr = """Donne le coefficient $\gamma_s$ dans le critere de la hauteur de
deferlement. Ne pas confondre avec le coefficient $\gamma$ qui
intervient dans la formule du spectre de JONSWAP.""",
                ang = """Fixes the coefficient $\gamma_s$ used in the criterion of the critical
breaking wave height. Do not confuse with coefficient $\gamma$
used in the JONSWAP spectrum.""",
            ),
#           -----------------------------------
            GDALLY = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.4],
                fr = """Donne le coefficient $\Gamma$ dans la formulation de la dissipation par
Dally et al., 1984. Ne pas confondre avec $\gamma$ (Formule de JONSWAP)
et Gammas (Critere de deferlement).""",
                ang = """Fixes the $\Gamma$ coefficient used in the formulation of Dally et al.,
1984, for the dissipation coefficient in surf-breaking. Do not confuse
with the coefficient $\gamma$ used in the JONSWAP formulae and
coefficient gammas used to determine the breaking wave height
criterion.""",
            ),
        ),
#       -----------------------------------
        BOTTOM_FRICTION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            FRICTION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Oui, si on veut prendre en compte le frottement sur le fond dans
la simulation.""",
                ang = """Yes, if one wants to include dissipation through bottom friction
in the computation.""",
            ),
#           -----------------------------------
            FLUID_KINEMATIC_VISCOSITY = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', max='**',
                defaut = [1.0E-6],
                fr = """Viscosite cinematique du fluide (eau) en m$^2$/s.
\begin{itemize}
\item 1.793E-6 : Pour une temperature de 0 C ;
\item 1.567E-6 : Pour une temperature de 4 C ;
\item 1.237E-6 : Pour une temperature de 12 C ;
\item 1.112E-6 : Pour une temperature de 16 C ;
\item 1.011E-6 : Pour une temperature de 20 C ;
\item 0.802E-6 : Pour une temperature de 30 C ;
\item 0.661E-6 : Pour une temperature de 40 C.
\item 1.0E-6   : Valeur par defaut.
\end{itemize}""",
                ang = """Kinematic viscosity of the fluid (water) in m$^2$/s.
\begin{itemize}
\item 1.793E-6: for a temperature of 0 C,
\item 1.567E-6: for a temperature of 4 C,
\item 1.237E-6: for a temperature of 12 C,
\item 1.112E-6: for a temperature of 16 C,
\item 1.011E-6: for a temperature of 20 C,
\item 0.802E-6: for a temperature of 30 C,
\item 0.661E-6: for a temperature of 40 C,
\item 1.0E-6: default value.
\end{itemize}""",
            ),
#           -----------------------------------
            DIAMETER90 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', max='**',
                defaut = [0.15E-3],
                fr = """\telfile{DIAM90} represente le diametre maximum, en m, de 90 \%
en poids des sediments.
\begin{itemize}
\item 1.0E-3   : Pour des sables tres grossiers ;
\item 0.5E-3   : Pour des sables grossiers ;
\item 0.25E-3  : Pour des sables moyens ;
\item 0.125E-3 : Pour des sables fins ;
\item 0.062E-3 : Pour des sables tres fins ;
\item 0.15E-3  : Valeur par defaut.
\end{itemize}""",
                ang = """\telfile{DIAM90} is the maximum grain diameter, in m, which defines
90 \% of the total weight of sediment.
\begin{itemize}
\item 1.0E-3: for very coarse sand,
\item 0.5E-3: for coarse sand,
\item 0.25E-3: for medium sand,
\item 0.125E-3: for fine sand,
\item 0.062E-3: for very fine sand,
\item 0.15E-3: default value.
\end{itemize}""",
            ),
#           -----------------------------------
            DIAMETER50 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', max='**',
                defaut = [0.10E-3],
                fr = """\telfile{DIAM50} represente le diametre maximum de 50 \% en poids des
sediments. En general, on a \telfile{DIAM90} = 1.5 $\times$
\telfile{DIAM50}
\telfile{DIAM50} est plus souvent donne dans des tables
\begin{itemize}
\item 0.66E-3  : Pour des sables tres grossiers ;
\item 0.33E-3  : Pour des sables grossiers ;
\item 0.17E-3  : Pour des sables moyens ;
\item 0.083E-3 : Pour des sables fins ;
\item 0.040E-3 : Pour des sables tres fins ;
\item 0.10E-3  : Valeur par defaut.
\end{itemize}""",
                ang = """\telfile{DIAM50} is the maximum grain diameter, in m, which defines
50 \% of the total weight of sediment. Usually, we have
\telfile{DIAM90} = 1.5 $\times$ \telfile{DIAM50}.
\telfile{DIAM50} is a more common value used.
\begin{itemize}
\item 0.66E-3: for very coarse sand,
\item 0.33E-3: for coarse sand,
\item 0.17E-3: for medium sand,
\item 0.083E-3: for fine sand,
\item 0.040E-3: for very fine sand,
\item 0.10E-3: default value.
\end{itemize}""",
            ),
#           -----------------------------------
            SEDIMENT_SPECIFIC_WEIGHT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [2650.0],
                fr = """Masse volumique du sediment en kg/m$^3$.""",
                ang = """Sediment specific weight in kg/m$^3$.""",
            ),
#           -----------------------------------
            FLUID_SPECIFIC_MASS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [1000.0],
                fr = """Masse volumique du fluide (eau) en kg/m$^3$.""",
                ang = """Fluid specific weight (water) in kg/m$^3$.""",
            ),
#           -----------------------------------
            HYDRAULIC_REGIME_IMPOSED = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Utilise avec l option \telkey{FROTTEMENT} = OUI.
Permet de choisir d imposer le regime hydraulique dans le cas
d un calcul automatique du facteur de frottement sur fonds sableux.""",
                ang = """Used with the option \telkey{FRICTION} = YES.
Enables to impose the hydraulic regime in the case of an automatic
calculation of the friction factor for sandy beds.""",
            ),
#           -----------------------------------
            SKIN_ROUGHNESS_ONLY = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Utilise avec l option \telkey{FROTTEMENT} = OUI.
Permet de choisir de ne prendre en compte
que la rugosite de peau dans le cas d un calcul automatique
du facteur de frottement sur fonds sableux.""",
                ang = """Used with the option \telkey{FRICTION} = YES.
Enables to restrict the total roughness to the skin roughness
in the case of an automatic calculation of the friction
factor for sandy beds.""",
            ),
#           -----------------------------------
            HYDRAULIC_REGIME_TYPE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', max='**',
                into = ["laminar regime","smooth-turbulent regime","rough-turbulent regime","transient regime"],
                defaut = ["laminar regime"],
                fr = """Utilise si le mot-cle \telkey{REGIME HYDRAULIQUE IMPOSE} = OUI.
Determine le regime hydraulique.
\begin{itemize}
\item regime laminaire ;
\item turbulent lisse ;
\item turbulent rugueux ;
\item transitoire.
\end{itemize}""",
                ang = """Used with option \telkey{HYDRAULIC REGIME IMPOSED} = YES.
Determines the type of the hydraulic regime (laminar,
smooth-turbulent, rough-turbulent, transient).
\begin{itemize}
\item laminar regime,
\item smooth-turbulent regime,
\item rough-turbulent regime,
\item transient regime.
\end{itemize}""",
            ),
#           -----------------------------------
            BOTTOM_FRICTION_LAW = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', max='**',
                into = ["Kostense formulation integrating Ue calculation (1986)","Putnam \& Johnson formulation (1949)"],
                defaut = ["Kostense formulation integrating Ue calculation (1986)"],
                fr = """Utilise avec l option \telkey{FROTTEMENT} = OUI.
Fixe le choix de la formulation du frottement :
\begin{itemize}
\item 1 : Kostense et al., 1986 ;
\item 2 : Putnam \& Johnson, 1949.
\end{itemize}""",
                ang = """Used with the option \telkey{FRICTION} = YES.
Fixes the formulation used for bottom friction law:
\begin{itemize}
\item 1: Kostense et al., 1986,
\item 2: Putnam \& Johnson, 1949.
\end{itemize}""",
            ),
#           -----------------------------------
            FRICTION_FACTOR_IMPOSED = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = [False],
                fr = """Utilise avec l option \telkey{FROTTEMENT} = OUI.
OUI permet de choisir d imposer un facteur de frottement, par un
mot-cle s il est uniforme (voir le reel d index 29) ou en
programmant dans le sous-programme \telfile{FWSPEC}.
Si NON, \artemis considere par defaut que les fonds sont sableux,
et calcule automatiquement le facteur de frottement avec les
caracteristiques du sediment et de l ecoulement.""",
                ang = """Used with the option \telkey{FRICTION} = YES.
If YES, enables the user to impose a friction factor, by a keyword
for a constant value (see real of index 29) or by programming in
the \telfile{FWSPEC} subroutine for non-uniform value.
If NO, \artemis automatically computes the friction factor assuming
that the bottom is sandy and uses the characteristics of sediment
and of motion.""",
            ),
#           -----------------------------------
            FRICTION_FACTOR = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Utilise si le mot-cle \telkey{FACTEUR DE FROTTEMENT IMPOSE} = OUI.
Fixe le facteur de frottement choisi uniforme sur le domaine.""",
                ang = """Used with the option \telkey{FRICTION FACTOR IMPOSED} = YES.
Fixes the value of the friction factor uniform over the domain.""",
            ),
#           -----------------------------------
            FRICTION_COEFFICIENT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """A ne pas confondre avec le \telkey{FACTEUR DE FROTTEMENT}.
Non utilise dans \artemis.
On le laisse par coherence avec \telemac{2d}.""",
                ang = """Do not confuse with the \telkey{FRICTION FACTOR}.
Not used in \artemis.
It is let here for consistence with \telemac{2d}.""",
            ),
#           -----------------------------------
            RIPPLES_COEFFICIENT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', max='**',
                defaut = [0.7],
                fr = """Specifie le coefficient de rides utilise dans la formule de
Van Rijn pour calculer le facteur de frottement :
\begin{itemize}
\item 1.0 : pour des rides seules,
\item 0.7 : pour des rides superposees a des vagues de sable.
\end{itemize}""",
                ang = """Fixes the ripples coefficient used in the formulae of Van Rijn
to calculate the friction factor:
\begin{itemize}
\item 1.0: for single ripples,
\item 0.7: for ripples superimposed to sand waves.
\end{itemize}""",
            ),
        ),
    ),
#   -----------------------------------
    CURRENT_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SUB_ITERATIONS_ACCURACY_FOR_CURRENT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1.E-2],
            fr = """Donne la precision requise pour les sous-iterations du calcul
du nombre d onde en presence de courant (vecteur d onde).""",
            ang = """Fixes the accuracy requested for sub-iterations necessary to
determine the wave vector.""",
        ),
#       -----------------------------------
        CURRENT = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Si = OUI : Prise en compte de la refraction de la houle par le courant.
Modele retenu : Kostense et Al. (1988).""",
            ang = """If = TRUE: Wave refraction due to current is described using
 Kostense model (1988).""",
        ),
    ),
#   -----------------------------------
    VARYING_TOPOGRAPHY = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        RAPIDLY_VARYING_TOPOGRAPHY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["mild-slope equation","gradient second order term","curvature second order term","gradient + curvature second order terms"],
            defaut = ["mild-slope equation"],
            fr = """Prise en compte des fortes pentes et courbures dans Berkhoff
\begin{itemize}
\item 0 : Berkhoff simple ;
\item 1 : prise en compte terme pente en grad($h$)$^2$ ;
\item 2 : prise en compte terme courbure en laplacien($h$) ;
\item 3 : prise en compte des termes de pente et courbure.
\end{itemize}
Modele retenu pour les fonctions E1 et E2 : Chamberlain
et Porter (1995).""",
            ang = """Extension of mild-slope equation with second
          order bottom effects
\begin{itemize}
\item 0: mild-slope equation,
\item 1: gradient second order term: grad($h$)$^2$,
\item 2: curvature second order term: laplacian($h$),
\item 3: gradient + curvature second order terms.
\end{itemize}
Model used for functions E1 and E2 expression: Chamberlain
and Porter (1995).""",
        ),
    ),
)
# -----------------------------------------------------------------------
INTERNAL = PROC(nom= "INTERNAL",op = None,
# -----------------------------------------------------------------------
    UIinfo = {"groupes": ("CACHE")},
#   -----------------------------------
    STEERING_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'), max='**',
        defaut = '',
        fr = """Nom du fichier contenant les parametres du calcul a realiser.""",
        ang = """Name of the steering file used for the computation.""",
    ),
#   -----------------------------------
    PARTITIONING_TOOL = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ['METIS','SCOTCH','PARMETIS','PTSCOTCH'],
        defaut = 'METIS',
        fr = """Choix du partitionneur :
\begin{itemize}
\item 1 : METIS ;
\item 2 : SCOTCH ;
\item 3 : PARMETIS ;
\item 4 : PTSCOTCH.
\end{itemize}""",
        ang = """Partitioning tool selection:
\begin{itemize}
\item 1: METIS,
\item 2: SCOTCH,
\item 3: PARMETIS,
\item 4: PTSCOTCH.
\end{itemize}""",
    ),
#   -----------------------------------
    DICTIONARY = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = 'artemis.dico',
        fr = """Dictionnaire des mots cles.""",
        ang = """Key word dictionary.""",
    ),
#   -----------------------------------
    CONCATENATE_PARTEL_OUTPUT = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Avec cette option, partel ne generera plus un fichier (GEO/CLI/PAR)
par processeur mais une concatenation de ceux-ci, ainsi qu un fichier
d index associe. Ainsi plutot que d avoir 3P fichiers, il n y en a
plus que 6.""",
        ang = """With this option, partel no more generates a file (GEO/CLI/PAR) per
process but a single concatenate file of them, associated to an index
file. Then instead of having partel generating 3P files, it only
generates 6 files.""",
    ),
)
# -----------------------------------------------------------------------
NESTING_WITH_TOMAWAC = PROC(nom= "NESTING_WITH_TOMAWAC",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NESTING_WITHIN_TOMAWAC_OUTER_MODEL = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["No nesting","Uses a TOMAWAC spectrum (1 point) as input to ARTEMIS","Uses several TOMAWAC spectra located along ARTEMIS boundary"],
        defaut = ["No nesting"],
        fr = """A besoin de NOMBRE DE DIRECTION DANS LE SPECTRE TOMAWAC et
NOMBRE DE FREQUENCES DANS LE SPECTRE TOMAWAC pour l option 1.
A besoin des noms de \telkey{FICHIER DE RESULTATS GLOBAL TOMAWAC},
\telkey{FICHIER DE SPECTRE GLOBAL TOMAWAC} et
\telkey{COORDONNEES POUR LE SPECTRE F DE REFERENCE} pour l option 2.
Donner \telkey{INSTANT DE LECTURE DU SPECTRE TOMAWAC} pour 1 et 2""",
        ang = """Remember to give NUMBER OF DIRECTIONS IN TOMAWAC SPECTRUM
and NUMBER OF FREQUENCIES IN TOMAWAC SPECTRUM with option 1.
Remember to give \telkey{TOMAWAC OUTER RESULT FILE},
\telkey{TOMAWAC OUTER SPECTRAL FILE} and
\telkey{COORDINATES OF THE REFERENCE F SPECTRUM} if using option 2.
\telkey{INSTANT FOR READING TOMAWAC SPECTRUM} required
in cases 1 and 2.""",
    ),
#   -----------------------------------
    INSTANT_FOR_READING_TOMAWAC_SPECTRUM = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [0.],
        fr = """Indique l instant de calcul \tomawac associe au spectre
qui doit etre importe dans \artemis.""",
        ang = """Give the instant of the \tomawac computation at which we
want to import the spectrum for \artemis.""",
    ),
#   -----------------------------------
    COORDINATES_OF_THE_REFERENCE_F_SPECTRUM = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R', min= 2, max= 2,
        defaut = [-99999.9,-99999.9],
        fr = """Coordonnees du point de reference pour le spectre en
frequence.""",
        ang = """Coordinates for the frequency spectrum reference point.""",
    ),
#   -----------------------------------
    INPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DATA = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            TOMAWAC_OUTER_SPECTRAL_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Nom du fichier de spectre provenant d un modele global
\tomawac.""",
                ang = """Name of the spectral file from an outer \tomawac wave
model.""",
            ),
#           -----------------------------------
            TOMAWAC_OUTER_SPECTRAL_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE SPECTRE GLOBAL TOMAWAC}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{TOMAWAC OUTER SPECTRAL FILE}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            TOMAWAC_LIQUID_BOUNDARY_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Nom du fichier de frontieres liquides provenant d un modele
global \tomawac.""",
                ang = """Name of a file containing the liquid boundaries derived
from an outer \tomawac wave model.""",
            ),
#           -----------------------------------
            TOMAWAC_LIQUID_BOUNDARY_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', max='**',
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE FRONTIERES LIQUIDES TOMAWAC}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{TOMAWAC LIQUID BOUNDARY FILE}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
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
            TOMAWAC_OUTER_RESULT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Nom du fichier de resultats 2D provenant d un modele global
\tomawac.""",
                ang = """Name of a file containing the results of an outer \tomawac
wave model.""",
            ),
#           -----------------------------------
            TOMAWAC_OUTER_RESULT_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DE RESULTATS GLOBAL TOMAWAC}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{TOMAWAC OUTER RESULT FILE}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
        ),
    ),
)
TEXTE_NEW_JDC = "\
COMPUTATION_ENVIRONMENT();\
GENERAL_PARAMETERS();\
NUMERICAL_PARAMETERS();\
BOUNDARY_CONDITIONS();\
PHYSICAL_PARAMETERS();\
NESTING_WITH_TOMAWAC();\
"
Ordre_Des_Commandes = (
'COMPUTATION_ENVIRONMENT',
'GENERAL_PARAMETERS',
'NUMERICAL_PARAMETERS',
'BOUNDARY_CONDITIONS',
'PHYSICAL_PARAMETERS',
'INTERNAL',
'NESTING_WITH_TOMAWAC')
try:
    import TelApy
    source = "eficas"
except Exception as excpt:
    source = "Telemac"
enum = source+'.artemis_enum_auto'
dicoCasEn = source+'.artemis_dicoCasEnToCata'
dicoCasFr = source+'.artemis_dicoCasFrToCata'
