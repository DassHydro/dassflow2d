
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



JdC = JDC_CATA (code = 'TELEMAC2D',
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
    INITIALIZATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TITLE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Titre du cas etudie.""",
            ang = """Title of the case being considered.""",
        ),
#       -----------------------------------
        INITIAL_CONDITIONS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ['ZERO ELEVATION','CONSTANT ELEVATION','ZERO DEPTH','CONSTANT DEPTH','SPECIAL','PARTICULIERES','PARTICULAR','TPXO SATELLITE ALTIMETRY'],
            defaut = 'ZERO ELEVATION',
            fr = """Permet de definir les conditions initiales sur
les hauteurs d''eau notamment.
Les valeurs possibles sont :
\begin{itemize}
\item COTE NULLE : Initialise la cote de surface libre a 0.
           Les hauteurs d''eau initiales sont alors retrouvees en
           faisant la difference entre les cotes de surface libre
           et du fond ;
\item COTE CONSTANTE : Initialise la cote de surface libre a la
valeur donnee par le mot-cle \telkey{COTE INITIALE}. Les hauteurs
d''eau initiales sont calculees comme precedemment ;
\item HAUTEUR NULLE : Initialise les hauteurs d''eau a 0 ;
\item HAUTEUR CONSTANTE : Initialise les hauteurs d''eau a la valeur
donnee par le mot-cle \telkey{HAUTEUR INITIALE} ;
\item ALTIMETRIE SATELLITE TPXO : Les conditions initiales sur la
hauteur d''eau et les vitesses sont etablies sur la base des donnees
           satellite TPXO dont les 8 premiers constistuents
ont ete extraites et
sauvees dans les fichiers
\telkey{BASE BINAIRE 1/2 DE DONNEES DE MAREE} ;
\item PARTICULIERES : Les conditions initiales sur la hauteur d''eau
doivent etre precisees dans le sous-programme \telkey{USER\_CONDIN\_H}.
\end{itemize}""",
            ang = """Makes it possible to define the initial conditions of
the water depth.
The possible values are as follows:
\begin{itemize}
\item ZERO ELEVATION: Initializes the free surface elevation to 0.
The initial water depths are then found by computing the difference
between the free surface and the bottom,
\item CONSTANT ELEVATION: Initializes the water elevation to the value
given by the keyword \telkey{INITIAL ELEVATION}.
The initial water depths are computed as in the previous case,
\item ZERO DEPTH: Initializes the water depths to 0.
\item CONSTANT DEPTH: Initializes the water depths to the value given
by the keyword \telkey{INITIAL DEPTH},
\item TPXO SATELITE ALTIMETRY: The initial conditions on the free
surface and velocities are established from the satellite program
data given by the harmonic constants database coming from OSU
(e.g. TPXO) and stored in the \telkey{BINARY DATABASE 1/2 FOR TIDE},
\item SPECIAL or PARTICULAR: The initial conditions with the water depth
should be stated in the \telfile{USER\_CONDIN\_H} subroutine.
\end{itemize}""",
        ),
#       -----------------------------------
        b_INITIAL_CONDITIONSG = BLOC(condition="INITIAL_CONDITIONS == 'CONSTANT ELEVATION'",
#       -----------------------------------
#           -----------------------------------
            INITIAL_ELEVATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Valeur utilisee avec l''option :
\telkey{CONDITIONS INITIALES} : ''COTE CONSTANTE''.""",
                ang = """Value to be used with the option:
\telkey{INITIAL CONDITIONS} : ''CONSTANT ELEVATION''.""",
            ),
        ),
#       -----------------------------------
        b_INITIAL_CONDITIONSH = BLOC(condition="INITIAL_CONDITIONS == 'CONSTANT DEPTH'",
#       -----------------------------------
#           -----------------------------------
            INITIAL_DEPTH = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Valeur utilisee avec l''option :
\telkey{CONDITIONS INITIALES} : ''HAUTEUR CONSTANTE''.""",
                ang = """Value to be used along with the option:
\telkey{INITIAL CONDITIONS} : ''CONSTANT DEPTH''.""",
            ),
        ),
#       -----------------------------------
        BINARY_DATA_FILE_1_FORMAT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['BIN','SERAFIN','SERAFIND','MED'],
            defaut = 'BIN',
            fr = """Format du \telkey{FICHIER DE DONNEES BINAIRE 1}.
Les valeurs possibles sont :
\begin{itemize}
\item BIN     : format binaire standard ;
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
            ang = """Format of the \telkey{BINARY DATA FILE 1}.
Possible values are:
\begin{itemize}
\item BIN     : standard binary format,
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
        ),
#       -----------------------------------
        BINARY_DATA_FILE_1 = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier de donnees code en binaire mis a la disposition
de l''utilisateur.""",
            ang = """Binary-coded data file available to the user.""",
        ),
#       -----------------------------------
        BINARY_DATA_FILE_2_FORMAT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['BIN','SERAFIN','SERAFIND','MED'],
            defaut = 'BIN',
            fr = """Format du \telkey{FICHIER DE DONNEES BINAIRE 2}.
Les valeurs possibles sont :
\begin{itemize}
\item BIN     : format binaire standard ;
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
            ang = """Format of the \telkey{BINARY DATA FILE 2}.
Possible values are:
\begin{itemize}
\item BIN     : standard binary format,
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
        ),
#       -----------------------------------
        BINARY_DATA_FILE_2 = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier de donnees code en binaire mis a la disposition
de l''utilisateur.""",
            ang = """Binary-coded data file available to the user.""",
        ),
#       -----------------------------------
        FORMATTED_DATA_FILE_1 = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier de donnees formate mis a la disposition de
l''utilisateur.""",
            ang = """Formatted data file available to the user.""",
        ),
#       -----------------------------------
        FORMATTED_DATA_FILE_2 = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier de donnees formate mis a la disposition de
l''utilisateur.""",
            ang = """Formatted data file available to the user.""",
        ),
#       -----------------------------------
        INPUT_FILES = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            GEOMETRY_FILE_FORMAT = SIMP(statut ='o',
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
            GEOMETRY_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                fr = """Nom du fichier contenant le maillage du calcul a realiser.""",
                ang = """Name of the file containing the mesh. This file may also
contain the topography and the friction coefficients.""",
            ),
#           -----------------------------------
            FORTRAN_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'FichierOuRepertoire',
                defaut = '',
                fr = """Nom du fichier ou repertoire FORTRAN a soumettre, contenant les
sous-programmes specifiques au modele.""",
                ang = """Name of the FORTRAN file or directory to be submitted,
including specific subroutines of the model.""",
            ),
#           -----------------------------------
            BOTTOM_TOPOGRAPHY_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier eventuel contenant la bathymetrie associee au
maillage.
Si ce mot-cle est utilise; c''est cette bathymetrie qui sera utilisee
pour le calcul.""",
                ang = """Name of the possible file containing the bathymetric data.
Where this keyword is used, these bathymetric data shall be used in
the computation.""",
            ),
#           -----------------------------------
            BOTTOM_SMOOTHINGS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Nombre de lissages effectues sur la topographie.
Chaque lissage, effectue a l''aide d''une matrice de masse,
est conservatif.
Utilise lorsque les donnees de bathymetrie donnent des resultats
trop irreguliers apres interpolation.""",
                ang = """Number of smoothings on bottom topography.
Each smoothing is mass conservative.
to be used when interpolation of bathymetry on the mesh gives
very rough results.""",
            ),
#           -----------------------------------
            BOTTOM_SMOOTHINGS_AFTER_USER_MODIFICATIONS = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Indique si le nombre d eventuels lissages sur la topographie est
effectue apres (ou sinon avant) les modifications de topographie
apportees par l utilisateur.""",
                ang = """Indicates if the number of potential smoothings on bottom topography
is done after (or before otherwise) the topography modifications
implemented by the user.""",
            ),
#           -----------------------------------
            BOUNDARY_CONDITIONS_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                fr = """Nom du fichier contenant les types de conditions aux limites.
Ce fichier est rempli de facon automatique par le mailleur au moyen de
couleurs affectees aux noeuds des frontieres du domaine de calcul.""",
                ang = """Name of the file containing the types of boundary conditions.
This file is filled automatically by the mesh generator through
through colours that are assigned to the boundary nodes.""",
            ),
#           -----------------------------------
            VALIDATION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Option utilisee principalement pour le dossier de validation. Le
\telkey{FICHIER DE REFERENCE} est alors considere comme une
reference a laquelle on va comparer le calcul. La comparaison est
effectuee par le sous-programme \telfile{BIEF\_VALIDA}
qui peut etre une comparaison avec une solution exacte par exemple.""",
                ang = """This option is primarily used for the validation documents.
The \telkey{REFERENCE FILE} is then considered as a reference
which the computation is going to be compared with.
The comparison is done by the subroutine \telfile{BIEF\_VALIDA},
which can be modified so as to include, for example,
a comparison with an exact solution.""",
            ),
#           -----------------------------------
            REFERENCE_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
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
                ang = """Format of the \telkey{REFERENCE FILE}.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            REFERENCE_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de resultats de reference pour la validation.""",
                ang = """Binary-coded result file for validation.""",
            ),
        ),
#       -----------------------------------
        GLOBAL = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            PARALLEL_PROCESSORS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Nombre de processeurs pour la decomposition en parallele:
\begin{itemize}
\item 0 : 1 machine, compilation sans bibliotheque de parallelisme ;
\item 1 : 1 machine, compilation avec bibliotheque de parallelisme ;
\item 2 : 2 processeurs ou machines en parallele etc...
\end{itemize}""",
                ang = """Number of processors for domain partition.
\begin{itemize}
\item 0: 1 machine, compiling without parallel library,
\item 1: 1 machine, compiling with a parallel library,
\item 2: 2 processors or machines in parallel etc...
\end{itemize}""",
            ),
#           -----------------------------------
            CHECKING_THE_MESH = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Si OUI on appelle le sous-programme \telfile{CHECKMESH} qui verifie
la coherence du maillage, points superposes, etc.""",
                ang = """If this keyword is equal to YES, a call to subroutine
\telfile{CHECKMESH} will look for errors in the mesh,
superimposed points, etc.""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_BOUNDARIES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 30,
                fr = """Nombre maximal de frontieres differentes dans le maillage.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
                ang = """Maximal number of boundaries in the mesh.
Used for dimensioning arrays. Can be increased if needed.""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 20,
                fr = """Nombre maximal de points sources dans le maillage.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
                ang = """Maximal number of punctual sources in the mesh.
Used for dimensioning arrays. Can be increased if needed.""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_TRACERS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 20,
                fr = """Nombre maximal de traceurs.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
                ang = """Maximal number of tracers.
Used for dimensioning arrays. Can be increased if needed.""",
            ),
#           -----------------------------------
            VECTOR_LENGTH = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Longueur du vecteur pour les machines vectorielles.""",
                ang = """Vector length on vector machines.""",
            ),
#           -----------------------------------
            SECURITY_COEFFICIENT_FOR_SCARACT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Coefficient de securite pour l''allocation de memoire pour
\telfile{SCARACT}.""",
                ang = """Security coefficient for memory allocation for \telfile{SCARACT}.""",
            ),
        ),
    ),
#   -----------------------------------
    RESTART = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        COMPUTATION_CONTINUED = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Determine si le calcul en cours est independant de tout autre
resultat ou est une reprise effectuee a partir du resultat d''un calcul
precedent.
\begin{itemize}
\item NON : Il s''agit du premier passage pour ce calcul et il est
necessaire de definir un jeu complet de conditions initiales
\item OUI : Il s''agit d''une reprise de calcul :
les conditions initiales sont constituees par le dernier pas de
temps du \telkey{FICHIER DU CALCUL PRECEDENT} du fichier des parametres
utilise pour soumettre le calcul.
\end{itemize}
Par contre, l''ensemble des donnees du fichier des parametres
peuvent etre redefinies, ce qui offre la possibilite de changer
par exemple, le pas de temps, le modele de turbulence, le
frottement, d''ajouter ou retirer un traceur\ldots\\
De meme, il est necessaire de definir des conditions aux limites
(sous-programme \telfile{BORD} ou valeurs placees dans le fichier des
parametres), qui peuvent egalement etre modifiees.""",
            ang = """Determines whether the computation under way is independent
result or is following an earlier result.
\begin{itemize}
\item NO: It is the first run for this computation and a whole set of
initial conditions should be defined,
\item YES: It follows a former computation:
the initial conditions consist in the last time step of the
\telkey{PREVIOUS COMPUTATION FILE} defined in the steering file
used for submitting the computation.
\end{itemize}
All the data from the steering file may be defined once again, which
provides an opportunity to change, for example, the time step,
the turbulence model, the friction, to add or remove a tracer\ldots\\
It is also possible to define new boundary conditions
(in the subroutine \telfile{BORD} or values defined
in the steering file).""",
        ),
#       -----------------------------------
        b_COMPUTATION_CONTINUEDG = BLOC(condition="COMPUTATION_CONTINUED == True",
#       -----------------------------------
#           -----------------------------------
            PREVIOUS_COMPUTATION_FILE_FORMAT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DU CALCUL PRECEDENT}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{PREVIOUS COMPUTATION FILE}.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            PREVIOUS_COMPUTATION_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom d''un fichier contenant les resultats d''un calcul precedent
realise sur le meme maillage et dont le dernier pas de temps enregistre
va fournir les conditions initiales pour une suite de calcul.""",
                ang = """Name of a file containing the results of an earlier computation
which was made on the same mesh. The last recorded time step will
provide the initial conditions for the new computation.""",
            ),
#           -----------------------------------
            RECORD_NUMBER_FOR_RESTART = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """En cas de suite de calcul, numero de l''enregistrement
de depart dans le fichier du calcul precedent. 0 signifie
que l''on prend le dernier enregistrement.""",
                ang = """In case of \telkey{COMPUTATION CONTINUED} = YES, record number
to start from in the \telkey{PREVIOUS COMPUTATION FILE}.
0 means the last record is taken.""",
            ),
        ),
#       -----------------------------------
        INITIAL_TIME_SET_TO_ZERO = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Remet le temps a zero en cas de suite de calcul.""",
            ang = """Initial time set to zero in case of restart.""",
        ),
    ),
#   -----------------------------------
    OUTPUT_FILES = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        RESULTS_FILES = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            NUMBER_OF_FIRST_TIME_STEP_FOR_GRAPHIC_PRINTOUTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Determine le nombre de pas de temps a partir duquel debute
l''ecriture des resultats dans le \telkey{FICHIER DES RESULTATS}.""",
                ang = """Determines the number of time steps after which the results
are first written into the \telkey{RESULTS FILE}.""",
            ),
#           -----------------------------------
            GRAPHIC_PRINTOUT_PERIOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode en nombre de pas de temps d''impression des
\telkey{VARIABLES POUR LES SORTIES GRAPHIQUES}
(voir ce mot-cle) dans le \telkey{FICHIER DES RESULTATS}.""",
                ang = """Determines, in number of time steps, the printout period for the
\telkey{VARIABLES FOR GRAPHIC PRINTOUTS}
in the \telkey{RESULTS FILE}.""",
            ),
#           -----------------------------------
            VARIABLES_FOR_GRAPHIC_PRINTOUTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["velocity along x axis (m/s)","velocity along y axis (m/s)","wave celerity (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","Froude number","scalar flowrate of fluid (m2/s)","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","tracer* etc.","tracer1 etc.","tracer2 etc.","tracer3 etc.","tracer4 etc.","tracer5 etc.","tracer6 etc.","tracer7 etc.","tracer8 etc.","tracer9 etc.","tracer10 etc.","tracer11 etc.","tracer12 etc.","tracer13 etc.","tracer14 etc.","tracer15 etc.","tracer16 etc.","tracer17 etc.","tracer18 etc.","tracer19 etc.","tracer** etc.","All the tracers 1 to 9.","All the tracers 10 to 19.","turbulent kinetic energy in k-epsilon model (J/kg)","dissipation of turbulent energy (W/kg)","turbulent viscosity (m2/s)","flowrate along x axis (m2/s)","flowrate along y axis (m2/s)","scalar velocity (m/s)","wind along x axis (m/s)","wind along y axis (m/s)","air pressure (Pa)","friction coefficient","drift along x (m)","drift along y (m)","Courant number","supplementary variable N","supplementary variable O","supplementary variable R","supplementary variable Z","maximum elevation (m)","time of maximum elevation (s)","maximum velocity (m/s)","time of maximum velocity (s)","friction velocity (m/s)","gradient 1, etc.","TAU_S","1/R (1/m)","OMEGA","distance to the closest wall (m)","reference level for Nestor (m)"],
                defaut = ["velocity along x axis (m/s)","velocity along y axis (m/s)","water depth (m)","bottom elevation (m)"],
                fr = """Noms des variables ecrites dans FICHIER DES RESULTATS.
Chaque variable est representee par une lettre (choix des
separateurs libre). Les possibilites sont :
\begin{itemize}
\item U : vitesse suivant $x$ (m/s) ;
\item V : vitesse suivant $y$ (m/s) ;
\item C : celerite des ondes (m/s) ;
\item H : hauteur d eau (m) ;
\item S : cote de surface libre (m) ;
\item B : cote du fond (m) ;
\item F : nombre de Froude ;
\item Q : debit scalaire du fluide (m$^2$/s) ;
\item Tn: traceur, avec n le numero du traceur ;
\item K : energie turbulente du modele $k-\epsilon$ (J/kg) ;
\item E : dissipation de l energie turbulente (W/kg) ;
\item D : viscosite turbulente (m$^2$/s) ;
\item I : debit suivant $x$ (m$^2$/s) ;
\item J : debit suivant $y$ (m$^2$/s) ;
\item M : vitesse scalaire (m/s) ;
\item X : vent suivant $x$ (m/s) ;
\item Y : vent suivant $y$ (m/s) ;
\item P : pression atmospherique (Pa) ;
\item W : coefficient de frottement sur le fond ;
\item A : derive en $x$ (m) ;
\item G : derive en $y$ (m) ;
\item L : nombre de Courant ;
\item MAXZ : cote maximum (m) ;
\item TMXZ : temps de la cote maximum (s) ;
\item MAXV : vitesse maximum (m/s) ;
\item TMXV : temps de la vitesse maximum (s) ;
\item US : vitesse de frottement (m/s) ;
\item Gn: gradient differencie, avec n le numero de reference du
gradient ;
\item TAU\_S : TAU\_S ;
\item 1/R : 1/R (1/m) ;
\item OMEGA : OMEGA ;
\item WDIST : distance au mur le plus proche (m) ;
\item ZRL : niveau de reference pour Nestor.
\end{itemize}
4 champs libres sont utilisables pour ecrire des variables creees
dans le fichier des resultats.
Ces variables doivent etre calculees dans le PRERES\_TELEMAC2D et leur
nom ecrit dans NOMVAR\_TELEMAC2D.
Ces 7 champs sont :
N, O, R, Z qui correspondent aux tableaux PRIVE(1,1) a PRIVE(1,4).
Contrairement aux variables precedentes, celles-ci sont conservees
dans tout le programme et peuvent donc etre reutilisees. Dans ce
dernier cas ne pas oublier de donner une taille suffisante au tableau
PRIVE.
Il est ainsi possible de limiter la taille des fichiers de resultats
pour de gros calculs.
Cependant, en cas de reprise de calcul, le code doit disposer dans le
fichier des resultats des informations necessaires a sa poursuite :
\begin{itemize}
\item vitesses U et V ;
\item hauteur d eau H ;
\item cote du fond B.
\end{itemize}
Toutefois, TELEMAC-2D peut recalculer certaines des variables a
partir d autres.""",
                ang = """Names of variables which will be written in the results file. Each
variable is represented by a letter (free separators).
The possible choices are:
\begin{itemize}
\item U: velocity along $x$ axis (m/s),
\item V: velocity along $y$ axis (m/s),
\item C: wave celerity (m/s),
\item H: water depth (m),
\item S: free surface elevation (m),
\item B: bottom elevation (m),
\item F: Froude number,
\item Q: scalar flowrate of fluid (m$^2$/s),
\item Tn: tracer, with n the tracer number,
\item K: turbulent kinetic energy in $k-\epsilon$ model (J/kg),
\item E: dissipation of turbulent energy (W/kg),
\item D: turbulent viscosity (m$^2$/s),
\item I: flowrate along $x$ axis (m$^2$/s),
\item J: flowrate along $y$ axis (m$^2$/s),
\item M: scalar velocity (m/s),
\item X: wind along $x$ axis (m/s),
\item Y: wind along $y$ axis (m/s),
\item P: air pressure (Pa),
\item W: friction coefficient,
\item A: drift along $x$ (m),
\item G: drift along $y$ (m),
\item L: Courant number,
\item MAXZ : maximum elevation (m),
\item TMXZ : time of maximum elevation (s),
\item MAXV : maximum velocity (m/s),
\item TMXV : time of maximum velocity (s),
\item US : friction velocity (m/s),
\item Gn: differentiated gradient, with n the gradient reference number,
\item TAU\_S : TAU\_S,
\item 1/R : 1/R (1/m),
\item OMEGA : OMEGA,
\item WDIST : distance to the closest wall (m),
\item ZRL : reference level for Nestor (m).
\end{itemize}
4 other variables are also available to the user
to write created variables results.
These user variables should be computed in
\telfile{PRERES\_TELEMAC2D} subroutine and their name
should be written in \telfile{NOMVAR\_TELEMAC2D} subroutine.
These seven variables are as follows:
N, O, R, Z which correspond to arrays \telfile{PRIVE(1,1)} up to
\telfile{PRIVE(1,4)}.
Unlike the previous variables, they are kept throughout the program,
so that they can be used again.
In the latter case, do not forget to provide the
array \telkey{PRIVE} with sufficiently large dimensions in
the FORTRAN FILE.
The size of the \telkey{RESULTS FILE} can be limited with this keyword.
However, if a computation must be continued, the
\telkey{RESULTS FILE} should contain the appropriate information
for running the code,i.e.:
\begin{itemize}
\item velocities $U$ and $V$,
\item water depth $H$,
\item bottom elevation $B$.
\end{itemize}
\telemac{2d} can compute some of these variables from others.""",
            ),
#           -----------------------------------
            NUMBER_OF_PRIVATE_ARRAYS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Nombre de tableaux mis a disposition de l utilisateur.""",
                ang = """Number of arrays for own user programming.""",
            ),
#           -----------------------------------
            NAMES_OF_PRIVATE_VARIABLES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                fr = """Noms des variables privees en 32 caracteres, 16 pour le nom
16 pour l''unite. Elles correspondent au bloc \telfile{PRIVE}
et peuvent etre lues dans le \telkey{FICHIER DE GEOMETRIE} si elles
y sont presentes avec leur nom.""",
                ang = """Name of private variables in 32 characters, 16 for the name,
16 for the unit.
They are stored in the block \telfile{PRIVE} and can be read
in the \telkey{GEOMETRY FILE} if they are here with their name.""",
            ),
#           -----------------------------------
            RESULTS_FILE_FORMAT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIERS DE RESULTATS}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{RESULTS FILE}. Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            RESULTS_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les resultats du
calcul avec la periodicite donnee par le mot cle \telkey{PERIODE POUR
LES SORTIES GRAPHIQUES}.""",
                ang = """Name of the file into which the computation results
are written with a periodicity given by the keyword
\telkey{GRAPHIC PRINTOUT PERIOD}.""",
            ),
#           -----------------------------------
            RESULT_FILE_IN_LONGITUDE_LATITUDE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Donne les coordonnees dans le fichier resultats en longitude-latitude
si le fichier geo est aussi donnee en longitude-latitude.""",
                ang = """Gives the coordinates of the result file in longitude-latitude
if the geometry file is also given in longitude-latitude.""",
            ),
#           -----------------------------------
            BINARY_RESULTS_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['BIN','SERAFIN','SERAFIND','MED'],
                defaut = 'BIN',
                fr = """Format du \telkey{FICHIER DE RESULTATS BINAIRE}.
Les valeurs possibles sont :
\begin{itemize}
\item BIN     : format binaire standard ;
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{BINARY RESULTS FILE}.
Possible values are:
\begin{itemize}
\item BIN     : standard binary format,
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            BINARY_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Fichier de resultats code en binaire mis a la disposition
de l''utilisateur.""",
                ang = """Additional binary-coded result file available to the user.""",
            ),
#           -----------------------------------
            FORMATTED_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Fichier de resultats formate mis a la disposition de
l''utilisateur.""",
                ang = """Formatted file of results available to the user.""",
            ),
        ),
#       -----------------------------------
        CONTROL_SECTION = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            CONTROL_SECTIONS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I', min=0, max='**',
                fr = """Couples de points (numeros globaux dans le maillage) entre
lesquels les debits instantanes et cumules seront donnes.""",
                ang = """Couples of points (global numbers
in the mesh) defining sections
where the instantaneous and cumulated discharges will be given.""",
            ),
#           -----------------------------------
            PRINTING_CUMULATED_FLOWRATES = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Impression du flux cumule a travers les sections de controle.""",
                ang = """Printing the cumulated flowrates through control sections.""",
            ),
#           -----------------------------------
            COMPATIBLE_COMPUTATION_OF_FLUXES = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Flux a travers les sections de controle, calcul compatible
avec l''impermeabilite sous forme faible.""",
                ang = """Flowrates through control sections, computation compatible
with the weak formulation of no-flux boundary condition.""",
            ),
#           -----------------------------------
            SECTIONS_INPUT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier des sections de controle, partitione.""",
                ang = """Sections input file, partitioned.""",
            ),
#           -----------------------------------
            SECTIONS_OUTPUT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Fichier de sortie des sections de controle, ecrit par le maitre.""",
                ang = """Sections output file, written by the master.""",
            ),
#           -----------------------------------
            FLUXLINE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Utiliser fluxline pour le calcul des flux sur les lignes.""",
                ang = """Use Fluxline to compute flux over lines.""",
            ),
#           -----------------------------------
            FLUXLINE_INPUT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de fluxline, avec des donnees sur les sections.""",
                ang = """Name of the fluxline file, with data on cross-sections.""",
            ),
        ),
#       -----------------------------------
        LISTING = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            NUMBER_OF_FIRST_TIME_STEP_FOR_LISTING_PRINTOUTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Determine le numero de pas de temps a partir duquel debute
l''ecriture des resultats dans le listing.""",
                ang = """Determines the number of time steps after which the results
are first written into the listing.""",
            ),
#           -----------------------------------
            LISTING_PRINTOUT_PERIOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode en nombre de pas de temps d''impression des
\telkey{VARIABLES A IMPRIMER} (voir ce mot-cle).
Pour la mise au point, il faut
savoir que la sortie des resultats est effectuee systematiquement sur le
listing (CAS.SORTIE sur station de travail).""",
                ang = """Determines, in number of time steps, the printout period of the
\telkey{VARIABLES TO BE PRINTED}.
The results are systematically printed out on
the listing file (file CAS.SORTIE at the workstation).""",
            ),
#           -----------------------------------
            LISTING_FOR_PRINTOUT_PERIOD = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode en nombre de pas de temps d''impression
des \telkey{VARIABLES A IMPRIMER} (voir ce mot-cle).
Pour la mise au point,
il faut savoir que la sortie des resultats est effectuee
systematiquement sur le listing
(CAS.SORTIE sur station de travail).
A la priorite sur \telkey{PERIODE DE SORTIE LISTING}.""",
                ang = """Determines, in number of time steps, the printout period of
the \telkey{VARIABLES TO BE PRINTED}.
The results are systematically printed out on the listing file
(file CAS.SORTIE at the workstation),
Has priority before \telkey{LISTING PRINTOUT PERIOD}.""",
            ),
#           -----------------------------------
            LISTING_PRINTOUT = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Sortie des resultats sur support papier ou a l ecran.
Si l''on met NON le listing ne contient que l''entete et la mention
FIN NORMALE DU PROGRAMME.
Commande a eviter.""",
                ang = """Result printout on hard copy.
When NO is selected, the listing only includes the heading and the
phrase "NORMAL END OF PROGRAM".
In addition, the options \telkey{MASS-BALANCE} and
\telkey{VALIDATION} are inhibited. Not recommended for use.""",
            ),
#           -----------------------------------
            VARIABLES_TO_BE_PRINTED = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                intoSug = ["velocity along x axis (m/s)","velocity along y axis (m/s)","wave celerity (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","Froude number","scalar flowrate of fluid (m2/s)","tracer 1, etc.","turbulent kinetic energy in k-epsilon model (J/kg)","dissipation of turbulent energy (W/kg)","turbulent viscosity (m2/s)","flowrate along x axis (m2/s)","flowrate along y axis (m2/s)","scalar velocity (m/s)","wind along x axis (m/s)","wind along y axis (m/s)","air pressure (Pa)","friction coefficient","drift along x (m)","drift along y (m)","Courant number","supplementary variable N","supplementary variable O","supplementary variable R","supplementary variable Z","maximum elevation (m)","time of maximum elevation (s)","maximum velocity (m/s)","time of maximum velocity (s)","friction velocity (m/s)","gradient 1, etc.","TAU_S","1/R (1/m)","OMEGA","distance to the closest wall (m)","reference level for Nestor (m)"],
                defaut = '',
                fr = """Nom des variables que l''utilisateur desire ecrire a l''ecran. Memes
possibilites que pour les sorties graphiques.""",
                ang = """Name of the variables that the user wants printed on screen.
Same values available as graphical outputs.""",
            ),
#           -----------------------------------
            MASS_BALANCE = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Determine si l''on effectue ou non le bilan de masse
sur le domaine.
Cette procedure calcule a chaque pas de temps :
\begin{itemize}
\item les flux aux entrees et sorties du domaine ;
\item le flux global a travers l''ensemble des parois du domaine
(liquides ou solides) ;
\item l''erreur relative sur la masse pour ce pas de temps.
\end{itemize}
En fin de listing, on trouve l''erreur relative sur la masse pour
l''ensemble du calcul.
Il ne s''agit que d''un calcul indicatif car il n''existe pas
d''expression compatible du debit en formulation c,u,v.""",
                ang = """Determines whether a check of the mass-balance over
the domain is mader or not.
This procedures computes the following at each time step:
\begin{itemize}
\item the domain inflows and outflows,
\item the overall flow across all the boundaries,
\item the relative error in the mass for that time step.
\end{itemize}
The relative error in the mass over the whole computation can be found
at the end of the listing.""",
            ),
#           -----------------------------------
            INFORMATION_ABOUT_SOLVER = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Donne a chaque pas de temps le nombre d''iterations necessaires
a la convergence du solveur de l''etape de propagation.""",
                ang = """If YES, prints the number of iterations
that have been necessary
to get the solution of the linear system.""",
            ),
#           -----------------------------------
            LIST_OF_POINTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I', min=0, max='**',
                fr = """Liste de points remarquables pour les impressions.""",
                ang = """List of remarkable points for printouts.""",
            ),
#           -----------------------------------
            NAMES_OF_POINTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                fr = """Noms des points remarquables pour les impressions.""",
                ang = """Names of remarkable points for printouts.""",
            ),
        ),
#       -----------------------------------
        FOURIER = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            FOURIER_ANALYSIS_PERIODS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Liste des periodes que l''on veut analyser.""",
                ang = """List of periods to be analysed.""",
            ),
#           -----------------------------------
            TIME_RANGE_FOR_FOURIER_ANALYSIS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min= 2, max= 2,
                defaut = [0.,0.],
                fr = """Pour le calcul du marnage et de la phase de la maree.""",
                ang = """For computing tidal range and phase of tide.""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
HYDRO = PROC(nom= "HYDRO",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    BOUNDARY_CONDITIONS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        PRESCRIBED_ELEVATIONS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', max='**',
            fr = """Valeurs des cotes imposees aux frontieres liquides.
Lire la partie du mode d''emploi consacree aux conditions aux limites.""",
            ang = """Values of the elevations prescribed at open boundaries.
The section about boundary conditions is to be read in the manual.""",
        ),
#       -----------------------------------
        PRESCRIBED_FLOWRATES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', max='**',
            fr = """Valeurs des debits imposes aux frontieres liquides.
Lire la partie du mode d''emploi consacree aux conditions aux limites.""",
            ang = """Values of the flowrates prescribed at open boundaries.
The section about boundary conditions is to be read in the manual.""",
        ),
#       -----------------------------------
        PRESCRIBED_VELOCITIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', max='**',
            fr = """Valeurs des vitesses imposees aux frontieres liquides.
Lire la partie du mode d''emploi consacree aux conditions aux limites.""",
            ang = """Values of the magnitudes of velocity prescribed at open boundaries.
Refer to the section dealing with the boundary conditions.""",
        ),
    ),
#   -----------------------------------
    BOUNDARY_CONDITIONS_OTHERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        STAGE_DISCHARGE_CURVES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', max='**',
            into = ["no","Z(Q)","Q(Z)"],
            fr = """Indique si une courbe de tarage doit etre utilisee pour une frontiere
(une valeur par frontiere liquide) :
\begin{itemize}
\item 0 : non ;
\item 1 : Z(Q) ;
\item 2 : Q(Z).
\end{itemize}""",
            ang = """Says if a discharge-elevation curve must be used for a given boundary
(one value per open boundary):
\begin{itemize}
\item 0: no,
\item 1: Z(Q),
\item 2: Q(Z).
\end{itemize}""",
        ),
#       -----------------------------------
        b_STAGE_DISCHARGE_CURVESG = BLOC(condition="STAGE_DISCHARGE_CURVES != 'no'",
#       -----------------------------------
#           -----------------------------------
            STAGE_DISCHARGE_CURVES_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les courbes de tarage.""",
                ang = """Name of the file containing stage-discharge curves.""",
            ),
        ),
#       -----------------------------------
        VELOCITY_PROFILES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["constant normal profile","u and v given in the conlim file","normal velocity given in ubor in the conlim file","velocity proportional to square root of depth","velocity proportional to square root of depth, variant"],
            fr = """Permet de specifier le type de profil horizontal de vitesse.
Les choix possibles sont :
\begin{itemize}
\item 1 : profil normal constant ;
\item 2 : $u$ et $v$ donnes dans le
\telkey{FICHIER DES CONDITIONS AUX LIMITES} ;
\item 3 : vitesse normale donnee dans \telfile{UBOR} dans le
\telkey{FICHIER DES CONDITIONS AUX LIMITES} ;
\item 4 : profil en $\sqrt{h}$ ;
\item 5 : comme 4 mais hauteur virtuelle calculee avec
la surface libre la plus basse de la frontiere.
\end{itemize}""",
            ang = """Specifies the type of horizontal profile of velocities.
Possible choices are:
\begin{itemize}
\item 1: constant normal profile,
\item 2: $u$ and $v$ given in the
\telkey{BOUNDARY CONDITION FILE},
\item 3: normal velocity given in \telfile{UBOR} in the
\telkey{BOUNDARY CONDITION FILE},
\item 4: $\sqrt{h}$ profile,
\item 5: like 4 but virtual depth based on
the lowest elevation of the boundary.
\end{itemize}""",
        ),
#       -----------------------------------
        OPTION_FOR_LIQUID_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', max='**',
            into = ["classical","Thompson method based on characteristics"],
            fr = """On donne un entier par frontiere liquide.
Les choix possibles sont :
\begin{itemize}
\item 1 : conditions aux limites classiques ;
\item 2 : methode de Thompson avec calcul de caracteristiques.
\end{itemize}""",
            ang = """One integer per liquid boundary is given.
Possible choices are:
\begin{itemize}
\item 1: classical boundary conditions,
\item 2: Thompson method based on characteristics.
\end{itemize}""",
        ),
#       -----------------------------------
        LIQUID_BOUNDARIES_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier de variations en temps des conditions aux limites.""",
            ang = """File containing the variations in time of boundary conditions.""",
        ),
#       -----------------------------------
        ELEMENTS_MASKED_BY_USER = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si OUI, remplir le sous-programme \telfile{USER\_MASKOB}.""",
            ang = """If YES, fill in the subroutine \telfile{USER\_MASKOB}.""",
        ),
#       -----------------------------------
        b_ELEMENTS_MASKED_BY_USERG = BLOC(condition="ELEMENTS_MASKED_BY_USER == True",
#       -----------------------------------
#           -----------------------------------
            Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#           -----------------------------------
                defaut = "Rewrite subroutine maskob"),
        ),
    ),
#   -----------------------------------
    PHYSICAL_PARAMETERS_HYDRO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        FRICTION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            MAXIMUM_NUMBER_OF_FRICTION_DOMAINS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 10,
                fr = """Nombre maximal de zones pouvant etre definies pour le
frottement. Peut etre augmente si necessaire.""",
                ang = """Maximal number of zones defined for the friction.
Could be increased if needed.""",
            ),
#           -----------------------------------
            FRICTION_DATA = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Logique qui dit si les lois de frottements sont definies par zone.""",
                ang = """Logical to say if friction laws are defined by area.""",
            ),
#           -----------------------------------
            FRICTION_DATA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de donnees pour le frottement.
Voir le manuel utilisateur de \telemac{2d} pour la description.""",
                ang = """Friction data file name.
See the \telemac{2d} user manual for its description.""",
            ),
#           -----------------------------------
            LAW_OF_BOTTOM_FRICTION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO FRICTION","HAALAND","CHEZY","STRICKLER","MANNING","NIKURADSE"],
                fr = """Selectionne le type de formulation utilisee pour le calcul
du frottement sur le fond.
Les lois possibles sont les suivantes (cf. Note de principe) :
\begin{itemize}
\item 0 : pas de frottement sur le fond ;
\item 1 : formule de Haaland ;
\item 2 : formule de Chezy ;
\item 3 : formule de Strickler ;
\item 4 : formule de Manning ;
\item 5 : formule de Nikuradse.
\end{itemize}""",
                ang = """Selects the type of formulation used for the bottom friction.
The possible laws are as follows (refer to the Principle note):
\begin{itemize}
\item 0: no friction against bottom,
\item 1: Haaland''s formula,
\item 2: Chezy''s formula,
\item 3: Strickler''s formula,
\item 4: Manning''s formula,
\item 5: Nikuradse''s formula.
\end{itemize}""",
            ),
#           -----------------------------------
            b_LAW_OF_BOTTOM_FRICTIONG = BLOC(condition="LAW_OF_BOTTOM_FRICTION != 'NO FRICTION'",
#           -----------------------------------
#               -----------------------------------
                FRICTION_COEFFICIENT = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 50.,
                    fr = """Fixe la valeur du coefficient de frottement pour la
formulation choisie.
Attention : la signification de ce chiffre varie suivant la formule
choisie :
\begin{itemize}
\item 1 : coefficient lineaire ;
\item 2 : coefficient de Chezy ;
\item 3 : coefficient de Strickler ;
\item 4 : coefficient de Manning ;
\item 5 : hauteur de rugosite de Nikuradse.
\end{itemize}""",
                    ang = """Sets the value of the friction coefficient for the selected
formulation. It is noteworthy that the meaning of this figure changes
according to the selected formula (Chezy, Strickler, etc.):
\begin{itemize}
\item 1: linear coefficient,
\item 2: Chezy coefficient,
\item 3: Strickler coefficient,
\item 4: Manning coefficient,
\item 5: Nikuradse grain size.
\end{itemize}""",
                ),
            ),
#           -----------------------------------
            MANNING_DEFAULT_VALUE_FOR_COLEBROOK_WHITE_LAW = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.02,
                fr = """Valeur par defaut du coefficient de Manning pour la loi de frottement
de Colebrook-White (loi numero 7).""",
                ang = """Manning default value for the friction law of Colebrook-White
(law number 7).""",
            ),
#           -----------------------------------
            DEPTH_IN_FRICTION_TERMS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["nodal","average"],
                defaut = "nodal",
                fr = """Valeurs possibles :
\begin{itemize}
\item 1 : nodale ;
\item 2 : moyenne.
\end{itemize}""",
                ang = """Possible choices:
\begin{itemize}
\item 1: nodal,
\item 2: average.
\end{itemize}""",
            ),
#           -----------------------------------
            VEGETATION_FRICTION = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Calcul du frottement du a la vegetation non submergee.""",
                ang = """Friction calculation of the non-submerged vegetation.""",
            ),
#           -----------------------------------
            LAW_OF_FRICTION_ON_LATERAL_BOUNDARIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO FRICTION","HAALAND","CHEZY","STRICKLER","MANNING","NIKURADSE","LOG LAW","COLEBROOK-WHITE"],
                defaut = "NO FRICTION",
                fr = """Selectionne le type de formulation utilisee pour le calcul
du frottement sur les parois laterales.
Les lois possibles sont les suivantes (cf. Note de principe) :
\begin{itemize}
\item 0 : pas de frottement ;
\item 1 : lineaire ;
\item 2 : Chezy ;
\item 3 : Strickler ;
\item 4 : Manning ;
\item 5 : formule de Nikuradse ;
\item 6 : loi en log ;
\item 7 : Colebrook-White.
\end{itemize}""",
                ang = """Selects the type of formulation used
for the friction on lateral boundaries.
The possible laws are as follows (refer to the Principle note):
\begin{itemize}
\item 0: no friction,
\item 1 : linear,
\item 2 : Chezy,
\item 3 : Strickler,
\item 4 : Manning,
\item 5: Nikuradse''s formula,
\item 6 : log law,
\item 7 : Colebrook-White.
\end{itemize}""",
            ),
#           -----------------------------------
            b_LAW_OF_FRICTION_ON_LATERAL_BOUNDARIESG = BLOC(condition="LAW_OF_FRICTION_ON_LATERAL_BOUNDARIES != 'NO FRICTION'",
#           -----------------------------------
#               -----------------------------------
                ROUGHNESS_COEFFICIENT_OF_BOUNDARIES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 100.,
                    fr = """Fixe la valeur du coefficient de frottement sur les frontieres
solides avec un regime turbulent rugueux sur les bords du domaine.
Meme convention que pour le coefficient de frottement:
\begin{itemize}
\item 1 : non programme ;
\item 2 : coefficient de Chezy ;
\item 3 : coefficient de Strickler ;
\item 4 : coefficient de Manning ;
\item 5 : hauteur de rugosite de Nikuradse.
\end{itemize}""",
                    ang = """Sets the value of the friction coefficient of the solid
boundary with the bed roughness option. Same meaning than friction
coefficient:
\begin{itemize}
\item 1: not implemented,
\item 2: Chezy coefficient,
\item 3: Strickler coefficient,
\item 4: Manning coefficient,
\item 5: Nikuradse grain size.
\end{itemize}""",
                ),
            ),
#           -----------------------------------
            DEFINITION_OF_ZONES = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Declenche l''appel au sous-programme \telfile{USER\_DEF\_ZONES},
pour donner un numero de zone a chaque point.""",
                ang = """Triggers the call to \telfile{USER\_DEF\_ZONES} subroutine
to give a zone number to every point.""",
            ),
#           -----------------------------------
            b_DEFINITION_OF_ZONESG = BLOC(condition="DEFINITION_OF_ZONES == True",
#           -----------------------------------
#               -----------------------------------
                Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#               -----------------------------------
                    defaut = "Rewrite subroutine def_zones"),
            ),
#           -----------------------------------
            ZONES_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier des zones avec sur chaque ligne :\\
numero de point  numero de zone.""",
                ang = """Zones file, with on every line:\\
point number   zone number.""",
            ),
        ),
#       -----------------------------------
        METEOROLOGY = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            WIND = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Prise en compte ou non des effets du vent.""",
                ang = """Determines whether the wind effects are to be taken into
account or not.""",
            ),
#           -----------------------------------
            b_WINDG = BLOC(condition="WIND == True",
#           -----------------------------------
#               -----------------------------------
                WIND_VELOCITY_ALONG_X = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 0.,
                    fr = """Composante de la vitesse du vent suivant
l''axe des $x$ (m/s), si constante.""",
                    ang = """Wind velocity, component along $x$ axis (m/s), if constant.""",
                ),
#               -----------------------------------
                WIND_VELOCITY_ALONG_Y = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 0.,
                    fr = """Composante de la vitesse du vent suivant
l''axe des $y$ (m/s), si constante.""",
                    ang = """Wind velocity, component along $y$ axis (m/s), if constant.""",
                ),
#               -----------------------------------
                THRESHOLD_DEPTH_FOR_WIND = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.,
                    fr = """Retire la force due au vent dans les petites profondeurs
plus petites que cette valeur.""",
                    ang = """Wind is not taken into account for depths smaller
than this value.""",
                ),
#               -----------------------------------
                COEFFICIENT_OF_WIND_INFLUENCE = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.55E-6,
                    fr = """Fixe la valeur du coefficient d''entrainement du vent.
Voir le manuel utilisateur ou la note de principe
pour la valeur a donner.""",
                    ang = """Sets the value of the wind driving coefficient.
See the User Manual or the principle note
for the value to give.""",
                ),
#               -----------------------------------
                OPTION_FOR_WIND = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["constant in time and space","variable in time given by formated file","variable in time and space given by formated file"],
                    defaut = "constant in time and space",
                    fr = """Donne les options pour introduire le vent :
\begin{itemize}
\item 1 : constant en temps et en espace (donne par les mots cle
\telkey{VITESSE ET DIRECTION DU VENT}
) ;
\item 2 : variable en temps donne par fichier formate ;
\item 3 : variable en temps et en espace donne par fichier formate
ou un fichier binaire SERAFIN.
\end{itemize}""",
                    ang = """Gives the option for managing the wind:
\begin{itemize}
\item 1: constant in time and space, given by the keyword
\telkey{SPEED AND DIRECTION OF WIND},
\item 2: variable in time and constant in space, given by
\telkey{ASCII ATMOSPHERIC DATA FILE},
\item 3: variable in time and space, given by formatted file or by
a binary SERAFIN file.
\end{itemize}""",
                ),
#               -----------------------------------
                b_OPTION_FOR_WINDG = BLOC(condition="OPTION_FOR_WIND == 'constant in time and space'",
#               -----------------------------------
#                   -----------------------------------
                    SPEED_AND_DIRECTION_OF_WIND = SIMP(statut ='o',
#                   -----------------------------------
                        typ = 'R', min= 2, max= 2,
                        defaut = [0.,0.],
                        fr = """Donne la vitesse et la direction (en degres de 0 a 360,
0 etant $y$ = 0 et $x$ = +inf) du vent lorsqu ils sont constants
en temps et en espace (mot cle \telkey{OPTION DU VENT} = 1).""",
                        ang = """Gives the speed and direction (in degrees (from 0 to 360),
0 given $y$ = 0 anx $x$ = +infinity) when they are constant
in time and space (keyword \telkey{OPTION FOR WIND} = 1).""",
                    ),
                ),
#               -----------------------------------
                b_OPTION_FOR_WINDH = BLOC(condition="OPTION_FOR_WIND == 'variable in time given by formated file' or OPTION_FOR_WIND == 'variable in time and space given by formated file'",
#               -----------------------------------
#                   -----------------------------------
                    Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#                   -----------------------------------
                        defaut = "Give the formatted file 3"),
                ),
            ),
#           -----------------------------------
            COEFFICIENT_OF_WIND_INFLUENCE_VARYING_WITH_WIND_SPEED = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Si OUI, la valeur du coefficient d''entrainement du vent est
calculee en fonction de la vitesse du vent.
La valeur de \telkey{COEFFICIENT D''INFLUENCE DU VENT} est ecrasee.
Ancienne valeur par defaut = NON jusqu a la version V8P1.""",
                ang = """If YES, the value of the wind driving coefficient is computed
with respect to the wind velocity.
The value of \telkey{COEFFICIENT OF WIND INFLUENCE} is overwritten.
Old default value = NO until V8P1.""",
            ),
#           -----------------------------------
            AIR_PRESSURE = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Permet de decider si l''on prend ou non en compte l''influence
d''un champ de pression.""",
                ang = """Provided to decide whether the influence of an atmosphere
pressure field is taken into account or not.""",
            ),
#           -----------------------------------
            b_AIR_PRESSUREG = BLOC(condition="AIR_PRESSURE == True",
#           -----------------------------------
#               -----------------------------------
                VALUE_OF_ATMOSPHERIC_PRESSURE = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 100000.,
                    fr = """Donne la valeur de la pression atmospherique lorsqu''elle est
constante en temps et en espace. En Pa.""",
                    ang = """Gives the value of atmospheric pressure when it is constant
in time and space. In Pa.""",
                ),
            ),
#           -----------------------------------
            RAIN_OR_EVAPORATION = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Pour ajouter un apport ou une perte d''eau en surface.
Voir le mot-cle \telkey{PLUIE OU EVAPORATION EN MM PAR JOUR}.""",
                ang = """Enables to add or remove water at the free surface.
See the keyword \telkey{RAIN OR EVAPORATION IN MM PER DAY}.""",
            ),
#           -----------------------------------
            b_RAIN_OR_EVAPORATIONG = BLOC(condition="RAIN_OR_EVAPORATION == True",
#           -----------------------------------
#               -----------------------------------
                RAIN_OR_EVAPORATION_IN_MM_PER_DAY = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 0.0,
                    fr = """Pour ajouter un apport ou une perte d''eau en surface.""",
                    ang = """To add or remove water at the free surface.""",
                ),
            ),
#           -----------------------------------
            RAINFALL_RUNOFF_MODEL = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["No infiltration","CN runoff model"],
                defaut = "No infiltration",
                fr = """Option pour le modele pluie-debit. Les options disponibles sont :
\begin{itemize}
\item 0 : Pas d infiltration (fonction de base) ;
\item 1 : Modele CN (Curve Number du SCS).
\end{itemize}""",
                ang = """Option for the rainfall-runoff model. Available options are:
\begin{itemize}
\item 0: No infiltration (basic function),
\item 1: CN runoff model (Curve Number method of the SCS).
\end{itemize}""",
            ),
#           -----------------------------------
            ANTECEDENT_MOISTURE_CONDITIONS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 2,
                fr = """Donne les conditions d humidite precedant un episode de pluie pour
le modele CN du SCS. Les options disponibles sont :
\begin{itemize}
\item 1 : conditions precedentes seches ;
\item 2 : conditions precedentes normales ;
\item 3 : conditions prcedentes mouillees.
\end{itemize}
Ce mot cle est uniquement utile pour le modele pluie-debit 1 (CN).""",
                ang = """Gives the antecedent moisture conditions before a rainfall
 event for the SCS CN runoff model. Available options are:
\begin{itemize}
\item 1: dry antecedent conditions,
\item 2: normal antecedent conditions,
\item 3: wet antecedent conditions.
\end{itemize}
This keyword is only usefull for runoff model 1 (SCS CN model).""",
            ),
#           -----------------------------------
            DURATION_OF_RAIN_OR_EVAPORATION_IN_HOURS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 1.E6,
                fr = """Donne la duree de la pluie en heures, par defaut pluie infinie.""",
                ang = """Gives the duration of the rain in hours,
default value is infinite.""",
            ),
#           -----------------------------------
            ASCII_ATMOSPHERIC_DATA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees ASCII contenant les informations
atmospheriques variables en temps.""",
                ang = """ASCII data file containing the atmospheric data varying in
time.""",
            ),
#           -----------------------------------
            BINARY_ATMOSPHERIC_DATA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees code en binaire contenant les informations
atmospheriques variables en temps et en espace sur le maillage.""",
                ang = """Binary-coded data file containing the atmospheric data varying in
time and space on the mesh.""",
            ),
#           -----------------------------------
            BINARY_ATMOSPHERIC_DATA_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER BINAIRE DE DONNEES ATMOSPHERIQUES}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{BINARY ATMOSPHERIC DATA FILE}.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            OPTION_FOR_INITIAL_ABSTRACTION_RATIO = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Donne le ratio entre pertes initiales IA et la retention potentielle
maximale S pour le modele pluie-debit SCS CN.
Les options disponibles sont :
\begin{itemize}
\item 1 : IA/S = 0.2 (methode standard) ;
\item 2 : IA/S = 0.05 (methode revisee, cf. Woodward, Hawkins et al.
2003.
\end{itemize}
Avec cette option les coefficients CN fournis en entree sont alors
automatiquement corriges, cf. manuel utilisateur).
Ce mot cle est uniquement utile pour le modele pluie-debit 1 (CN).""",
                ang = """Gives the ratio for Initial Abstraction to Maximal Potential
Retention S for the SCS CN runoff model. Available options are:
\begin{itemize}
\item 1: IA/S = 0.2 (standard method),
\item 2: IA/S = 0.05 (revised method, see Woodward, Hawkins et al. 2003.
\end{itemize}
With this option the CN values given in input are automatically
converted see user manual).
This keyword is only useful for runoff model 1 (SCS CN model).""",
            ),
#           -----------------------------------
            AIR_TEMPERATURE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 20.,
                fr = """Donne la valeur de la temperature de l air lorsqu elle est
constante en temps et en espace. En $^{\circ}$C.""",
                ang = """Gives the value of air temperature when it is constant
in time and space. In $^{\circ}$C.""",
            ),
#           -----------------------------------
            CLOUD_COVER = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 5.,
                fr = """Donne la valeur de la nebulosite lorsqu elle est
constante en temps et en espace. En Octas ou dizieme.""",
                ang = """Gives the value of cloud cover when it is constant
in time and space. In Octas or tenth.""",
            ),
#           -----------------------------------
            SOLAR_RADIATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 160.,
                fr = """Donne la valeur du rayonnement solaire lorsqu il est
constant en temps et en espace. En W/m$^2$.""",
                ang = """Gives the value of solar radiation when it is constant
in time and space. In W/m$^2$.""",
            ),
#           -----------------------------------
            VAPOROUS_PRESSURE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1000.,
                fr = """Donne la valeur de la pression de vapeur saturante lorsqu elle est
constante en temps et en espace. En Pa.""",
                ang = """Gives the value of vaporous pressure when it is constant
in time and space. In Pa.""",
            ),
        ),
#       -----------------------------------
        WAVE = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            WAVE_DRIVEN_CURRENTS = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Active la prise en compte des courants de houle.""",
                ang = """Wave driven currents are taken into account.""",
            ),
#           -----------------------------------
            b_WAVE_DRIVEN_CURRENTSG = BLOC(condition="WAVE_DRIVEN_CURRENTS == True",
#           -----------------------------------
#               -----------------------------------
                RECORD_NUMBER_IN_WAVE_FILE = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 1,
                    fr = """Numero d enregistrement dans le fichier des courants de houle.""",
                    ang = """Record number to be read in the wave driven currents file.""",
                ),
            ),
#           -----------------------------------
            WAVE_ENHANCED_FRICTION_FACTOR = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Active la prise en compte des interactions non-lineaires entre la
houle et les courants pour le calcul du courant de houle (cf OConnor and
Yoo, 1988, Coast Eng.12.).""",
                ang = """Wave friction enhancement for the calculation of the wave generated
longshore current (cf OConnor and Yoo, 1988, Coast Eng.12.).""",
            ),
        ),
#       -----------------------------------
        ESTIMATION = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            PARAMETER_ESTIMATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ['','FRICTION','FRICTION, STEADY'],
                defaut = '',
                fr = """Liste des parametres a estimer, choix :
\begin{itemize}
\item FROTTEMENT ;
\item FROTTEMENT, PERMANENT.
\end{itemize}""",
                ang = """List of parameter to be estimated, choices:
\begin{itemize}
\item FRICTION,
\item FRICTION, STEADY.
\end{itemize}""",
            ),
#           -----------------------------------
            COST_FUNCTION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["Computed with h,u,v","Computed with c,u,v"],
                defaut = "Computed with h,u,v",
                fr = """Choix possibles :
\begin{itemize}
\item calculee sur $h$, $u$, $v$ ;
\item calculee avec $c$, $u$ , $v$.
\end{itemize}""",
                ang = """Possible choices:
\begin{itemize}
\item computed with $h$, $u$, $v$,
\item computed with $c$, $u$, $v$.
\end{itemize}""",
            ),
#           -----------------------------------
            IDENTIFICATION_METHOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["list of tests","gradient simple","conjugate gradient","Lagrange interpolation"],
                defaut = "gradient simple",
                fr = """Choix possibles :
\begin{itemize}
\item 0 : plan d''experience ;
\item 1 : gradient simple ;
\item 2 : gradient conjugue ;
\item 3 : interpolation de Lagrange.
\end{itemize}""",
                ang = """Possible choices:
\begin{itemize}
\item 0: list of tests,
\item 1: gradient,
\item 2: conjugate gradient,
\item 3: Lagrangian interpolation.
\end{itemize}""",
            ),
#           -----------------------------------
            TOLERANCES_FOR_IDENTIFICATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min= 4, max= 4,
                defaut = [1.E-3,1.E-3,1.E-3,1.E-4],
                fr = """4 nombres : precision absolue sur $H$, $U$, $V$,
et precision relative sur la fonction cout.""",
                ang = """4 numbers: absolute precision on $H$, $U$, $V$,
and relative precision on the cost function.""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_ITERATIONS_FOR_IDENTIFICATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 20,
                fr = """Chaque iteration comprend au moins un
calcul direct et un calcul adjoint.""",
                ang = """Every iteration implies at least a direct and
an adjoint computation.""",
            ),
        ),
#       -----------------------------------
        SOURCES = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            ABSCISSAE_OF_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Valeurs des abscisses des sources de debit et de traceur.
La source sera placee au noeud du maillage le plus proche.""",
                ang = """Abscissae of sources of flowrate and/or tracer.
The source will be located at the nearest node in the mesh.""",
            ),
#           -----------------------------------
            ORDINATES_OF_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Valeurs des ordonnees des sources de debit et de traceur.
La source sera placee au noeud du maillage le plus proche.""",
                ang = """Ordinates of sources of flowrate and/or tracer.
The source will be located at the nearest node in the mesh.""",
            ),
#           -----------------------------------
            WATER_DISCHARGE_OF_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Specifie le debit de chaque source.
Un debit positif signifie qu''il s''agit d''un apport de fluide.""",
                ang = """Specifies the discharge for every source.
A positive discharge means that fluid is added.""",
            ),
#           -----------------------------------
            VELOCITIES_OF_THE_SOURCES_ALONG_X = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Vitesses du courant selon $x$ a chacune des sources.
Si elles ne sont pas donnees, on considere que la vitesse est
celle du courant.""",
                ang = """Velocities along $x$ at the sources.
If they are not given, the velocity of the flow at this location
is taken.""",
            ),
#           -----------------------------------
            VELOCITIES_OF_THE_SOURCES_ALONG_Y = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Vitesses du courant selon $y$ a chacune des sources.
Si elles ne sont pas donnees, on considere que la vitesse est
celle du courant.""",
                ang = """Velocities along $y$ at the sources.
If they are not given, the velocity of the flow at this location
is taken.""",
            ),
#           -----------------------------------
            TYPE_OF_SOURCES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Normal","Dirac"],
                defaut = "Normal",
                fr = """Definit comment les sources sont calculees :
\begin{itemize}
\item 1 : Source portee par une base elements finis ;
\item 2 : Source portee par une fonction de Dirac
(recommande quand il y a beaucoup de sources).
\end{itemize}""",
                ang = """Defines how the sources are computed:
\begin{itemize}
\item 1: Source term multiplied by a finite element basis,
\item 2: Source term multiplied by a Dirac function
(recommended with high numbers of sources).
\end{itemize}""",
            ),
#           -----------------------------------
            SOURCES_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les informations variables
en temps des sources.""",
                ang = """Name of the file containing time-dependent
information on sources.""",
            ),
#           -----------------------------------
            GLOBAL_NUMBERS_OF_SOURCE_NODES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I', min=0, max='**',
                fr = """Numeros globaux des noeuds du maillage sur lequels sont affectes des
points source.""",
                ang = """Global numbers of nodes in the mesh that correspond to source point
locations.""",
            ),
#           -----------------------------------
            SOURCE_REGIONS_DATA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees en ASCII contenant les coordonnees des polygones
qui contiennent les sources.""",
                ang = """ASCII data file containing sources informations: coordinates of the
polygons containing sources.""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_POINTS_FOR_SOURCES_REGIONS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 10,
                fr = """Nombre maximal de points pour definir les regions contenant les
sources.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
                ang = """Maximal number of points to define regions containing sources.
Used for dimensioning arrays. It can be increased if needed.""",
            ),
        ),
#       -----------------------------------
        WATER_QUALITY_INFO = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            WATER_QUALITY_PROCESS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Donne le numero du processus de qualite d''eau, defini
comme une combinaison multiplicative de nombres premiers
(2,3,5,7,11 et 13) avec les cas particuliers 0 et 1 :
\begin{itemize}
\item 0 : tous ;
\item 1 : rien ;
\item 2 : O2 ;
\item 3 : BIOMASS ;
\item 5 : EUTRO ;
\item 7 : MICROPOL ;
\item 11 : THERMIC ;
\item 17 : Loi de degradation ;
\item 19 : Court-cicuit temporaire pour les glaces.
\end{itemize}
Exemple: 110 = 2x5x11 activera O2, EUTRO et THERMIC ensemble.
On notera que AED2, pour l instant, n est pas disponible en 2D.""",
                ang = """Gives the water quality process number, defined as
a multiplicative combination of prime numbers (2,3,5,7,11, 17
and 19) with 0 and 1 having a special role:
\begin{itemize}
\item 0: all,
\item 1: none,
\item 2: O2,
\item 3: BIOMASS,
\item 5: EUTRO,
\item 7: MICROPOL,
\item 11: THERMIC,
\item 17: Degradation law,
\item 19: Ghost process for ice modelling.
\end{itemize}
Example: 110 = 2x5x11 activate O2, EUTRO and THERMIC together.
It is noted that AED2 is not available in 2D, for the time being.""",
            ),
        ),
#       -----------------------------------
        ADVANCED_PHY = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            WATER_DENSITY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1000.,
                fr = """Fixe la valeur de la masse volumique de l''eau.""",
                ang = """Sets the value of water density.""",
            ),
#           -----------------------------------
            GRAVITY_ACCELERATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 9.81,
                fr = """Fixe la valeur de l''acceleration de la pesanteur en m/s$^2$.""",
                ang = """Sets the value of the acceleration due to gravity in m/s$^2$.""",
            ),
#           -----------------------------------
            VERTICAL_STRUCTURES = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Prise en compte de la force de trainee de structures verticales
(il faut alors remplir le sous-programme \telfile{DRAGFO}).""",
                ang = """Drag forces from vertical structures are taken into account.
(subroutine \telfile{DRAGFO} must then be implemented).""",
            ),
#           -----------------------------------
            b_VERTICAL_STRUCTURESG = BLOC(condition="VERTICAL_STRUCTURES == True",
#           -----------------------------------
#               -----------------------------------
                Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#               -----------------------------------
                    defaut = "Fill the subroutine DRAGFO"),
            ),
        ),
    ),
#   -----------------------------------
    NUMERICAL_PARAMETERS_HYDRO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        EQUATIONS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ['SAINT-VENANT FE','SAINT-VENANT FV','BOUSSINESQ'],
            defaut = 'SAINT-VENANT FE',
            fr = """Choix des equations a resoudre :
\begin{itemize}
\item Saint-Venant Elements Finis ;
\item Saint-Venant Volumes Finis ;
\item Boussinesq 20 caracteres.
\end{itemize}""",
            ang = """Choice of equations to solve:
\begin{itemize}
\item Shallow Water Finite Elements,
\item Shallow Water Finite Volumes,
\item Boussinesq 20 Characters.
\end{itemize}""",
        ),
#       -----------------------------------
        b_EQUATIONSG = BLOC(condition="EQUATIONS == 'SAINT-VENANT FV'",
#       -----------------------------------
#           -----------------------------------
            FINITE_VOLUME_SCHEME = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Roe scheme","kinetic","Zokagoa scheme","Tchamen scheme","HLLC scheme","WAF scheme"],
                defaut = "kinetic",
                fr = """Choix possibles :
\begin{itemize}
\item 0 : schema de Roe ;
\item 1 : cinetique ;
\item 3 : schema de Zokagoa ;
\item 4 : schema de Tchamen ;
\item 5 : HLLC ;
\item 6 : WAF.
\end{itemize}""",
                ang = """Possible choices:
\begin{itemize}
\item 0: Roe scheme,
\item 1: kinetic,
\item 3: Zokagoa scheme,
\item 4: Tchamen scheme,
\item 5: HLLC,
\item 6: WAF.
\end{itemize}""",
            ),
        ),
#       -----------------------------------
        TREATMENT_OF_THE_LINEAR_SYSTEM = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["Coupled","Wave equation"],
            defaut = "Wave equation",
            fr = """Choix possibles :
\begin{itemize}
\item 1 : Traitement couple ;
\item 2 : Equation d onde.
\end{itemize}
Ancienne valeur par defaut = 1 (couple) jusqu a la version
V8P1.""",
            ang = """Possible choices:
\begin{itemize}
\item 1: Coupled,
\item 2: Wave equation.
\end{itemize}
Old default value = 1 (coupled) until release V8P1.""",
        ),
    ),
#   -----------------------------------
    FLUID = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        CORIOLIS_EFFECT = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            CORIOLIS = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Prise en compte ou non de la force de Coriolis.""",
                ang = """The Coriolis force is taken into account or ignored.""",
            ),
#           -----------------------------------
            CORIOLIS_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Fixe la valeur du coefficient de la force de Coriolis.
Celui-ci doit etre calcule en fonction de la latitude $l$
par la formule
$FCOR = 2 \omega sin(l)$ ,
$\omega$ etant la vitesse de rotation de la terre.
$\omega$ = 7.2921 10-5 rad/s.\\
Les composantes de la force de Coriolis sont alors :\\
$FU =  FCOR \times V,$\\
$FV = -FCOR \times U.$
Lorsqu''on utilise les coordonnees spheriques, le coefficient de
Coriolis est calcule automatiquement.""",
                ang = """Sets the value of the Coriolis force coefficient,
in cartesian coordinates.
This coefficient, denoted \telfile{FCOR} in the code, should be equal to
$2 \omega \sin(l)$  where $\omega$ denotes the earth angular speed of
rotation and $l$ the latitude. $\omega$ = 7.2921 10-5 rad/s.\\
The Coriolis force components are then:\\
$FU =  FCOR \times V,$\\
$FV = -FCOR \times U.$\\
When using the spherical coordinates, the Coriolis coefficient is
automatically computed.""",
            ),
        ),
#       -----------------------------------
        TSUNAMI = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            OPTION_FOR_TSUNAMI_GENERATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["No Tsunami","Tsunami generated on the basis of the Okada model 1992"],
                defaut = "No Tsunami",
                fr = """Choix possibles :
\begin{itemize}
\item 0 : pas de tsunami ;
\item 1 : generation d un Tsunami sur la base du modele de Okada
(1992).
\end{itemize}""",
                ang = """Possible choices:
\begin{itemize}
\item 0: no tsunami,
\item 1: tsunami generated on the basis of the Okada model (1992).
\end{itemize}""",
            ),
#           -----------------------------------
            PHYSICAL_CHARACTERISTICS_OF_THE_TSUNAMI = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min=10, max=10,
                defaut = [100.,210000.,75000.,13.6,81.,41.,110.,0.,0.,3.],
                fr = """Caracteristiques physiques du tsunami au nombre de 10, dans l ordre :
\begin{itemize}
\item $HH$ profondeur focale (en m) ;
\item $L$ longueur de fault (en m) ;
\item $W$ largeur de fault (en m) ;
\item $D$ dislocation (en m) ;
\item $TH$ direction de strike (en degres decimaux) ;
\item $DL$ angle d immersion (en degres decimaux) ;
\item $RD$ angle de glissement (en degres decimaux) ;
\item $Y0$ latitude de l epicentre (en degres decimaux) ;
\item $X0$ longitude de l epicentre (en degres decimaux) ;
\item $C0$ taille de l ellipse d influence ($L$ $\times$ $W$).
\end{itemize}""",
                ang = """Physical characteristics of the tsunami.
There are 10 of them:
\begin{itemize}
\item $HH$ focal depth (in m),
\item $L$ fault length (in m),
\item $W$ fault width (in m),
\item $D$ dislocation (in m),
\item $TH$ strike direction (in decimal degrees),
\item $DL$ dip angle (in decimal degrees),
\item $RD$ slip angle (in decimal degrees),
\item $Y0$ epicentre latitude (in decimal degrees),
\item $X0$ epicentre longitude (in decimal degrees),
\item $C0$ size of the ellipse of influence ($L$ $\times$ $W$).
\end{itemize}""",
            ),
        ),
#       -----------------------------------
        SECONDARY_CURRENTS_INFO = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            SECONDARY_CURRENTS = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Pour prendre en compte les courants secondaires.""",
                ang = """Using the parametrisation for secondary currents.""",
            ),
#           -----------------------------------
            b_SECONDARY_CURRENTSG = BLOC(condition="SECONDARY_CURRENTS == True",
#           -----------------------------------
#               -----------------------------------
                PRODUCTION_COEFFICIENT_FOR_SECONDARY_CURRENTS = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 7.071,
                    fr = """Une constante dans les termes de creation de $\Omega$.""",
                    ang = """A constant in the production terms of $\Omega$.""",
                ),
#               -----------------------------------
                DISSIPATION_COEFFICIENT_FOR_SECONDARY_CURRENTS = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 5.E-1,
                    fr = """Coefficient de dissipation de $\Omega$.""",
                    ang = """Coefficient of dissipation term of $\Omega$.""",
                ),
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
NUMERICAL_PARAMETERS = PROC(nom= "NUMERICAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    ADVANCED = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        FINITE_VOLUME_SCHEME_SPACE_ORDER = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Choix possibles :
\begin{itemize}
\item 1 : ordre 1 en espace ;
\item 2 : ordre 2 en espace.
\end{itemize}""",
            ang = """Possible choices:
\begin{itemize}
\item 1: first order in space,
\item 2: second order in space.
\end{itemize}""",
        ),
#       -----------------------------------
        FINITE_VOLUME_SCHEME_TIME_ORDER = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Choix possibles :
\begin{itemize}
\item 1 : ordre 1 en temps ;
\item 2 : ordre 2 en temps.
\end{itemize}""",
            ang = """Possible choices:
\begin{itemize}
\item 1: first order in time,
\item 2: second order in time.
\end{itemize}""",
        ),
#       -----------------------------------
        FLUX_LIMITOR_FOR_H_PLUS_Z = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Choix possibles :
\begin{itemize}
\item 1 : Minmod ;
\item 2 : Van Albada ;
\item 3 : MC (Monotonized Central-difference) ;
\item 4 : GenMinmod.
\end{itemize}""",
            ang = """Possible choices:
\begin{itemize}
\item 1 : Minmod,
\item 2 : Van Albada,
\item 3 : MC (Monotonized Central-difference),
\item 4 : GenMinmod.
\end{itemize}""",
        ),
#       -----------------------------------
        FLUX_LIMITOR_FOR_U_AND_V = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 2,
            fr = """Choix possibles :
\begin{itemize}
\item 1 : Minmod ;
\item 2 : Van Albada ;
\item 3 : MC (Monotonized Central-difference) ;
\item 4 : GenMinmod.
\end{itemize}""",
            ang = """Possible choices:
\begin{itemize}
\item 1 : Minmod,
\item 2 : Van Albada,
\item 3 : MC (Monotonized Central-difference),
\item 4 : GenMinmod.
\end{itemize}""",
        ),
#       -----------------------------------
        FLUX_LIMITOR_FOR_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 2,
            fr = """Choix possibles :
\begin{itemize}
\item 1 : Minmod ;
\item 2 : Van Albada ;
\item 3 : MC (Monotonized Central-difference) ;
\item 4 : GenMinmod.
\end{itemize}""",
            ang = """Possible choices:
\begin{itemize}
\item 1 : Minmod,
\item 2 : Van Albada,
\item 3 : MC (Monotonized Central-difference),
\item 4 : GenMinmod.
\end{itemize}""",
        ),
#       -----------------------------------
        TYPE_OF_BOUNDARY_CONDITION_FOR_KINETIC_SCHEME = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Choix possibles :
\begin{itemize}
\item 1 : Imposition faible ;
\item 2 : Imposition forte.
\end{itemize}""",
            ang = """Possible choices:
\begin{itemize}
\item 1: Weak imposition,
\item 2: Strong imposition.
\end{itemize}""",
        ),
#       -----------------------------------
        MATRIX_STORAGE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["classical EBE","Edge-based storage"],
            defaut = "Edge-based storage",
            fr = """Permet de definir la methode de stockage des matrices.
Les choix possibles sont :
\begin{itemize}
\item 1 : EBE classique ;
\item 3 : stockage par segments.
\end{itemize}
L option 3 est obligatoire avec un schema distributif pour la convection
(= 3, 4, 5, 13, 14 ou 15).""",
            ang = """Defines the method to store matrices. The possible choices are:
\begin{itemize}
\item 1: classical EBE,
\item 3: edge-based storage.
\end{itemize}
Option 3 is mandatory with a distributive scheme for advection
(= 3, 4, 5, 13, 14 or 15).""",
        ),
#       -----------------------------------
        MATRIX_VECTOR_PRODUCT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Choix possibles :
\begin{itemize}
\item 1 : classique ;
\item 2 : frontal.
Attention, avec 2, il faut une numerotation speciale des points.
\end{itemize}""",
            ang = """Possible choices are:
\begin{itemize}
\item 1: classic,
\item 2: frontal.
Beware, with option 2, a special numbering of points is required.
\end{itemize}""",
        ),
#       -----------------------------------
        NEWMARK_TIME_INTEGRATION_COEFFICIENT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.5,
            fr = """Choix possibles :
\begin{itemize}
\item 1. : Euler explicite ;
\item 0.5 : ordre 2 en temps.
\end{itemize}
Seulement pour Volumes Finis.""",
            ang = """Possible choices are:
\begin{itemize}
\item 1.: Euler explicit,
\item 0.5: order 2 in time.
\end{itemize}
Only for Finite Volumes.""",
        ),
#       -----------------------------------
        ZERO = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-10,
            fr = """Non active pour l''instant.""",
            ang = """Not yet implemented""",
        ),
#       -----------------------------------
        OPTION_OF_THE_HYDROSTATIC_RECONSTRUCTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Donne l option de la reconstruction hydrostatique
(option utile uniquement pour les Volumes Finis,
pour les schemas cinetique, HLLC et WAF) :
\begin{itemize}
\item 1 : option d Audusse et al. ;
\item 2 : option de Chen et Noelle.
\end{itemize}""",
            ang = """Gives the option for hydrostatic reconstruction
(only used for Finite Volumes with kinetic, HLLC and WAF schemes):
\begin{itemize}
\item 1: Audusse et al. option ;
\item 2: Chen and Noelle option.
\end{itemize}""",
        ),
#       -----------------------------------
        CONVERGENCE_STUDY = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Active une etude de convergence par rapport a une
solution analytique sur un maillage fin.""",
            ang = """Activates a convergence study compared
to an analytical solution on a fine mesh.""",
        ),
#       -----------------------------------
        REFINEMENT_LEVELS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Donne le nombre de raffinements que l''utilisateur
veut utiliser pour l''etude de convergence
(en activant \telkey{ETUDE DE CONVERGENCE} = OUI).
Chaque niveau multiplie par 4 le nombre d''elements.""",
            ang = """Gives the number of refinement levels that the
user wants to use in the convergence study (when activating
\telkey{CONVERGENCE STUDY} = YES).
Each level multiplies the number of elements by 4.""",
        ),
    ),
#   -----------------------------------
    SOLVER_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SOLVER = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","cgstab","gmres","direct"],
            defaut = "conjugate gradient on a normal equation",
            fr = """Permet de choisir le solveur utilise pour la resolution de
l''etape de propagation.
Les choix possibles sont :
\begin{itemize}
\item 1 : gradient conjugue,
\item 2 : residu conjugue,
\item 3 : gradient conjugue sur equation normale,
\item 4 : erreur minimale,
\item 5 : gradient conjugue carre (non programme),
\item 6 : gradient conjugue carre stabilise (cgstab),
\item 7 : GMRES (voir aussi \telkey{OPTION DU SOLVEUR}),
\item 8 : direct.
\end{itemize}""",
            ang = """Makes it possible to select the solver used for solving the
propagation step.
Possible choices are:
\begin{itemize}
\item 1: conjugate gradient,
\item 2: conjugate residual,
\item 3: conjugate gradient on a normal equation,
\item 4: minimum error,
\item 5: squared conjugate gradient (not implemented),
\item 6: conjugate gradient squared stabilised (cgstab),
\item 7: GMRES (see \telkey{SOLVER OPTION}),
\item 8: direct.
\end{itemize}""",
        ),
#       -----------------------------------
        b_SOLVERG = BLOC(condition="SOLVER == 'gmres'",
#       -----------------------------------
#           -----------------------------------
            SOLVER_OPTION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 2,
                fr = """Si le solveur est GMRES (7), le mot cle est la dimension de
l''espace de Krylov (valeurs conseillees entre 2 et 15).""",
                ang = """When GMRES (7) is chosen for solver, dimension of the Krylov space.
Try values between 2 and 15.""",
            ),
        ),
#       -----------------------------------
        SOLVER_ACCURACY = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-4,
            fr = """Precision demandee pour la resolution de l''etape de propagation
(cf. Note de principe).""",
            ang = """Required accuracy for solving the propagation step
(refer to Principle note).""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_ITERATIONS_FOR_SOLVER = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 100,
            fr = """Les algorithmes utilises pour la resolution de l''etape de
propagation etant iteratifs. Il est necessaire de limiter le nombre
d''iterations autorisees.
Remarque : un maximum de 40 iterations par pas de temps semble
raisonnable.""",
            ang = """Since the algorithms used for solving the propagation step are
iterative, the allowed number of iterations should be limited.
NOTE: a maximum number of 40 iterations per time step seems to be
reasonable.""",
        ),
#       -----------------------------------
        CONTINUITY_CORRECTION = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Corrige les vitesses sur les points avec hauteur imposee ou
l''equation de continuite n''a pas ete resolue.
Doit etre activee avec des bancs decouvrants et
\telkey{TRAITEMENT DES HAUTEURS NEGATIVES} = 2 ou 3.""",
            ang = """Correction of the velocities on points with a prescribed
elevation, where the continuity equation has not been solved.
It has to be activated with tidal flats and
\telkey{TREATMENT OF NEGATIVE DEPTHS} = 2 or 3.""",
        ),
#       -----------------------------------
        PRECONDITIONING = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["diagonal","no preconditioning","block-diagonal (4-9 matrices)","absolute value of diagonal","Crout","Gauss-Seidel","Supplied by the user","diagonal and Crout"],
            defaut = "diagonal",
            fr = """Permet de preconditionner le systeme de l''etape de propagation
afin d''accelerer la convergence lors de sa resolution.
\begin{itemize}
\item 0 : pas de preconditionnement ;
\item 2 : preconditionnement diagonal ;
\item 3 : preconditionnement diagonal-bloc ;
\item 5 : preconditionnement valeur absolue diagonale ;
\item 7 : preconditionnement de Crout par element ou segment
(ne marche pas en parallele) ;
\item 11 : preconditionnement de Gauss-Seidel par element ou segment ;
\item 13 : preconditionnement fourni par l''utilisateur.
\end{itemize}
Certains preconditionnements sont cumulables
(les diagonaux 2 ou 3 avec les autres)
Pour cette raison on ne retient que les nombres premiers pour
designer les preconditionnements. Si l''on souhaite en cumuler
plusieurs on formera le produit des options correspondantes.""",
            ang = """Choice of the preconditioning in the propagation step linear
system that the convergence is speeded up when it is being solved.
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal preconditioning,
\item 3: block-diagonal preconditioning (systemes a 4 ou 9 matrices),
\item 5: diagonal preconditioning with absolute value,
\item 7: Crout''s preconditioning per element or segment
(does not work in parallel),
\item 11: Gauss-Seidel''s preconditioning per element or segment,
\item 13: preconditioning supplied by the user.
\end{itemize}
Some operations (either 2 or 3 diagonal preconditioning) can be
performed concurrently with the others.
Only prime numbers are therefore kept to denote the preconditioning
operations. When several of them are to be performed concurrently,
the product of relevant options shall be made.""",
        ),
#       -----------------------------------
        C_U_PRECONDITIONING = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Changement de variable de $H$ en $C$ dans le systeme lineaire final.
Cette option est desactivee avec l equation d onde i.e.
\telkey{TRAITEMENT DU SYSTEME LINEAIRE} = 2.""",
            ang = """Change of variable from $H$ to $C$ in the final linear system.
This option is deactivated with wave equation i.e.
\telkey{TREATMENT OF THE LINEAR SYSTEM} = 2.""",
        ),
#       -----------------------------------
        FINITE_ELEMENT_ASSEMBLY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["Normal","Integer I8","Compensated"],
            defaut = "Normal",
            fr = """Les choix possibles sont :
\begin{itemize}
\item 1 : normal ;
\item 2 : avec des entiers I8 ;
\item 3 : compense (pour reproductibilite).
\end{itemize}""",
            ang = """Possible choices are:
\begin{itemize}
\item 1: normal,
\item 2: with I8 integers,
\item 3: compensation (for reproducibility).
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    DISCRETISATIONS_IMPLICITATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        IMPLICITATION_FOR_DEPTH = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.55,
            fr = """Fixe la valeur du coefficient d''implicitation sur C dans l''etape de
propagation (cf. Note de principe). Les valeurs inferieures a 0.5
donnent un schema instable (et sont donc interdites).""",
            ang = """Sets the value of the implicitation coefficient for C (the celerity of
waves) in the propagation step (refer to principle note). Values below
0.5 result in an unstable scheme (and are then forbidden).""",
        ),
#       -----------------------------------
        IMPLICITATION_FOR_VELOCITY = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.55,
            fr = """Fixe la valeur du coefficient d''implicitation sur la vitesse dans
l''etape de propagation (cf. Note de principe). Les valeurs inferieures
a 0.5 donnent un schema instable (et sont donc interdites).""",
            ang = """Sets the value of the implicitation coefficient for velocity in the
propagation step (refer to principle note). Values below 0.5 result in
an unstable condition (and are then forbidden).""",
        ),
#       -----------------------------------
        DISCRETIZATIONS_IN_SPACE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', max='**',
            into = ["linear","quasi-bubble","quadratic"],
            defaut = ["linear","linear","linear","linear","linear"],
            fr = """Choix de la discretisation pour chaque variable.
Ces coefficients sont respectivement appliques a :
\begin{itemize}
\item 1) $U$ et $V$ ;
\item 2) $H$ ;
\item 3) $T$ ;
\item 4) $k$ et $\epsilon$ (modele $k-\epsilon$) ;
\item 5) $\tilde{\nu}$ (modele Spalart-Allmaras).
\end{itemize}
Les choix possibles sont :
\begin{itemize}
\item 11 : lineaire ;
\item 12 : quasi-bulle ;
\item 13 : quadratique.
\end{itemize}
Les elements quadratiques (13) ne sont pas implementes pour l equation
d onde (\telkey{TRAITEMENT DU SYSTEME LINEAIRE} = 2).""",
            ang = """Choice of space discretisation for every variable.
These coefficients are respectively applied to:
\begin{itemize}
\item 1) $U$ and $V$,
\item 2) $H$,
\item 3) $T$,
\item 4) $k$ and $\epsilon$ ($k-\epsilon$ model),
\item 5) $\tilde{\nu}$ (Spalart-Allmaras model).
\end{itemize}
Possible choices are:
\begin{itemize}
\item 11: linear,
\item 12: quasi-bubble,
\item 13: quadratic.
\end{itemize}
Quadratic elements (13) have not been implemented for wave equation
(\telkey{TREATMENT OF THE LINEAR SYSTEM} = 2).""",
        ),
#       -----------------------------------
        b_DISCRETIZATIONS_IN_SPACEG = BLOC(condition="DISCRETIZATIONS_IN_SPACE != None",
#       -----------------------------------
#           -----------------------------------
            Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#           -----------------------------------
                defaut = "Choice of space discretisation for every variable. These coefficients are applied respectively to 1) U and V 2) H 3) T 4) K and EPSILON"),
        ),
    ),
#   -----------------------------------
    PROPAGATION_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        PROPAGATION = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non de la propagation de la vitesse et de la
hauteur d''eau.
La diffusion etant contenue dans cette etape sera supprimee aussi.""",
            ang = """Determines whether the propagation step is taken into account
or not.
The diffusion being included in that step will be deleted as well.""",
        ),
#       -----------------------------------
        b_PROPAGATIONG = BLOC(condition="PROPAGATION == True",
#       -----------------------------------
#           -----------------------------------
            MEAN_DEPTH_FOR_LINEARIZATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Fixe la hauteur d''eau autour de laquelle s''effectue la linearisation
lorsque l''option \telkey{PROPAGATION LINEARISEE} est choisie.""",
                ang = """Sets the water depth around which the linearization is done
when the \telkey{LINEARIZED PROPAGATION} option is selected.""",
            ),
#           -----------------------------------
            INITIAL_GUESS_FOR_U = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["zero","previous","extrapolation"],
                defaut = "previous",
                fr = """Tir initial du solveur de l''etape de propagation.
Offre la possibilite de modifier la valeur initiale de U, a chaque
iteration, dans l''etape de propagation en utilisant les valeurs
finales de cette variable aux pas de temps precedents. Ceci peut
permettre d''accelerer la vitesse de convergence lors de la resolution
du systeme. Trois possibilites sont offertes :
\begin{itemize}
\item 0 : U = 0 ;
\item 1 : U = U(n) ;
\item 2 : U = 2 U(n)- U(n-1) (extrapolation).
\end{itemize}""",
                ang = """Initial guess for the solver in the propagation step.
Makes it possible to modify the initial value of U, upon each
iteration in the propagation step, by using the ultimate values this
variable had in the earlier time steps. Thus, the convergence can be
speeded up when the system is being solved. 3 options are available:
\begin{itemize}
\item 0 : U = 0,
\item 1 : U = U(n),
\item 2 : U = 2 U(n)- U(n-1) (extrapolation).
\end{itemize}""",
            ),
        ),
#       -----------------------------------
        INITIAL_GUESS_FOR_H = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["zero","previous","extrapolation"],
            defaut = "previous",
            fr = """Tir initial du solveur de l''etape de propagation. Offre la
possibilite de modifier la valeur initiale de DH, accroissement de H,
a chaque iteration, dans l''etape de propagation en utilisant les
valeurs finales de cette variable aux pas de temps precedents. Ceci peut
permettre d''accelerer la vitesse de convergence lors de la resolution
du systeme. Trois possibilites sont offertes :
\begin{itemize}
\item 0 : DH = 0 ;
\item 1 : DH = DHn (valeur finale de DH au pas de temps precedent) ;
\item 2 : DH = 2.DHn - DHn-1 (extrapolation).
\end{itemize}""",
            ang = """Initial guess for the solver in the propagation step.
Makes it possible to modify the initial value of C, upon each
iteration in the propagation step, by using the ultimate values this
variable had in the earlier time steps. Thus, the convergence can be
speeded up when the system is being solved. 3 options are available:
\begin{itemize}
\item 0: DH = 0,
\item 1: DH = DHn (ultimate DH value in the next previous time step),
\item 2: DH = 2.DHn - DHn-1 (extrapolation).
\end{itemize}""",
        ),
#       -----------------------------------
        LINEARIZED_PROPAGATION = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Permet de lineariser l''etape de propagation;
par exemple lors de la realisation de cas tests
pour lesquels on dispose
d''une solution analytique dans le cas linearise.""",
            ang = """Provided for linearizing the propagation step, e.g. when
performing test-cases for which an analytical
solution in the linearized
case is available.""",
        ),
    ),
#   -----------------------------------
    ADVECTION_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ADVECTION = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non des termes de convection.
En cas de reponse positive,
on peut encore supprimer certains termes de convection avec
les mots-cles \telkey{CONVECTION DE}...""",
            ang = """Are the advection terms taken into account or not?
If YES, some advection terms can still be deleted using the keywords
\telkey{ADVECTION OF}...""",
        ),
#       -----------------------------------
        ADVECTION_OF_H = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non de la convection de $H$.""",
            ang = """The advection of $H$ is taken into account or ignored.""",
        ),
#       -----------------------------------
        ADVECTION_OF_U_AND_V = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non de la convection de $U$ et $V$.""",
            ang = """The advection of $U$ and $V$ is taken into account or ignored.""",
        ),
#       -----------------------------------
        b_ADVECTION_OF_U_AND_VG = BLOC(condition="ADVECTION_OF_U_AND_V == True",
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_VELOCITIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO ADVECTION","CHARACTERISTICS","EXPLICIT + SUPG","EXPLICIT LEO POSTMA","EXPLICIT + MURD SCHEME N","EXPLICIT + MURD SCHEME PSI","N-SCHEME FOR TIDAL FLATS","N-SCHEME FOR TIDAL FLATS","ERIA SCHEME"],
                defaut = "CHARACTERISTICS",
                fr = """Choix du schema de convection pour les vitesses,
remplace \telkey{FORME DE LA CONVECTION}.""",
                ang = """Choice of the advection scheme for the velocities,
replaces \telkey{TYPE OF ADVECTION}.""",
            ),
        ),
#       -----------------------------------
        TYPE_OF_ADVECTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', max='**',
            into = ["CHARACTERISTICS","SUPG","CONSERVATIVE N-SCHEME LP","CONSERVATIVE N-SCHEME","CONSERVATIVE PSI-SCHEME","EDGE-BASED N-SCHEME LP","EDGE-BASED N-SCHEME","ERIA SCHEME"],
            defaut = ["CHARACTERISTICS","CONSERVATIVE PSI-SCHEME","CHARACTERISTICS","CHARACTERISTICS"],
            fr = """Choix du schema de convection pour chaque variable.
Ces coefficients sont respectivement appliques a\\
\begin{itemize}
\item 1) $U$ et $V$ ;
\item 2) $H$ ;
\item 3) $T$ ;
\item 4) $k$ et $\epsilon$.
\end{itemize}
Les choix possibles sont :
\begin{itemize}
\item 1 : caracteristiques sur $h$ ;
\item 2 : SUPG ;
\item 3 : Schema N conservatif ;
\item 4 : Schema N conservatif ;
\item 5 : Schema PSI conservatif ;
\item 13 : Schema N par segment ;
\item 14 : Schema N par segment ;
\item 15 : Schema ERIA.
\end{itemize}
Le 2e entier doit etre 5.""",
            ang = """Choice of advection schemes for every variable.
These coefficients are applied respectively to
\begin{itemize}
\item 1) $U$ and $V$,
\item 2) $H$,
\item 3) $T$,
\item 4) $k$ and $\epsilon$.
\end{itemize}
Possible choices are:
\begin{itemize}
\item 1: characteristics,
\item 2: SUPG,
\item 3: Conservative N-scheme,
\item 4: Conservative N-scheme,
\item 5: Conservative PSI-scheme,
\item 13: Edge-based N-scheme,
\item 14: Edge-based N-scheme,
\item 15: ERIA scheme.
\end{itemize}
The second integer must be 5.""",
        ),
#       -----------------------------------
        b_TYPE_OF_ADVECTIONG = BLOC(condition="TYPE_OF_ADVECTION != None",
#       -----------------------------------
#           -----------------------------------
            Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#           -----------------------------------
                defaut = "Choice of space discretisation for every variable. These coefficients are applied respectively to 1) U and V 2) H 3) T 4) K and EPSILON"),
        ),
#       -----------------------------------
        OPTION_FOR_CHARACTERISTICS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["strong","weak"],
            defaut = "strong",
            fr = """Les choix possibles sont :
\begin{itemize}
\item 1 : forme forte ;
\item 2 : forme faible.
\end{itemize}
Si une des composantes du tableau  \telkey{FORME DE LA CONVECTION} = 1
ou \telkey{SCHEMA POUR LA CONVECTION...} = 1,
ainsi que le mot-cle correspondant
\telkey{OPTION DU SCHEMA POUR LA CONVECTION...} = 2,
\telkey{OPTION POUR LES CARACTERISTIQUES} est automatiquement mis a 2.""",
            ang = """Possible choices are:
\begin{itemize}
\item 1: strong form,
\item 2: weak form.
\end{itemize}
If one component of array \telkey{TYPE OF ADVECTION} = 1 or
\telkey{SCHEME FOR ADVECTION OF...} = 1,
and also the corresponding keyword
\telkey{SCHEME OPTION FOR ADVECTION OF...} = 2,
\telkey{OPTION FOR CHARACTERISTICS} is automatically set to 2.""",
        ),
#       -----------------------------------
        SUPG_OPTION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I', min=0, max='**',
            defaut = [2,2,2,2],
            fr = """Les choix possibles sont :
\begin{itemize}
\item 0:pas de decentrement SUPG ;
\item 1:SUPG classique ;
\item 2:SUPG modifiee.
\end{itemize}
Ces coefficients sont respectivement appliques a :
\begin{itemize}
\item 1) $U$ et $V$ ;
\item 2) $H$ ;
\item 3) $T$ ;
\item 4) $k$ et $\epsilon$.
\end{itemize}
En cas d utilisation d un schema distributif (3, 4, 5, 13, 14 ou 15),
le coefficient applique a $H$ est automatiquement mis a 0.
De plus, si \telkey{TRAITEMENT DES HAUTEURS NEGATIVES} = 2 ou 3
avec \telkey{OPTION DE TRAITEMENT DES BANCS DECOUVRANTS} = 1,
il est obligatoire de choisir 0 pour la 2e composante du mot-cle
\telkey{OPTION DE SUPG} (hauteur d eau).""",
            ang = """Possible choices are:
\begin{itemize}
\item 0: no upwinding,
\item 1: classical SUPG,
\item 2: modified SUPG.
\end{itemize}
These coefficients are applied respectively to:
\begin{itemize}
\item 1) $U$ and $V$,
\item 2) $H$,
\item 3) $T$,
\item 4) $k$ and $\epsilon$.
\end{itemize}
If using a distributive scheme (3, 4, 5, 13, 14, 15),
the coefficient applied to $H$ is automatically set to 0.
Moreover, if using \telkey{TREATMENT OF NEGATIVE DEPTHS} = 2 or 3
with \telkey{OPTION FOR THE TREATMENT OF TIDAL FLATS} = 1,
it is mandatory to choose 0 for the 2nd component of
\telkey{SUPG OPTION} (water depth).""",
        ),
#       -----------------------------------
        NUMBER_OF_GAUSS_POINTS_FOR_WEAK_CHARACTERISTICS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 3,
            fr = """Voir les release notes 6.3.
Nombre de points de Gauss utilises pour le calcul des caracteristiques
faibles.
Les choix possibles sont :
\begin{itemize}
\item 1 point ;
\item 3 points ;
\item 4 points ;
\item 6 points ;
\item 7 points ;
\item 12 points.
\end{itemize}
Plus le nombre est grand, plus le schema sera conservatif,
mais plus les couts de calculs seront eleves.""",
            ang = """See release notes 6.3.
Number of Gauss points used to compute the weak characteristics.
Possible choices are:
\begin{itemize}
\item 1 point,
\item 3 points,
\item 4 points,
\item 6 points,
\item 7 points,
\item 12 points.
\end{itemize}
The bigger the number is, the more conservative the scheme is,
but the higher the computational costs are.""",
        ),
#       -----------------------------------
        MASS_LUMPING_FOR_WEAK_CHARACTERISTICS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Fixe le taux de mass-lumping qui est applique a la matrice de
masse lors de l''utilisation des caracteristiques faibles.""",
            ang = """Sets the amount of mass-lumping that is applied to the mass
matrix when using weak characteristics.""",
        ),
#       -----------------------------------
        b_MAXIMUM_NUMBER_OF_ITERATIONS_FOR_ADVECTION_SCHEMESF = BLOC(condition="(ADVECTION_OF_TRACERS == True and SCHEME_FOR_ADVECTION_OF_TRACERS == 'EDGE-BASED N-SCHEME') or (ADVECTION_OF_K_AND_EPSILON == True and SCHEME_FOR_ADVECTION_OF_K_EPSILON == 'EDGE-BASED N-SCHEME') or (ADVECTION_OF_U_AND_V == True and SCHEME_FOR_ADVECTION_OF_VELOCITIES == 'EDGE-BASED N-SCHEME')",
#       -----------------------------------
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_ITERATIONS_FOR_ADVECTION_SCHEMES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 50,
            fr = """Limite le nombre d''iterations pour les schemas de convection,
seulement pour schemes 13, 14 et 15.
Ancienne valeur par defaut = 10 jusqu en version 8.1.""",
            ang = """Limits the number of solver iterations for the advection
schemes, only for schemes 13, 14 and 15.
Old default value = 10 until release 8.1.""",
        ),
#       -----------------------------------
        MASS_LUMPING_ON_H = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """\telemac{2d} offre la possibilite d''effectuer du mass-lumping
sur $H$ ou sur la vitesse.
Ceci revient a ramener tout ou partie (suivant la valeur de ce
coefficient) des matrices AM1 (h) ou AM2 (U) et AM3 (V) sur
leur diagonale.
Cette technique permet d''accelerer le code dans des proportions tres
importantes et de le rendre egalement beaucoup plus stable. Cependant
les solutions obtenues se trouvent lissees.
Ce parametre fixe le taux de mass-lumping effectue sur $H$.""",
            ang = """\telemac{2d} provides an opportunity to carry out mass-lumping
either on $H$ or on the velocity.
This is equivalent to bringing the matrices AM1 (h) or AM2 (U) and
AM3 (V) wholly or partly, back onto their diagonal.
Thanks to that technique, the code can be speeded up to a quite
significant extent and it can also be made much more stable. The
resulting solutions, however, become artificially smoothed.
This parameter sets the extent of mass-lumping that is performed
on $h$.""",
        ),
#       -----------------------------------
        MASS_LUMPING_ON_VELOCITY = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Fixe le taux de mass-lumping effectue sur la vitesse.
Le mot-cle \telkey{TRAITEMENT DU SYSTEME LINEAIRE} change la valeur
utilisee a 1.""",
            ang = """Sets the amount of mass-lumping that is performed on the velocity.
The keyword \telkey{TREATMENT OF THE LINEAR SYSTEM} changes the
used value to 1.""",
        ),
#       -----------------------------------
        SCHEME_OPTION_FOR_ADVECTION_OF_VELOCITIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Si present remplace et a priorite sur :
\telkey{OPTION POUR LES CARACTERISTIQUES}
\telkey{OPTION DE SUPG}.
Si caracteristiques :
\begin{itemize}
\item 1 = forme forte,
\item 2 = forme faible.
\end{itemize}
Si schema PSI ou N :
\begin{itemize}
\item 1 = explicite ;
\item 2 = predicteur-correcteur ;
\item 3 = predicteur-correcteur 2e ordre en temps ;
\item 4 = implicite.
\end{itemize}""",
            ang = """If present replaces and has priority over:
\telkey{OPTION FOR CHARACTERISTICS}
\telkey{SUPG OPTION}.
If characteristics:
\begin{itemize}
\item 1 = strong form,
\item 2 = weak form.
\end{itemize}
If N or PSI scheme:
\begin{itemize}
\item 1 = explicit,
\item 2 = predictor-corrector,
\item 3 = predictor-corrector second-order in time,
\item 4 = implicit.
\end{itemize}""",
        ),
#       -----------------------------------
        FREE_SURFACE_GRADIENT_COMPATIBILITY = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """Des valeurs comprises entre 0 et 1 peuvent supprimer les
oscillations parasites.""",
            ang = """Values between 0 and 1 may suppress spurious oscillations.""",
        ),
#       -----------------------------------
        NUMBER_OF_SUB_ITERATIONS_FOR_NON_LINEARITIES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Permet de reactualiser, pour un meme pas de temps, les champs
convecteur et propagateur au cours de plusieurs sous-iterations. A la
premiere sous-iteration, ces champs sont donnes par C et le champ de
vitesses au pas de temps precedent. Aux iterations suivantes, ils sont
pris egaux au champ de vitesse obtenu a la fin de la sous-iteration
precedente. Cette technique permet d''ameliorer la prise en compte des
non linearites.""",
            ang = """Used for updating, within one time step, the advection and
propagation field.
upon the first sub-iteration, these fields are given by
C and the velocity field in the previous time step. At subsequent
iterations, the results of the previous sub-iteration is used to
update the advection and propagation field.
The non-linearities can be taken into account through this technique.""",
        ),
#       -----------------------------------
        b_TREATMENT_OF_FLUXES_AT_THE_BOUNDARIESF = BLOC(condition="(ADVECTION_OF_TRACERS == True and SCHEME_FOR_ADVECTION_OF_TRACERS in ['EDGE-BASED N-SCHEME','SUPG','CONSERVATIVE N-SCHEME','CONSERVATIVE PSI-SCHEME']) or (ADVECTION_OF_K_AND_EPSILON == True and SCHEME_FOR_ADVECTION_OF_K_EPSILON in ['EDGE-BASED N-SCHEME','SUPG','CONSERVATIVE N-SCHEME','CONSERVATIVE PSI-SCHEME']) or (ADVECTION_OF_U_AND_V == True and SCHEME_FOR_ADVECTION_OF_VELOCITIES in ['EDGE-BASED N-SCHEME','SUPG','CONSERVATIVE N-SCHEME','CONSERVATIVE PSI-SCHEME'])",
#       -----------------------------------
        ),
#       -----------------------------------
        TREATMENT_OF_FLUXES_AT_THE_BOUNDARIES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["Priority to prescribed values","Priority to fluxes"],
            defaut = ["Priority to prescribed values"],
            fr = """Utilise pour les schemas SUPG, PSI et N, avec option 2, on ne retrouve
pas exactement les valeurs imposees des traceurs, mais le flux est
correct. Une seule meme valeur pour toutes les frontieres liquides.""",
            ang = """Used so far only with the SUPG, PSI and N schemes. With option 2,
Dirichlet prescribed values are not obeyed, but the fluxes are correct.
One single and same value for every liquid boundary.""",
        ),
#       -----------------------------------
        NUMBER_OF_CORRECTIONS_OF_DISTRIBUTIVE_SCHEMES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """Pour les options avec predicteur-correcteur et schema de convection
de type (3, 4, 5, LIPS ou pas et ERIA).
Nombre d iterations pour tous les pas de temps (ou sous pas) pour
converger vers la solution.
Il est utile pour les cas non stationnaires.
Pour les ecoulements quasi-stationnaires, ce mot-cle n a pas d impact
sur la solution, il peut donc etre fixe a 0.
D un autre cote, pour les ecoulements instationnaires, il est suggere
de fixer ce mot-cle a 2 (au moins), ce qui est un bon compromis
entre precision et temps CPU.
En effet, en augmentant le nombre de corrections, le schema est plus
precis mais le temps CPU augmente rapidement.""",
            ang = """For predictor-corrector options with advection scheme of type
3, 4, 5, LIPS or not, and ERIA).
Number of iterations for every time step (or sub-time step) to converge
to the solution.
It is useful for unsteady cases.
For quasi-steady flows, this keyword does not have a large impact
on the solution, so it can be set to 0.
On the other hand, for unsteady flows, it is suggested to set this
keyword to 2 (at least), which is a good compromise between accuracy
and computational time. Indeed, increasing the number of corrections
the scheme is more accurate but the CPU time rapidly increases.""",
        ),
#       -----------------------------------
        NUMBER_OF_SUB_STEPS_OF_DISTRIBUTIVE_SCHEMES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Pour les options predicteur-correcteur avec schema localement
implicite (3, 4 or 5).
Ce mot-cle permet de diviser le pas de temps donne par l utilisateur
dans le \telkey{FICHIER DES PARAMETRES} en plusieurs sous-pas.
A nouveau, il produit un effet sur la precision du schema et
il est pratique d ajuster ce mot-cle afin d avoir des nombres de
Courant pas trop grands (autour de 1).""",
            ang = """Only for implicit scheme with predictor-corrector (3, 4 or 5).
This keyword allows to subdivide the time step given by the user in the
\telkey{STEERING FILE}, into several sub-steps.
Again, it produces an effect on the precision of the scheme and
it is convenient to set this keyword in order to have Courant numbers
not too large (around 1).""",
        ),
    ),
#   -----------------------------------
    DIFFUSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DIFFUSION_OF_VELOCITY = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Permet de decider si l''on prend ou non en compte la diffusion
des vitesses.""",
            ang = """Makes it possible to decide whether the diffusion of velocity
(i.e. viscosity) is taken into account or not.""",
        ),
#       -----------------------------------
        b_DIFFUSION_OF_VELOCITYG = BLOC(condition="DIFFUSION_OF_VELOCITY == True",
#       -----------------------------------
#           -----------------------------------
            IMPLICITATION_FOR_DIFFUSION_OF_VELOCITY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Fixe la valeur du coefficient d''implicitation sur les termes de
diffusion des vitesses.""",
                ang = """Sets the value of the implicitation coefficient for the diffusion of
velocity.""",
            ),
#           -----------------------------------
            OPTION_FOR_THE_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["div( nu grad(U) )","1/h div ( h nu grad(U)"],
                defaut = "div( nu grad(U) )",
                fr = """Les choix possibles sont :
\begin{itemize}
\item 1: Diffusion de la forme div( $\nu$ grad($U$) ) ;
\item 2: Diffusion de la forme 1/$h$ div ( $h$ $\nu$ grad($U$) ).
\end{itemize}""",
                ang = """Possible choices are:
\begin{itemize}
\item 1: Diffusion in the form div( $\nu$ grad($U$) ),
\item 2: Diffusion in the form 1/$h$ div ( $h$ $\nu$ grad($U$) ).
\end{itemize}""",
            ),
        ),
    ),
#   -----------------------------------
    AUTOMATIC_DIFFERENTIATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        AD_NUMBER_OF_DERIVATIVES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Definit le nombre de derivees utilisateurs, dans le cadre
de la differentiation algorithmique.""",
            ang = """Defines the number of user derivatives, within the framework
of the algorithmic differentiation.""",
        ),
#       -----------------------------------
        AD_NAMES_OF_DERIVATIVES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min= 2, max= 2,
            fr = """Noms des differentiateurs utilisateurs en 32 caracteres,
         16 pour le nom, 16 pour l''unite.""",
            ang = """Name of user differentiators in 32 characters,
         16 for the name, 16 for the unit.""",
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
GENERAL_PARAMETERS = PROC(nom= "GENERAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    DEBUGGER = SIMP(statut ='o',
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
        TIME_STEP = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """Definit le pas de temps en secondes. Remarque : Pour une bonne
precision; il est souhaitable de choisir le pas de temps de telle sorte
que le nombre de Courant de propagation soit inferieur a 2 ; voir 3.
Ceci peut etre realisable en hydraulique fluviale ; mais ne l''est
pratiquement jamais en hydraulique maritime ou l''on peut atteindre des
valeurs de 50.""",
            ang = """Specifies the time step in seconds.""",
        ),
#       -----------------------------------
        NUMBER_OF_TIME_STEPS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """ Definit le nombre de pas de temps effectues lors de l''execution du
code.""",
            ang = """ Specifies the number of time steps performed when running the code.""",
        ),
#       -----------------------------------
        DURATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Duree de la simulation. Alternative au parametre
\telkey{NOMBRE DE PAS DE TEMPS}.
On en deduit le nombre de pas de temps en prenant l''entier le
plus proche de (duree du calcul/pas de temps).
Si le \telkey{NOMBRE DE PAS DE TEMPS} est aussi donne,
on prend la plus grande valeur.""",
            ang = """Sets the duration of simulation in seconds.
May be used instead of the parameter \telkey{NUMBER OF TIME STEPS}.
The nearest integer to (duration/time step) is taken.
If \telkey{NUMBER OF TIME STEPS} is also given,
the greater value is taken.""",
        ),
#       -----------------------------------
        ORIGINAL_DATE_OF_TIME = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I', min= 3, max= 3,
            defaut = [1900,1,1],
            fr = """Permet de fixer la date d''origine des temps du modele lorsque
la maree est prise en compte (force generatrice de la maree et/ou les
conditions aux limites de maritimes).
Egalement utilise pour les flotteurs,
les echanges thermiques avec l atmosphere, le frazil,
le chainage avec DELWAQ.""",
            ang = """Enables to set the date of the time origin of the model when
taking into account of the tide (tide generator force and/or the tidal
boundary conditions).
Also used with drogues, heat exchange with atmosphere, frazil,
chaining with DELWAQ.""",
        ),
#       -----------------------------------
        ORIGINAL_HOUR_OF_TIME = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I', min= 3, max= 3,
            defaut = [0,0,0],
            fr = """Permet de fixer l''heure d''origine des temps du modele lorsque
la maree est prise en compte (force generatrice de la maree et/ou les
conditions aux limites de maritimes.
Egalement utilise pour les flotteurs,
les echanges thermiques avec l atmosphere, le frazil,
le chainage avec DELWAQ.""",
            ang = """Enables to set the time of the time origin of the model when
taking into account of the tide (tide generator force and/or the tidal
boundary conditions).
Also used with drogues, heat exchange with atmosphere, frazil,
chaining with DELWAQ.""",
        ),
#       -----------------------------------
        STOP_IF_A_STEADY_STATE_IS_REACHED = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """A utiliser avec le mot-cle : \telkey{CRITERES D''ARRET}.""",
            ang = """To be used with the keyword: \telkey{STOP CRITERIA}.""",
        ),
#       -----------------------------------
        b_STOP_IF_A_STEADY_STATE_IS_REACHEDG = BLOC(condition="STOP_IF_A_STEADY_STATE_IS_REACHED == True",
#       -----------------------------------
#           -----------------------------------
            STOP_CRITERIA = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min= 3, max= 3,
                defaut = [1.E-4,1.E-4,1.E-4],
                fr = """Criteres d''arret pour un ecoulement permanent.
Ces coefficients sont respectivement appliques a:
\begin{itemize}
\item $U$ et $V$ ;
\item $H$ ;
\item Tracers.
\end{itemize}
A utiliser avec le mot-cle
\telkey{ARRET SI UN ETAT PERMANENT EST ATTEINT}.""",
                ang = """Stop criteria for a steady state.
These coefficients are applied respectively to:
\begin{itemize}
\item $U$ and $V$,
\item $H$,
\item Tracers.
\end{itemize}
To be used with the keyword
\telkey{STOP IF A STEADY STATE IS REACHED}.""",
            ),
        ),
#       -----------------------------------
        CONTROL_OF_LIMITS = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Utiliser avec le mot-cle : \telkey{VALEURS LIMITES}.
Le programme s''arrete si les limites sur
$U$, $V$, $H$ ou traceurs sont depassees.""",
            ang = """Use with the key-word: \telkey{LIMIT VALUES}.
The program is stopped if
the limits on $U$, $V$, $H$ or tracers are trepassed.""",
        ),
#       -----------------------------------
        b_CONTROL_OF_LIMITSG = BLOC(condition="CONTROL_OF_LIMITS == True",
#       -----------------------------------
#           -----------------------------------
            LIMIT_VALUES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min= 8, max= 8,
                defaut = [-1000.,9000.,-1000.,1000.,-1000.,1000.,-1000.,1000.],
                fr = """Utilise avec le mot-cle \telkey{CONTROLE DES LIMITES}.
Valeurs minimales et maximales acceptables pour
$H$, $U$ ,$V$ et traceurs dans l''ordre suivant :
min(H) max(H) min(U) max(U) min(V) max(V) min(T) max(T).""",
                ang = """To be used with the key-word \telkey{CONTROL OF LIMITS}.
Min and max acceptable values for
$H$, $U$, $V$ and tracers in the following order:
min(H) max(H) min(U) max(U) min(V) max(V) min(T) max(T).""",
            ),
        ),
#       -----------------------------------
        VARIABLE_TIME_STEP = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Pas de temps variable pour avoir un nombre de Courant souhaite.""",
            ang = """Variable time-step to get a given Courant number.""",
        ),
#       -----------------------------------
        b_VARIABLE_TIME_STEPG = BLOC(condition="VARIABLE_TIME_STEP == True",
#       -----------------------------------
#           -----------------------------------
            DESIRED_COURANT_NUMBER = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Nombre de Courant souhaite en cas de pas de temps variable.""",
                ang = """Desired Courant number when \telkey{VARIABLE TIME-STEP} is set to YES.""",
            ),
        ),
    ),
#   -----------------------------------
    LOCATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SPHERICAL_COORDINATES = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Choix des coordonnees spheriques pour la realisation du calcul
(pour les grands domaines de calcul).
Attention : cette option est etroitement liee au
maillage qui doit avoir
ete saisi sur une carte marine en projection de Mercator.
Il faut de plus relever sur la carte la
\telkey{LATITUDE DU POINT ORIGINE}
qui correspond dans le maillage a l''ordonnee $y$ = 0.""",
            ang = """Selection of spherical coordinates to perform the computation
(for large computation domains).
Warning: this option is closely related to the mesh that should have
been entered onto a nautical chart drawn as per Mercator projection
The \telkey{LATITUDE OF ORIGIN POINT}, which corresponds to
ordinate $y$ = 0 in the mesh, must moreover be given.""",
        ),
#       -----------------------------------
        LATITUDE_OF_ORIGIN_POINT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 48.,
            fr = """Determine l''origine utilisee pour le calcul de latitudes
lorsque l''on effectue un calcul en coordonnees spheriques.
Egalement utilise pour les echanges thermiques avec l atmosphere,
le modele d Okada pour les tsunamis, le frazil.""",
            ang = """Determines the origin used for computing latitudes when
a computation is made in spherical coordinates.
This latitude is in particular used to compute the Coriolis force.
In cartesian coordinates, Coriolis coefficient is considered constant.
Also used for heat exchange with atmosphere, Okada model for tsunamis,
frazil.""",
        ),
#       -----------------------------------
        LONGITUDE_OF_ORIGIN_POINT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Fixe la valeur de la longitude du point origine du modele,
lors de l''utilisation de la force generatrice de la maree.
Pour la projection de Mercator, voir le mot cle
\telkey{SYSTEME GEOGRAPHIQUE}.
Egalement utilise pour la force generatrice de la maree,
les echanges thermiques avec l atmosphere,
le modele d Okada pour les tsunamis.""",
            ang = """Give the value of the longitude of the origin point
of the model, when taking into account of the tide generator force.
For the Mercator projection, see the keyword
\telkey{GEOGRAPHIC SYSTEM}.
Also used for tide generating force, heat exchange with atmosphere,
Okada model for tsunamis.""",
        ),
#       -----------------------------------
        NORTH = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Angle que fait le nord en degres, par rapport a l''axe des $y$
dans le sens trigonometrique. Notation decimale 10.5 signifie 10 degres
et 30 minutes.
Lu mais non utilise.""",
            ang = """Angle of the North with the $y$ axis, counted counter-clockwise,
in degrees. 10.5 means 10 degrees and 30 minutes. Read but not used.""",
        ),
#       -----------------------------------
        SPATIAL_PROJECTION_TYPE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["CARTESIAN, NOT GEOREFERENCED","MERCATOR","LATITUDE LONGITUDE"],
            defaut = "CARTESIAN, NOT GEOREFERENCED",
            fr = """Permet de specifier le type de projection spatiale utilisee dans
le cas de l''utilisation des coordonnees spheriques par exemple.
Les choix possibles sont :
\begin{itemize}
\item 1 : Lambert Cartesien non georeference ;
\item 2 : Mercator ;
\item 3 : Latitude/longitude (exprimees en degres).
\end{itemize}
Option 2 ou 3 obligatoire pour les coordonnees spheriques.
Option 3 : latitude et longitude en degres !
Dans le cas de l''option 3, \telemac{2d} convertit les informations
latitude/longitude a l''aide de la projection de Mercator.""",
            ang = """Specifies the type of spatial projection used
(for example when using spherical coordinates).
Possible choices are:
\begin{itemize}
\item 1: Cartesian, not georeferenced,
\item 2: Mercator,
\item 3: latitude/longitude (in degrees).
\end{itemize}
Option 2 or 3 mandatory for spherical coordinates. Option 3: latitude
and longitude in degrees! When using option 3, the coordinates are
automatically
treated by \telemac{2d} using Mercator projection.""",
        ),
    ),
)
# -----------------------------------------------------------------------
TURBULENCE = PROC(nom= "TURBULENCE",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    VELOCITY_DIFFUSIVITY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-6,
        fr = """Fixe de facon uniforme pour l''ensemble du domaine,
la valeur du coefficient de diffusion de viscosite globale (dynamique +
turbulente). Cette valeur peut avoir une influence non negligeable sur
la forme et la taille des recirculations.""",
        ang = """Sets, in an even way for the whole domain, the value of the
coefficient of global (dynamic+turbulent) viscosity. This value may
have a significant effect both on the shapes and sizes of
recirculation zones.""",
    ),
#   -----------------------------------
    TURBULENCE_MODEL = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        into = ["CONSTANT VISCOSITY","ELDER","K-EPSILON MODEL","SMAGORINSKI","MIXING LENGTH","SPALART-ALLMARAS"],
        defaut = "CONSTANT VISCOSITY",
        fr = """6 choix sont possibles actuellement :
\begin{itemize}
\item 1 : viscosite constante ;
\item 2 : modele de Elder ;
\item 3 : modele $k$-$\epsilon$ ;
\item 4 : modele de Smagorinski ;
\item 5 : modele de longueur de melange ;
\item 6 : modele de Spalart-Allmaras.
\end{itemize}
Attention : si on choisit l''option 1
il ne faut pas oublier d''ajuster la valeur du mot-cle
\telkey{COEFFICIENT DE DIFFUSION DES VITESSES}.
Si on choisit l''option 2,
il ne faut pas oublier d''ajuster les deux valeurs du mot-cle :
\telkey{COEFFICIENTS ADIMENSIONNELS DE DISPERSION}.
Si on choisit l''option 3,
ce meme parametre doit retrouver sa vraie valeur physique car elle est
utilisee comme telle dans le modele de turbulence.""",
        ang = """The current alternatives are as follows:
\begin{itemize}
\item 1: constant viscosity,
\item 2: elder''s model,
\item 3: $k$-$\epsilon$ model,
\item 4: Smagorinski model,
\item 5: mixing length model,
\item 6: Spalart-Allmaras model.
\end{itemize}
NOTE: when option 1 is chosen, it should be kept in mind that the
value of the keyword \telkey{VELOCITY DIFFUSIVITY} has to be ajusted.
When option 2 is chosen, the two values of keyword
\telkey{NON-DIMENSIONAL DISPERSION COEFFICIENTS} are used.
When option 3 is chosen, this parameter should recover its true
physical value, since it is used as such in the turbulence model.""",
    ),
#   -----------------------------------
    b_TURBULENCE_MODELG = BLOC(condition="TURBULENCE_MODEL == 'CONSTANT VISCOSITY'",
#   -----------------------------------
    ),
#   -----------------------------------
    b_TURBULENCE_MODELH = BLOC(condition="TURBULENCE_MODEL == 'Elder'",
#   -----------------------------------
#       -----------------------------------
        NON_DIMENSIONAL_DISPERSION_COEFFICIENTS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min= 2, max= 2,
            defaut = [6.,0.6],
            fr = """Coefficients longitudinal et transversal dans la formule de Elder.
Utilises uniquement avec \telkey{MODELE DE TURBULENCE} = 2.""",
            ang = """Longitudinal and transversal coefficients in Elder s formula.
Used only with \telkey{TURBULENCE MODEL} = 2.""",
        ),
    ),
#   -----------------------------------
    ACCURACY_OF_SPALART_ALLMARAS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-9,
        fr = """Fixe la precision demandee sur le modele de Spalart-Allmaras pour
le test d''arret dans l''etape de diffusion et termes sources de
$\tilde{\nu}$.""",
        ang = """Sets the required accuracy for the Spalart-Allmaras model in
the diffusion and source terms step of the $\tilde{\nu}$.""",
    ),
#   -----------------------------------
    INFORMATION_ABOUT_SPALART_ALLMARAS_MODEL = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Si oui, les informations du solveur du modele de Spalart-Allmaras
sont imprimees.""",
        ang = """If yes, informations about solver of Spalart-Allmaras model
are printed to the listing.""",
    ),
#   -----------------------------------
    ADVANCED = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TURBULENCE_REGIME_FOR_SOLID_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["smooth","rough"],
            defaut = "rough",
            fr = """ Permet de choisir le regime de turbulence aux parois.
Les choix possibles sont :
\begin{itemize}
\item 1 : regime turbulent lisse ;
\item 2 : regime turbulent rugueux.
\end{itemize}""",
            ang = """ Provided for selecting the type of friction on the walls.
Possible choices are:
\begin{itemize}
\item 1: smooth,
\item 2: rough.
\end{itemize}""",
        ),
#       -----------------------------------
        INFORMATION_ABOUT_K_EPSILON_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Donne le nombre d''iterations du solveur de l''etape de
diffusion et termes sources du modele $k-\epsilon$.""",
            ang = """Gives the number of iterations of the solver in the diffusion
and source terms step of the $k-\epsilon$ model.""",
        ),
#       -----------------------------------
        ADVECTION_OF_K_AND_EPSILON = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non de la convection de $k$ et $\epsilon$
(pour modele $k-\epsilon$) ou $\tilde{\nu}$
(pour modele de Spalart-Allmaras).""",
            ang = """The $k$ and $\epsilon$ advection is taken into account or ignored
(for $k-\epsilon$ model) or $\tilde{\nu}$ advection
(for Spalart-Allmaras model).""",
        ),
#       -----------------------------------
        b_ADVECTION_OF_K_AND_EPSILONG = BLOC(condition="ADVECTION_OF_K_AND_EPSILON == True",
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_K_EPSILON = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO ADVECTION","CHARACTERISTICS","EXPLICIT + SUPG","EXPLICIT LEO POSTMA","EXPLICIT + MURD SCHEME N","EXPLICIT + MURD SCHEME PSI","LEO POSTMA FOR TIDAL FLATS","N-SCHEME FOR TIDAL FLATS","ERIA SCHEME FOR TIDAL FLATS"],
                defaut = "CHARACTERISTICS",
                fr = """Choix du schema de convection pour $k$ et $\epsilon$ (pour modele
$k-\epsilon$) ou $\tilde{\nu}$ (pour modele de Spalart-Allmaras),
remplace \telkey{FORME DE LA CONVECTION}.""",
                ang = """Choice of the advection scheme for $k$ and $\epsilon$ (for
$k-\epsilon$ model) or $\tilde{\nu}$ (for Spalart-Allmaras model),
replaces \telkey{TYPE OF ADVECTION}.""",
            ),
        ),
#       -----------------------------------
        SCHEME_OPTION_FOR_ADVECTION_OF_K_EPSILON = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Si present remplace et a priorite sur :
\telkey{OPTION POUR LES CARACTERISTIQUES}
\telkey{OPTION DE SUPG}.
Si caracteristiques :
\begin{itemize}
\item 1 = forme forte,
\item 2 = forme faible.
\end{itemize}
Si schema PSI ou N :
\begin{itemize}
\item 1 = explicite ;
\item 2 = predicteur-correcteur ;
\item 3 = predicteur-correcteur 2e ordre en temps ;
\item 4 = implicite.
\end{itemize}
Mot-cle commun pour les variables $k$, $\epsilon$
(pour modele $k-\epsilon$)
et $\tilde{\nu}$ (pour modele de Spalart-Allmaras).""",
            ang = """If present replaces and has priority over:
\telkey{OPTION FOR CHARACTERISTICS}
\telkey{SUPG OPTION}.
If characteristics:
\begin{itemize}
\item 1 = strong form,
\item 2 = weak form.
\end{itemize}
If N or PSI scheme:
\begin{itemize}
\item 1 = explicit,
\item 2 = predictor-corrector,
\item 3 = predictor-corrector second-order in time,
\item 4 = implicit.
\end{itemize}
Common keyword for variables $k$, $\epsilon$ (for $k-\epsilon$ model)
and $\tilde{\nu}$ (for Spalart-Allmaras model).""",
        ),
#       -----------------------------------
        MIXING_LENGTH_MODEL_COEFFICIENTS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min= 2, max= 2,
            defaut = [0.1066667,0.0666667],
            fr = """Coefficients de calage $C_l$ et $\alpha_t$ dans la formule de longueur
de melange.
Utilises uniquement avec \telkey{MODELE DE TURBULENCE} = 5.""",
            ang = """Calibration coefficients $C_l$ and $\alpha_t$ in mixing length formula.
Only used with \telkey{TURBULENCE MODEL} = 5.""",
        ),
    ),
#   -----------------------------------
    SOLVER_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SOLVER_FOR_K_EPSILON_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["conjugate gradient","conjugate residual","conjugate gradient on normal equation","minimum error","squared conjugate gradient","conjugate gradient squared stabilised (cgstab)","gmres (see option for the solver for k-epsilon model)","direct"],
            defaut = "conjugate gradient",
            fr = """Permet de choisir le solveur utilise pour la resolution
du systeme de diffusion de $k$, $\epsilon$ (pour modele $k-\epsilon$)
ou $\tilde{\nu}$ (pour modele de Spalart-Allmaras).
Les choix possibles sont :
\begin{itemize}
\item 1 : gradient conjugue,
\item 2 : residu conjugue,
\item 3 : gradient conjugue sur equation normale,
\item 4 : erreur minimale,
\item 5 : gradient conjugue carre (non programme),
\item 6 : gradient conjugue carre stabilise (cgstab),
\item 7 : GMRES (voir aussi \telkey{OPTION DU SOLVEUR}
pour le modele $k$-$\epsilon$),
\item 8 : direct.
\end{itemize}""",
            ang = """Makes it possible to select the solver used for solving
the system of the diffusion of $k$, $\epsilon$ (for $k-\epsilon$ model)
or $\tilde{\nu}$ (for Spalart-Allmaras model).
Possible choices are:
\begin{itemize}
\item 1: conjugate gradient,
\item 2: conjugate residual,
\item 3: conjugate gradient on a normal equation,
\item 4: minimum error,
\item 5: squared conjugate gradient (not implemented),
\item 6: conjugate gradient squared stabilised (cgstab),
\item 7: GMRES (see \telkey{OPTION FOR THE SOLVER FOR K-EPSILON MODEL}),
\item 8: direct.
\end{itemize}""",
        ),
#       -----------------------------------
        OPTION_FOR_THE_SOLVER_FOR_K_EPSILON_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 2,
            fr = """Si le solveur est GMRES (7), le mot cle est la dimension de
l''espace de Krylov (valeurs conseillees entre 2 et 15).
Mot-cle commun pour les variables $k$, $\epsilon$
(pour modele $k-\epsilon$)
et $\tilde{\nu}$ (pour modele de Spalart-Allmaras).""",
            ang = """When GMRES (7) is chosen for solver, dimension of the Krylov space.
Try values between 2 and 15.
Common keyword for variables $k$, $\epsilon$ (for $k-\epsilon$ model)
and $\tilde{\nu}$ (for Spalart-Allmaras model).""",
        ),
#       -----------------------------------
        PRECONDITIONING_FOR_K_EPSILON_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["diagonal","no preconditioning","Crout","diagonal and Crout"],
            defaut = "diagonal",
            fr = """Permet de preconditionner le systeme relatif a la diffusion de
$k$, $\epsilon$ (pour modele $k-\epsilon$) ou $\tilde{\nu}$
(pour modele de Spalart-Allmaras)
afin d''accelerer la convergence lors de sa resolution.
\begin{itemize}
\item 0 : pas de preconditionnement ;
\item 2 : preconditionnement diagonal ;
\item 7 : preconditionnement de Crout par element ou segment
(ne marche pas en parallele).
\end{itemize}
Certains preconditionnements sont cumulables
(les diagonaux 2 ou 3 avec les autres)
Pour cette raison on ne retient que les nombres premiers pour
designer les preconditionnements. Si l''on souhaite en cumuler
plusieurs on formera le produit des options correspondantes.""",
            ang = """Choice of the preconditioning of the linear system in the diffusion
step of $k$, $\epsilon$ (for $k-\epsilon$ model) or $\tilde{\nu}$
(for Spalart-Allmaras model) so that the convergence is speeded up
when it is being solved.
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal preconditioning,
\item 7: Crout''s preconditioning per element or segment
(does not work in parallel).
\end{itemize}
Some operations (either 2 or 3 diagonal preconditioning) can be
performed concurrently with the others.
Only prime numbers are therefore kept to denote the preconditioning
operations. When several of them are to be performed concurrently,
the product of relevant options shall be made.""",
        ),
    ),
#   -----------------------------------
    ACCURACY = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ACCURACY_OF_K = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-9,
            fr = """Fixe la precision demandee sur $k$ pour le test d''arret dans
l''etape de diffusion et termes sources de l equation sur $k$.""",
            ang = """Sets the required accuracy for computing $k$ in the diffusion
and source terms step of the $k$ transport equation.""",
        ),
#       -----------------------------------
        ACCURACY_OF_EPSILON = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-9,
            fr = """Fixe la precision demandee sur $\epsilon$ pour le test d''arret dans
l''etape de diffusion et termes sources de l equation sur $\epsilon$.""",
            ang = """Sets the required accuracy for computing $\epsilon$ in the diffusion
and source terms step of the $\epsilon$ transport equation.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_ITERATIONS_FOR_K_AND_EPSILON = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 50,
            fr = """Fixe le nombre maximum d''iterations accepte lors de la resolution
du systeme diffusion-termes sources de $k$, $\epsilon$ (pour modele
$k-\epsilon$) ou $\tilde{\nu}$ (pour modele de Spalart-Allmaras).""",
            ang = """Sets the maximum number of iterations that are acceptable when
solving the diffusion source-terms step of $k$, $\epsilon$ (for
$k-\epsilon$ model) or $\tilde{\nu}$ (for Spalart-Allmaras model).""",
        ),
    ),
)
# -----------------------------------------------------------------------
TIDAL_FLATS_INFO = PROC(nom= "TIDAL_FLATS_INFO",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    TIDAL_FLATS = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Permet de supprimer les tests sur les bancs decouvrants, dans
les cas ou l''on est certain qu''il n''y en aura pas.
En cas de doute, utiliser OUI.""",
        ang = """When NO, the specific treatments for tidal flats are by-passed.
This spares time, but of course you must be sure that you have no
tidal flats.""",
    ),
#   -----------------------------------
    b_TIDAL_FLATSG = BLOC(condition="TIDAL_FLATS == True",
#   -----------------------------------
#       -----------------------------------
        OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["EQUATIONS SOLVED EVERYWHERE WITH CORRECTION ON TIDAL FLATS","DRY ELEMENTS FROZEN","LIKE 1 BUT WITH POROSITY (DEFINA METHOD)"],
            defaut = "EQUATIONS SOLVED EVERYWHERE WITH CORRECTION ON TIDAL FLATS",
            fr = """Utilise si \telkey{BANCS DECOUVRANTS} est OUI.
Les choix possibles sont :
\begin{itemize}
\item 1 : equations resolues partout avec correction
sur les bancs decouvrants (correction du gradient de surface libre) ;
\item 2 : gel des elements decouvrants (zones de bancs decouvrants
sont masquees).
Attention : la conservation de la masse peut etre alteree ;
\item 3 : comme 1 mais avec porosite (methode defina).
\end{itemize}""",
            ang = """Used if \telkey{TIDAL FLATS} is YES.
Possible choices are:
\begin{itemize}
\item 1: equations solved everywhere with correction on tidal flats
(corrected free surface gradient),
\item 2: dry elements are frozen (tidal flats area are masked).
Warning: mass-conservation may be altered,
\item 3: like 1 but with porosity (defina method).
\end{itemize}""",
        ),
#       -----------------------------------
        b_OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATSG = BLOC(condition="OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATS == 'EQUATIONS SOLVED EVERYWHERE WITH CORRECTION ON TIDAL FLATS'",
#       -----------------------------------
#           -----------------------------------
            TREATMENT_OF_NEGATIVE_DEPTHS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO TREATMENT","SMOOTHING","FLUX CONTROL","FLUX CONTROL (ERIA)"],
                defaut = "SMOOTHING",
                fr = """Seulement avec
\telkey{OPTION DE TRAITEMENT DES BANCS DECOUVRANTS} = 1.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de traitement ;
\item 1 : lissage ;
\item 2 : limitation des flux, approche par segment ;
\item 3 : limitation des flux, approche par triangle.
\end{itemize}
Si les options 2 ou 3 avec bancs decouvrants sont utilisees,
il est obligatoire d avoir \telkey{MASS-LUMPING SUR H} = 1.
+ \telkey{CORRECTION DE CONTINUITE} = OUI
+ \telkey{OPTION DE SUPG} pour hauteur d eau = 0
(pas de decentrement SUPG sur la hauteur d eau).""",
                ang = """Only with
\telkey{OPTION FOR THE TREATMENT OF TIDAL FLATS} = 1.
Possible choices are:
\begin{itemize}
\item 0: no treatment,
\item 1: smoothing,
\item 2: flux control, by segment,
\item 3: flux control, by element.
\end{itemize}
If using options 2 or 3 with tidal flats, it is mandatory to set
\telkey{MASS-LUMPING ON H} = 1. + \telkey{CONTINUITY CORRECTION} = YES +
\telkey{SUPG OPTION} for water depth = 0 (no SUPG upwinding on depth).""",
            ),
        ),
#       -----------------------------------
        THRESHOLD_FOR_NEGATIVE_DEPTHS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """En dessous du seuil, les hauteurs negatives sont lissees.
Seulement utilise avec
\telkey{TRAITEMENT DES HAUTEURS NEGATIVES} = 1.""",
            ang = """Below the threshold the negative depths are smoothed.
Only used with \telkey{TREATMENT OF NEGATIVE DEPTHS} = 1.""",
        ),
#       -----------------------------------
        THRESHOLD_DEPTH_FOR_RECEDING_PROCEDURE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Si > 0., declenche la procedure de ressuyage qui evite le
franchissement parasite des digues mal discretisees.""",
            ang = """If > 0., will trigger the receeding procedure that avoids overwhelming
of dykes which are too loosely discretised.""",
        ),
    ),
#   -----------------------------------
    H_CLIPPING = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Determine si l''on desire ou non limiter par valeur inferieure
la hauteur d''eau $H$ (dans le cas des bancs decouvrants par exemple).""",
        ang = """Determines whether limiting the water depth $H$ by a lower value
desirable or not (for instance in the case of tidal flats).
This keyword may have an influence on mass conservation since
the truncation of depth is equivalent to adding mass.""",
    ),
#   -----------------------------------
    b_H_CLIPPINGG = BLOC(condition="H_CLIPPING == True",
#   -----------------------------------
#       -----------------------------------
        MINIMUM_VALUE_OF_DEPTH = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Fixe la valeur minimale de $H$ lorsque l''option \telkey{CLIPPING DE H}
est activee.""",
            ang = """Sets the minimum water depth $H$ value when option \telkey{H CLIPPING}
is implemented. Not fully implemented.""",
        ),
    ),
)
# -----------------------------------------------------------------------
TRACERS = PROC(nom= "TRACERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    BOUNDARY_CONDITIONS_FOR_TRACERS = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        PRESCRIBED_TRACERS_VALUES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', max='**',
            fr = """Valeurs du traceur imposees aux frontieres liquides entrantes.
Lire la partie du mode d''emploi consacree aux conditions aux limites.""",
            ang = """Tracer values prescribed at the inflow boundaries.
Read the user manual section dealing with the boundary conditions.""",
        ),
    ),
#   -----------------------------------
    SETTING = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Definit le nombre de traceurs.""",
            ang = """Defines the number of tracers""",
        ),
#       -----------------------------------
        NAMES_OF_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            fr = """Noms des traceurs en 32 caracteres, 16 pour le nom 16 pour l''unite.""",
            ang = """Name of tracers in 32 characters, 16 for the name, 16 for the unit.""",
        ),
#       -----------------------------------
        INITIAL_VALUES_OF_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [0.,0.],
            fr = """Fixe la valeur initiale du/des traceur(s).
Les valeurs requises sont separees par un point virgule ; si plus d une.
Le nombre de valeurs fournies doit etre egal au nombre de traceurs
declares.""",
            ang = """Sets the initial value of the tracer(s).
Required value(s) separated with a semicolumn ; if more than one.
The number of supplied values must be equal to the number of declared
tracers.""",
        ),
#       -----------------------------------
        DENSITY_EFFECTS = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Prise en compte du gradient horizontal de densite.
Le premier traceur est alors la salinite.
$\rho_{\rm{eau}} = 999.972.(1-7.10^{-6}(T_{\rm{moy}}-4)^2)$.""",
            ang = """The horizontal gradient of density is taken into account.
The 1st tracer is then the salinity.
$\rho_{\rm{water}} = 999.972.(1-7.10^{-6}(T_{\rm{mean}}-4)^2)$.""",
        ),
#       -----------------------------------
        b_DENSITY_EFFECTSG = BLOC(condition="DENSITY_EFFECTS == True",
#       -----------------------------------
#           -----------------------------------
            MEAN_TEMPERATURE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 20.,
                fr = """Temperature de reference pour le calcul des effets de densite.
A utiliser avec le mot-cle \telkey{EFFETS DE DENSITE}.""",
                ang = """Reference temperature for density effects.
To be used with the keyword \telkey{DENSITY EFFECTS}.""",
            ),
        ),
#       -----------------------------------
        b_DENSITY_EFFECTSH = BLOC(condition="DENSITY_EFFECTS == True",
#       -----------------------------------
#           -----------------------------------
            Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#           -----------------------------------
                defaut = "The first tracer must be the salinity in kg/m3"),
        ),
    ),
#   -----------------------------------
    SOLVER_TRA = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SOLVER_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres (see solver option for tracers diffusion)","direct"],
            defaut = ["conjugate gradient","conjugate gradient"],
            fr = """Permet de choisir le solveur utilise pour la resolution
du systeme de diffusion de traceur(s) :
Les choix possibles sont :
\begin{itemize}
\item 1 : gradient conjugue,
\item 2 : residu conjugue,
\item 3 : gradient conjugue sur equation normale,
\item 4 : erreur minimale,
\item 5 : gradient conjugue carre (non programme),
\item 6 : gradient conjugue carre stabilise (cgstab),
\item 7 : GMRES (voir aussi
\telkey{OPTION DU SOLVEUR POUR LA DIFFUSION DES TRACEURS}),
\item 8 : direct.
\end{itemize}""",
            ang = """Makes it possible to select the solver used for solving
the system of tracer(s) diffusion.
Possible choices are:
\begin{itemize}
\item 1: conjugate gradient,
\item 2: conjugate residual,
\item 3: conjugate gradient on a normal equation,
\item 4: minimum error,
\item 5: squared conjugate gradient (not implemented),
\item 6: cgstab,
\item 7: GMRES (see \telkey{SOLVER OPTION FOR TRACERS DIFFUSION},
\item 8: direct.
\end{itemize}""",
        ),
#       -----------------------------------
        SOLVER_OPTION_FOR_TRACERS_DIFFUSION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I', min=0, max='**',
            defaut = [2],
            fr = """Si le solveur est GMRES (7), le mot cle est la dimension de
l''espace de Krylov (valeurs conseillees entre 2 et 15).""",
            ang = """When GMRES (7) is chosen for solver, dimension of the Krylov space.
Try values between 2 and 15.""",
        ),
#       -----------------------------------
        PRECONDITIONING_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["no preconditioning","diagonal","Crout","diagonal and Crout"],
            defaut = ["diagonal"],
            fr = """Permet de preconditionner le systeme relatif au traceur.
afin d''accelerer la convergence lors de sa resolution.
\begin{itemize}
\item 0 : pas de preconditionnement ;
\item 2 : preconditionnement diagonal ;
\item 7 : preconditionnement de Crout par element ou segment
(ne marche pas en parallele).
\end{itemize}
Certains preconditionnements sont cumulables
(les diagonaux 2 ou 3 avec les autres)
Pour cette raison on ne retient que les nombres premiers pour
designer les preconditionnements. Si l''on souhaite en cumuler
plusieurs on formera le produit des options correspondantes.""",
            ang = """Choice of the preconditioning of the linear system of the tracer
diffusion so that the convergence is speeded up when it is being solved.
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal preconditioning,
\item 7: Crout''s preconditioning per element or segment
(does not work in parallel).
\end{itemize}
Some operations (either 2 or 3 diagonal preconditioning) can be
performed concurrently with the others.
Only prime numbers are therefore kept to denote the preconditioning
operations. When several of them are to be performed concurrently,
the product of relevant options shall be made.""",
        ),
    ),
#   -----------------------------------
    ACCURACY_TRA = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ACCURACY_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-6,
            fr = """Fixe la precision demandee pour le calcul de la diffusion
du traceur.""",
            ang = """Sets the required accuracy for computing the tracer
diffusion.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_ITERATIONS_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 60,
            fr = """Limite le nombre d''iterations du solveur a chaque pas de temps pour
le calcul de la diffusion du ou des traceur(s).""",
            ang = """Limits the number of solver iterations at each time step for
the diffusion of tracer(s).""",
        ),
    ),
#   -----------------------------------
    SOURCES_TRA = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        VALUES_OF_THE_TRACERS_AT_THE_SOURCES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Valeurs des traceurs a chacune des sources.
toutes les sources pour le premier traceur
puis toutes les sources du deuxieme traceur, etc.
(cf. manuel utilisateur).
Par exemple, s''il y a 3 traceurs (T1, T2 et T3) et 2 sources
(S1 et S2), la syntaxe suivante est utilisee :\\
S1\_T1;S1\_T2;S1\_T3;S2\_T1;S2\_T2;S2\_T3\\
10.0; 10.0; 0.0;  0.0; 10.0; 10.0""",
            ang = """Values of the tracers at the sources.
All sources for the first tracer, then
all sources for the second tracer, etc.
(see user manual).
For example, if there are 3 tracers (T1, T2 and T3)
and 2 sources (S1 and S2), the following syntax is used:\\
S1\_T1;S1\_T2;S1\_T3;S2\_T1;S2\_T2;S2\_T3\\
10.0; 10.0; 0.0;  0.0; 10.0; 10.0""",
        ),
    ),
#   -----------------------------------
    METEOROLOGY_TRA = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        VALUES_OF_TRACERS_IN_THE_RAIN = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Generalement, ce traceur est la temperature, dans ce cas
cette valeur est a modifier, sinon la valeur 0 est raisonnable.""",
            ang = """Most often, this tracer is temperature, in this case
this value should be modified, otherwise, default value of 0 seems
reasonable.""",
        ),
    ),
#   -----------------------------------
    NUMERICAL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ADVECTION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non de la convection du traceur passif.""",
            ang = """The advection of the passive tracer is taken into account
or ignored.""",
        ),
#       -----------------------------------
        b_ADVECTION_OF_TRACERSG = BLOC(condition="ADVECTION_OF_TRACERS == True",
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_TRACERS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["NO ADVECTION","CHARACTERISTICS","EXPLICIT + SUPG","EXPLICIT LEO POSTMA","EXPLICIT + MURD SCHEME N","EXPLICIT + MURD SCHEME PSI","LEO POSTMA FOR TIDAL FLATS","N-SCHEME FOR TIDAL FLATS","ERIA SCHEME FOR TIDAL FLATS"],
                defaut = ["CHARACTERISTICS"],
                fr = """Choix du schema de convection pour les traceurs,
remplace \telkey{FORME DE LA CONVECTION}.""",
                ang = """Choice of the advection scheme for the tracers,
replaces \telkey{TYPE OF ADVECTION}.""",
            ),
        ),
#       -----------------------------------
        IMPLICITATION_COEFFICIENT_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.6,
            fr = """Fixe la valeur du coefficient d''implicitation du traceur.
Si un schema de convection pour les traceurs est un schema distributif
(par ex : 3, 4, 5, 13, 14 ou 15),
\telkey{COEFFICIENT D''IMPLICITATION DES TRACEURS} est impose a 0.
(explicite).""",
            ang = """Sets the value of the implicitation coefficient
for the tracer.
If an advection scheme for tracers is a distributive scheme
(e.g.: 3, 4, 5, 13, 14 or 15),
\telkey{IMPLICITATION COEFFICIENT OF TRACERS} is prescribed at 0.
(explicit).""",
        ),
#       -----------------------------------
        DIFFUSION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non de la diffusion du traceur passif.""",
            ang = """The diffusion of the passive tracer is taken into account
or ignored.""",
        ),
#       -----------------------------------
        b_DIFFUSION_OF_TRACERSG = BLOC(condition="DIFFUSION_OF_TRACERS == True",
#       -----------------------------------
#           -----------------------------------
            COEFFICIENT_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                defaut = [1.E-6],
                fr = """Fixe la valeur du coefficient de diffusion du traceur.
L''influence de ce parametre sur l''evolution des traceurs dans
le temps est importante.
C est un tableau depuis la version 8.2, avec une valeur par traceur,
separation par un point virgule.""",
                ang = """Sets the value of the tracer diffusivity.
These values may have a significant effect on the evolution of
tracers in time.
Since release 8.2, it has been an array, with one value per tracer,
separated by semicolons.""",
            ),
        ),
#       -----------------------------------
        OPTION_FOR_THE_DIFFUSION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["div( nu grad(T) )","1/h div ( h nu grad(T)"],
            defaut = ["div( nu grad(T) )"],
            fr = """Choix possibles :
\begin{itemize}
\item 1: Diffusion de la forme div( nu grad(T) ),
\item 2: Diffusion de la forme 1/h div ( h nu grad(T) ).
\end{itemize}""",
            ang = """Possible choices:
\begin{itemize}
\item 1: Diffusion in the form div( nu grad(T) ),
\item 2: Diffusion in the form 1/h div ( h nu grad(T) ).
\end{itemize}""",
        ),
#       -----------------------------------
        SCHEME_OPTION_FOR_ADVECTION_OF_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I', min=0, max='**',
            defaut = [1],
            fr = """Si present remplace et a priorite sur :
\telkey{OPTION POUR LES CARACTERISTIQUES}
\telkey{OPTION DE SUPG}.
Si caracteristiques :
\begin{itemize}
\item 1 = forme forte,
\item 2 = forme faible.
\end{itemize}
Si schema PSI ou N :
\begin{itemize}
\item 1 = explicite ;
\item 2 = predicteur-correcteur ;
\item 3 = predicteur-correcteur 2e ordre en temps ;
\item 4 = implicite.
\end{itemize}""",
            ang = """If present replaces and has priority over:
\telkey{OPTION FOR CHARACTERISTICS}
\telkey{SUPG OPTION}.
If characteristics:
\begin{itemize}
\item 1 = strong form,
\item 2 = weak form.
\end{itemize}
If N or PSI scheme:
\begin{itemize}
\item 1 = explicit,
\item 2 = predictor-corrector,
\item 3 = predictor-corrector second-order in time,
\item 4 = implicit.
\end{itemize}""",
        ),
#       -----------------------------------
        MASS_LUMPING_ON_TRACERS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Fixe le taux de mass-lumping effectue sur le traceur.
Lu mais remplace par la valeur de \telkey{MASS-LUMPING SUR H}
pour assurer la conservation de la masse de traceur.""",
            ang = """Sets the amount of mass-lumping that is performed on
the tracer.
Read but replaced by the value of \telkey{MASS-LUMPING ON H}
to ensure tracer mass conservation.""",
        ),
    ),
)
# -----------------------------------------------------------------------
PARTICLE_TRANSPORT = PROC(nom= "PARTICLE_TRANSPORT",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    DROGUES = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        MAXIMUM_NUMBER_OF_DROGUES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Permet d''effectuer un suivi de flotteurs.""",
            ang = """Maximum number of drogues in the computation.""",
        ),
#       -----------------------------------
        b_MAXIMUM_NUMBER_OF_DROGUESG = BLOC(condition="MAXIMUM_NUMBER_OF_DROGUES != 0",
#       -----------------------------------
#           -----------------------------------
            ASCII_DROGUES_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Fichier de resultats ASCII avec les positions des flotteurs.""",
                ang = """ASCII results file with positions of drogues.""",
            ),
#           -----------------------------------
            PRINTOUT_PERIOD_FOR_DROGUES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Nombre de pas de temps entre 2 sorties de positions de
flotteurs dans le fichier des resultats supplementaire.
N affecte pas la qualite du calcul de la trajectoire.""",
                ang = """Number of time steps between 2 outputs of drogues
positions in the output file.""",
            ),
        ),
#       -----------------------------------
        INITIAL_DROGUES_SAMPLING_DENSITY = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [1000,1000],
            fr = """Densite initiale utilisee a la repartition des flotteurs
au debut (ou en cours) de la simulation.""",
            ang = """Initial density of drogues, or number of drogues per m2,
used to spatially place the initial drogues in a simulation.""",
        ),
#       -----------------------------------
        BINARY_DROGUES_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats binaire avec les positions des flotteurs.""",
            ang = """Binary results file with positions of drogues.""",
        ),
#       -----------------------------------
        DROGUES_FILE_FORMAT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['BKBINPCL','TECPLOT'],
            defaut = 'TECPLOT',
            fr = """Format du \telkey{FICHIER BINAIRE DES FLOTTEURS}.
Les valeurs possibles sont :
\begin{itemize}
\item BKBINPCL: format PCL binaire natif de BlueKenue;
\item TECPLOT: format Tecplot original (ASCII).
\end{itemize}""",
            ang = """Format of the \telkey{BINARY DROGUES FILE}.
Possible choices are:
\begin{itemize}
\item BKBINPCL: format binary PCL native to BlueKenue,
\item TECPLOT: format Tecplot original (ASCII).
\end{itemize}""",
        ),
#       -----------------------------------
        PREVIOUS_DROGUES_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom d''un fichier contenant les resultats d''un calcul
precedent avec flotteurs.""",
            ang = """Name of a file containing the results of an earlier
computation with drogues.""",
        ),
#       -----------------------------------
        PREVIOUS_DROGUES_FILE_FORMAT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['BKBINPCL','BKASCPCL'],
            defaut = 'BKBINPCL',
            fr = """Format du \telkey{FICHIER DES FLOTTEURS PRECEDENT}.
Les valeurs possibles sont :
\begin{itemize}
\item BKBINPCL: format PCL binaire natif de BlueKenue;
\item BKASCPCL: format PCL ASCII natif de BlueKenue.
\end{itemize}""",
            ang = """Format of the \telkey{PREVIOUS DROGUES FILE}.
Possible choices are:
\begin{itemize}
\item BKBINPCL: format binary PCL native of BlueKenue,
\item BKASCPCL: format SCII PCL native of BlueKenue.
\end{itemize}""",
        ),
#       -----------------------------------
        DROGUES_INITIAL_POSITIONING_DATA_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier de donnees en ASCII contenant les coordonnees des polygones
ou des points qui vont permettre un positionnement des particles
au depart de la simulation.""",
            ang = """ASCII data file containing polygons or points defining the
 initial positioning of drogues at the start of the simulation.""",
        ),
#       -----------------------------------
        FORMAT_OF_THE_DROGUES_POSITIONING_DATA_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['BKASCI2S'],
            defaut = 'BKASCI2S',
            fr = """Format du
 \telkey{FICHIER POSITIONNANT LES DROGUES INITIALES}.
La seule valeur possible est :
\begin{itemize}
\item BKASCI2S: format I2S ASCII natif de BlueKenue.
\end{itemize}""",
            ang = """Format of the
 \telkey{DROGUES INITIAL POSITIONING DATA FILE}.
Single possible choice is:
\begin{itemize}
\item BKASCI2S: format AXCII I2S native of BlueKenue.
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    ALGAE = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        ALGAE_TRANSPORT_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si OUI, une partie ou tous les flotteurs seront des algues.""",
            ang = """If YES, some or all the floats or particles will be algae.""",
        ),
#       -----------------------------------
        b_ALGAE_TRANSPORT_MODELG = BLOC(condition="ALGAE_TRANSPORT_MODEL == True",
#       -----------------------------------
#           -----------------------------------
            ALGAE_TYPE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["SPHERE","IRIDAEA FLACCIDA (CLOSE TO ULVA)","PELVETIOPSIS LIMITATA","GIGARTINA LEPTORHYNCHOS"],
                defaut = ["SPHERE","SPHERE"],
                fr = """Type des algues.
Les choix possibles sont :
\begin{itemize}
\item 1 : Sphere ;
\item 2 : Iridaeca Flaccida ;
\item 3 : Pelvetiopsis Limitata ;
\item 4 : Gigartina Leptorhynchos.
\end{itemize}
Pour le choix 1 les algues seront
modelisees comme des spheres, pour les autres choix voir Gaylord
et al. (1994).""",
                ang = """Algae type.
Possible choices are:
\begin{itemize}
\item 1: Sphere,
\item 2: Iridaeca Flaccida,
\item 3: Pelvetiopsis Limitata,
\item 4: Gigartina Leptorhynchos.
\end{itemize}
For choice 1 the algae particles will be
modeled as spheres, and for the other choices see Gaylord et
al. (1994).""",
            ),
#           -----------------------------------
            DIAMETER_OF_ALGAE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                defaut = [0.1,0.1],
                fr = """Diametre des algues en m.""",
                ang = """Diameter of algae in m.""",
            ),
#           -----------------------------------
            DENSITY_OF_ALGAE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                defaut = [1050.,1050.],
                fr = """Masse volumique des algues en kg/m$^3$.""",
                ang = """Density of algae in kg/m$^3$.""",
            ),
#           -----------------------------------
            THICKNESS_OF_ALGAE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                defaut = [0.01,0.01],
                fr = """Epaisseur des algues en m.""",
                ang = """Thickness of algae in m.""",
            ),
        ),
#       -----------------------------------
        NUMBER_OF_ALGAE_CLASSES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Nombre de classes d''algues. Chaque classe sera associee
 a une propriete particuliere.""",
            ang = """Number of algae classes. Each class will be associated
 with a particular property.""",
        ),
#       -----------------------------------
        DURATION_BEFORE_ALGAE_RELEASE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [0.0,0.0],
            fr = """Duree en secondes avant le relachemet des algues depuis le
depart de la simulation.""",
            ang = """Duration in seconds before the release of the algae from
the start of the simulation.""",
        ),
#       -----------------------------------
        ALGAE_RELEASE_TYPE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["TIMED","DISLODGEMENT"],
            defaut = ["TIMED","TIMED"],
            fr = """Les valeurs possibles sont:
\begin{itemize}
\item 1 : PERIODE : Les algues bougent apres une certaine periode de
temps,
\item 2 : DEPLACEMENT : Les algues bougent lorsque la vitesse orbitale
de vague critique depasse un certain seuil.
\end{itemize}""",
            ang = """Possible values are:
\begin{itemize}
\item 1: TIMED       : Algae move after a specified time has elapsed,
\item 2: DISLODGEMENT: Algae move after a critical wave orbital velocity
is exceeded.
\end{itemize}""",
        ),
#       -----------------------------------
        WAVE_ORBITAL_VELOCITY_THRESHOLD_FOR_ALGAE_1 = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [2.,2.],
            fr = """Vitesse orbitale de vague 1 pour le deplacement des algues en m/s.""",
            ang = """Wave orbital velocity 1 for algae dislodgement in m/s.""",
        ),
#       -----------------------------------
        WAVE_ORBITAL_VELOCITY_THRESHOLD_FOR_ALGAE_2 = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [0.,0.],
            fr = """Vitesse orbitale de vague 2 pour le deplacement des algues en m/s.""",
            ang = """Wave orbital velocity 2 for algae dislodgement in m/s.""",
        ),
#       -----------------------------------
        RATE_OF_DEGRADATION_FOR_ALGAE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [0.,0.],
            fr = """Taux de degradation pour les algues.""",
            ang = """Rate of degradation for algae.""",
        ),
    ),
#   -----------------------------------
    OIL_SPILL = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        OIL_SPILL_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Pour declencher le modele de derive de nappes, dans ce cas
le fichier de commandes migrhycar est necessaire.""",
            ang = """Will trigger the oil spill model, in this case the
\telkey{OIL SPILL STEERING FILE} is needed.""",
        ),
#       -----------------------------------
        b_OIL_SPILL_MODELG = BLOC(condition="OIL_SPILL_MODEL == True",
#       -----------------------------------
#           -----------------------------------
            OIL_SPILL_STEERING_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Contient les donnees pour le modele de derive de nappes.""",
                ang = """Contains data for the \telkey{OIL SPILL MODEL}.""",
            ),
        ),
    ),
#   -----------------------------------
    BROWNIAN_MOTION = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        STOCHASTIC_DIFFUSION_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["No model","brownian movement"],
            defaut = "No model",
            fr = """Pour les particules : flotteurs, hydrocarbures.
Si aucune turbulence n est activee, la diffusion stochastique n est pas
prise en compte pendant le transport de particules.""",
            ang = """Meant for particles: drogues, oil spills.
If no turbulence is activated, this stochastic diffusion is not
considered during the particle transport.""",
        ),
    ),
#   -----------------------------------
    LAGRANGIAN_DRIFTS = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_LAGRANGIAN_DRIFTS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Permet d''effectuer simultanement plusieurs calculs de derives
lagrangiennes initiees a des pas differents.
Ajouter A et G au mot-cle
\telkey{VARIABLES POUR LES SORTIES GRAPHIQUES}.""",
            ang = """Provided for performing several computations of Lagrangian
drifts starting at different times.
Add A and G in the \telkey{VARIABLES FOR GRAPHIC PRINTOUTS} keyword.""",
        ),
#       -----------------------------------
        b_NUMBER_OF_LAGRANGIAN_DRIFTSG = BLOC(condition="NUMBER_OF_LAGRANGIAN_DRIFS != 0",
#       -----------------------------------
#           -----------------------------------
            Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#           -----------------------------------
                defaut = "Add 'drift along x (m)' and 'drift along y (m)' in VARIABLES FOR GRAPHIC PRINTOUTS"),
        ),
    ),
)
# -----------------------------------------------------------------------
HYDRAULIC_STRUCTURES = PROC(nom= "HYDRAULIC_STRUCTURES",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    WEIRS = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_WEIRS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Nombre de seuils qui seront traites par des conditions aux
limites. Ces seuils doivent etre decrits comme des frontieres du
domaine de calcul, et leurs caracteristiques sont donnees dans le
\telkey{FICHIER DE DONNEES DES SEUILS} (voir la documentation ecrite).""",
            ang = """Number of weirs that will be treated by boundary conditions.
They must be described as boundaries of the domain and their features
are given in the \telkey{WEIRS DATA FILE} (see written documentation)""",
        ),
#       -----------------------------------
        b_NUMBER_OF_WEIRSG = BLOC(condition="NUMBER_OF_WEIRS != 0",
#       -----------------------------------
#           -----------------------------------
            WEIRS_DATA_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de description des seuils presents dans le modele.""",
                ang = """Description of weirs existing in the model.""",
            ),
#           -----------------------------------
            TYPE_OF_WEIRS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["HORIZONTAL WITH SAME NUMBER OF NODES UPSTREAM/DOWNSTREAM","GENERAL"],
                defaut = "HORIZONTAL WITH SAME NUMBER OF NODES UPSTREAM/DOWNSTREAM",
                fr = """Methode de traitement des seuils. Deux Solutions:
\begin{itemize}
\item horizontal avec meme nombre de noeuds amont/aval
(Solution historique avec le sous-programme \telfile{BORD}) ;
\item generale (nouvelle solution avec pts sources).
\end{itemize}""",
                ang = """Method for treatment of weirs. Two options:
\begin{itemize}
\item horizontal with same number of nodes upstream/downstream
(Historical solution with the \telfile{BORD} subroutine),
\item general (new solution with sources points).
\end{itemize}""",
            ),
#           -----------------------------------
            WEIRS_DISCHARGE_OUTPUT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Fichier de sortie des debits sur les seuils presents dans le modele.""",
                ang = """Output file of discharge of weirs existing in the model.""",
            ),
        ),
    ),
#   -----------------------------------
    CULVERTS = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_CULVERTS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Nombre de buses ou ponts traites comme des termes sources ou
puits. Ces buses doivent etre decrits comme des sources dans le
fichier des parametres .
Leurs caracteristiques sont donnees dans le
\telkey{FICHIER DE DONNEES DES BUSES}
(voir la documentation ecrite).""",
            ang = """Number of culverts, tubes or bridges treated as source terms.
They must be described as sources in the domain and their features
are given in the \telfile{CULVERTS DATA FILE}
(see written documentation).""",
        ),
#       -----------------------------------
        b_NUMBER_OF_CULVERTSG = BLOC(condition="NUMBER_OF_CULVERTS != 0",
#       -----------------------------------
#           -----------------------------------
            CULVERTS_DATA_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de description des buses/ponts presents dans le modele.""",
                ang = """Description of culverts/tubes/bridges existing in the model.""",
            ),
        ),
#       -----------------------------------
        OPTION_FOR_CULVERTS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Option pour le traitement des buses. Il existe deux formulations
dans \telemac{2d}, basees sur les formules de Bodhaine (1968)
et Carlier (1976).
Lire le \telemac{3d} theory guide plus plus d informations.""",
            ang = """Option for the treatment of culverts. There are two options in
\telemac{2d} based on Bodhaine (1968) and Carlier (1976) formulae.
Read the \telemac{3d} theory guide for more informations.""",
        ),
    ),
#   -----------------------------------
    BREACHES = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        BREACH = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Prise en compte de breches dans le calcul par
modification altimetrique dans le maillage. La description
des breches se fait avec le \telkey{FICHIER DE DONNEES DES BRECHES}.""",
            ang = """Take in account some breaches during the computation
by modifying the bottom level of the mesh. Breach description
is done with the \telkey{BREACHES DATA FILE}.""",
        ),
#       -----------------------------------
        b_BREACHG = BLOC(condition="BREACH == True",
#       -----------------------------------
#           -----------------------------------
            BREACHES_DATA_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de description des breches.""",
                ang = """Description of breaches.""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
TIDES = PROC(nom= "TIDES",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    BINARY_DATABASE_1_FOR_TIDE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Base de donnees binaire 1 de constantes harmoniques.
 Dans le cas des donnees satellitaires de TPXO, ce fichier correspond
 aux donnees de niveau d''eau, par exemple h\_tpxo7.2""",
        ang = """Binary database 1 of harmonic constants.
 In the case of the TPXO satellite altimetry model, this file should
 be for free surface level, for instance h\_tpxo7.2""",
    ),
#   -----------------------------------
    BINARY_DATABASE_2_FOR_TIDE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Base de donnees binaire 2 de constantes harmoniques.
 Dans le cas des donnees satellitaires de TPXO, ce fichier correspond
 aux donnees de vitesse de marrees, par exemple u\_tpxo7.2""",
        ang = """Binary database 2 of harmonic constants.
 In the case of the TPXO satellite altimetry model, this file should
 be for tidal velocities, for instance u\_tpxo7.2""",
    ),
#   -----------------------------------
    GEOGRAPHIC_SYSTEM = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        into = ["NO DEFAULT VALUE","DEFINED BY USER","WGS84 LONGITUDE/LATITUDE IN REAL DEGREES","WGS84 NORTHERN UTM","WGS84 SOUTHERN UTM","LAMBERT","MERCATOR FOR TELEMAC"],
        defaut = "NO DEFAULT VALUE",
        fr = """Systeme de coordonnees geographiques dans lequel est construit
le modele numerique.
Indiquer la zone correspondante avec le mot-cle.
Indique le systeme de coordonnees geographiques dans lequel est
construit le modele numerique. Les choix possibles sont :
\begin{itemize}
\item 0 : defini par l''utilisateur ;
\item 1 : WGS84 longitude/latitude en degres reels ;
\item 2 : WGS84 nord UTM ;
\item 3 : WGS84 sud UTM ;
\item 4 : Lambert ;
\item 5 : projection Mercator.
\end{itemize}""",
        ang = """Geographic coordinates system in which the numerical model is
built. Indicate the corresponding zone with the keyword.
The possible choices are:
\begin{itemize}
\item 0: defined by the user,
\item 1: WGS84 longitude/latitude in real degrees,
\item 2: WGS84 Northern UTM,
\item 3: WGS84 Southern UTM,
\item 4: Lambert,
\item 5: Mercator projection.
\end{itemize}""",
    ),
#   -----------------------------------
    b_GEOGRAPHIC_SYSTEMG = BLOC(condition="GEOGRAPHIC_SYSTEM in ['WGS84 NOTHERN UTM','WGS84 SOUTHERN UTM','LAMBERT']",
#   -----------------------------------
#       -----------------------------------
        ZONE_NUMBER_IN_GEOGRAPHIC_SYSTEM = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = -1,
            fr = """Numero de zone (fuseau ou type de projection) lors de
l''utilisation d''une projection plane. Indiquer le systeme
geographique dans lequel est construit le modele numerique avec le
mot-cle \telkey{SYSTEME GEOGRAPHIQUE}.
Les choix possibles sont :
\begin{itemize}
\item 1 : Lambert 1 nord ;
\item 2 : Lambert 2 centre ;
\item 3 : Lambert 3 sud ;
\item 4 : Lambert 4 Corse ;
\item 22 : Lambert 2 etendu ;
\item 93 : Lambert 93 ;
\item X : Valeur UTM de la zone WGS84 (X est le numero de la zone).
\end{itemize}""",
            ang = """Number of zone when using a plane projection.
Indicate the geographic system in which the numerical model is built
with the keyword \telkey{GEOGRAPHIC SYSTEM}.
Possible choices are:
\begin{itemize}
\item 1: Lambert 1 north,
\item 2: Lambert 2 center,
\item 3: Lambert 3 south,
\item 4: Lambert 4 Corsica,
\item 22: Lambert 22 extended,
\item 93: Lambert 93,
\item X: UTM zone with WGS84 (X is the number of the zone).
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    LAMBERT_93_CONVERSION_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Nom du fichier gr3df97a.txt, grille de conversion pour Lambert 93.""",
        ang = """Name of file gr3df97a.txt, conversion grid for Lambert 93.""",
    ),
#   -----------------------------------
    COEFFICIENT_TO_CALIBRATE_SEA_LEVEL = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Coefficient pour ajuster le niveau de mer.
Il peut dependre de la reference altimetrique utilisee dans le modele,
par exemple Cartes Marines, Niveau Moyen de la Mer, NGF...""",
        ang = """Coefficient to calibrate the sea level.
It may depend on the altimetric reference used in the model,
for example Chart Datum, Mean Sea Level...""",
    ),
#   -----------------------------------
    GLOBAL_NUMBER_OF_THE_POINT_TO_CALIBRATE_HIGH_WATER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 0,
        fr = """Numero global du point
(entre 1 et le nombre de points frontieres du maillage 2D)
par rapport auquel les ondes de maree sont dephasees
pour debuter le calcul par une pleine mer
(en marees schematiques seulement).
Ce point doit etre un point de frontiere maritime.
Ne concerne que les bases de constantes harmoniques de type TPXO.""",
        ang = """Global number of the point
(between 1 and the number of boundary nodes in the 2D mesh)
with respect to which the tidal constituents have their phase shifted
to start the calculation with a high water
(for schematic tides only).
This point has to be a maritime boundary node.
Only harmonic constants databases like TPXO are concerned.""",
    ),
#   -----------------------------------
    MINOR_CONSTITUENTS_INFERENCE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Pour la base de donnees TPXO uniquement.
Interpolation de composantes harmoniques mineures
a partir de celles lues dans les fichiers d''entree
lies aux mots-cles \telkey{BASE BINAIRE 1 DE DONNEES DE MAREE}
et \telkey{BASE BINAIRE 2 DE DONNEES DE MAREE}.""",
        ang = """For TPXO tidal data base only.
Inference of minor constituents from the one read in input files
linked to keywords \telkey{BINARY DATABASE 1 FOR TIDE}
and \telkey{BINARY DATABASE 2 FOR TIDE}.""",
    ),
#   -----------------------------------
    INITIAL_VELOCITIES_COMPUTED_BY_TPXO = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Composantes de vitesses initiales calculees a partir des
solutions de maree de l OSU (ex TPXO).
Prendre NON pour eviter une initialisation avec des vitesses de maree
trop grandes.
Pour les solutions de marees provenant de l OSU uniquement (ex TPXO).""",
        ang = """Initial velocity components computed from a tidal solution from
OSU (e.g. TPXO).
NO to prevent from an initialisation with too big tidal velocities.
For tidal solutions coming from OSU only (e.g. TPXO).""",
    ),
#   -----------------------------------
    MINIMUM_DEPTH_TO_COMPUTE_TIDAL_VELOCITIES_INITIAL_CONDITIONS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.1,
        fr = """Valeur minimale de hauteur d eau au dessus de laquelle les conditions
initiales de courants de maree sont calculees.
Autrement, les vitesses sont annulees.
Pour les solutions de maree provenant de l OSU uniquement (ex TPXO).""",
        ang = """Minimum value of water depth above which initial conditions for tidal
velocities are computed.
Otherwise, the velocity components are equal to 0.
For tidal solutions coming from OSU only (e.g. TPXO).""",
    ),
#   -----------------------------------
    MINIMUM_DEPTH_TO_COMPUTE_TIDAL_VELOCITIES_BOUNDARY_CONDITIONS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.1,
        fr = """Valeur minimale de hauteur d eau utilisee pour le calcul des conditions
aux limites de courants de maree si les profondeurs sont trop petites.
Pour les solutions de marees provenant de l OSU uniquement (ex TPXO).""",
        ang = """Minimum value of water depth used to compute tidal boundary conditions
for velocities if the water depths are too small.
For tidal solutions coming from OSU only (e.g. TPXO).""",
    ),
#   -----------------------------------
    BOUNDARY_CONDITIONS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        OPTION_FOR_TIDAL_BOUNDARY_CONDITIONS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', max='**',
            into = ["No tide","Real tide (recommended methodology)","Astronomical tide","Mean spring tide","Mean tide","Mean neap tide","Astronomical neap tide","Real tide (methodology before 2010)"],
            fr = """Option pour les conditions aux limites de maree. Pour des marees
reelles, l option 1 est recommandee. Depuis la version 7.1, ce mot-cle
est un tableau avec une valeur donnee par frontiere liquide, separee par
point-virgules. Ceci permet d''avoir des conditions de maree (ou pas)
calculees sur des frontieres liquides avec vitesses ou hauteur d eau
imposees. Ca evite un conflit lors de l utilisation de seuils dans le
domaine. 0 est le code pour des conditions autres que des conditions de
maree. ATTENTION depuis la version 7.1 ! Les anciens modeles doivent
etre changes si la frontiere de maree n a pas le numero 1. Dans ce cas,
le mot-cle doit etre change et plus de valeurs doivent etre donnees.
Calage possible par les mots-cles
\telkey{COEFFICIENT DE CALAGE DU MARNAGE},
\telkey{COEFFICIENT DE CALAGE DES VITESSES DE COURANT} et
\telkey{COEFFICIENT DE CALAGE DU NIVEAU DE MER}.""",
            ang = """Option for tidal boundary conditions. For real tides, option 1 is
recommended. This keyword has been an array with a value given per
liquid boundary, separated by semicolons, since release 7.1. This
enables to have tidal conditions (or not) computed on liquid boundaries
with prescribed velocities or depths, avoiding a clash when using weirs
in the domain. 0 codes for conditions other than tidal. BEWARE since
release 7.1! Old models must be changed if their tidal boundary is not
number 1. In that case this keyword must be changed and more values
given. Possible calibration with the keywords
\telkey{COEFFICIENT TO CALIBRATE TIDAL RANGE},
\telkey{COEFFICIENT TO CALIBRATE TIDAL VELOCITIES} and
\telkey{COEFFICIENT TO CALIBRATE SEA LEVEL}.""",
        ),
#       -----------------------------------
        TIDAL_DATA_BASE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["NO DEFAULT VALUE","JMJ","TPXO","MISCELLANEOUS (LEGOS-NEA, FES20XX, PREVIMER...)"],
            defaut = "NO DEFAULT VALUE",
            fr = """Pour JMJ, renseigner la localisation du fichier bdd\_jmj et geofin
dans les mots-cles \telkey{BASE ASCII DE DONNEES DE MAREE} et
\telkey{FICHIER DU MODELE DE MAREE}.
Pour TPXO, LEGOS-NEA, FES20XX et PREVIMER, l''utilisateur doit
telecharger les fichiers de constantes harmoniques sur internet.""",
            ang = """For JMJ, indicate the location of the files bdd\_jmj and geofin with
keywords \telkey{ASCII DATABASE FOR TIDE} and \telkey{TIDAL MODEL FILE}.
For TPXO, LEGOS-NEA,
FES20XX and PREVIMER, the user has to download files of harmonic
constituents on the internet.""",
        ),
#       -----------------------------------
        b_TIDAL_DATA_BASEG = BLOC(condition="TIDAL_DATA_BASE == 'TPXO'",
#       -----------------------------------
        ),
#       -----------------------------------
        HARMONIC_CONSTANTS_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier contenant les constantes harmoniques pour le calcul des
conditions aux limites de maree.""",
            ang = """File containing the harmonic constants to compute the
tidal boundary conditions.""",
        ),
#       -----------------------------------
        TIDAL_MODEL_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier de geometrie du modele dont sont extraites
les constantes harmoniques.""",
            ang = """Geometry file of the model from which harmonic constituents
are extracted.""",
        ),
#       -----------------------------------
        TIDAL_MODEL_FILE_FORMAT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ['SERAFIN','SERAFIND','MED'],
            defaut = 'SERAFIN',
            fr = """Format du \telkey{FICHIER DU MODELE DE MAREE}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
            ang = """Format of the \telkey{TIDAL MODEL FILE}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
        ),
#       -----------------------------------
        ASCII_DATABASE_FOR_TIDE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Base de donnees de constantes harmoniques
tirees du \telkey{FICHIER DU MODELE DE MAREE}.
Ancien nom en version 6.1 : \telkey{BASE DE DONNEES DE MAREE}.""",
            ang = """Tide data base of harmonic constituents
extracted from the \telkey{TIDAL MODEL FILE}.
Old name in release 6.1: \telkey{TIDE DATA BASE}.""",
        ),
#       -----------------------------------
        COEFFICIENT_TO_CALIBRATE_TIDAL_RANGE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """Coefficient pour ajuster le marnage de l''onde de maree
aux frontieres maritimes.""",
            ang = """Coefficient to calibrate the tidal range of tidal wave
at tidal open boundary conditions.""",
        ),
#       -----------------------------------
        COEFFICIENT_TO_CALIBRATE_TIDAL_VELOCITIES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 999999.,
            fr = """Coefficient pour ajuster les composantes de vitesse
de l''onde de maree aux frontieres maritimes.
La valeur par defaut 999~999. signifie que c''est la racine carree
du \telkey{COEFFICIENT DE CALAGE DU MARNAGE} qui est prise.""",
            ang = """Coefficient to calibrate the tidal velocities of tidal wave
at tidal open boundary conditions.
Default value 999,999. means that the square root of
\telkey{COEFFICIENT TO CALIBRATE TIDAL RANGE} is taken.""",
        ),
#       -----------------------------------
        LOCAL_NUMBER_OF_THE_POINT_TO_CALIBRATE_HIGH_WATER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Numero local du point entre 1 et le nombre de points de frontiere
maritime (du \telkey{FICHIER DES CONSTANTES HARMONIQUES}) ou les
conditions aux limites de maree sont calculees avec les bases de donnees
JMJ, NEA, FES, PREVIMER (sauf les bases de type TPXO).
Les ondes de maree sont
dephasees par rapport a ce point pour debuter le calcul par une pleine
mer (en marees schematiques seulement).""",
            ang = """Local number between 1 and the number of tidal boundary points (of the
\telkey{HARMONIC CONSTANTS FILE}) where the tidal boundary conditions
are computed with JMJ, NEA, FES, PREVIMER databases (except TPXO-type
databases). The tidal constituents have their phase shifted with respect
to this point to start the simulation with a high water (for schematic
tides only).""",
        ),
    ),
#   -----------------------------------
    PHYSICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TIDE_GENERATING_FORCE = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Active la prise en compte de la force generatrice de la maree.
Il faut utiliser \telkey{COORDONNEES SPHERIQUES} = OUI,
il est impossible de prendre en compte la force generatrice de la maree
en coordonnees cartesiennes.""",
            ang = """The tide generating force is taken into account.
The keyword \telkey{SPHERICAL COORDINATES} has to be activated,
it is impossible to account tide generating force in cartesian
coordinates.""",
        ),
#       -----------------------------------
        b_TIDE_GENERATING_FORCEG = BLOC(condition="TIDE_GENERATING_FORCE == True",
#       -----------------------------------
        ),
    ),
)
# -----------------------------------------------------------------------
COUPLING = PROC(nom= "COUPLING",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    COUPLING_WITH = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM', min=0, max='**',
        into = ['SISYPHE','TOMAWAC','TOMAWAC2','WAQTEL','KHIONE','GAIA','DELWAQ'],
        defaut = '',
        fr = """Liste des codes avec lesquels on couple \telemac{2d} :
\begin{itemize}
\item \sisyphe : couplage interne avec \sisyphe ;
\item \tomawac : couplage interne avec \tomawac ;
\item \waqtel : couplage interne avec \waqtel ;
\item \khione : couplage interne avec \khione ;
\item DELWAQ : sortie de fichiers de resultats pour Delwaq ;
\item \gaia : couplage interne avec \gaia.
\end{itemize}""",
        ang = """List of codes to be coupled with \telemac{2d}:
\begin{itemize}
\item \sisyphe: internal coupling with \sisyphe,
\item \tomawac: internal coupling with \tomawac,
\item \waqtel: internal coupling with \waqtel,
\item \khione: internal coupling with \khione,
\item DELWAQ: will yield results file for DELWAQ,
\item \gaia: internal coupling with \gaia.
\end{itemize}""",
    ),
#   -----------------------------------
    NAMES_OF_CLANDESTINE_VARIABLES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM', min= 2, max= 2,
        fr = """Noms de variables qui ne sont pas utilisees par \telemac{2d},
mais qui doivent etre conservees lors de son execution.
Ceci peut etre utilise entre autres lors du couplage de \telemac{2d}
avec un autre code.
Les variables clandestines sont alors des variables propres a l''autre
code et sont rendues dans le fichier de resultats.""",
        ang = """Names of variables that are not used by \telemac{2d}, but should be
preserved when it is being run. This keyword may be used, for instance
when \telemac{2d} is coupled with another code. Thus, the clandestine
variables belong to the other code and are given back in the results
file.""",
    ),
#   -----------------------------------
    DELWAQ = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DELWAQ_STEERING_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de commande pour le chainage avec DELWAQ.""",
            ang = """Steering file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        DELWAQ_PRINTOUT_PERIOD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Periode de sortie des resultats pour DELWAQ.""",
            ang = """Printout period for DELWAQ files.""",
        ),
#       -----------------------------------
        EXCHANGES_BETWEEN_NODES_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        NODES_DISTANCES_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        BOTTOM_SURFACES_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        VOLUMES_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        EXCHANGE_AREAS_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        VERTICAL_FLUXES_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        VELOCITY_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        DIFFUSIVITY_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        TEMPERATURE_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        SALINITY_DELWAQ_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'),
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        VELOCITY_FOR_DELWAQ = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie de la vitesse pour DELWAQ.""",
            ang = """Triggers the output of velocity for DELWAQ.""",
        ),
#       -----------------------------------
        DIFFUSIVITY_FOR_DELWAQ = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie du coefficient de diffusion pour DELWAQ.""",
            ang = """Triggers the output of diffusion for DELWAQ.""",
        ),
#       -----------------------------------
        TEMPERATURE_FOR_DELWAQ = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie de la temperature pour DELWAQ.""",
            ang = """Triggers the output of temperature for DELWAQ.""",
        ),
#       -----------------------------------
        SALINITY_FOR_DELWAQ = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie de la salinite pour DELWAQ.""",
            ang = """Triggers the output of salinity for DELWAQ.""",
        ),
    ),
#   -----------------------------------
    SISYPHE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SISYPHE_STEERING_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres de \sisyphe en cas de couplage
interne.""",
            ang = """\sisyphe parameter file in case of internal coupling.""",
        ),
#       -----------------------------------
        COUPLING_PERIOD_FOR_SISYPHE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Fixe la periode de couplage avec le module \sisyphe,
en nombre de pas de temps.
Par defaut, il est couple a chaque pas de temps.""",
            ang = """Sets the coupling period with the \sisyphe module, in number
of time steps. By default, it is coupled at every time step.""",
        ),
    ),
#   -----------------------------------
    GAIA = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        GAIA_STEERING_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres de \gaia en cas de couplage
interne.""",
            ang = """\gaia parameter file in case of internal coupling.""",
        ),
    ),
#   -----------------------------------
    TOMAWAC = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TOMAWAC_STEERING_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres de \tomawac en cas de couplage
interne.""",
            ang = """\tomawac parameter file in case of internal coupling.""",
        ),
#       -----------------------------------
        COUPLING_PERIOD_FOR_TOMAWAC = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Fixe la periode de couplage avec le module \tomawac,
en nombre de pas de temps.
Par defaut, il est couple a chaque pas de temps.""",
            ang = """Sets the coupling period with the \tomawac module, in number
of time steps. By default, it is coupled at every time step.""",
        ),
    ),
#   -----------------------------------
    WAQTEL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WAQTEL_STEERING_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres physiques pour les processus de qualite d eau
(internes pas ceux de DELWAQ).""",
            ang = """File for physical parameters of water quality processes
(local ones of \telemac{2d}-\waqtel not those of DELWAQ).""",
        ),
    ),
#   -----------------------------------
    KHIONE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        KHIONE_STEERING_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres physiques pour les processus lies aux glaces.""",
            ang = """Steering file for physical parameters of ice processes.""",
        ),
#       -----------------------------------
        ICE_PROCESSES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Donne le numero du processus des glaces, avec le numero defini sur la
base multiplicative de nombres premiers (2, 3, 5, 7, 11, 13\ldots).
Par exemple, 14 (= 2 $\times$ 7) rend actif les processus 2 et 7.
Exception a la regle :
\begin{itemize}
\item si l entier est 1, il n y a pas de couplage ;
\item si l entier est 0, tous les effets suivants sont inclus :
\begin{itemize}
\item 2 : bilan thermique ;
\item 3 : \ldots
\end{itemize}
\end{itemize}""",
            ang = """Provides the ice process number with the number being defined on the
basis of a multiplication of primary numbers (2, 3, 5, 7, 11, 13\ldots).
For instance, 14 (= 2 $\times$ 7) activates processes 2 and 7.
Exception is for:
\begin{itemize}
\item if the number is 1, there will be no ice processes included,
\item if the number is 0, all processes are included, as follows:
\begin{itemize}
\item 2: thermical budget,
\item 3: \ldots
\end{itemize}
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    NESTOR = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        NESTOR = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Active l utilisation du module \nestor pour changer les fonds.""",
            ang = """Activates the use of the \nestor module to change the bottom.""",
        ),
#       -----------------------------------
        NESTOR_ACTION_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom du fichier des parametres de \nestor.""",
            ang = """Name of the \nestor steering file.""",
        ),
#       -----------------------------------
        NESTOR_POLYGON_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom du fichier de polygones de \nestor indiquant la localisation.""",
            ang = """Name of the \nestor polygon file which indicates the location.""",
        ),
#       -----------------------------------
        NESTOR_SURFACE_REFERENCE_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'), max='**',
            defaut = '',
            fr = """Nom du fichier de \nestor qui contient la surface de reference.""",
            ang = """Name of the \nestor file which contains the reference water surface.""",
        ),
#       -----------------------------------
        NESTOR_RESTART_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom du fichier de reprise de \nestor.""",
            ang = """Name of the \nestor restart file.""",
        ),
    ),
)
# -----------------------------------------------------------------------
INTERNAL = PROC(nom= "INTERNAL",op = None,
# -----------------------------------------------------------------------
    UIinfo = {"groupes": ("CACHE")},
#   -----------------------------------
    LANGUAGE = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["FRANCAIS","ANGLAIS"],
        defaut = "ANGLAIS",
        fr = """1 : FRANCAIS   2 : ANGLAIS""",
        ang = """1: FRENCH   2: ENGLISH""",
    ),
#   -----------------------------------
    STEERING_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Nom du fichier contenant les parametres du calcul a realiser.""",
        ang = """Name of the file containing the parameters of the computation
Written by the user.""",
    ),
#   -----------------------------------
    DICTIONARY = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = 'telemac2d.dico',
        fr = """Dictionnaire des mots cles.""",
        ang = """Key word dictionary.""",
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
    CONCATENATE_PARTEL_OUTPUT = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Avec cette option partel ne generera non plus un fichier (GEO/CLI/PAR)
par processeur mais une concatenation de ceux-ci, ainsi qu''un fichier
d''index associe. Ainsi plutot que d''avoir 3P fichiers, il n''y en a
plus que 6.""",
        ang = """With this option partel no more generates a file (GEO/CLI/PAR) per
process but a single concatenate file of them, associated to an index
file. Then instead of having partel generating 3P files, it only
generates 6 files.""",
    ),
)
TEXTE_NEW_JDC = "\
COMPUTATION_ENVIRONMENT();\
HYDRO();\
NUMERICAL_PARAMETERS();\
GENERAL_PARAMETERS();\
TURBULENCE();\
"
Ordre_Des_Commandes = (
'COMPUTATION_ENVIRONMENT',
'HYDRO',
'NUMERICAL_PARAMETERS',
'GENERAL_PARAMETERS',
'TURBULENCE',
'TIDAL_FLATS_INFO',
'TRACERS',
'PARTICLE_TRANSPORT',
'HYDRAULIC_STRUCTURES',
'TIDES',
'COUPLING',
'INTERNAL')
try:
    import TelApy
    source = "eficas"
except Exception as excpt:
    source = "Telemac"
enum = source+'.telemac2d_enum_auto'
dicoCasEn = source+'.telemac2d_dicoCasEnToCata'
dicoCasFr = source+'.telemac2d_dicoCasFrToCata'
