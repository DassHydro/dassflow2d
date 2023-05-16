
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



JdC = JDC_CATA (code = 'TELEMAC3D',
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
        TITLE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Titre du cas etudie.""",
            ang = """Title of the case being considered.""",
        ),
#       -----------------------------------
        PARALLEL_PROCESSORS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Nombre de processeurs pour la decomposition en parallele.
La valeur 0 correspond a un calcul scalaire.""",
            ang = """Number of processors for domain partition.
Value 0 corresponds to a scalar computation.""",
        ),
#       -----------------------------------
        CHECKING_THE_MESH = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si OUI on appelle le sous-programme \telfile{CHECKMESH} qui verifie
la coherence du maillage, points superposes, etc.""",
            ang = """If this keyword is equal to YES, a call to subroutine
\telfile{CHECKMESH} will look for errors in the mesh,
superimposed points, etc.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 30,
            fr = """Nombre maximal de frontieres differentes dans le maillage.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
            ang = """Maximal number of boundaries in the mesh.
Used for dimensioning arrays. Can be increased if needed.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 20,
            fr = """Nombre maximal de traceurs.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
            ang = """Maximal number of tracers.
Used for dimensioning arrays. Can be increased if needed.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_SOURCES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 20,
            fr = """Nombre maximal de points sources dans le maillage,
incluant les sources ponctuelles et 2 fois le nombre de buses.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
            ang = """Maximal number of source points in the mesh, including
punctual sources and twice the number of culverts.
Used for dimensioning arrays. Can be increased if needed.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_BOUNDARIES_ON_THE_BED = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 30,
            fr = """Nombre maximal de frontieres liquides sur le fond.
Sert au dimensionnement de la memoire, a augmenter si necessaire.""",
            ang = """Maximal number of liquid boundaries on the bed.
Used for dimensioning arrays. Can be increased if needed.""",
        ),
#       -----------------------------------
        VECTOR_LENGTH = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Longueur du vecteur pour les machines vectorielles.""",
            ang = """Vector length on vector machines.""",
        ),
    ),
#   -----------------------------------
    INPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DATA = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            GEOMETRY_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                fr = """Nom du fichier contenant le maillage du calcul a realiser.""",
                ang = """Name of the file containing the mesh. This file may also
contain the topography and the friction coefficients.""",
            ),
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
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
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
colours that are assigned to the boundary nodes.""",
            ),
#           -----------------------------------
            BINARY_BOUNDARY_DATA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees code en binaire contenant les informations
de conditions aux limites variables en temps et en espace
provenant de jeux de donnees externes par exemple.""",
                ang = """Binary-coded data file containing the boundary conditions data
varying in time and space.""",
            ),
#           -----------------------------------
            BINARY_BOUNDARY_DATA_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER BINAIRE DE DONNEES DE FRONTIERE}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{BINARY BOUNDARY DATA FILE}.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            FORTRAN_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'FichierOuRepertoire',
                defaut = '',
                fr = """Nom du fichier ou repertoire FORTRAN a soumettre,
contenant les sous-programmes specifiques au modele.""",
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
Si ce mot-cle est utilise, c''est cette bathymetrie qui sera utilisee
pour le calcul.""",
                ang = """Name of the possible file containing the bathymetric data.
Where this keyword is used, these bathymetric data shall be used in
the computation.""",
            ),
#           -----------------------------------
            NUMBER_OF_BOTTOM_SMOOTHINGS = SIMP(statut ='o',
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
To be used when interpolation of bathymetry on the mesh gives
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
            FORMATTED_DATA_FILE_1 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees formate mis a la disposition de
l''utilisateur.""",
                ang = """Formatted data file available to the user.""",
            ),
#           -----------------------------------
            FORMATTED_DATA_FILE_2 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees formate mis a la disposition de
l''utilisateur.""",
                ang = """Formatted data file available to the user.""",
            ),
#           -----------------------------------
            BINARY_DATA_FILE_1 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees code en binaire mis a la disposition
de l''utilisateur.""",
                ang = """Data file in binary mode available to the user.""",
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
                ang = """Format of the \telkey{BINARY DATA FILE 1}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            BINARY_DATA_FILE_2 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees code en binaire mis a la disposition
de l''utilisateur.""",
                ang = """Data file in binary mode available to the user.""",
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
            b_VALIDATIONG = BLOC(condition="VALIDATION == True",
#           -----------------------------------
#               -----------------------------------
                REFERENCE_FILE = SIMP(statut ='f',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)'), max='**',
                    defaut = '',
                    fr = """Fichier de resultats de reference pour la validation.""",
                    ang = """Binary-coded result file for validation.""",
                ),
#               -----------------------------------
                REFERENCE_FILE_FORMAT = SIMP(statut ='o',
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
                    ang = """Format of the \telkey{REFERENCE FILE}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
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
            RD_RESULT_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les resultats 3D du
calcul avec la periodicite donnee par le mot cle \telkey{PERIODE POUR
LES SORTIES GRAPHIQUES}.""",
                ang = """Name of the file into which the 3D results of the computation
are written with a periodicity given by the keyword
\telkey{GRAPHIC PRINTOUT PERIOD}.""",
            ),
#           -----------------------------------
            RD_RESULT_FILE_FORMAT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DES RESULTATS 3D}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{3D RESULT FILE}. Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            ED_RESULT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les resultats 2D du
calcul avec la periodicite donnee par le mot cle \telkey{PERIODE POUR
LES SORTIES GRAPHIQUES}.""",
                ang = """Name of the file into which the 2D results of the computation
are written with a periodicity given by the keyword
\telkey{GRAPHIC PRINTOUT PERIOD}.""",
            ),
#           -----------------------------------
            ED_RESULT_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER DES RESULTATS 2D}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{2D RESULT FILE}. Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            RESULT_FILE_IN_LONGITUDE_LATITUDE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Donne les coordonnees dans le fichier resultats en longitude-latitude
si le fichier de geometrie est aussi donne en longitude-latitude.""",
                ang = """Gives the coordinates of the result file in longitude-latitude
if the geometry file is also given in longitude-latitude.""",
            ),
#           -----------------------------------
            VARIABLES_FOR_3D_GRAPHIC_PRINTOUTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["velocity along x axis (m/s)","velocity along y axis (m/s)","velocity along z axis (m/s)","elevation z (m)","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracers from 1 to 9","concentrations for tracers from 10 to 99","viscosity for U and V along x axis (m2/s)","viscosity for U and V along y axis (m2/s)","viscosity for U and V along z axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along x axis (m2/s)","viscosity for tracer1 along x axis (m2/s)","viscosity for tracer2 along x axis (m2/s)","viscosity for tracer3 along x axis (m2/s)","viscosity for tracer4 along x axis (m2/s)","viscosity for tracer5 along x axis (m2/s)","viscosity for tracer6 along x axis (m2/s)","viscosity for tracer7 along x axis (m2/s)","viscosity for tracer8 along x axis (m2/s)","viscosity for tracer9 along x axis (m2/s)","viscosity for tracer10 along x axis (m2/s)","viscosity for tracer11 along x axis (m2/s)","viscosity for tracer12 along x axis (m2/s)","viscosity for tracer13 along x axis (m2/s)","viscosity for tracer14 along x axis (m2/s)","viscosity for tracer15 along x axis (m2/s)","viscosity for tracer16 along x axis (m2/s)","viscosity for tracer17 along x axis (m2/s)","viscosity for tracer18 along x axis (m2/s)","viscosity for tracer19 along x axis (m2/s)","viscosity for tracer** along x axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along y axis (m2/s)","viscosity for tracer1 along y axis (m2/s)","viscosity for tracer2 along y axis (m2/s)","viscosity for tracer3 along y axis (m2/s)","viscosity for tracer4 along y axis (m2/s)","viscosity for tracer5 along y axis (m2/s)","viscosity for tracer6 along y axis (m2/s)","viscosity for tracer7 along y axis (m2/s)","viscosity for tracer8 along y axis (m2/s)","viscosity for tracer9 along y axis (m2/s)","viscosity for tracer10 along y axis (m2/s)","viscosity for tracer11 along y axis (m2/s)","viscosity for tracer12 along y axis (m2/s)","viscosity for tracer13 along y axis (m2/s)","viscosity for tracer14 along y axis (m2/s)","viscosity for tracer15 along y axis (m2/s)","viscosity for tracer16 along y axis (m2/s)","viscosity for tracer17 along y axis (m2/s)","viscosity for tracer18 along y axis (m2/s)","viscosity for tracer19 along y axis (m2/s)","viscosity for tracer** along y axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","viscosity for tracer* along z axis (m2/s)","viscosity for tracer1 along z axis (m2/s)","viscosity for tracer2 along z axis (m2/s)","viscosity for tracer3 along z axis (m2/s)","viscosity for tracer4 along z axis (m2/s)","viscosity for tracer5 along z axis (m2/s)","viscosity for tracer6 along z axis (m2/s)","viscosity for tracer7 along z axis (m2/s)","viscosity for tracer8 along z axis (m2/s)","viscosity for tracer9 along z axis (m2/s)","viscosity for tracer10 along z axis (m2/s)","viscosity for tracer11 along z axis (m2/s)","viscosity for tracer12 along z axis (m2/s)","viscosity for tracer13 along z axis (m2/s)","viscosity for tracer14 along z axis (m2/s)","viscosity for tracer15 along z axis (m2/s)","viscosity for tracer16 along z axis (m2/s)","viscosity for tracer17 along z axis (m2/s)","viscosity for tracer18 along z axis (m2/s)","viscosity for tracer19 along z axis (m2/s)","viscosity for tracer** along z axis (m2/s)","Richardson number in case of mixing length model","turbulent energie for k-epsilon model (J/kg)","dissipation of turbulent energie (W/kg)","dynamic pressure (multiplied by DT/RHO)","hydrostatic pressure (in Pascals)","relative density","wall distance","private variable 1","private variable 2","private variable 3","private variable 4","Stokes velocity along x axis (m/s)","Stokes velocity along y axis (m/s)","Stokes velocity along z axis (m/s)"],
                defaut = [],
                fr = """Noms des variables que l''utilisateur veut ecrire dans
le \telkey{FICHIER DES RESULTATS 3D}.
Le choix des separateurs est libre.
Les possibilites offertes sont les suivantes :
\begin{itemize}
\item U   : vitesse suivant l''axe des $x$ (m/s) ;
\item V   : vitesse suivant l''axe des $y$ (m/s) ;
\item W   : vitesse suivant l''axe des $z$ (m/s) ;
\item Z   : cote $z$ (m) ;
\item TAx : concentrations des traceurs ;
\item NUX : viscosite pour $U$ et $V$ suivant l''axe des $x$ (m$^2$/s) ;
\item NUY : viscosite pour $U$ et $V$ suivant l''axe des $y$ (m$^2$/s) ;
\item NUZ : viscosite pour $U$ et $V$ suivant l''axe des $z$ (m$^2$/s) ;
\item NAX : viscosites pour les traceurs suivant l''axe des $x$
(m$^2$/s) ;
\item NAY : viscosites pour les traceurs suivant l''axe des $y$
(m$^2$/s) ;
\item NAZ : viscosites pour les traceurs suivant l''axe des $z$
(m$^2$/s) ;
\item RI  : nombre de Richardson en cas de modele de longueur de
melange ;
\item K   : energie turbulente du modele k-epsilon (J/kg) ;
\item EPS : dissipation de l''energie turbulente (W/kg) ;
\item DP  : pression dynamique (multipliee par DT/RHO) ;
\item PH  : pression hydrostatique (en Pascals) ;
\item RHO : densite relative ;
\item P1  : variable privee 1 ;
\item P2  : variable privee 2 ;
\item P3  : variable privee 3 ;
\item P4  : variable privee 4 ;
\item US  : vitesse de Stokes suivant l''axe des $x$ (m/s) ;
\item VS  : vitesse de Stokes suivant l''axe des $y$ (m/s) ;
\item WS  : vitesse de Stokes suivant l''axe des $z$ (m/s).
\end{itemize}""",
                ang = """Names of variables to be written in the
\telkey{3D RESULT FILE}. Free choice of separator. You can ask for:
\begin{itemize}
\item U  : velocity along $x$ (m/s),
\item V  : velocity along $y$ (m/s),
\item W  : velocity along $z$ (m/s),
\item Z  : elevation $z$ (m),
\item TAx: concentration of tracers,
\item NUX: viscosity for $U$ and $V$ along $x$ (m$^2$/s),
\item NUY: viscosity for $U$ and $V$ along $y$ (m$^2$/s),
\item NUZ: viscosity for $U$ and $V$ along $z$ (m$^2$/s),
\item NAX: viscosity for tracers along $x$ (m$^2$/s),
\item NAY: viscosity for tracers along $y$ (m$^2$/s),
\item NAZ: viscosity for tracers along $z$ (m$^2$/s),
\item RI : Richardson number for mixing length model,
\item K  : turbulent kinetic energy for $k$-$\epsilon$ model (J/kg),
\item EPS: dissipation of turbulent kinetic energy (W/kg),
\item DP : dynamic pressure (multiplied by DT/RHO),
\item PH : hydrostatic pressure (Pa),
\item RHO: relative density,
\item P1 : private variable 1,
\item P2 : private variable 2,
\item P3 : private variable 3,
\item P4 : private variable 4,
\item US : Stokes velocity along $x$ axis (m/s),
\item VS : Stokes velocity along $y$ axis (m/s),
\item WS : Stokes velocity along $z$ axis (m/s).
\end{itemize}""",
            ),
#           -----------------------------------
            VARIABLES_FOR_2D_GRAPHIC_PRINTOUTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["depth averaged velocity along x axis (m/s)","depth averaged velocity along y axis (m/s)","celerity (m/s)","water depth (m)","free surface elevation (m)","bottom elevation (m)","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracer*","concentrations for tracer1","concentrations for tracer2","concentrations for tracer3","concentrations for tracer4","concentrations for tracer5","concentrations for tracer6","concentrations for tracer7","concentrations for tracer8","concentrations for tracer9","concentrations for tracer10","concentrations for tracer11","concentrations for tracer12","concentrations for tracer13","concentrations for tracer14","concentrations for tracer15","concentrations for tracer16","concentrations for tracer17","concentrations for tracer18","concentrations for tracer19","concentrations for tracer**","concentrations for tracers from 1 to 9","concentrations for tracers from 10 to 99","Froude number","scalar discharge (m2/s)","discharge along x (m2/s)","discharge along y (m2/s)","norm of velocity (m/s)","wind along x axis (m/s)","wind along y axis (m/s)","atmospheric pressure (Pa)","friction coefficient","non erodible bottom elevation (m)","thickness of the sediment bed layer (m)","erosion rate (kg/m2/s)","deposition flux (kg/m2/s)","bed evolution","work array PRIVE 1","work array PRIVE 2","work array PRIVE 3","work array PRIVE 4","solid discharge (m2/s)","solid discharge along x (m2/s)","solid discharge along y (m2/s)","friction velocity (m/s)","maximum value of the free surface elevation (m)","time corresponding to this maximum elevation (s)","air temperature (degree C)","surface velocity along x axis (m/s)","surface velocity along y axis (m/s)","surface velocity along z axis (m/s)","magnitude of velocity at the surface (m/s)","TAi conc for tracers at the surface, i is the tracer number","TA concentrations for tracers at the surface from 1 to 9","TA concentrations for tracers at the surface from 10 to 99"],
                defaut = [],
                fr = """Noms des variables que l''utilisateur veut ecrire dans
le \telkey{FICHIER DES RESULTATS 2D}.
Chaque variable est representee par une lettre.
Le choix des separateurs est libre.
Les possibilites offertes sont les suivantes :
\begin{itemize}
\item U : vitesse moyenne suivant l''axe des x (m/s) ;
\item V : vitesse moyenne suivant l''axe des y (m/s) ;
\item C : celerite (m/s) ;
\item H : hauteur d''eau (m) ;
\item S : cote de surface libre (m) ;
\item B : cote du fond (m) ;
\item TAi : TAi concentrations des traceurs, i numero du traceur ;
\item TA* : TA concentrations des traceurs de 1 a 9 ;
\item TA** : TA concentrations des traceurs de 10 a 99 ;
\item F : nombre de Froude ;
\item Q : debit scalaire (m$^2$/s) ;
\item I : debit suivant x (m$^2$/s) ;
\item J : debit suivant y (m$^2$/s) ;
\item M : norme de la vitesse (m/s) ;
\item X : vent suivant l''axe des x (m/s) ;
\item Y : vent suivant l''axe des y (m/s) ;
\item P : pression atmospherique (Pa) ;
\item W : coefficient de frottement ;
\item RB : cote des fonds non erodables (m) ;
\item HD : epaisseur des depots frais (m) ;
\item EF : flux d''erosion (kg/m$^2$/s) ;
\item DF : probabilite de depot (kg/m$^2$/s) ;
\item DZF : evolution du lit ;
\item PRIVE1 : tableau de travail PRIVE 1 ;
\item PRIVE2 : tableau de travail PRIVE 2 ;
\item PRIVE3 : tableau de travail PRIVE 3 ;
\item PRIVE4 : tableau de travail PRIVE 4 ;
\item QS : debit solide (m$^2$/s) ;
\item QSX : debit solide suivant x (m$^2$/s) ;
\item QSY : debit solide suivant y (m$^2$/s) ;
\item US : vitesse de frottement (m/s) ;
\item MAXZ : valeur maximum de la cote de l eau au cours du calcul (m) ;
\item TMXZ : temps correspondant a ce niveau maximum (s) ;
\item TAIR : temperature de l air ($^{\circ}$C) ;
\item USURF : vitesse en surface suivant l''axe des x (m/s) ;
\item VSURF : vitesse en surface suivant l''axe des y (m/s) ;
\item WSURF : vitesse en surface suivant l''axe des z (m/s) ;
\item MSURF : norme de la vitesse en surface (m/s) ;
\item TASURFi : TAi conc des traceurs en surface, i numero du traceur ;
\item TASURF* : TA concentrations des traceurs en surface de 1 a 9 ;
\item TASURF** : TA concentrations des traceurs en surface de 10 a 99.
\end{itemize}""",
                ang = """Names of variables that may be written in the
\telkey{2D RESULT FILE}.
Every variable is represented by a group of letters with
any separator between them , ; or blank
possibilities are the following:
\begin{itemize}
\item U: depth averaged velocity along x axis (m/s),
\item V: depth averaged velocity along y axis (m/s),
\item C: celerity (m/s),
\item H: water depth (m),
\item S: free surface elevation (m),
\item B: bottom elevation (m),
\item TAi: TAi concentrations for tracers, i is the tracer number,
\item TA*: TA concentrations for tracers from 1 to 9,
\item TA**: TA concentrations for tracers from 10 to 99,
\item F: Froude number,
\item Q: scalar discharge (m$^2$/s),
\item I: discharge along x (m$^2$/s),
\item J: discharge along y (m$^2$/s),
\item M: norm of velocity (m/s),
\item X: wind along x axis (m/s),
\item Y: wind along y axis (m/s),
\item P: atmospheric pressure (Pa),
\item W: friction coefficient,
\item RB: non erodible bottom elevation (m),
\item HD: thickness of the fresh deposits (m),
\item EF: erosion rate (kg/m$^2$/s),
\item DF: probability of deposition (kg/m$^2$/s),
\item DZF : bed evolution,
\item PRIVE1: work array PRIVE 1,
\item PRIVE2: work array PRIVE 2,
\item PRIVE3: work array PRIVE 3,
\item PRIVE4: work array PRIVE 4,
\item QS: solid discharge (m$^2$/s),
\item QSX: solid discharge along x (m$^2$/s),
\item QSY: solid discharge along y (m$^2$/s),
\item US: friction velocity (m/s),
\item MAXZ: maximum value of the free surface
elevation during the computation (m),
\item TMXZ: time corresponding to this maximum elevation (s),
\item TAIR: air temperature ($^{\circ}$C),
\item USURF: velocity along x axis at the surface (m/s),
\item VSURF: velocity along y axis at the surface (m/s),
\item WSURF: velocity along z axis at the surface (m/s),
\item MSURF: magnitude of velocity at the surface (m/s),
\item TASURFi: TAi conc for tracers at the surface, i is the tracer
number,
\item TASURF*: TA conc for tracers at the surface from 1 to 9,
\item TASURF**: TA conc for tracers at the surface from 10 to 99.
\end{itemize}""",
            ),
#           -----------------------------------
            GRAPHIC_PRINTOUT_PERIOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode en nombre de pas de temps d''impression des
\telkey{VARIABLES POUR LES SORTIES GRAPHIQUES 2D ou 3D}
(voir ces mot-cles) dans le \telkey{FICHIER DES RESULTATS 2D ou 3D}.""",
                ang = """Determines, in number of time steps, the printout period for the
\telkey{VARIABLES FOR 2D (or 3D) GRAPHIC PRINTOUTS}
in the \telkey{2D or 3D RESULT FILE}.""",
            ),
#           -----------------------------------
            NUMBER_OF_FIRST_TIME_STEP_FOR_GRAPHIC_PRINTOUTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Determine le numero de pas de temps a partir duquel debute
l''ecriture des resultats dans le \telkey{FICHIER DES RESULTATS 2D}
ou \telkey{3D}''.""",
                ang = """Determines the number of time steps after which the results
are first written into the \telkey{2D} or \telkey{3D RESULT FILE}.""",
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
            NUMBER_OF_2D_PRIVATE_ARRAYS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Nombre de tableaux 2D mis a disposition de l utilisateur
dans le bloc \telfile{PRIVE2D}. Il doit etre inferieur ou egal a 4.""",
                ang = """Number of 2D arrays for own user programming
in block \telfile{PRIVE2D}. It has to be lower or equal to 4.""",
            ),
#           -----------------------------------
            NAMES_OF_2D_PRIVATE_VARIABLES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min= 4, max= 4,
                fr = """Noms des variables dans les tableaux prives 2D en 32
caracteres, 16 pour le nom 16 pour l''unite. Elles seront lues dans le
\telkey{FICHIER DE GEOMETRIE} si elles y sont.
Nombre maximum de 4 noms.""",
                ang = """Name of variables in 2D private arrays in 32 characters,
16 for the name, 16 for the unit. If present, will be read
in the \telkey{GEOMETRY FILE}. Maximum number of 4 names.""",
            ),
#           -----------------------------------
            FORMATTED_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier de resultats formate mis a la disposition de
l''utilisateur.""",
                ang = """Formatted file of results available to the user.""",
            ),
#           -----------------------------------
            BINARY_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier de resultats code en binaire mis a la disposition
de l''utilisateur.""",
                ang = """Additional binary-coded result file available to the user.""",
            ),
        ),
#       -----------------------------------
        LISTING = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            LISTING_PRINTOUT = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Sortie des resultats sur support papier ou a l ecran.
Si l''on met NON le listing ne contient que l''entete et la mention
FIN NORMALE DU PROGRAMME. Commande a eviter.""",
                ang = """Result printout on hard copy.
When NO is selected, the listing only includes the heading and the
phrase "NORMAL END OF PROGRAM".
In addition, the options \telkey{MASS-BALANCE} and
\telkey{VALIDATION} are inhibited. Not recommended for use.""",
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
the domain is done or not.
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
            INFORMATION_ABOUT_MASS_BALANCE_FOR_EACH_LISTING_PRINTOUT = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = True,
                fr = """Donne a \telkey{PERIODE POUR LES SORTIES LISTING}
une information sur le bilan de masse.""",
                ang = """Gives the information about mass-balance
at every \telkey{LISTING PRINTOUT PERIOD}.""",
            ),
#           -----------------------------------
            LISTING_PRINTOUT_PERIOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode en nombre de pas de temps d''impression des
VARIABLES A IMPRIMER. Pour la mise au point,
il faut savoir que la sortie des resultats est effectuee
systematiquement sur le listing.""",
                ang = """Determines, in number of time steps, the printout period of the
VARIABLES TO BE PRINTED.
The results are systematically printed out on the listing file.""",
            ),
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
(sous-programme \telfile{BORD3D} ou valeurs placees dans le fichier des
parametres), qui peuvent egalement etre modifiees.\\
Afin d''obtenir une suite de calcul parfaite, l''utilisateur doit
activer le \telkey{MODE SUITE} dans un calcul precedent afin de generer
le fichier a partir duquel le calcul suivant commence
(\telkey{FICHIER POUR SUITE}).""",
            ang = """Determines whether the computation under way is independent
or is following an earlier result.
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
(in the subroutine \telfile{BORD3D} or values defined
in the steering file).\\
In order to get a perfect continued computation, the user has to
activate the \telkey{RESTART MODE} in a previous computation to generate
the file from which the following computation starts
(\telkey{RESTART FILE}).""",
        ),
#       -----------------------------------
        b_COMPUTATION_CONTINUEDG = BLOC(condition="COMPUTATION_CONTINUED == True",
#       -----------------------------------
#           -----------------------------------
            PREVIOUS_COMPUTATION_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom d''un fichier contenant les resultats d''un calcul precedent
realise sur le meme maillage et dont le dernier pas de temps enregistre
va fournir les conditions initiales pour une suite de calcul.
Dans le cas d''une suite de calcul que l''on souhaite parfaite,
le \telkey{FICHIER DU CALCUL PRECEDENT} doit etre le
\telkey{FICHIER POUR SUITE} du dernier calcul, ce dernier fichier
etant alors un fichier de sortie du dernier calcul.
Le \telkey{FORMAT DU FICHIER DU CALCUL PRECEDENT} et le
\telkey{FORMAT DU FICHIER POUR SUITE} doivent alors etre mis a
 ''SERAFIND'' ou ''MED''.""",
                ang = """Name of a file containing the results of an earlier computation
which was made on the same mesh. The last recorded time step will
provide the initial conditions for the new computation.
In case of a perfect continued computation, the
\telkey{PREVIOUS COMPUTATION FILE} has to be the \telkey{RESTART FILE}
of the last computation.
This last file is then an output file of the last computation.
The \telkey{PREVIOUS COMPUTATION FILE FORMAT} and the
\telkey{RESTART FILE FORMAT} have to be set with ''SERAFIND''
or ''MED''.""",
            ),
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
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            RECORD_NUMBER_FOR_RESTART = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """En cas de suite de calcul, numero de l''enregistrement
de depart dans le fichier du calcul precedent. 0 signifie
que l''on prend le dernier enregistrement.""",
                ang = """In case of \telkey{COMPUTATION CONTINUED} = YES, record number
to start from in the \telkey{PREVIOUS COMPUTATION FILE}.
0 means that the last record is taken.""",
            ),
        ),
#       -----------------------------------
        INITIAL_TIME_SET_TO_ZERO = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Remet le temps a zero en cas de suite de calcul.""",
            ang = """Initial time set to zero in case of restart.""",
        ),
#       -----------------------------------
        RESTART_MODE = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Declenche le remplissage du
\telkey{FICHIER POUR SUITE}, qui permet une suite de calcul
parfaite, contrairement au \telkey{FICHIER DES RESULTATS 3D}.""",
            ang = """Triggers the filling of the \telkey{RESTART FILE},
which ensures a perfect restart of a computation,
unlike using the \telkey{3D RESULT FILE}.""",
        ),
#       -----------------------------------
        b_RESTART_MODEG = BLOC(condition="RESTART_MODE == True",
#       -----------------------------------
#           -----------------------------------
            RESTART_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les resultats du
dernier calcul pour obtenir une suite de calcul parfaite.
C''est donc un fichier de sortie pour le calcul en cours,
qui servira de fichier d''entree lors de la suite de calcul que l''on
souhaite parfaite (le mot-cle \telkey{FICHIER DU CALCUL PRECEDENT}
est alors utilise).
Le \telkey{FORMAT DU FICHIER POUR SUITE} et le
\telkey{FORMAT DU FICHIER DU CALCUL PRECEDENT} doivent alors etre mis a
 ''SERAFIND'' ou ''MED''.""",
                ang = """Name of the file into which the last computation results shall
be written in order to get a perfect continued computation.
It is then an output file for the current computation,
which will be used as an input file when a continued computation
is expected to be perfect (the keyword
\telkey{PREVIOUS COMPUTATION FILE} is then used).
The \telkey{RESTART FILE FORMAT} and the
\telkey{PREVIOUS COMPUTATION FILE FORMAT} have to be set with
 ''SERAFIND'' or ''MED''.""",
            ),
#           -----------------------------------
            RESTART_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIND',
                fr = """Format du \telkey{FICHIER POUR SUITE}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}
Seul les formats double precision assurent une suite parfaite.""",
                ang = """Format of the \telkey{RESTART FILE}.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}
Only double precision formats ensure a perfect restart.""",
            ),
        ),
#       -----------------------------------
        ED_CONTINUATION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Permet d''utiliser un \telkey{FICHIER DES RESULTATS 2D}
stocke dans le \telkey{FICHIER POUR SUITE 2D} comme fichier de
conditions initiales.""",
            ang = """Enables to use a \telkey{2D RESULT FILE} in
\telkey{FILE FOR 2D CONTINUATION} as initial conditions file.""",
        ),
#       -----------------------------------
        b_ED_CONTINUATIONG = BLOC(condition="ED_CONTINUATION == True",
#       -----------------------------------
#           -----------------------------------
            FILE_FOR_2D_CONTINUATION = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier utilise en cas de suite 2D.""",
                ang = """File to be used in case of 2D continuation.""",
            ),
#           -----------------------------------
            FILE_FOR_2D_CONTINUATION_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du \telkey{FICHIER POUR SUITE 2D}.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Format of the \telkey{FILE FOR 2D CONTINUATION}.
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
# -----------------------------------------------------------------------
GENERAL_PARAMETERS = PROC(nom= "GENERAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    DEBUGGER = SIMP(statut ='o',
#   -----------------------------------
        typ = 'I',
        defaut = 0,
        fr = """Pour imprimer la sequence des appels, mettre 1.""",
        ang = """If 1, additional writings will be printed in the listing,
in particular the calls of subroutines.""",
    ),
#   -----------------------------------
    TIME = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TIME_STEP = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """Definit le pas de temps en secondes.
Remarque : Pour une bonne precision, il est souhaitable de choisir
le pas de temps de telle sorte que le nombre de Courant de propagation
soit inferieur a 2, voire 3.
Ceci peut etre realisable en hydraulique fluviale, mais ne l''est
pratiquement jamais en hydraulique maritime ou l''on peut atteindre
des valeurs de 50.""",
            ang = """Specifies the time step in seconds.""",
        ),
#       -----------------------------------
        NUMBER_OF_TIME_STEPS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Definit le nombre de pas de temps effectues lors de
l''execution du code.""",
            ang = """Specifies the number of time steps performed when running
the code.""",
        ),
#       -----------------------------------
        DURATION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Duree de la simulation en secondes. Alternative au parametre
\telkey{NOMBRE DE PAS DE TEMPS}.
On en deduit le nombre de pas de temps en prenant l''entier le
plus proche de (duree du calcul/pas de temps).
Si le \telkey{NOMBRE DE PAS DE TEMPS} est aussi donne,
on prend la plus grande valeur.""",
            ang = """Sets the duration of the simulation in seconds.
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
conditions aux limites de maritimes.
Egalement utilise pour les flotteurs,
les echanges thermiques avec atmosphere, le chainage avec DELWAQ.""",
            ang = """Enables to set the date of the time origin of the model when
taking into account of the tide (tide generator force and/or the tidal
boundary conditions).
Also used with drogues, heat exchange with atmosphere,
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
les echanges thermiques avec atmosphere, en chainage avec DELWAQ.""",
            ang = """Enables to set the time of the time origin of the model when
taking into account of the tide (tide generator force and/or the tidal
boundary conditions).
Also used with drogues, heat exchange with atmosphere,
chaining with DELWAQ.""",
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
Attention : cette option est etroitement liee au maillage qui doit avoir
ete saisi sur une carte marine en projection de Mercator. Il faut de
plus relever sur la carte la \telkey{LATITUDE DU POINT ORIGINE}
qui correspond dans le maillage a l''ordonnee $y$ = 0.""",
            ang = """Selection of spherical coordinates to perform the computation
(for large computation domains).
Warning: this option is closely related to the mesh that should have
been entered onto a nautical chart drawn as per Mercator projection
The \telkey{LATITUDE OF ORIGIN POINT}, which corresponds to
ordinate $y$ = 0 in the mesh, must moreover be given.""",
        ),
#       -----------------------------------
        SPATIAL_PROJECTION_TYPE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["CARTESIAN, NOT GEOREFERENCED","MERCATOR","LATITUDE LONGITUDE"],
            defaut = "MERCATOR",
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
Dans le cas de l''option 3, \telemac{3d} convertit les informations
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
treated by \telemac{3d} using Mercator projection.""",
        ),
#       -----------------------------------
        LATITUDE_OF_ORIGIN_POINT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Donne la valeur de la latitude du point origine du maillage
(pour la projection de Mercator, voir le mot cle
\telkey{SYSTEME GEOGRAPHIQUE}).
Egalement utilise pour les echanges thermiques avec l atmosphere.""",
            ang = """Gives the value of the latitude of the origin point of the
mesh (for the Mercator projection, see the keyword
\telkey{GEOGRAPHIC SYSTEM}).
Also used for heat exchange with atmosphere.""",
        ),
#       -----------------------------------
        LONGITUDE_OF_ORIGIN_POINT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Donne la valeur de la longitude du point origine du maillage
(pour la projection de Mercator, voir le mot cle
\telkey{SYSTEME GEOGRAPHIQUE}).
Egalement utilise pour la force generatrice de la maree,
echanges thermiques avec l atmosphere.""",
            ang = """Gives the value of the longitude of the origin point of the
mesh (for the Mercator projection, see the keyword
\telkey{GEOGRAPHIC SYSTEM}).
Also used for tide generating force, heat exchange with atmosphere.""",
        ),
#       -----------------------------------
        NORTH = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Angle que fait le nord,
dans le sens trigonometrique, avec L, axe Oy. Lu mais non utilise.""",
            ang = """Angle of North, counted counter-clockwise, with Oy. Read but not used.""",
        ),
    ),
)
# -----------------------------------------------------------------------
VERTICAL = PROC(nom= "VERTICAL",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NUMBER_OF_HORIZONTAL_LEVELS = SIMP(statut ='o',
#   -----------------------------------
        typ = 'I',
        defaut = 2,
        fr = """Definit le nombre de plans du maillage entre le fond et la
surface. Vaut au moins 2.""",
        ang = """Gives the number of planes from bottom to free surface. Must
be at least 2.""",
    ),
#   -----------------------------------
    MESH_TRANSFORMATION = SIMP(statut ='o',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Permet que specifier la methode de repartition des plans
verticaux du maillage. Les choix possibles sont :
\begin{itemize}
\item 0 : utilisateur (sous-programme \telfile{CALCOT} a programmer) ;
\item 1 : sigma ;
\item 2 : zstar ;
\item 3 : plans fixes ;
\item 5 : adaptatif.
\end{itemize}
Ce mot-clef doit etre coherent avec le sous-programme
\telfile{CONDIM}.""",
        ang = """Specifies the distribution of vertical planes of the mesh.
Possible choices are:
\begin{itemize}
\item 0: user defined (then subroutine \telfile{CALCOT} to be
implemented),
\item 1: sigma,
\item 2: zstar,
\item 3: horizontal fixed planes,
\item 5: adaptive mesh.
\end{itemize}
This keyword must comply with what is done in \telkey{CONDIM}
subroutine.""",
    ),
#   -----------------------------------
    MINIMUM_VOLUME_OF_3D_ELEMENTS = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-6,
        fr = """Volume minimal des elements 3D exprime en m$^3$.""",
        ang = """Minimum volume of 3D elements in m$^3$.""",
    ),
#   -----------------------------------
    MINIMUM_DISTANCE_BETWEEN_PLANES_CLOSE_TO_THE_BOTTOM = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.2,
        fr = """Distance minimale entre plans pres du fond.
Seulement pour les plans de type \telfile{TRANSF\_PLANE\%I(...)} = 3.""",
        ang = """Minimum distance between planes close to the bottom.
Only for planes of type \telfile{TRANSF\_PLANE\%I(...)} = 3.""",
    ),
#   -----------------------------------
    MINIMUM_DISTANCE_BETWEEN_PLANES_CLOSE_TO_THE_FREE_SURFACE = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.2,
        fr = """Distance minimale entre plans pres de la surface libre.
Seulement pour les plans de type \telfile{TRANSF\_PLANE\%I(...)} = 3.""",
        ang = """Minimum distance between planes close to the free surface.
Only for planes of type \telfile{TRANSF\_PLANE\%I(...)} = 3.""",
    ),
#   -----------------------------------
    THRESHOLD_HEIGHT_BEFORE_CRUSHED_ELEMENTS = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Hauteur minimale sous laquelle les elements 3D sont traites comme
ecrases.
Ce n est pas fait pour le plan de la surface libre.""",
        ang = """Minimum height below which 3D elements are treated as crushed.
This is not done for the free surface plane.""",
    ),
#   -----------------------------------
    NUMBER_OF_TRACER_FOR_AMR = SIMP(statut ='o',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Permet de donner le numero de traceur utilise par l algorithme de
maillage adaptatif (\telkey{TRANSFORMATION DU MAILLAGE} = 5).""",
        ang = """Specifies the number of tracer used by the adaptive mesh algorithm
(\telkey{MESH TRANSFORMATION} = 5).""",
    ),
)
# -----------------------------------------------------------------------
NUMERICAL_PARAMETERS = PROC(nom= "NUMERICAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NUMBER_OF_SUB_ITERATIONS_FOR_NON_LINEARITIES = SIMP(statut ='o',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Permet de reactualiser, pour un meme pas de temps, les champs
convecteur et propagateur au cours de plusieurs sous-iterations. A la
premiere sous-iteration, ces champs sont donnes par $C$ et le champ de
vitesses au pas de temps precedent. Aux iterations suivantes, ils sont
pris egaux au champ de vitesse obtenu a la fin de la sous-iteration
precedente. Cette technique permet d''ameliorer la prise en compte des
non linearites.""",
        ang = """Used for updating, within one time step, the advection and
propagation fields.
Upon the first sub-iteration, these fields are given by
$C$ and the velocity field in the previous time step. At subsequent
iterations, the results of the previous sub-iteration is used to
update the advection and propagation field.
The non-linearities can be taken into account through this technique.""",
    ),
#   -----------------------------------
    ZERO = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.E-10,
        fr = """Non active pour l''instant.""",
        ang = """Not used so far.""",
    ),
#   -----------------------------------
    ADVECTION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ADVECTION_STEP = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Prise en compte ou non des termes de convection.
En cas de reponse positive,
on peut encore supprimer certains termes de convection avec
les mots-cles \telkey{SCHEMA POUR LA CONVECTION...}""",
            ang = """Takes into account the advection terms or not.
If YES, some advection terms can still be ignored with the keywords
\telkey{SCHEME FOR ADVECTION OF...}""",
        ),
#       -----------------------------------
        b_TREATMENT_OF_FLUXES_AT_THE_BOUNDARIESF = BLOC(condition="(ADVECTION_STEP == True and ((SCHEME_FOR_ADVECTION_OF_TRACERS in ['N-SCHEME FOR TIDAL FLATS','LEO POSTMA FOR TIDAL FLATS','EXPLICIT + SUPG','EXPLICIT + MURD SCHEME N','EXPLICIT LEO POSTMA','EXPLICIT + MURD SCHEME PSI']) or (SCHEME_FOR_ADVECTION_OF_K_EPSILON in ['N-SCHEME FOR TIDAL FLATS','LEO POSTMA FOR TIDAL FLATS','EXPLICIT + SUPG','EXPLICIT + MURD SCHEME N','EXPLICIT LEO POSTMA','EXPLICIT + MURD SCHEME PSI']) or (SCHEME_FOR_ADVECTION_OF_VELOCITIES in ['N-SCHEME FOR TIDAL FLATS','LEO POSTMA FOR TIDAL FLATS','EXPLICIT + SUPG','EXPLICIT + MURD SCHEME N','EXPLICIT LEO POSTMA','EXPLICIT + MURD SCHEME PSI'])))",
#       -----------------------------------
        ),
#       -----------------------------------
        TREATMENT_OF_FLUXES_AT_THE_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["Priority to prescribed values","Priority to fluxes"],
            fr = """Option s''utilisant uniquement pour les schemas SUPG, PSI et N
(un entier par frontiere liquide).
Les choix possibles sont :
\begin{itemize}
\item 1 : priorite aux valeurs imposees ;
\item 2 : priorite aux flux.
\end{itemize}
Avec l''option 2, on ne retrouve pas exactement les valeurs imposees
des traceurs, mais le flux est correct.""",
            ang = """Used so far only with the SUPG, PSI and N schemes
(one integer per open boundary).
Possible choices are:
\begin{itemize}
\item 1: priority to prescribed values,
\item 2: priority to fluxes.
\end{itemize}
With option 2, Dirichlet prescribed values are not obeyed,
but the fluxes are correct.""",
        ),
#       -----------------------------------
        SUPG_OPTION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I', min=0, max='**',
            defaut = [1,0,1,1],
            fr = """Permet de specifier le type de decentrement utilise.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de decentrement SUPG ;
\item 1 : SUPG classique ;
\item 2 : SUPG modifiee.
\end{itemize}
Seul le premier coefficient (applique a $U$, $V$ et $W$) est utilise.
C est aussi ce coefficient qui est applique aux traceurs, $k$
et $\epsilon$ eventuellement.""",
            ang = """Specifies the type of upwinding used.
Possible choices are:
\begin{itemize}
\item 0: no upwinding,
\item 1: classical SUPG,
\item 2: modified SUPG.
\end{itemize}
Only the 1st coefficient is used (and applied to $U$, $V$ and $W$).
It is also this coefficient which is applied to tracer(s), $k$
and $\epsilon$ if needed.""",
        ),
#       -----------------------------------
        b_MAXIMUM_NUMBER_OF_ITERATIONS_FOR_ADVECTION_SCHEMESF = BLOC(condition="(ADVECTION_STEP == True and ((SCHEME_FOR_ADVECTION_OF_TRACERS == 'N-SCHEME FOR TIDAL FLATS') or (SCHEME_FOR_ADVECTION_OF_K_EPSILON == 'N-SCHEME FOR TIDAL FLATS') or (SCHEME_FOR_ADVECTION_OF_VELOCITIES == 'N-SCHEME FOR TIDAL FLATS') or (SCHEME_FOR_ADVECTION_OF_TRACERS == 'LEO POSTMA FOR TIDAL FLATS') or (SCHEME_FOR_ADVECTION_OF_K_EPSILON == 'LEO POSTMA FOR TIDAL FLATS') or (SCHEME_FOR_ADVECTION_OF_VELOCITIES == 'LEO POSTMA FOR TIDAL FLATS')))",
#       -----------------------------------
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_ITERATIONS_FOR_ADVECTION_SCHEMES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 50,
            fr = """Limite le nombre d''iterations pour les schemas de convection,
seulement pour schemes 13 et 14.
Ancienne valeur par defaut = 10 jusqu en version 8.1.""",
            ang = """Limits the number of solver iterations for the advection
schemes, only for schemes 13 and 14.
Old default value = 10 until release 8.1.""",
        ),
#       -----------------------------------
        NUMBER_OF_SUB_STEPS_OF_DISTRIBUTIVE_SCHEMES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Pour les options predicteur-correcteur avec schema localement
implicite (4 ou 5).
Ce mot-cle permet de diviser le pas de temps donne par l utilisateur
dans le \telkey{FICHIER DES PARAMETRES} en plusieurs sous-pas.
A nouveau, il produit un effet sur la precision du schema et
il est pratique d ajuster ce mot-cle afin d avoir des nombres de
Courant pas trop grands (autour de 1).""",
            ang = """Only for implicit scheme with predictor-corrector (4 or 5).
This keyword allows to subdivide the time step given by the user in the
\telkey{STEERING FILE}, into several sub-steps.
Again, it produces an effect on the precision of the scheme and
it is convenient to set this keyword in order to have Courant numbers
not too large (around 1).""",
        ),
#       -----------------------------------
        NUMBER_OF_CORRECTIONS_OF_DISTRIBUTIVE_SCHEMES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Pour les options avec predicteur-correcteur et schema de convection
de type (3, 4 et 5, LIPS ou pas).
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
3, 4, 5, LIPS or not).
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
    ),
#   -----------------------------------
    DIFFUSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        MASS_LUMPING_FOR_DIFFUSION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Mass-lumping de la matrice de masse dans la diffusion.
Si la diffusion est explicite (\telkey{IMPLICITATION FOR DIFFUSION}
< 0.001) ou qu un schema de convection est de type 3, 4, 5, 13 ou 14,
il est automatiquement mis a 1. dans \telfile{DIFF3D}.""",
            ang = """Mass-lumping of the mass-matrix in the diffusion step.
If diffusion is explicit (\telkey{IMPLICITATION FOR DIFFUSION} < 0.001)
or if the advection scheme is of type 3, 4, 5, 13 or 14,
it is automatically set to 1. in \telfile{DIFF3D}.""",
        ),
    ),
)
# -----------------------------------------------------------------------
HYDRODYNAMICS = PROC(nom= "HYDRODYNAMICS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NON_HYDROSTATIC_VERSION = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Permet de specifier s''il y a utilisation ou non de la version
non-hydrostatique.""",
        ang = """Specifies the use of the non-hydrostatic code version or not.""",
    ),
#   -----------------------------------
    b_NON_HYDROSTATIC_VERSIONG = BLOC(condition="NON_HYDROSTATIC_VERSION == True",
#   -----------------------------------
#       -----------------------------------
        SOLVER_FOR_PPE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres","direct solver"],
            defaut = "gmres",
            fr = """Permet de choisir le solveur utilise pour la resolution de
l''equation de Poisson.
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
\end{itemize}
Ancienne valeur par defaut = 1 (gradient conjugue) jusqu a la version
V8P0.""",
            ang = """Choice of the solver for the Poisson Pressure Equation.
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
\end{itemize}
Old default value = 1 (conjugate gradient) until version V8P0.""",
        ),
#       -----------------------------------
        b_SOLVER_FOR_PPEG = BLOC(condition="SOLVER_FOR_PPE == 'gmres'",
#       -----------------------------------
#           -----------------------------------
            OPTION_OF_SOLVER_FOR_PPE = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 5,
                fr = """Dimension de l''espace de Krylov pour la methode GMRES (7).
Ancienne valeur par defaut = 3 jusqu a la version V8P0.""",
                ang = """Dimension of Krylov space for the GMRES method (7).
Old default value = 3 until version V8P0.""",
            ),
        ),
#       -----------------------------------
        ACCURACY_FOR_PPE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-8,
            fr = """Fixe la precision pour l''equation de Poisson.
Ancienne valeur par defaut = 1.E-4 jusqu a la version V8P0.""",
            ang = """Sets the precision needed for the computation of the Poisson
Pressure Equation.
Old default value = 1.E-4 until version V8P0.""",
        ),
#       -----------------------------------
        MAXIMUM_NUMBER_OF_ITERATIONS_FOR_PPE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 100,
            fr = """Limite le nombre d iterations pour l''equation de Poisson.""",
            ang = """Limits the number of solver iterations for the Poisson
Pressure Equation.""",
        ),
#       -----------------------------------
        PRECONDITIONING_FOR_PPE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["no preconditioning","diagonal","diagonal condensed","diagonal with absolute values","Crout","Gauss-Seidel EBE","Matrix defined by the user","diagonal and Crout","direct solver on the vertical","diagonal condensed and Crout","diagonal and direct solver on the vertical"],
            defaut = "diagonal",
            fr = """Preconditionnement pour l''equation de Poisson.
Les choix possibles sont :
\begin{itemize}
\item 0 : aucun ;
\item 2 : diagonal ;
\item 3 : diagonal avec matrice condensee ;
\item 5 : diagonal avec valeurs absolues ;
\item 7 : Crout ;
\item 11 : Gauss-Seidel EBE ;
\item 13 : matrice fournie par l''utilisateur ;
\item 14 : diagonal et Crout ;
\item 17 : solveur direct sur la verticale ;
\item 21 : diagonal condensee et Crout ;
\item 34 : diagonal et solveur direct sur la verticale.
\end{itemize}""",
            ang = """Preconditioning for the Poisson Pressure Equation.
Possible choices are:
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal,
\item 3: diagonal with the condensed matrix,
\item 5: diagonal with absolute values,
\item 7: Crout,
\item 11: Gauss-Seidel EBE,
\item 13: matrix defined by the user,
\item 14: diagonal and Crout,
\item 17: direct solver on the vertical,
\item 21: diagonal condensed and Crout,
\item 34: diagonal and direct solver on the vertical.
\end{itemize}""",
        ),
#       -----------------------------------
        DYNAMIC_PRESSURE_IN_WAVE_EQUATION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Definit si une estimation du gradient de pression dynamique
est prise en compte dans l''equation d''onde.""",
            ang = """Defines if an estimated pressure gradient is taken into
account in the wave equation.""",
        ),
#       -----------------------------------
        DYNAMIC_BOUNDARY_CONDITION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si OUI, permet d imposer une vitesse en surface
selon la condition a la limite dynamique.""",
            ang = """If YES, it enables to prescribe a velocity at the free surface
coherent with the dynamic boundary condition.""",
        ),
#       -----------------------------------
        CONTINUITY_CORRECTION_ON_OPEN_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Modifie les vitesses libres sur les frontieres ouvertes
pour avoir un meilleur champ a divergence nulle.""",
            ang = """Changes the free velocities on open boundaries to get
a better divergence-free field.""",
        ),
    ),
#   -----------------------------------
    ELEMENTS_MASKED_BY_USER = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Si OUI, remplir le sous-programme \telfile{MASKOB}.""",
        ang = """If YES, fill in the subroutine \telfile{MASKOB}.""",
    ),
#   -----------------------------------
    b_ELEMENTS_MASKED_BY_USERG = BLOC(condition="ELEMENTS_MASKED_BY_USER == True",
#   -----------------------------------
#       -----------------------------------
        Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#       -----------------------------------
            defaut = "Rewrite subroutine maskob"),
    ),
#   -----------------------------------
    PHYSICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        GRAVITY_ACCELERATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 9.81,
            fr = """Fixe la valeur de l''acceleration de la pesanteur en m/s$^2$.""",
            ang = """Sets the value of the acceleration due to gravity in m/s$^2$.""",
        ),
#       -----------------------------------
        AVERAGE_WATER_DENSITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1025.,
            fr = """Valeur de la densite moyenne dans le domaine, voir
\telfile{DRSURR}.""",
            ang = """Average water density in the domain, see subroutine
\telfile{DRSURR}.""",
        ),
#       -----------------------------------
        FRICTION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            LAW_OF_BOTTOM_FRICTION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO FRICTION","HAALAND","CHEZY","STRICKLER","MANNING","NIKURADSE"],
                defaut = "NIKURADSE",
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
\end{itemize}
Ancienne valeur par defaut = 2 (Chezy) jusqu a la version
V7P3 et 0 (pas de frottement) en V8P0.""",
                ang = """Selects the type of formulation used for the bottom friction.
The possible laws are as follows (refer to the Principle note):
\begin{itemize}
\item 0: no friction against bottom,
\item 1: Haaland''s formula,
\item 2: Chezy''s formula,
\item 3: Strickler''s formula,
\item 4: Manning''s formula,
\item 5: Nikuradse''s formula.
\end{itemize}
Old default value = 2 (Chezy) until version V7P3
and 0 (no friction) in V8P0.""",
            ),
#           -----------------------------------
            b_LAW_OF_BOTTOM_FRICTIONG = BLOC(condition="LAW_OF_BOTTOM_FRICTION != 'NO FRICTION'",
#           -----------------------------------
#               -----------------------------------
                FRICTION_COEFFICIENT_FOR_THE_BOTTOM = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 0.01,
                    fr = """Fixe la valeur du coefficient de frottement au fond,
si constant.
Ancienne valeur par defaut = 60. jusqu a la version V8P0.""",
                    ang = """Friction coefficient on the bottom, if constant.
Old default value = 60. until version V8P0.""",
                ),
            ),
#           -----------------------------------
            LAW_OF_FRICTION_ON_LATERAL_BOUNDARIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO FRICTION","COEFFICIENT TAKEN IN BOUNDARY CONDITIONS FILE","IDEM","IDEM","IDEM","NIKURADSE"],
                defaut = "NO FRICTION",
                fr = """Selectionne le type de formulation utilisee pour le calcul
du frottement sur les parois laterales.
Les lois possibles sont les suivantes (cf. Note de principe) :
\begin{itemize}
\item 0 : pas de frottement, ou \telfile{AUBOR} donne par le
\telkey{FICHIER DES CONDITIONS AUX LIMITES} ;
\item 5 : formule de Nikuradse.
\end{itemize}""",
                ang = """Selects the type of formulation used for the friction on
lateral boundaries. The possible laws are as follows (refer to the
Principle note):
\begin{itemize}
\item 0: no friction, or \telfile{AUBOR} given by the
\telkey{BOUNDARY CONDITION FILE},
\item 5: Nikuradse''s formula.
\end{itemize}""",
            ),
#           -----------------------------------
            b_LAW_OF_FRICTION_ON_LATERAL_BOUNDARIESG = BLOC(condition="LAW_OF_FRICTION_ON_LATERAL_BOUNDARIES != 'NO FRICTION'",
#           -----------------------------------
#               -----------------------------------
                FRICTION_COEFFICIENT_FOR_LATERAL_SOLID_BOUNDARIES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 0.01,
                    fr = """Fixe la valeur du coefficient de frottement sur les parois,
si constant.
Ancienne valeur par defaut = 60. jusqu a la version V8P0.""",
                    ang = """Friction coefficient on the lateral boundaries, if constant.
Old default value = 60. until version V8P0.""",
                ),
            ),
        ),
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
            b_CORIOLISG = BLOC(condition="CORIOLIS == True",
#           -----------------------------------
#               -----------------------------------
                CORIOLIS_COEFFICIENT = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 0.,
                    fr = """Fixe la valeur du coefficient de la force de Coriolis.
Celui-ci doit etre calcule en fonction de la latitude $l$
par la formule
$FCOR = 2 \omega sin(l)$ ,
$\omega$ etant la vitesse de rotation de la terre.
$\omega$ = 7.29 10-5 rad/s.\\
Les composantes de la force de Coriolis sont alors :\\
$FU =  FCOR \times V,$\\
$FV = -FCOR \times U.$
Lorsqu''on utilise les coordonnees spheriques, le coefficient de
Coriolis est calcule automatiquement.""",
                    ang = """Sets the value of the Coriolis force coefficient,
in cartesian coordinates.
This coefficient, denoted \telfile{FCOR} in the code, should be equal to
$2 \omega \sin(l)$  where $\omega$ denotes the earth angular speed of
rotation and $l$ the latitude. $\omega$ = 7.29 10-5 rad/s.\\
The Coriolis force components are then:\\
$FU =  FCOR \times V,$\\
$FV = -FCOR \times U.$\\
When using the spherical coordinates, the Coriolis coefficient is
automatically computed.""",
                ),
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
                OPTION_FOR_WIND = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["constant in time and space","variable in time given by formatted file","variable in time and space given by formatted file"],
                    defaut = "constant in time and space",
                    fr = """Donne les options pour introduire le vent :
\begin{itemize}
\item 1 : constant en temps et en espace (donne par les mots cles
\telkey{VITESSE DU VENT SUIVANT X} et \telkey{VITESSE DU VENT SUIVANT Y}
) ;
\item 2 : variable en temps donne par fichier formate ;
\item 3 : variable en temps et en espace donne par fichier formate
ou un fichier binaire.
\end{itemize}""",
                    ang = """Gives the option for managing the wind:
\begin{itemize}
\item 1: constant in time and space, given by the keywords
\telkey{WIND VELOCITY ALONG X} and \telkey{WIND VELOCITY ALONG Y},
\item 2: variable in time and constant in space, given by formatted
file,
\item 3: variable in time and space, given by formatted file or by
a binary file.
\end{itemize}""",
                ),
#               -----------------------------------
                b_OPTION_FOR_WINDG = BLOC(condition="OPTION_FOR_WIND == 'variable in time given by formatted file' or OPTION_FOR_WIND == 'variable in time and space given by formatted file'",
#               -----------------------------------
#                   -----------------------------------
                    Consigne = SIMP(statut ="o", homo="information", typ="TXM",
#                   -----------------------------------
                        defaut = "Give the ascii atmospheric data file"),
                ),
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
                COEFFICIENT_OF_WIND_INFLUENCE = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.55E-6,
                    fr = """Fixe la valeur du coefficient d''entrainement du vent.
Voir le manuel utilisateur pour la valeur a donner.""",
                    ang = """Sets the value of the wind driving coefficient.
See the User Manual for the value to give.""",
                ),
#               -----------------------------------
                COEFFICIENT_OF_WIND_INFLUENCE_VARYING_WITH_WIND_SPEED = SIMP(statut ='o',
#               -----------------------------------
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
            ),
#           -----------------------------------
            AIR_PRESSURE = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Permet de decider si l''on prend ou non en compte l''influence
d''un champ de pression.""",
                ang = """Sets whether the influence of an atmosphere
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
                    defaut = 0.,
                    fr = """Pour ajouter un apport ou une perte d''eau en surface.""",
                    ang = """Specifies the amount of water to add or remove at the
free surface.""",
                ),
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
                ang = """Binary-coded data file containing the atmospheric data varying
in time and space on the mesh.""",
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
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
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
            RELATIVE_HUMIDITY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 50.,
                fr = """Donne la valeur de l humidite relative lorsqu elle est
constante en temps et en espace. En \%.""",
                ang = """Gives the value of relative humidity when it is constant
in time and space. In \%.""",
            ),
        ),
#       -----------------------------------
        SOURCES = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            SOURCES_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Nom du fichier contenant les informations variables
en temps des sources.""",
                ang = """Name of the file containing time-dependent
information on sources.""",
            ),
#           -----------------------------------
            GLOBAL_NUMBERS_OF_SOURCE_NODES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I', min= 2, max= 2,
                fr = """Numeros globaux des noeuds du maillage 2D sur lequels sont affectes des
points source.""",
                ang = """Global numbers of nodes in the 2D mesh that correspond to source point
locations.""",
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
            ABSCISSAE_OF_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Nombres reels donnant les abscisses d eventuelles sources de
debit (en metres). La source sera placee au noeud du maillage le plus
proche.""",
                ang = """Floats giving the abscissae of potential sources of flow rates
(in meters). The source will be located at the nearest node in the
mesh.""",
            ),
#           -----------------------------------
            ORDINATES_OF_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Nombres reels donnant les ordonnees d''eventuelles sources de
debit (en metres). La source sera placee au noeud du maillage le plus
proche.""",
                ang = """Floats giving the ordinates of potential sources of flow rates
(in meters). The source will be located at the nearest node in the
mesh.""",
            ),
#           -----------------------------------
            ELEVATIONS_OF_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Fixe la hauteur des sources.
Les sources sont automatiquement recalees sur le plan le plus proche.
L''utilisation d''un plan fixe est alors conseillee afin d''eviter que
le plan le plus proche ne change en cas de variation de la hauteur
d''eau locale.""",
                ang = """Sets the height of the sources.
The source will be located at the nearest plane in the mesh.
The use of a fixed plane is then recommended to avoid the change
of the nearest plane in case of variation of local water height.""",
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
                fr = """Permet de specifier la composante selon $x$ de la vitesse aux
sources. Si rien n''est specifie, les sources diffusent sans vitesse
dans toutes les directions (cf. cas de validation source).""",
                ang = """Specifies the compoment along $x$ of the velocities of the
sources. If nothing is specified, the sources diffuse without any
velocity in every direction (cf. validation case source).""",
            ),
#           -----------------------------------
            VELOCITIES_OF_THE_SOURCES_ALONG_Y = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Permet de specifier la composante selon $y$ de la vitesse aux
sources. Si rien n''est specifie, les sources diffusent sans vitesse
dans toutes les directions (cf. cas de validation source).""",
                ang = """Specifies the compoment along y of the velocities of the
sources.  If nothing is specified, the sources diffuse without any
velocity in every direction (cf. validation case source).""",
            ),
#           -----------------------------------
            VELOCITIES_OF_THE_SOURCES_ALONG_Z = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Permet de specifier la composante selon $z$ de la vitesse aux
sources. Si rien n''est specifie, les sources diffusent sans vitesse
dans toutes les directions (cf. cas de validation source).""",
                ang = """Specifies the compoment along $z$ of the velocities of the
sources. If nothing is specified, the sources diffuse without any
velocity in every direction (cf. validation case source).""",
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
                fr = """Active la prise en compte des courants de houle
(voir le sous-programme \telfile{TRISOU}).""",
                ang = """Wave driven currents are taken into account,
see subroutine \telfile{TRISOU}.""",
            ),
#           -----------------------------------
            b_WAVE_DRIVEN_CURRENTSG = BLOC(condition="WAVE_DRIVEN_CURRENTS == True",
#           -----------------------------------
#               -----------------------------------
                RECORD_NUMBER_IN_WAVE_FILE = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 1,
                    fr = """Numero d''enregistrement a lire par \telemac{3d} dans le
fichier des courants de houle.""",
                    ang = """Record number to be read by \telemac{3d} in the wave driven
currents file.""",
                ),
            ),
        ),
    ),
#   -----------------------------------
    BOUNDARY_CONDITIONS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        PRESCRIBED_ELEVATIONS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Valeurs des cotes imposees aux frontieres liquides.""",
            ang = """Values of the elevations prescribed at open boundaries.""",
        ),
#       -----------------------------------
        PRESCRIBED_FLOWRATES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Valeurs des debits imposes aux frontieres liquides.""",
            ang = """Values of the flowrates prescribed at open boundaries.""",
        ),
#       -----------------------------------
        PRESCRIBED_VELOCITIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Valeurs des vitesses imposees aux frontieres liquides.""",
            ang = """Values of the magnitudes of velocity prescribed at open boundaries.""",
        ),
#       -----------------------------------
        LIQUID_BOUNDARIES_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'), max='**',
            defaut = '',
            fr = """Fichier de variations en temps des conditions aux limites.""",
            ang = """File containing the variations in time of boundary conditions.""",
        ),
#       -----------------------------------
        VELOCITY_PROFILES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min= 2, max= 2,
            into = ["constant normal profile","ubor and vbor given in the conlim file","normal velocity given in ubor in the conlim file","Velocity = square root elevation","like 4 with virtual depth, see help"],
            fr = """Permet de specifier le type de profil horizontal de vitesse
(un entier par frontiere liquide).
Les choix possibles sont :
\begin{itemize}
\item 1 : profil normal constant ;
\item 2 : $u$ et $v$ donnes dans le
\telkey{FICHIER DES CONDITIONS AUX LIMITES} ;
\item 3 : vitesse normale donnee dans \telfile{UBOR} dans le
\telkey{FICHIER DES CONDITIONS AUX LIMITES} ;
\item 4 : vitesse normale en $\sqrt{h}$ ;
\item 5 : comme 4 mais hauteur virtuelle calculee avec
la surface libre la plus basse de la frontiere.
\end{itemize}""",
            ang = """Specifies the type of horizontal profile of velocities
(one integer per open boundary).
Possible choices are:
\begin{itemize}
\item 1: constant normal profile,
\item 2: $u$ and $v$ given in the
\telkey{BOUNDARY CONDITION FILE},
\item 3: normal velocity given in \telfile{UBOR} in the
\telkey{BOUNDARY CONDITION FILE},
\item 4: normal velocity in $\sqrt{h}$,
\item 5: like 4 but virtual depth based on
the lowest elevation of the boundary.
\end{itemize}""",
        ),
#       -----------------------------------
        VELOCITY_VERTICAL_PROFILES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min= 2, max= 2,
            into = ["User defined","Constant","Logarithmic"],
            fr = """Permet de specifier le type de profil vertical de vitesse
(un entier par frontiere liquide).
Les choix possibles sont :
\begin{itemize}
\item 0 : programmation utilisateur ;
\item 1 : constant ;
\item 2 : logarithmique.
\end{itemize}""",
            ang = """Specifies the type of vertical profile of velocity
(one integer per open boundary).
Possible choices are:
\begin{itemize}
\item 0: defined by user,
\item 1: constant,
\item 2: logarithmic.
\end{itemize}""",
        ),
#       -----------------------------------
        STAGE_DISCHARGE_CURVES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I', min=10, max=10,
            fr = """Indique si une courbe de tarage doit etre utilisee
pour une frontiere (une valeur par frontiere liquide) :
\begin{itemize}
\item 0 : non ;
\item 1 : Z(Q) ;
\item 2 : Q(Z). Pas encore programme.
\end{itemize}""",
            ang = """Specifies if a discharge-elevation curve must be used
for a given boundary (one value per open boundary):
\begin{itemize}
\item 0: no,
\item 1: Z(Q),
\item 2: Q(Z). Not yet implemented.
\end{itemize}""",
        ),
#       -----------------------------------
        b_STAGE_DISCHARGE_CURVESG = BLOC(condition="STAGE_DISCHARGE_CURVES != 'no'",
#       -----------------------------------
#           -----------------------------------
            STAGE_DISCHARGE_CURVES_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Nom du fichier contenant les courbes de tarage.""",
                ang = """Name of the file containing stage-discharge curves.""",
            ),
        ),
#       -----------------------------------
        OPTION_FOR_LIQUID_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
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
        TURBULENCE_REGIME_FOR_THE_BOTTOM = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["smooth","rough","rough compatibility with old versions"],
            defaut = "rough",
            fr = """Permet de definir le regime de turbulence pour le fond dans le
cas du modele de longueur de melange ou du modele $k$-$\epsilon$ :
\begin{itemize}
\item 1 : lisse ;
\item 2 : rugueux ;
\item 3 : rugueux (compatibilite avec anciennes versions).
\end{itemize}""",
            ang = """Defines the turbulence regime for the bottom in the case of a
$k$-$\epsilon$ or mixing-length model:
\begin{itemize}
\item 1: smooth,
\item 2: rough,
\item 3: rough also (for compatibility with old versions).
\end{itemize}""",
        ),
#       -----------------------------------
        TURBULENCE_REGIME_FOR_LATERAL_SOLID_BOUNDARIES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["smooth","rough"],
            defaut = "rough",
            fr = """Definit le regime de turbulence pour les parois laterales :
\begin{itemize}
\item 1 : lisse ;
\item 2 : rugueux.
\end{itemize}""",
            ang = """Defines the turbulence regime for the lateral boundaries:
\begin{itemize}
\item 1: smooth,
\item 2: rough.
\end{itemize}""",
        ),
#       -----------------------------------
        BOUNDARY_CONDITION_ON_THE_BOTTOM = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["LOG LAW FOR VELOCITIES ON BOTTOM","NO SLIP FOR VELOCITIES ON BOTTOM"],
            defaut = "LOG LAW FOR VELOCITIES ON BOTTOM",
            fr = """Specifie le type de conditions aux limites au fond.
Les choix possibles sont :
\begin{itemize}
\item 1 : conditions de Neumann pour les vitesses au fond ;
\item 2 : vitesses nulles au fond. Va de pair logiquement avec
 un bon raffinement du maillage au fond.
\end{itemize}""",
            ang = """Specifies the type of boundary conditions on the bottom
layer. Possible choices are:
\begin{itemize}
\item 1: Neumann conditions on velocity on bottom,
\item 2: velocities will be set to 0. Should be linked to
a refined mesh near the bottom.
\end{itemize}""",
        ),
#       -----------------------------------
        VELOCITY_PROJECTED_ON_SOLID_LATERAL_BOUNDARIES = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """$\vec{U}.\vec{n} = 0$ sur les parois laterales solides est force
en fin de boucle en temps.""",
            ang = """Will ensure $\vec{U}.\vec{n} = 0$ on solid lateral boundaries
by a projection at the end of time loop.""",
        ),
#       -----------------------------------
        VELOCITY_PROJECTED_ON_BOTTOM = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """$\vec{U}.\vec{n} = 0$ sur le fond est force en fin de boucle
en temps.""",
            ang = """Will ensure $\vec{U}.\vec{n} = 0$ on bottom by a projection
at the end of time loop.""",
        ),
#       -----------------------------------
        OPEN_BOUNDARY_CONDITIONS_ON_THE_BED = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Determine s''il y a des conditions ouvertes sur le fond.""",
            ang = """Defines if there are open boundary conditions
on the bed.""",
        ),
#       -----------------------------------
        PRESCRIBED_FLOWRATES_ON_THE_BED = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            defaut = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
            fr = """Fixe le debit sur les frontieres a debit impose du fond.""",
            ang = """Sets the value for flow rate on flow
rate-imposed bed boundaries.""",
        ),
    ),
#   -----------------------------------
    INITIALIZATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        INITIAL_CONDITIONS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ['ZERO ELEVATION','CONSTANT ELEVATION','ZERO DEPTH','CONSTANT DEPTH','SPECIAL','PARTICULAR','TPXO SATELLITE ALTIMETRY'],
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
satellite TPXO dont les 8 premieres composantes ont ete extraites et
sauvees dans les fichiers
\telkey{BASE BINAIRE 1/2 DE DONNEES DE MAREE} ;
\item PARTICULIERES : Les conditions initiales sur la hauteur d''eau
doivent etre precisees dans le sous-programme \telkey{USER\_CONDI3D\_H}.
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
should be stated in the \telfile{USER\_CONDI3D\_H} subroutine.
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
    ),
#   -----------------------------------
    NUMERICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
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
\end{itemize}""",
            ang = """Defines the method to store matrices. The possible choices are:
\begin{itemize}
\item 1: classical EBE,
\item 3: edge-based storage.
\end{itemize}""",
        ),
#       -----------------------------------
        MASS_LUMPING_FOR_DEPTH = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """\telemac{3d} offre la possibilite d''effectuer du
mass-lumping sur $H$.
Ceci revient a ramener tout ou partie (suivant la valeur de ce
coefficient) de la matrice \telfile{AM1 (H)} sur sa diagonale.
Cette technique permet d''accelerer le code dans des proportions tres
importantes et de le rendre egalement beaucoup plus stable. Cependant
les solutions obtenues se trouvent lissees.
Ce parametre fixe le taux de mass-lumping effectue sur $H$.""",
            ang = """\telemac{3d} offers the possibility to perform mass-lumping
on $H$.
This gathers all or part (given the value of the coefficient)
of the \telfile{AM1(H)} matrices on their diagonal.
This technique can speed-up the code a lot and also render it
more stable.
Yet, the solutions are smoothened.
This parameter sets the mass-lumping amount done for $H$.""",
        ),
#       -----------------------------------
        HYDROSTATIC_INCONSISTENCY_FILTER = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Permet de filtrer les inconsistances hydrostatiques.""",
            ang = """Allows to filter hydrostatic inconsistencies.""",
        ),
#       -----------------------------------
        DISCRETISATION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            ELEMENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                defaut = 'PRISM',
                fr = """Permet de specifier le type d''element utilise pour le calcul.
Les choix possibles sont :
\begin{itemize}
\item PRISME : maillages de triangles empiles ;
\item TETRAEDRE : decoupage en tetraedres des prismes.
\end{itemize}""",
                ang = """Specifies the type of elements used in the computation.
The possible choices are:
\begin{itemize}
\item PRISM: superimposed meshes of triangles,
\item TETRAHEDRON: the same but prisms are split into tetrahedrons.
\end{itemize}""",
            ),
        ),
#       -----------------------------------
        PROPAGATION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SOLVER_FOR_PROPAGATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres","direct solver"],
                defaut = "gmres",
                fr = """Permet de choisir le solveur utilise pour la resolution de
l''etape de propagation.
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
\end{itemize}
Ancienne valeur par defaut = 1 (gradient conjugue) jusqu a la version
V8P0.""",
                ang = """Choice of the solver for the propagation equation.
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
\end{itemize}
Old default value = 1 (conjugate gradient) until version V8P0.""",
            ),
#           -----------------------------------
            b_SOLVER_FOR_PROPAGATIONG = BLOC(condition="SOLVER_FOR_PROPAGATION == 'gmres'",
#           -----------------------------------
#               -----------------------------------
                OPTION_OF_SOLVER_FOR_PROPAGATION = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 5,
                    fr = """Dimension de l''espace de Krylov pour la methode GMRES (7).
Ancienne valeur par defaut = 3 jusqu a la version V8P0.""",
                    ang = """Dimension of Krylov space for the GMRES method (7).
Old default value = 3 until version V8P0.""",
                ),
            ),
#           -----------------------------------
            ACCURACY_FOR_PROPAGATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.E-8,
                fr = """Fixe la precision demandee pour l''etape de propagation.
Ancienne valeur par defaut = 1.E-6 jusqu a la version V8P0.""",
                ang = """Sets the accuracy needed for the computation
of the propagation step.
Old default value = 1.E-6 until version V8P0.""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_ITERATIONS_FOR_PROPAGATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 200,
                fr = """Les algorithmes utilises pour la resolution de l''etape de
propagation etant iteratifs; il est necessaire de limiter le nombre
d''iterations autorisees.
Remarque : un maximum de 40 iterations par pas de temps semble
raisonnable.""",
                ang = """Since the algorithms used for solving the propagation step are
iterative, the allowed number of iterations should be limited.
NOTE: a maximum number of 40 iterations per time step seems to be
reasonable.""",
            ),
#           -----------------------------------
            PRECONDITIONING_FOR_PROPAGATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["no preconditioning","diagonal","diagonal condensed","diagonal with absolute values","Crout","Gauss-Seidel EBE","Matrix defined by the user","diagonal and Crout","direct solver on the vertical","diagonal condensed and Crout","diagonal and direct solver on the vertical"],
                defaut = "diagonal",
                fr = """Permet de preconditionner le systeme de l''etape de propagation
afin d''accelerer la convergence lors de sa resolution. Les choix
possibles sont :
\begin{itemize}
\item 0 : aucun ;
\item 2 : diagonal ;
\item 3 : diagonal avec matrice condensee ;
\item 5 : diagonal avec valeurs absolues ;
\item 7 : Crout ;
\item 11 : Gauss-Seidel EBE ;
\item 13 : matrice fournie par l''utilisateur ;
\item 14 : diagonal et Crout ;
\item 17 : solveur direct sur la verticale ;
\item 21 : diagonal condensee et Crout ;
\item 34 : diagonal et solveur direct sur la verticale.
\end{itemize}
Certains preconditionnements sont cumulables
(les diagonaux 2 ou 3 avec les autres).
Pour cette raison on ne retient que les nombres premiers pour
designer les preconditionnements. Si l''on souhaite en cumuler
plusieurs on formera le produit des options correspondantes.""",
                ang = """Choice of the preconditioning in the propagation step linear
system that the convergence is speeded up when it is being solved.
Possible choices are:
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal,
\item 3: diagonal with the condensed matrix,
\item 5: diagonal with absolute values,
\item 7: Crout,
\item 11: Gauss-Seidel EBE,
\item 13: matrix defined by the user,
\item 14: diagonal and Crout,
\item 17: direct solver on the vertical,
\item 21: diagonal condensed and Crout,
\item 34: diagonal and direct solver on the vertical.
\end{itemize}
Some operations (either 2 or 3 diagonal preconditioning) can be
performed concurrently with the others.
Only prime numbers are therefore kept to denote the preconditioning
operations. When several of them are to be performed concurrently,
the product of relevant options shall be done.""",
            ),
#           -----------------------------------
            INITIAL_GUESS_FOR_DEPTH = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["zero","previous","extrapolation"],
                defaut = "previous",
                fr = """Tir initial du solveur de l''etape de propagation.
Offre la possibilite de modifier la valeur initiale de $\delta h$,
accroissement de $h$, a chaque iteration,
dans l''etape de propagation en utilisant les valeurs
finales de cette variable aux pas de temps precedents. Ceci peut
permettre d''accelerer la vitesse de convergence lors de la resolution
du systeme. Trois possibilites sont offertes :
\begin{itemize}
\item 0 : $\delta h$ = 0,
\item 1 : $\delta h$ = $\delta h_n$ (valeur finale de $\delta h$
 au pas de temps precedent),
\item 2 : $\delta h$ = 2 $\delta h_n$ - $\delta h_{n-1}$
(extrapolation).
\end{itemize}
Si l option 2 est utilisee avec la version non-hydrostatique,
\telkey{ORDRE DU TIR INITIAL POUR LA HAUTEUR} est automatiquement
mis a 1.""",
                ang = """Initial guess for the solver in the propagation step.
Makes it possible to modify the initial value of $\delta h$, upon each
iteration in the propagation step, by using the ultimate values this
variable had in the earlier time steps. Thus, the convergence can be
speeded up when the system is being solved. 3 options are available:
\begin{itemize}
\item 0: $\delta h$ = 0,
\item 1: $\delta h$ = $\delta h_n$  (ultimate $\delta h$ value
in the next previous time step),
\item 2: $\delta h$ = 2 $\delta h_n$ - $\delta h_{n-1}$ (extrapolation).
\end{itemize}
If option 2 with the non-hydrostatic version,
\telkey{INITIAL GUESS FOR DEPTH} is automatically set to 1.""",
            ),
#           -----------------------------------
            LINEARIZED_PROPAGATION = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Permet de lineariser l''etape de propagation,
par exemple lors de la realisation de cas tests pour lesquels on dispose
d''une solution analytique dans le cas linearise.
Il faut alors preciser la valeur de
\telkey{PROFONDEUR MOYENNE POUR LA LINEARISATION}.""",
                ang = """Provided for linearizing the propagation step, e.g. when
performing test-cases for which an analytical solution in the linearized
case is available.
Thus, the value of \telkey{MEAN DEPTH FOR LINEARIZATION} has
to be given.""",
            ),
#           -----------------------------------
            MEAN_DEPTH_FOR_LINEARIZATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Fixe la hauteur d''eau autour de laquelle s''effectue la linearisation
lorsque l''option \telkey{PROPAGATION LINEARISEE} est choisie.""",
                ang = """Sets the water depth about which the linearization is done
when the \telkey{LINEARIZED PROPAGATION} option is selected.""",
            ),
        ),
#       -----------------------------------
        ADVECTION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_DEPTH = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 5,
                fr = """Le schema conservatif (5) est desormais impose.""",
                ang = """The conservative scheme (5) is now mandatory.""",
            ),
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_VELOCITIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO ADVECTION","CHARACTERISTICS","EXPLICIT + SUPG","EXPLICIT LEO POSTMA","EXPLICIT + MURD SCHEME N","EXPLICIT + MURD SCHEME PSI","LEO POSTMA FOR TIDAL FLATS","N-SCHEME FOR TIDAL FLATS"],
                defaut = "EXPLICIT + MURD SCHEME PSI",
                fr = """Fixe le schema utilise pour la convection des vitesses.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de convection ;
\item 1 : caracteristiques ;
\item 2 : explicite + SUPG ;
\item 3 : explicite Leo Postma ;
\item 4 : explicite + MURD schema N ;
\item 5 : explicite + MURD schema PSI ;
\item 13 : Leo Postma pour bancs decouvrants ;
\item 14 : schema N pour bancs decouvrants.
\end{itemize}
Valeurs par defaut de
\telkey{SCHEMA POUR LA CONVECTION DES VITESSES} = 5
+ \telkey{OPTION DU SCHEMA POUR LA CONVECTION DES VITESSES} = 4 donnent
le schema LIPS.
Ancienne valeur par defaut = 1 (caracteristiques) jusqu a la version
V8P0.""",
                ang = """Sets the advection scheme for the velocities.
Possible choices are:
\begin{itemize}
\item 0: no convection,
\item 1: characteristics,
\item 2: explicit + SUPG,
\item 3: explicit Leo Postma,
\item 4: explicit + MURD scheme N,
\item 5: explicit + MURD scheme PSI,
\item 13: Leo Postma for tidal flats,
\item 14: N-scheme for tidal flats.
\end{itemize}
Default values for \telkey{SCHEME FOR ADVECTION OF VELOCITIES} = 5
+ \telkey{SCHEME OPTION FOR ADVECTION OF VELOCITIES} = 4
give LIPS scheme.
Old default value = 1 (characteristics) until version V8P0.""",
            ),
#           -----------------------------------
            FREE_SURFACE_GRADIENT_COMPATIBILITY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Des valeurs comprises entre 0 et 1 peuvent supprimer les
oscillations parasites.""",
                ang = """Values between 0 and 1 may suppress spurious oscillations.""",
            ),
#           -----------------------------------
            BYPASS_VOID_VOLUMES = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Accelere les schemas de convection distributifs et volumes finis
en cas de bancs decouvrants ou de transformation sigma generalisee.""",
                ang = """Will speed-up distributive and finite volumes advection
schemes in case of tidal flats or generalised sigma transformation.""",
            ),
#           -----------------------------------
            MASS_LUMPING_FOR_VELOCITIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Fixe le taux de mass-lumping effectue sur la vitesse.
Lu mais non utilise.""",
                ang = """Sets the amount of mass-lumping that is performed on
the velocity. Read but not used.""",
            ),
#           -----------------------------------
            SCHEME_OPTION_FOR_ADVECTION_OF_VELOCITIES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 4,
                fr = """Si present remplace et a priorite sur :
\telkey{OPTION POUR LES CARACTERISTIQUES} et
\telkey{OPTION DE SUPG}.
Si schema PSI ou N :
\begin{itemize}
\item 1 : explicite ;
\item 2 : predicteur-correcteur ;
\item 3 : predicteur-correcteur deuxieme ordre en temps ;
\item 4 : implicite (compatible avec bancs decouvrants).
\end{itemize}
Si pas de bancs decouvrants, l option 2 est plus rapide.
Ancienne valeur par defaut = 1 (explicite) jusqu a la version V8P0.""",
                ang = """If present replaces and has priority over:
\telkey{OPTION FOR CHARACTERISTICS} and
\telkey{SUPG OPTION}.
If N or PSI scheme:
\begin{itemize}
\item 1: explicit,
\item 2: predictor-corrector,
\item 3: predictor-corrector second-order in time,
\item 4: implicit (compatible with tidal flats).
\end{itemize}
If no tidal flats, option 2 is faster.
Old default value = 1 (explicit) until version V8P0.""",
            ),
#           -----------------------------------
            OPTION_FOR_CHARACTERISTICS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["strong","weak"],
                defaut = "strong",
                fr = """Les choix possibles sont :
\begin{itemize}
\item 1: forme forte ;
\item 2: forme faible.
\end{itemize}
Si \telkey{SCHEMA POUR LA CONVECTION...} = 1
et le mot-cle correspondant
\telkey{OPTION DU SCHEMA POUR LA CONVECTION...} = 2,
\telkey{OPTION POUR LES CARACTERISTIQUES} est automatiquement mis a 2.""",
                ang = """Possible choices are:
\begin{itemize}
\item 1: strong form,
\item 2: weak form.
\end{itemize}
If \telkey{SCHEME FOR ADVECTION OF...} = 1
and also the corresponding keyword
\telkey{SCHEME OPTION FOR ADVECTION OF...} = 2,
\telkey{OPTION FOR CHARACTERISTICS} is automatically set to 2.""",
            ),
#           -----------------------------------
            b_OPTION_FOR_CHARACTERISTICSG = BLOC(condition="OPTION_FOR_CHARACTERISTICS == 2",
#           -----------------------------------
#               -----------------------------------
                NUMBER_OF_GAUSS_POINTS_FOR_WEAK_CHARACTERISTICS = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'I',
                    defaut = 6,
                    fr = """Voir les release notes v6.3.
Nombre de points de Gauss utilises pour le calcul des caracteristiques
faibles.
6 (points) est le seul choix pour \telemac{3d}.""",
                    ang = """See release notes v6.3.
Number of Gauss points used to compute the weak characteristics.
6 (points) is the only choice for \telemac{3d}.""",
                ),
#               -----------------------------------
                MASS_LUMPING_FOR_WEAK_CHARACTERISTICS = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 0.,
                    fr = """Fixe le taux de mass-lumping qui est applique a la matrice de
masse lors de l''utilisation des caracteristiques faibles.""",
                    ang = """Sets the amount of mass-lumping that is applied to the mass
matrix when using weak characteristics.""",
                ),
            ),
        ),
#       -----------------------------------
        DIFFUSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO DIFFUSION","IMPLICIT"],
                defaut = "IMPLICIT",
                fr = """Permet de specifier si l''on utilise ou non la diffusion
des vitesses horizontales $U$ et $V$.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de diffusion,
\item 1 : implicite.
\end{itemize}""",
                ang = """Monitors the choice of the diffusion scheme
for velocities.
Possible choices are:
\begin{itemize}
\item 0: no diffusion,
\item 1: implicit.
\end{itemize}""",
            ),
#           -----------------------------------
            b_SCHEME_FOR_DIFFUSION_OF_VELOCITIESG = BLOC(condition="SCHEME_FOR_DIFFUSION_OF_VELOCITIES == 'IMPLICIT'",
#           -----------------------------------
#               -----------------------------------
                SOLVER_FOR_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres","direct solver"],
                    defaut = "conjugate gradient",
                    fr = """Permet de choisir le solveur utilise pour la resolution
de la diffusion des vitesses $U$ et $V$.
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
                    ang = """Choice of the solver for the diffusion of velocities
$U$ and $V$.
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
#               -----------------------------------
                b_SOLVER_FOR_DIFFUSION_OF_VELOCITIESG = BLOC(condition="SOLVER_FOR_DIFFUSION_OF_VELOCITIES == 'gmres'",
#               -----------------------------------
#                   -----------------------------------
                    OPTION_OF_SOLVER_FOR_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#                   -----------------------------------
                        typ = 'I',
                        defaut = 5,
                        fr = """Dimension de l''espace de Krylov pour la methode GMRES (7).
Ancienne valeur par defaut = 3 jusqu a la version V8P0.""",
                        ang = """Dimension of Krylov space for the GMRES method (7).
Old default value = 3 until version V8P0.""",
                    ),
                ),
#               -----------------------------------
                ACCURACY_FOR_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.E-8,
                    fr = """Fixe la precision demandee pour le calcul de la diffusion
de la vitesse.
Ancienne valeur par defaut = 1.E-5 jusqu a la version V8P0.""",
                    ang = """Sets the accuracy needed for the computation of the
diffusion of the velocities.
Old default value = 1.E-5 until version V8P0.""",
                ),
#               -----------------------------------
                MAXIMUM_NUMBER_OF_ITERATIONS_FOR_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 60,
                    fr = """Limite le nombre d''iterations du solveur a chaque pas
de temps pour le calcul de la diffusion de la vitesse.""",
                    ang = """Limits the number of solver iterations for the diffusion of
velocities.""",
                ),
#               -----------------------------------
                PRECONDITIONING_FOR_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["no preconditioning","diagonal","diagonal condensed","diagonal with absolute values","Crout","Gauss-Seidel EBE","Matrix defined by the user","diagonal and Crout","direct solver on the vertical","diagonal condensed and Crout","diagonal and direct solver on the vertical"],
                    defaut = "diagonal",
                    fr = """Permet de preconditionner le systeme relatif
a la diffusion des vitesses. Les choix possibles sont :
\begin{itemize}
\item 0 : aucun ;
\item 2 : diagonal ;
\item 3 : diagonal avec matrice condensee ;
\item 5 : diagonal avec valeurs absolues ;
\item 7 : Crout ;
\item 11 : Gauss-Seidel EBE ;
\item 13 : matrice fournie par l''utilisateur ;
\item 14 : diagonal et Crout ;
\item 17 : solveur direct sur la verticale ;
\item 21 : diagonal condensee et Crout ;
\item 34 : diagonal et solveur direct sur la verticale.
\end{itemize}""",
                    ang = """Choice of preconditioning for the diffusion of
velocities. Possible choices are:
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal,
\item 3: diagonal with the condensed matrix,
\item 5: diagonal with absolute values,
\item 7: Crout,
\item 11: Gauss-Seidel EBE,
\item 13: matrix defined by the user,
\item 14: diagonal and Crout,
\item 17: direct solver on the vertical,
\item 21: diagonal condensed and Crout,
\item 34: diagonal and direct solver on the vertical.
\end{itemize}""",
                ),
#               -----------------------------------
                IMPLICITATION_FOR_DIFFUSION = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.,
                    fr = """Fixe la valeur du coefficient d''implication pour l''etape de
diffusion.
Dans le cas de \telkey{OPTION POUR LA DIFFUSION} = 2, cette valeur est
ecrasee a 0 et un traitement particulier est fait pour la diffusion.""",
                    ang = """Sets the value of the implicitation coefficient for the
diffusion step.
When \telkey{OPTION FOR THE DIFFUSION} = 2, this value is changed at 0
and a specific treatment is done for the diffusion.""",
                ),
            ),
        ),
#       -----------------------------------
        NON_HYDROSTATIC = FACT(statut='o',
#       -----------------------------------
        ),
#       -----------------------------------
        IMPLICITATION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            IMPLICITATION_FOR_DEPTH = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.55,
                fr = """Fixe la valeur du coefficient d''implicitation sur la hauteur
d''eau dans l''etape de propagation (cf. Note de principe).
Les valeurs inferieures a 0.5 donnent un schema instable.""",
                ang = """Sets the value of the implicitation coefficient for water
depth in the propagation step (cf. Principe note).
The values lower than 0.5 give an instable scheme.""",
            ),
#           -----------------------------------
            IMPLICITATION_FOR_VELOCITIES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.55,
                fr = """Fixe la valeur du coefficient d''implicitation sur la vitesse
dans l''etape de propagation (cf.  Note de principe).
Les valeurs inferieures a 0.5 donnent un schema instable.
Ancienne valeur par defaut = 1. jusqu a la version V8P0.""",
                ang = """Sets the value of the implicitation coefficient
for the velocity
in the propagation step (cf. Principe note).
The values lower than 0.5 give an instable scheme.
Old default value = 1. until V8P0.""",
            ),
        ),
    ),
#   -----------------------------------
    TIDAL_FLATS_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TIDAL_FLATS = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = True,
            fr = """Permet de supprimer les tests sur les bancs decouvrants, dans
les cas ou l''on est certain qu''il n''y en aura pas.
En cas de doute, utiliser OUI.""",
            ang = """When NO, the specific treatments for tidal flats are by-passed.
This spares time, but of course you must be sure that you have no
tidal flats.""",
        ),
#       -----------------------------------
        b_TIDAL_FLATSG = BLOC(condition="TIDAL_FLATS == True",
#       -----------------------------------
#           -----------------------------------
            OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATS = SIMP(statut ='o',
#           -----------------------------------
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
#           -----------------------------------
            b_OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATSG = BLOC(condition="OPTION_FOR_THE_TREATMENT_OF_TIDAL_FLATS == 'EQUATIONS SOLVED EVERYWHERE WITH CORRECTION ON TIDAL FLATS'",
#           -----------------------------------
#               -----------------------------------
                TREATMENT_OF_NEGATIVE_DEPTHS = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["NO TREATMENT","SMOOTHING","FLUX CONTROL"],
                    defaut = "SMOOTHING",
                    fr = """Seulement avec \telkey{OPTION DE TRAITEMENT DES BANCS DECOUVRANTS}
= 1. Les choix possibles sont :
\begin{itemize}
\item 0 : pas de traitement ;
\item 1 : lissage ;
\item 2 : limitation des flux.
\end{itemize}
Si l option 2 avec bancs decouvrants est utilisee,
il est obligatoire d avoir \telkey{MASS-LUMPING SUR H} = 1.""",
                    ang = """Only with \telkey{OPTION FOR THE TREATMENT OF TIDAL FLATS}
= 1. Possible choices are:
\begin{itemize}
\item 0: no treatment,
\item 1: smoothing,
\item 2: flux control.
\end{itemize}
If using option 2 with tidal flats, it is mandatory to set
\telkey{MASS-LUMPING ON H} = 1.""",
                ),
            ),
#           -----------------------------------
            TREATMENT_ON_TIDAL_FLATS_FOR_VELOCITIES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["FORCED TO ZERO","VALUE BEFORE MASKED"],
                defaut = "FORCED TO ZERO",
                fr = """Traitement sur les bancs decouvrants a l''etape de diffusion.
\begin{itemize}
\item 0 : forcage a zero ;
\item 1 : valeur avant masquage.
\end{itemize}""",
                ang = """Treatment of tidal flats at the diffusion step for velocities.
\begin{itemize}
\item 0: forced to zero,
\item 1: value before masked.
\end{itemize}""",
            ),
#           -----------------------------------
            THRESHOLD_FOR_VISCOSITY_CORRECTION_ON_TIDAL_FLATS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.2,
                fr = """Pour les profondeurs inferieures, la viscosite sera
progressivement reduite. Voir le sous-programme
\telfile{CLIP}.""",
                ang = """Below the threshold, viscosity will be progressively
cancelled. See \telfile{CLIP} subroutine.""",
            ),
        ),
#       -----------------------------------
        MINIMAL_VALUE_FOR_DEPTH = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = -1000.,
            fr = """Fixe la valeur minimale de $H$.""",
            ang = """Sets the minimum water depth value $H$.""",
        ),
    ),
#   -----------------------------------
    TIDES = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        BINARY_DATABASE_1_FOR_TIDE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom du fichier de la base de donnees binaire 1. Dans le cas des
donnees satellitaires de l''OSU (type TPXO), ce fichier correspond aux
donnees de niveau d''eau, par exemple h\_tpxo7.2.""",
            ang = """File name for the binary database 1 of tidal harmonic
constants. In the case of the OSU satellite altimetry model (TPXO type),
this file should be for free surface level, for instance h\_tpxo7.2.""",
        ),
#       -----------------------------------
        BINARY_DATABASE_2_FOR_TIDE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom du fichier de la base de donnees binaire 2. Dans le cas des
donnees satellitaires de l''OSU (type TPXO), ce fichier correspond aux
donnees de vitesse de marees, par exemple u\_tpxo7.2.""",
            ang = """File name for the binary database 2 of tidal harmonic
constants. In the case of the OSU satellite altimetry model (TPXO type),
this file should be for tidal velocities, for instance u\_tpxo7.2.""",
        ),
#       -----------------------------------
        GEOGRAPHIC_SYSTEM = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["NO DEFAULT VALUE","DEFINED BY USER","WGS84 LONGITUDE/LATITUDE IN REAL DEGREES","WGS84 NORTHERN UTM","WGS84 SOUTHERN UTM","LAMBERT","MERCATOR PROJECTION"],
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
#       -----------------------------------
        b_GEOGRAPHIC_SYSTEMG = BLOC(condition="GEOGRAPHIC_SYSTEM in ['WGS84 NORTHERN UTM','WGS84 SOUTHERN UTM','LAMBERT']",
#       -----------------------------------
#           -----------------------------------
            ZONE_NUMBER_IN_GEOGRAPHIC_SYSTEM = SIMP(statut ='o',
#           -----------------------------------
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
#       -----------------------------------
        LAMBERT_93_CONVERSION_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom du fichier gr3df97a.txt, grille de conversion pour Lambert 93.""",
            ang = """Name of file gr3df97a.txt, conversion grid for Lambert 93.""",
        ),
#       -----------------------------------
        COEFFICIENT_TO_CALIBRATE_SEA_LEVEL = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Coefficient pour ajuster le niveau de mer.
Ce coefficient correspond d''habitude au niveau moyen de la mer
ou une valeur proche.""",
            ang = """Coefficient to calibrate the sea level.
This coefficient usually corresponds to the mean sea level
or a close value.""",
        ),
#       -----------------------------------
        b_GLOBAL_NUMBER_OF_THE_POINT_TO_CALIBRATE_HIGH_WATERF = BLOC(condition="(TIDAL_DATA_BASE == 'TPXO' and OPTION_FOR_TIDAL_BOUNDARY_CONDITIONS >=2 and OPTION_FOR_TIDAL_BOUNDARY_CONDITIONS <= 6)",
#       -----------------------------------
        ),
#       -----------------------------------
        GLOBAL_NUMBER_OF_THE_POINT_TO_CALIBRATE_HIGH_WATER = SIMP(statut ='o',
#       -----------------------------------
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
#       -----------------------------------
        MINOR_CONSTITUENTS_INFERENCE = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Pour les solutions developpees par OSU (ex. TPXO) uniquement.
Interpolation de composantes harmoniques mineures
a partir de celles lues dans les fichiers d''entree
lies aux mots-cles \telkey{BASE BINAIRE 1 DE DONNEES DE MAREE}
et \telkey{BASE BINAIRE 2 DE DONNEES DE MAREE}.""",
            ang = """For tidal solutions developed by OSU (e.g. TPXO) only.
Inference of minor constituents from the ones read in input files
linked to keywords \telkey{BINARY DATABASE 1 FOR TIDE}
and \telkey{BINARY DATABASE 2 FOR TIDE}.""",
        ),
#       -----------------------------------
        INITIAL_VELOCITIES_COMPUTED_BY_TPXO = SIMP(statut ='f',
#       -----------------------------------
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
#       -----------------------------------
        MINIMUM_DEPTH_TO_COMPUTE_TIDAL_VELOCITIES_INITIAL_CONDITIONS = SIMP(statut ='f',
#       -----------------------------------
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
#       -----------------------------------
        PHYSICAL_PARAMETERS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            TIDE_GENERATING_FORCE = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Active la prise en compte de la force generatrice de la maree""",
                ang = """The tide generating force is taken into account.""",
            ),
        ),
#       -----------------------------------
        BOUNDARY_CONDITIONS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            OPTION_FOR_TIDAL_BOUNDARY_CONDITIONS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["No tide","Real tide (recommended methodology)","Astronomical tide","Mean spring tide","Mean tide","Mean neap tide","Astronomical neap tide","Real tide (methodology before 2010)"],
                fr = """Option pour les conditions aux limites de maree.
Pour des marees reelles, l''option 1 est recommandee.
Depuis la version 7.1, ce mot-cle est un tableau avec une valeur
donnee par frontiere liquide, separee par point-virgules.
Ceci permet d''avoir des conditions de maree (ou pas) calculees
sur des frontieres liquides avec vitesses ou hauteur d''eau imposees.
Ca evite un conflit lors de l''utilisation de seuils dans le domaine.
0 est le code pour des conditions autres que des conditions de maree.
ATTENTION depuis la version 7.1 !
Les anciens modeles doivent etre changes si la frontiere de maree
n''a pas le numero 1. Dans ce cas, le mot-cle doit etre change et
plus de valeurs doivent etre donnees.
Calage possible par les mots-cles
\telkey{COEFFICIENT POUR CALAGE EN MARNAGE},
\telkey{COEFFICIENT DE CALAGE DES VITESSES DE COURANT},
et \telkey{COEFFICIENT POUR CALAGE EN NIVEAU}.
Les choix possibles sont :
\begin{itemize}
\item 0 : Pas de maree ;
\item 1 : Maree reelle (methodologie recommandee) ;
\item 2 : Maree de vive-eau exceptionnelle (coef. presque 120) ;
\item 3 : Maree de vive-eau moyenne (coef. presque 95) ;
\item 4 : Maree moyenne (coef. presque 70) ;
\item 5 : Maree de morte-eau moyenne (coef. presque 45) ;
\item 6 : Maree de morte-eau exceptionnelle (coef. presque 20) ;
\item 7 : Maree reelle (methodologie d avant 2010).
\end{itemize}""",
                ang = """Option for tidal boundary conditions.
For real tides, option 1 is recommended.
This keyword has been an array with a value given per liquid boundary,
separated by semicolons, since version 7.1.
This enables to have tidal conditions (or not) computed
on liquid boundaries with prescribed velocities or depths,
avoiding a clash when using weirs in the domain.
0 codes for conditions other than tidal.
BEWARE since version 7.1!
Old models must be changed if their tidal boundary is not number 1.
In that case this keyword must be changed and more values given.
Possible calibration with the keywords
\telkey{COEFFICIENT TO CALIBRATE TIDAL RANGE},
\telkey{COEFFICIENT TO CALIBRATE TIDAL VELOCITIES},
and \telkey{COEFFICIENT TO CALIBRATE SEA LEVEL}.
Possible choices are:
\begin{itemize}
\item 0: No tide,
\item 1: Real tide (recommended methodology),
\item 2: Astronomical tide,
\item 3: Mean spring tide,
\item 4: Mean tide,
\item 5: Mean neap tide,
\item 6: Astronomical neap tide,
\item 7: Real tide (methodology before 2010).
\end{itemize}""",
            ),
#           -----------------------------------
            b_OPTION_FOR_TIDAL_BOUNDARY_CONDITIONSG = BLOC(condition="OPTION_FOR_TIDAL_BOUNDARY_CONDITIONS != 0",
#           -----------------------------------
#               -----------------------------------
                TIDAL_DATA_BASE = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["NO DEFAULT VALUE","JMJ","TPXO","MISCELLANEOUS (LEGOS-NEA, FES20XX, PREVIMER...)"],
                    defaut = "NO DEFAULT VALUE",
                    fr = """Fournit le nom de la base de donnees utilisee pour la generation
automatique des conditions aux limites. Les choix possibles sont :
\begin{itemize}
\item 1 : JMJ ;
\item 2 : TPXO ;
\item 3 : divers (LEGOS-NEA, FES20XX, PREVIMER).
\end{itemize}
Pour JMJ, renseigner la localisation du fichier bdd\_jmj et
geofin dans les mots-cles \telkey{BASE ASCII DE DONNEES DE MAREE} et
\telkey{FICHIER DU MODELE DE MAREE}.
Pour TPXO, LEGOS-NEA, FES20XX et PREVIMER, l''utilisateur
doit telecharger les fichiers de constantes harmoniques sur internet.""",
                    ang = """Gives the name of the data base used to automatically generate
the boundary conditions. Possible choices are:
\begin{itemize}
\item 1: JMJ,
\item 2: TPXO,
\item 3: MISCELLANEOUS (LEGOS-NEA, FES20XX, PREVIMER...).
\end{itemize}
For JMJ, indicate the location of the files bdd\_jmj and geofin
with keywords \telkey{ASCII DATABASE FOR TIDE} and
\telkey{TIDAL MODEL FILE}. For TPXO, LEGOS-NEA,
FES20XX and PREVIMER, the user has to download files of harmonic
constituents on the internet.""",
                ),
#               -----------------------------------
                b_TIDAL_DATA_BASEG = BLOC(condition="TIDAL_DATA_BASE == 'TPXO'",
#               -----------------------------------
                ),
#               -----------------------------------
                b_TIDAL_DATA_BASEH = BLOC(condition="(TIDAL_DATA_BASE == 'JMJ') or (TIDAL_DATA_BASE == 'MISCELLANEOUS (LEGOS-NEA, FES20XX, PREVIMER...)')",
#               -----------------------------------
#                   -----------------------------------
                    HARMONIC_CONSTANTS_FILE = SIMP(statut ='o',
#                   -----------------------------------
                        typ = ('Fichier','All Files (*)'),
                        defaut = '',
                        fr = """Nom du fichier contenant les constantes harmoniques extraites
du fichier du modele de maree (JMJ)
ou autres atlas (FES, NEA, PREVIMER).""",
                        ang = """Name of the file containing the harmonic constants extracted
from the tidal model file (JMJ) or other atlases (FES, NEA, PREVIMER).""",
                    ),
                ),
#               -----------------------------------
                COEFFICIENT_TO_CALIBRATE_TIDAL_RANGE = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.,
                    fr = """Coefficient pour ajuster le marnage de l''onde de maree
aux frontieres maritimes.""",
                    ang = """Coefficient to calibrate the tidal range of tidal wave
at tidal open boundary conditions.""",
                ),
#               -----------------------------------
                COEFFICIENT_TO_CALIBRATE_TIDAL_VELOCITIES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 999999.,
                    fr = """Coefficient pour ajuster les composantes de vitesse
de l''onde de maree aux frontieres maritimes.
La valeur par defaut 999999. signifie que c''est la racine carree
du \telkey{COEFFICIENT DE CALAGE DU MARNAGE} qui est prise.""",
                    ang = """Coefficient to calibrate the tidal velocities of tidal wave
at tidal open boundary conditions.
Default value 999999. means that the square root of
\telkey{COEFFICIENT TO CALIBRATE TIDAL RANGE} is taken.""",
                ),
            ),
#           -----------------------------------
            TIDAL_MODEL_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de geometrie du modele dont sont extraites
les constantes harmoniques (JMJ seulement).""",
                ang = """Geometry file of the model from which harmonic constituents
are extracted (JMJ only).""",
            ),
#           -----------------------------------
            TIDAL_MODEL_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
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
#           -----------------------------------
            ASCII_DATABASE_FOR_TIDE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom de la base de donnees de constantes harmoniques
tirees du \telkey{FICHIER DU MODELE DE MAREE}.""",
                ang = """File name for the tide data base of harmonic constituents
extracted from the \telkey{TIDAL MODEL FILE}.""",
            ),
#           -----------------------------------
            LOCAL_NUMBER_OF_THE_POINT_TO_CALIBRATE_HIGH_WATER = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Numero local du point entre 1 et le nombre de points
de frontiere maritime (du FICHIER DES CONSTANTES HARMONIQUES)
ou les conditions aux limites de maree sont calculees
avec les bases de donnees JMJ, NEA, FES, PREVIMER
(sauf les bases de type TPXO).
Les ondes de maree sont dephasees par rapport a ce point
pour debuter le calcul par une pleine mer
(en marees schematiques seulement).""",
                ang = """Local number between 1 and the number of tidal boundary points
(of the \telkey{HARMONIC CONSTANTS FILE}) where the tidal boundary
conditions are computed with JMJ, NEA, FES, PREVIMER databases
(except TPXO-type databases).
The tidal constituents have their phase shifted with respect to
this point to start the simulation with a high water
(for schematic tides only).""",
            ),
#           -----------------------------------
            MINIMUM_DEPTH_TO_COMPUTE_TIDAL_VELOCITIES_BOUNDARY_CONDITIONS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.1,
                fr = """Valeur minimale de hauteur d eau utilisee pour le calcul des conditions
aux limites de courants de maree si les profondeurs sont trop petites.
Pour les solutions de marees provenant de l OSU uniquement (ex TPXO).""",
                ang = """Minimum value of water depth used to compute tidal boundary conditions
for velocities if the water depths are too small.
For tidal solutions coming from OSU only (e.g. TPXO).""",
            ),
        ),
    ),
#   -----------------------------------
    PARTICLES_TRANSPORT = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        DROGUES = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            MAXIMUM_NUMBER_OF_DROGUES = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Permet d''effectuer un suivi de flotteurs.
Fixe le nombre maximum de flotteurs a traiter lors du calcul.
Il est alors necessaire de mettre a jour le sous-programme
\telfile{USER\_FLOT3D} afin de fournir les informations sur les
positions de largage et les temps de suivi des flotteurs.
La trajectoire des flotteurs est enregistree dans le
\telkey{FICHIER ASCII DES FLOTTEURS} ou le
\telkey{FICHIER BINAIRE DES FLOTTEURS}
qui doivent etre donnes dans le fichier des parametres.""",
                ang = """Maximum number of drogues in the computation.
The user must then fill the subroutine \telfile{USER\_FLOT3D}
specifying the coordinates of the starting points,
their departure and arrival times.
The trajectory of drogues is recorded in the
\telkey{ASCII DROGUES FILE} or the \telkey{BINARY DROGUES FILE}
that must be given in the steering file.""",
            ),
#           -----------------------------------
            b_MAXIMUM_NUMBER_OF_DROGUESG = BLOC(condition="MAXIMUM NUMBER_OF_DROGUES != 0",
#           -----------------------------------
#               -----------------------------------
                ASCII_DROGUES_FILE = SIMP(statut ='o',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)','Sauvegarde'),
                    defaut = '',
                    fr = """Fichier de resultat ASCII avec les positions des flotteurs.""",
                    ang = """ASCII results file with positions of drogues.""",
                ),
#               -----------------------------------
                BINARY_DROGUES_FILE = SIMP(statut ='o',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)','Sauvegarde'),
                    defaut = '',
                    fr = """Fichier de resultat binaire avec les positions des flotteurs.""",
                    ang = """Binary results file with positions of drogues.""",
                ),
#               -----------------------------------
                DROGUES_FILE_FORMAT = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'TXM',
                    into = ['BKBINPCL','TECPLOT'],
                    defaut = 'TECPLOT',
                    fr = """Format du \telkey{FICHIER BINAIRE DES FLOTTEURS}.
Les valeurs possibles sont :
\begin{itemize}
\item BKBINPCL : format PCL binaire natif de Blue Kenue ;
\item TECPLOT : format TecPlot original (ASCII).
\end{itemize}""",
                    ang = """Format of the \telkey{BINARY DROGUES FILE}.
Possible choices are:
\begin{itemize}
\item BKBINPCL: binary PCL format native to Blue Kenue,
\item TECPLOT: original TecPlot format (ASCII).
\end{itemize}""",
                ),
#               -----------------------------------
                PRINTOUT_PERIOD_FOR_DROGUES = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 1,
                    fr = """Nombre de pas de temps entre 2 sorties de positions de
flotteurs dans le fichier des resultats binaire supplementaire
N affecte pas la qualite du calcul de la trajectoire.""",
                    ang = """Number of time steps between 2 outputs of drogues
positions in the binary file.
It does not disturb the quality of the computation of the trajectory.""",
                ),
            ),
        ),
#       -----------------------------------
        OIL_SPILL = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            OIL_SPILL_MODEL = SIMP(statut ='o',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Pour declencher le modele de nappes d''hydrocarbures,
dans ce cas le
\telkey{FICHIER DE COMMANDES HYDROCARBURES} est necessaire.""",
                ang = """Will trigger the oil spill model, in this case
the \telkey{OIL SPILL STEERING FILE} is needed.""",
            ),
#           -----------------------------------
            b_OIL_SPILL_MODELG = BLOC(condition="OIL_SPILL_MODEL == True",
#           -----------------------------------
#               -----------------------------------
                OIL_SPILL_STEERING_FILE = SIMP(statut ='o',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)'),
                    defaut = '',
                    fr = """Contient les donnees pour le
modele de nappes d''hydrocarbures.""",
                    ang = """Contains data for the oil spill model.""",
                ),
            ),
        ),
    ),
#   -----------------------------------
    HYDRAULIC_STRUCTURES = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        CULVERTS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            NUMBER_OF_CULVERTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Nombre de buses ou ponts traites comme des termes sources ou
puits. Ces buses doivent etre decrites comme des sources dans le
fichier cas. Leurs caracteristiques sont donnees dans le
\telkey{FICHIER DE DONNEES DES BUSES} (voir la documentation ecrite).""",
                ang = """Number of culverts, tubes or bridges treated as source terms.
They must be described as sources in the domain and their features
are given in the \telkey{CULVERTS DATA FILE} (see written
documentation).""",
            ),
#           -----------------------------------
            b_NUMBER_OF_CULVERTSG = BLOC(condition="NUMBER_OF_CULVERTS != 0",
#           -----------------------------------
#               -----------------------------------
                CULVERTS_DATA_FILE = SIMP(statut ='o',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)'), max='**',
                    defaut = '',
                    fr = """Fichier de description des buses/ponts presents dans le modele.""",
                    ang = """Description of culverts/bridges existing in the model.""",
                ),
#               -----------------------------------
                OPTION_FOR_CULVERTS = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 1,
                    fr = """Option pour le traitement des buses. Il existe deux formulations
dans \telemac{3d}, basees sur les formules de Bodhaine (1968)
et Carlier (1976).
Lire le \telemac{3d} theory guide plus plus d informations.""",
                    ang = """Option for the treatment of culverts. There are two options in
\telemac{3d} based on Bodhaine (1968) and Carlier (1976) formulae.
Read the \telemac{3d} theory guide for more informations.""",
                ),
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
TURBULENCE = PROC(nom= "TURBULENCE",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    PHYSICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        HORIZONTAL_TURBULENCE_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["CONSTANT VISCOSITY","K-EPSILON MODEL","SMAGORINSKI","SPALART-ALLMARAS","K-OMEGA MODEL","DES (DETACHED EDDY SIMULATION) MODEL"],
            defaut = "CONSTANT VISCOSITY",
            fr = """Permet de specifier le modele de turbulence horizontal.
Les choix possibles sont :
\begin{itemize}
\item 1 : viscosite constante ;
\item 3 : modele $k$-$\epsilon$ ;
\item 4 : Smagorinski ;
\item 5 : Spalart-Allmaras ;
\item 7 : modele $k$-$\omega$ ;
\item 9 : modele DES (Detached Eddy Simulation).
\end{itemize}
Attention : si on choisit l''option 1, il ne faut pas oublier d''ajuster
la valeur du mot cle \telkey{COEFFICIENT DE DIFFUSION}\ldots
Si on choisit les autres options, ce meme parametre doit retrouver sa
vraie valeur physique car elle est utilisee comme telle dans le modele
de turbulence.
Si on choisit l''option 3 ou 7, ce meme parametre doit retrouver sa
vraie valeur physique, en general environ $10^{-6}$
car elle est utilisee comme telle dans le modele de turbulence.""",
            ang = """Specifies the horizontal turbulence model.
The available choices are:
\begin{itemize}
\item 1: constant viscosity,
\item 3: $k$-$\epsilon$ model,
\item 4: Smagorinski,
\item 5: Spalart-Allmaras,
\item 7: $k$-$\omega$ model,
\item 9: DES model (Detached Eddy Simulation).
\end{itemize}
Caution: if option 1 is chosen, give the right
\telkey{COEFFICIENT FOR\ldots\ DIFFUSION OF VELOCITIES}\ldots\
If option 3 ou 7 is chosen, this parameter must get its real physical
value of molecular diffusivity, generally about $10^{-6}$
because it is used as well in the turbulence model.""",
        ),
#       -----------------------------------
        VERTICAL_TURBULENCE_MODEL = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["CONSTANT VISCOSITY","MIXING LENGTH","K-EPSILON MODEL","SMAGORINSKI","SPALART-ALLMARAS","GOTM","K-OMEGA MODEL","DES (DETACHED EDDY SIMULATION) MODEL"],
            defaut = "CONSTANT VISCOSITY",
            fr = """Permet de specifier le modele de turbulence vertical.
Les choix possibles sont :
\begin{itemize}
\item 1 : viscosite constante ;
\item 2 : longueur de melange ;
\item 3 : modele $k$-$\epsilon$ ;
\item 4 : Smagorinski ;
\item 5 : Spalart-Allmaras ;
\item 6 : GOTM ;
\item 7 : modele $k$-$\omega$ ;
\item 9 : modele DES (Detached Eddy Simulation).
\end{itemize}
Attention : si on choisit l''option 1, il ne faut pas oublier d''ajuster
la valeur du mot cle \telkey{COEFFICIENT DE DIFFUSION}\ldots
Si on choisit les autres options, ce meme parametre doit retrouver sa
vraie valeur physique car elle est utilisee comme telle dans le modele
de turbulence.
Si on choisit l''option 3 ou 7, ce meme parametre doit retrouver sa
vraie valeur physique, en general environ $10^{-6}$
car elle est utilisee comme telle dans le modele de turbulence.""",
            ang = """Specifies the vertical turbulence model.
The available choices are:
\begin{itemize}
\item 1: constant viscosity,
\item 2: mixing length,
\item 3: $k$-$\epsilon$ model,
\item 4: Smagorinski,
\item 5: Spalart-Allmaras,
\item 6: GOTM,
\item 7: $k$-$\omega$ model,
\item 9: DES model (Detached Eddy Simulation).
\end{itemize}
Caution: if option 1 is chosen, give the right
\telkey{COEFFICIENT FOR\ldots\ DIFFUSION OF VELOCITIES}\ldots\
If option 3 ou 7 is chosen, this parameter must get its real physical
value of molecular diffusivity, generally about $10^{-6}$
because it is used as well in the turbulence model.""",
        ),
#       -----------------------------------
        b_VERTICAL_TURBULENCE_MODELG = BLOC(condition="VERTICAL_TURBULENCE_MODEL == 'MIXING LENGTH'",
#       -----------------------------------
#           -----------------------------------
            MIXING_LENGTH_MODEL = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Permet de specifier le modele de longueur utilise pour la
turbulence verticale. Les choix possibles sont :
\begin{itemize}
\item 1: Prandtl (notamment pour ecoulements avec forte composante
barotropique comme les courants de maree) ;
\item 3: Nezu et Nakawaga ;
\item 5: Quetin (meilleure representation du vent) ;
\item 6: Tsanis (meilleure representation du vent).
\end{itemize}
4 (jet) a ete supprime.""",
                ang = """Specifies the mixing length model used for vertical turbulence.
Possible choices are:
\begin{itemize}
\item 1: Prandtl (suits such flows with a strong barotropic component
as tidal flows),
\item 3: Nezu and Nakawaga,
\item 5: Quetin (better representation of wind drift),
\item 6: Tsanis (better representation of wind drift).
\end{itemize}
4 (jet) has been suppressed.""",
            ),
#           -----------------------------------
            DAMPING_FUNCTION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Specifie le type de fonction d''amortissement utilise
(quand un modele de longueur de melange est utilise).
Les choix possibles sont :
\begin{itemize}
\item 0: rien ;
\item 1: fait par l''utilisateur (dans \telkey{USER\_DRIUTI}) ;
\item 2: Viollet ;
\item 3: Munk et Anderson.
\end{itemize}""",
                ang = """Specifies the type of damping function used (when using mixing
length turbulence model). The possible choices are:
\begin{itemize}
\item 0: nothing,
\item 1: user programmed (in \telkey{USER\_DRIUTI}),
\item 2: Viollet,
\item 3: Munk and Anderson.
\end{itemize}""",
            ),
        ),
#       -----------------------------------
        COEFFICIENT_FOR_HORIZONTAL_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-6,
            fr = """Fixe de facon uniforme pour l''ensemble du domaine;
la valeur du coefficient de diffusion de viscosite globale (dynamique +
turbulente). Cette valeur peut avoir une influence non negligeable sur
la forme et la taille des recirculations.""",
            ang = """Sets, in an even way for the whole domain, the value of the
coefficient of global (dynamic+turbulent) viscosity
for the horizontal direction. This value may
have a significant effect both on the shapes and sizes of
recirculation zones.""",
        ),
#       -----------------------------------
        COEFFICIENT_FOR_VERTICAL_DIFFUSION_OF_VELOCITIES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-6,
            fr = """Fixe de facon uniforme pour l''ensemble du domaine;
la valeur du coefficient de diffusion de viscosite globale (dynamique +
turbulente). Cette valeur peut avoir une influence non negligeable sur
la forme et la taille des recirculations.""",
            ang = """Sets, in an even way for the whole domain, the value of the
coefficient of global (dynamic+turbulent) viscosity
for the horizontal direction. This value may
have a significant effect both on the shapes and sizes of
recirculation zones.""",
        ),
#       -----------------------------------
        PRANDTL_NUMBER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.0,
            fr = """Rapport entre viscosite et diffusivite turbulente.""",
            ang = """Ratio between eddy viscosity and eddy diffusivity.""",
        ),
#       -----------------------------------
        KARMAN_CONSTANT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.4,
            fr = """Valeur de la constante de Von Karman.""",
            ang = """Value of Von Karman constant.""",
        ),
#       -----------------------------------
        GOTM_STEERING_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Fichier contenant les parametres de GOTM. Utilise avec
\telkey{MODELE DE TURBULENCE VERTICAL} = 6.""",
            ang = """File containing parameters of GOTM. Use in combination with
\telkey{VERTICAL TURBULENCE MODEL} = 6.""",
        ),
    ),
#   -----------------------------------
    BOUNDARY_CONDITIONS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        FICTITIOUS_BED_LEVEL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 2.0,
            fr = """Rapport entre le fond fictif et la hauteur de
la premiere  maille utilisee par le modele de turbulence
$k$-$\epsilon$ et pour le transport du sable.""",
            ang = """Ratio between the fictitious bed and
the grid size above the bed.""",
        ),
#       -----------------------------------
        OPTION_FOR_THE_BOUNDARY_CONDITIONS_OF_K_EPSILON = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["no turbulence","Hans Burchard"],
            defaut = "no turbulence",
            fr = """Calcul des conditions aux limites laterales sur $k$ et
$\epsilon$. Les choix possibles sont :
\begin{itemize}
\item 1: pas de turbulence = les valeurs minimales \telfile{KMIN}
et \telfile{EMIN} definies dans \telfile{CSTKEP} ;
\item 2: formule de Hans Burchard.
\end{itemize}""",
            ang = """Computation of the lateral boundary conditions of $k$
and $\epsilon$. Possible choices are:
\begin{itemize}
\item 1: no turbulence = the minimum values \telfile{KMIN} and
\telfile{EMIN} defined in \telfile{CSTKEP},
\item 2: Hans Burchard formula.
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    NUMERICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        VERTICAL_VELOCITY_DERIVATIVES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Mode de calcul des derivees des vitesses suivant $z$
uniquement pour un modele de longueur de melange sur la verticale :
\begin{itemize}
\item 1 : derivee lineaire (classique) ;
\item 2 : derivee logarithmique (mieux pour profils logarithmiques)
entre le fond et 0.2 fois la hauteur d eau.
\end{itemize}
L option 2 permet une meilleure representation des resultats
pour la modelisation d un profil de vitesse pres du fond.""",
            ang = """Way of computing the velocity derivatives along $z$
only for a mixing length model over the vertical:
\begin{itemize}
\item 1: linear derivative (classic),
\item 2: logarithmic derivative (better for logarithmic profiles)
between the bottom and 0.2 times the water depth.
\end{itemize}
Option 2 allows getting better results when modelling the velocity
profile near the bottom.""",
        ),
#       -----------------------------------
        ADVECTION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            b_SCHEME_FOR_ADVECTION_OF_K_EPSILONF = BLOC(condition="((VERTICAL_TURBULENCE_MODEL == 3) or (HORIZONTAL_TURBULENCE_MODEL == 3) or (VERTICAL_TURBULENCE_MODEL == 7) or (HORIZONTAL_TURBULENCE_MODEL == 7))",
#           -----------------------------------
            ),
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_K_EPSILON = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO ADVECTION","CHARACTERISTICS","SUPG","LEO POSTMA","MURD SCHEME N","MURD SCHEME PSI","LEO POSTMA FOR TIDAL FLATS","EXPLICIT N-SCHEME FOR TIDAL FLATS"],
                defaut = "MURD SCHEME PSI",
                fr = """Fixe le schema utilise pour la convection du modele
$k$-$\epsilon$ ou Spalart-Allmaras.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de convection ;
\item 1 : caracteristiques ;
\item 2 : SUPG ;
\item 3 : schema de Leo Postma ;
\item 4 : MURD schema N - options supplementaires possibles :
  \begin{itemize}
  \item OPTION DU SCHEMA POUR LA CONVECTION DES VITESSES/TRACEURS/
K-EPSILON) :
  \begin{itemize}
  \item 1: schema explicite classique ;
  \item 2: predicteur-correcteur au premier order ;
  \item 3: predicteur-correcteur au deuxieme order ;
  \item 4: schema implicite compatible avec les bancs decouvrants;
  \end{itemize}
  \item NOMBRE DE CORRECTIONS POUR LES SCHEMAS DISTRIBUTIFS :
  fixe le nombre de corrections pour les schemas distributifs
  predicteur-correcteur du premier et deuxieme ordre
  (une value autour de 4 est recommandee) ;
  \end{itemize}
\item 5 : MURD schema PSI ;
\item 13 : Leo Postma pour bancs decouvrants ;
\item 14 : schema N explicite pour bancs decouvrants.
\end{itemize}
Valeurs par defaut de
\telkey{SCHEMA POUR LA CONVECTION DU K-EPSILON} = 5
+ \telkey{OPTION DU SCHEMA POUR LA CONVECTION DU K-EPSILON} = 4 donnent
le schema LIPS.
Ancienne valeur par defaut = 1 (caracteristiques) jusqu a la version
V8P0.""",
                ang = """Sets the advection scheme for the $k$-$\epsilon$ model
or Spalart-Allmaras model.
Possible choices are:
\begin{itemize}
\item 0: no convection,
\item 1: characteristics,
\item 2: SUPG,
\item 3: Leo Postma scheme,
\item 4: MURD scheme N,
  \begin{itemize}
  \item SCHEME OPTION FOR ADVECTION OF VELOCITIES/TRACERS/K-EPSILON):
  \begin{itemize}
  \item 1: classical explicit scheme;
  \item 2: 1st order predictor-corrector;
  \item 3: 2nd order predictor-corrector;
  \item 4: implicit scheme compatible with tidal flats;
  \end{itemize}
  \item \telkey{NUMBER OF CORRECTIONS OF DISTRIBUTIVE SCHEMES}:
  sets the number of sub-iterations for the 1st and 2nd order
  predictor-corrector schemes (a value around 4 is recommended).
  \end{itemize}
\item 5: MURD scheme PSI,
  \begin{itemize}
  \item SCHEME OPTION FOR ADVECTION OF VELOCITIES/TRACERS/K-EPSILON):
  \begin{itemize}
  \item 1: classical explicit scheme;
  \item 2: 1st order predictor-corrector;
  \item 3: 2nd order predictor-corrector;
  \item 4: implicit scheme compatible with tidal flats;
  \end{itemize}
  \item \telkey{NUMBER OF CORRECTIONS OF DISTRIBUTIVE SCHEMES}:
  sets the number of sub-iterations for the 1st and 2nd order
  predictor-corrector schemes (a value around 4 is recommended).
  \end{itemize}
\item 13: Leo Postma for tidal flats,
\item 14: explicit N-scheme for tidal flats.
\end{itemize}
Default values for \telkey{SCHEME FOR ADVECTION OF K-EPSILON} = 5
+ \telkey{SCHEME OPTION FOR ADVECTION OF K-EPSILON} = 4
give LIPS scheme.
Old default value = 1 (characteristics) until version V8P0.""",
            ),
#           -----------------------------------
            b_SCHEME_OPTION_FOR_ADVECTION_OF_K_EPSILONF = BLOC(condition="((VERTICAL_TURBULENCE_MODEL == 3) or (HORIZONTAL_TURBULENCE_MODEL == 3) or (VERTICAL_TURBULENCE_MODEL == 7) or (HORIZONTAL_TURBULENCE_MODEL == 7))",
#           -----------------------------------
            ),
#           -----------------------------------
            SCHEME_OPTION_FOR_ADVECTION_OF_K_EPSILON = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 4,
                fr = """Si present remplace et a priorite sur :
\telkey{OPTION POUR LES CARACTERISTIQUES} et
\telkey{OPTION DE SUPG}.
Si schema PSI ou N :
\begin{itemize}
\item 1 : explicite ;
\item 2 : predicteur-correcteur ;
\item 3 : predicteur-correcteur deuxieme ordre en temps ;
\item 4 : implicite (compatible avec bancs decouvrants).
\end{itemize}
Si pas de bancs decouvrants, l option 2 est plus rapide.
Ancienne valeur par defaut = 1 (explicite) jusqu a la version V8P0.""",
                ang = """If present replaces and has priority over:
\telkey{OPTION FOR CHARACTERISTICS} and
\telkey{SUPG OPTION}.
If N or PSI scheme:
\begin{itemize}
\item 1: explicit,
\item 2: predictor-corrector,
\item 3: predictor-corrector second-order in time,
\item 4: implicit (compatible with tidal flats).
\end{itemize}
If no tidal flats, option 2 is faster.
Old default value = 1 (explicit) until version V8P0.""",
            ),
        ),
#       -----------------------------------
        DIFFUSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_DIFFUSION_OF_K_EPSILON = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO DIFFUSION","IMPLICIT"],
                defaut = "IMPLICIT",
                fr = """Permet de specifier si l''on utilise ou non la diffusion
de $k$ et $\epsilon$.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de diffusion,
\item 1 : implicite.
\end{itemize}""",
                ang = """Monitors the choice of the diffusion scheme
for $k$ and $\epsilon$.
Possible choices are:
\begin{itemize}
\item 0: no diffusion,
\item 1: implicit.
\end{itemize}""",
            ),
#           -----------------------------------
            b_SCHEME_FOR_DIFFUSION_OF_K_EPSILONG = BLOC(condition="SCHEME_FOR_DIFFUSION_OF_K_EPSILON == 'IMPLICIT'",
#           -----------------------------------
#               -----------------------------------
                SOLVER_FOR_DIFFUSION_OF_K_EPSILON = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres","direct solver"],
                    defaut = "conjugate gradient",
                    fr = """Permet de choisir le solveur utilise pour la resolution de
la diffusion du modele $k$-$\epsilon$ ou Spalart-Allmaras.
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
                    ang = """Choice of the solver for the diffusion of $k$ and $\epsilon$.
but also $\tilde{\nu}$.
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
#               -----------------------------------
                b_SOLVER_FOR_DIFFUSION_OF_K_EPSILONG = BLOC(condition="SOLVER_FOR_DIFFUSION_OF_K_EPSILON == 'gmres'",
#               -----------------------------------
#                   -----------------------------------
                    OPTION_OF_SOLVER_FOR_DIFFUSION_OF_K_EPSILON = SIMP(statut ='o',
#                   -----------------------------------
                        typ = 'I',
                        defaut = 5,
                        fr = """Dimension de l''espace de Krylov pour la methode GMRES (7).
Ancienne valeur par defaut = 3 jusqu a la version V8P0.""",
                        ang = """Dimension of Krylov space for the GMRES method (7).
Old default value = 3 until version V8P0.""",
                    ),
                ),
#               -----------------------------------
                ACCURACY_FOR_DIFFUSION_OF_K_EPSILON = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.E-8,
                    fr = """Fixe la precision demandee pour le calcul de la diffusion
du $k$-$\epsilon$ ou $\tilde{\nu}$.
Ancienne valeur par defaut = 1.E-6 jusqu a la version V8P0.""",
                    ang = """Sets the accuracy needed for the computation of the
diffusion of the $k$-$\epsilon$ or Spalart-Allmaras models.
Old default value = 1.E-6 until version V8P0.""",
                ),
#               -----------------------------------
                MAXIMUM_NUMBER_OF_ITERATIONS_FOR_DIFFUSION_OF_K_EPSILON = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 200,
                    fr = """Fixe le nombre maximum d''iterations accepte lors de la
resolution du systeme diffusion-termes sources du modele
$k$-$\epsilon$ ou du modele de Spalart-Allmaras.""",
                    ang = """Limits the number of solver iterations for the diffusion of
$k$-$\epsilon$ or $\tilde{\nu}$.""",
                ),
#               -----------------------------------
                PRECONDITIONING_FOR_DIFFUSION_OF_K_EPSILON = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ["no preconditioning","diagonal","diagonal condensed","diagonal with absolute values","Crout","Gauss-Seidel EBE","Matrix defined by the user","diagonal and Crout","direct solver on the vertical","diagonal condensed and Crout","diagonal and direct solver on the vertical"],
                    defaut = "diagonal",
                    fr = """Permet de preconditionner le systeme relatif
a la diffusion du modele $k$-$\epsilon$ ou du modele Spalart-Allmaras.
Les choix possibles sont :
\begin{itemize}
\item 0 : aucun ;
\item 2 : diagonal ;
\item 3 : diagonal avec matrice condensee ;
\item 5 : diagonal avec valeurs absolues ;
\item 7 : Crout ;
\item 11 : Gauss-Seidel EBE ;
\item 13 : matrice fournie par l''utilisateur ;
\item 14 : diagonal et Crout ;
\item 17 : solveur direct sur la verticale ;
\item 21 : diagonal condensee et Crout ;
\item 34 : diagonal et solveur direct sur la verticale.
\end{itemize}""",
                    ang = """Choice of preconditioning for the diffusion of
the $k$-$\epsilon$ model or Spalart-Allmaras model.
Possible choices are:
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal,
\item 3: diagonal with the condensed matrix,
\item 5: diagonal with absolute values,
\item 7: Crout,
\item 11: Gauss-Seidel EBE,
\item 13: matrix defined by the user,
\item 14: diagonal and Crout,
\item 17: direct solver on the vertical,
\item 21: diagonal condensed and Crout,
\item 34: diagonal and direct solver on the vertical.
\end{itemize}""",
                ),
            ),
        ),
    ),
#   -----------------------------------
    TIDAL_FLATS_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TREATMENT_ON_TIDAL_FLATS_FOR_K_EPSILON = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["FORCED TO ZERO","VALUE BEFORE MASKED"],
            defaut = "FORCED TO ZERO",
            fr = """Traitement sur les bancs decouvrants a l''etape de diffusion
pour $k$ et $\epsilon$.
\begin{itemize}
\item 0 : forcage a zero ;
\item 1 : valeur avant masquage.
\end{itemize}""",
            ang = """Treatment of tidal flats at the diffusion step for $k$ and
$\epsilon$.
\begin{itemize}
\item 0: forced to zero,
\item 1: value before masked.
\end{itemize}""",
        ),
    ),
)
# -----------------------------------------------------------------------
TRACERS = PROC(nom= "TRACERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    NUMBER_OF_TRACERS = SIMP(statut ='o',
#   -----------------------------------
        typ = 'I',
        defaut = 0,
        fr = """Definit le nombre de traceurs.""",
        ang = """Defines the number of tracers.""",
    ),
#   -----------------------------------
    b_NUMBER_OF_TRACERSG = BLOC(condition="NUMBER_OF_TRACERS != 0",
#   -----------------------------------
#       -----------------------------------
        NAMES_OF_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            fr = """Noms des traceurs en 32 caracteres, 16 pour le nom,
16 pour l''unite.""",
            ang = """Name of tracers in 32 characters, 16 for the name,
16 for the unit.""",
        ),
    ),
#   -----------------------------------
    PHYSICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        METEOROLOGY = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            VALUES_OF_TRACERS_IN_THE_RAIN = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Fixe la valeur des traceurs dans la pluie.""",
                ang = """Sets the value of the tracers in the rain.""",
            ),
        ),
#       -----------------------------------
        SOURCES = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            VALUE_OF_THE_TRACERS_AT_THE_SOURCES = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Fixe la valeur des traceurs aux sources.
Toutes les sources pour le premier traceur
puis toutes les sources du deuxieme traceur, etc.
(cf. manuel utilisateur).
Par exemple, s''il y a 3 traceurs (T1, T2 et T3) et 2 sources
(S1 et S2), la syntaxe suivante est utilisee :\\
S1\_T1;S1\_T2;S1\_T3;S2\_T1;S2\_T2;S2\_T3\\
10.0; 10.0; 0.0;  0.0; 10.0; 10.0""",
                ang = """Sets the value of the tracers at the sources.
All sources for the first tracer, then
all sources for the second tracer, etc.
(see user manual).
For example, if there are 3 tracers (T1, T2 and T3)
and 2 sources (S1 and S2), the following syntax is used:\\
S1\_T1;S1\_T2;S1\_T3;S2\_T1;S2\_T2;S2\_T3\\
10.0; 10.0; 0.0;  0.0; 10.0; 10.0""",
            ),
        ),
#       -----------------------------------
        DENSITY = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            DENSITY_LAW = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Donne le type de loi de densite utilisee dans le cas de
l''utilisation de traceurs actifs. Le sediment est considere
actif par defaut. Les choix possibles sont :
\begin{itemize}
\item 0 : rien (et du sediment si present) ;
\item 1 : fonction de temperature (et du sediment si present) ;
\item 2 : fonction de la salinite (et du sediment si present) ;
\item 3 : fonction de temperature et salinite
          (et du sediment si present) ;
\item 4 : BETA donnes par l''utilisateur pour chaque traceur
          (et effets du sediment avec un beta predefini si present);
\item 5 : le sediment et les autres traceurs sont forces a etre
passifs ;
\item 6 : Jackett et al. 2006.
\end{itemize}""",
                ang = """Gives the type of the law of density used in the case of
active tracers. The sediment is considered active with its own
law by default. The possible choices are:
\begin{itemize}
\item 0: nothing (and sediment if present),
\item 1: function of the temperature (and sediment if present),
\item 2: function of the salinity (and sediment if present),
\item 3: function of the temperature and salinity
(and sediment if present),
\item 4: user-defined BETA coefficients (and effect of the
         sediment with its own behaviour if present),
\item 5: the sediment and other tracers are forced to be passive,
\item 6: Jackett et al. 2006.
\end{itemize}""",
            ),
#           -----------------------------------
            b_DENSITY_LAWG = BLOC(condition="DENSITY_LAW == 4",
#           -----------------------------------
#               -----------------------------------
                BETA_EXPANSION_COEFFICIENT_FOR_TRACERS = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'R', min=0, max='**',
                    fr = """Unite : K$^{-1}$.
Ce coefficient permet de definir l''evolution de la densite de
l''eau en fonction de la concentration en traceur lors de
l''utilisation de la valeur 4 du mot cle
\telkey{LOI DE DENSITE}
(une valeur par traceur).""",
                    ang = """Unit: K$^{-1}$.
This coefficient is used to define the evolution of the water density
with respect to the tracer concentration when using
\telkey{DENSITY LAW} = 4 (one value per tracer).""",
                ),
#               -----------------------------------
                STANDARD_VALUES_FOR_TRACERS = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'R', min=0, max='**',
                    fr = """Valeur du traceur pour laquelle la densite est donnee lors de
l''utilisation de la valeur 4 du mot cle
\telkey{LOI DE DENSITE} (une valeur par traceur).""",
                    ang = """Reference value of tracers corresponding to the given density
when using \telkey{DENSITY LAW} = 4 (one value per tracer).""",
                ),
            ),
        ),
    ),
#   -----------------------------------
    BOUNDARY_CONDITIONS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        PRESCRIBED_TRACERS_VALUES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Determine la valeur imposee des traceurs a la premiere
frontiere, puis a la deuxieme etc. suivant la meme logique que pour les
\telkey{VALEURS DES TRACEURS DES SOURCES}.""",
            ang = """Determines the imposed value of tracers at the first boundary,
then at the second, and so on, with the same logic as
\telkey{VALUE OF THE TRACERS AT THE SOURCES}.""",
        ),
#       -----------------------------------
        TRACERS_VERTICAL_PROFILES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["User defined","Constant","Constant or Rouse if sediment", "Normalized Rouse profile and imposed conc", "Modified Rouse profile accounting for molecular viscosity"],
            fr = """Permet de specifier le type de profil de concentration des
traceurs sur la verticale.
Il y a autant de valeurs a donner que le produit du nombre de traceurs
par le nombre de frontieres liquides.
Les choix possibles sont :
\begin{itemize}
\item 0 : Programmation utilisateur ;
\item 1 : Constant ;
\item 2 : Rouse equilibrium concentration ;
\item 3 : Rouse (normalise) et concentration imposee.
\item 4 : Rouse modifie avec viscosite moleculaire.
\end{itemize}""",
            ang = """Specifies the type of profiles of tracer concentration on the
vertical.
There are as many values to be given as the product of the number of
tracers and the number of open boundaries.
Possible choices are:
\begin{itemize}
\item 0: user defined,
\item 1: constant,
\item 2: Rouse equilibrium, constant (diluted tracer)
or Rouse (sediment),
\item 3: Rouse (normalized) and imposed concentration.
\item 4: Rouse modified with molecular viscosity.
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    INITIALIZATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        INITIAL_VALUES_OF_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Fixe la valeur initiale du/des traceur(s).
Les valeurs requises sont separees par un point virgule ; si plus d une.
Le nombre de valeurs fournies doit etre egal au nombre de traceurs
declares.""",
            ang = """Sets the initial values of tracer(s).
Required value(s) separated with a semicolumn ; if more than one.
The number of supplied values must be equal to the number of declared
tracers.""",
        ),
    ),
#   -----------------------------------
    NUMERICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ADVECTION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_ADVECTION_OF_TRACERS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["NO ADVECTION","CHARACTERISTICS","EXPLICIT + SUPG","EXPLICIT LEO POSTMA","EXPLICIT + MURD SCHEME N","EXPLICIT + MURD SCHEME PSI","LEO POSTMA FOR TIDAL FLATS","N-SCHEME FOR TIDAL FLATS"],
                fr = """Fixe le schema utilise pour la convection des traceurs
(un entier par traceur).
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de convection ;
\item 1 : caracteristiques ;
\item 2 : explicite + SUPG ;
\item 3 : explicite Leo Postma ;
\item 4 : explicite + MURD schema N ;
\item 5 : explicite + MURD schema PSI ;
\item 13 : Leo Postma pour bancs decouvrants ;
\item 14 : schema N pour bancs decouvrants.
\end{itemize}
Valeurs par defaut de
\telkey{SCHEMA POUR LA CONVECTION DES TRACEURS} = 5
+ \telkey{OPTION DU SCHEMA POUR LA CONVECTION DES TRACEURS} = 4 donnent
le schema LIPS.""",
                ang = """Sets the advection scheme for the tracers (one integer per tracer).
Possible choices are:
\begin{itemize}
\item 0: no convection,
\item 1: characteristics,
\item 2: explicit + SUPG,
\item 3: explicit Leo Postma,
\item 4: explicit + MURD scheme N,
\item 5: explicit + MURD scheme PSI,
\item 13: Leo Postma for tidal flats,
\item 14: N-scheme for tidal flats.
\end{itemize}
Default values for \telkey{SCHEME FOR ADVECTION OF TRACERS} = 5
+ \telkey{SCHEME OPTION FOR ADVECTION OF TRACERS} = 4
give LIPS scheme.""",
            ),
#           -----------------------------------
            SCHEME_OPTION_FOR_ADVECTION_OF_TRACERS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I', min=0, max='**',
                fr = """Si present remplace et a priorite sur :
\telkey{OPTION POUR LES CARACTERISTIQUES} et
\telkey{OPTION DE SUPG}.
Si schema PSI ou N :
\begin{itemize}
\item 1 : explicite ;
\item 2 : predicteur-correcteur ;
\item 3 : predicteur-correcteur deuxieme ordre en temps ;
\item 4 : implicite (compatible avec bancs decouvrants).
\end{itemize}
Valeur par defaut = 4.
Si pas de bancs decouvrants, l option 2 est plus rapide.
Ancienne valeur par defaut = 1 (explicite) jusqu a la version V8P0.""",
                ang = """If present replaces and has priority over:
\telkey{OPTION FOR CHARACTERISTICS} and
\telkey{SUPG OPTION}.
If N or PSI scheme:
\begin{itemize}
\item 1: explicit,
\item 2: predictor-corrector,
\item 3: predictor-corrector second-order in time,
\item 4: implicit (compatible with tidal flats).
\end{itemize}
Default value = 4.
If no tidal flats, option 2 is faster.
Old default value = 1 (explicit) until version V8P0.""",
            ),
        ),
#       -----------------------------------
        DIFFUSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SCHEME_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["NO DIFFUSION","IMPLICIT","VERTICAL DIFFUSION ONLY"],
                defaut = "IMPLICIT",
                fr = """Permet de specifier si l''on utilise ou non la diffusion
des traceurs
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de diffusion,
\item 1 : implicite;
\item 2 : diffusion verticale seulement.
\end{itemize}""",
                ang = """Monitors the choice of the diffusion scheme
for tracers.
Possible choices are:
\begin{itemize}
\item 0: no diffusion,
\item 1: implicit,
\item 2: vertical diffusion only.
\end{itemize}""",
            ),
#           -----------------------------------
            b_SCHEME_FOR_DIFFUSION_OF_TRACERSG = BLOC(condition="SCHEME_FOR_DIFFUSION_OF_TRACERS != 'NO DIFFUSION'",
#           -----------------------------------
#               -----------------------------------
                SOLVER_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'TXM', min=0, max='**',
                    into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres","direct solver"],
                    fr = """Permet de choisir le solveur utilise pour la resolution de
la diffusion des traceurs (un entier par traceur).
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
                    ang = """Choice of the solver for the diffusion of tracers
(one integer per tracer).
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
#               -----------------------------------
                ACCURACY_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'R',
                    defaut = 1.E-8,
                    fr = """Fixe la precision demandee pour le calcul de la diffusion
des traceurs. Une seule valeur pour tous les traceurs.
Ancienne valeur par defaut = 1.E-6 jusqu a la version V7P3.""",
                    ang = """Sets the accuracy needed for the computation of
the diffusion of the tracers. One single value for every tracer.
Old default value = 1.E-6 until version V7P3.""",
                ),
#               -----------------------------------
                MAXIMUM_NUMBER_OF_ITERATIONS_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'I',
                    defaut = 60,
                    fr = """Limite le nombre d''iterations du solveur a chaque pas
de temps pour le calcul de la diffusion du ou des traceurs.""",
                    ang = """Limits the number of solver iterations for the diffusion of
tracer(s).""",
                ),
#               -----------------------------------
                PRECONDITIONING_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#               -----------------------------------
                    typ = 'TXM', min=0, max='**',
                    into = ["no preconditioning","diagonal","diagonal condensed","diagonal with absolute values","Crout","Gauss-Seidel EBE","Matrix defined by the user","diagonal and Crout","direct solver on the vertical","diagonal condensed and Crout","diagonal and direct solver on the vertical"],
                    fr = """Permet de preconditionner le systeme relatif
a la diffusion des traceurs. Les choix possibles sont :
\begin{itemize}
\item 0 : aucun ;
\item 2 : diagonal ;
\item 3 : diagonal avec matrice condensee ;
\item 5 : diagonal avec valeurs absolues ;
\item 7 : Crout ;
\item 11 : Gauss-Seidel EBE ;
\item 13 : matrice fournie par l''utilisateur ;
\item 14 : diagonal et Crout ;
\item 17 : solveur direct sur la verticale ;
\item 21 : diagonal condensee et Crout ;
\item 34 : diagonal et solveur direct sur la verticale.
\end{itemize}""",
                    ang = """Choice of preconditioning for the diffusion of tracers.
Possible choices are:
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal,
\item 3: diagonal with the condensed matrix,
\item 5: diagonal with absolute values,
\item 7: Crout,
\item 11: Gauss-Seidel EBE,
\item 13: matrix defined by the user,
\item 14: diagonal and Crout,
\item 17: direct solver on the vertical,
\item 21: diagonal condensed and Crout,
\item 34: diagonal and direct solver on the vertical.
\end{itemize}""",
                ),
            ),
#           -----------------------------------
            OPTION_OF_SOLVER_FOR_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 5,
                fr = """Dimension de l''espace de Krylov pour la methode GMRES (7).
Ancienne valeur par defaut = 3 jusqu a la version V8P0.""",
                ang = """Dimension of Krylov space for the GMRES method (7).
Old default value = 3 until version V8P0.""",
            ),
        ),
    ),
#   -----------------------------------
    TIDAL_FLATS_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TREATMENT_ON_TIDAL_FLATS_FOR_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["FORCED TO ZERO","VALUE BEFORE MASKED"],
            defaut = "FORCED TO ZERO",
            fr = """Traitement sur les bancs decouvrants a l''etape de diffusion.
\begin{itemize}
\item 0 : forcage a zero ;
\item 1 : valeur avant masquage.
\end{itemize}
Utiliser le choix 1 pour assurer la convervation du/des traceur(s)""",
            ang = """Treatment of tidal flats at the diffusion step for tracers.
\begin{itemize}
\item 0: forced to zero,
\item 1: value before masked.
\end{itemize}
Use choice 1 to ensure conservation of tracer(s).""",
        ),
    ),
#   -----------------------------------
    TURBULENCE = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        COEFFICIENT_FOR_HORIZONTAL_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Fixe les valeurs de coefficients de diffusion horizontal des
traceurs.  L''influence de ce parametre sur l''evolution des traceurs
dans le temps est importante.
C est un tableau depuis la version 7.1, avec une valeur par traceur,
separation par un point virgule.""",
            ang = """Sets the values of the horizontal diffusion of tracers.
These values may have a significant effect on the evolution of
tracers in time.
Since version 7.1, it has been an array, with one value per tracer,
separated by semicolons.""",
        ),
#       -----------------------------------
        COEFFICIENT_FOR_VERTICAL_DIFFUSION_OF_TRACERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Fixe les valeurs de coefficients de diffusion vertical des
traceurs.  L''influence de ce parametre sur l''evolution des traceurs
dans le temps est importante.
C est un tableau depuis la version 7.1, avec une valeur par traceur,
separation par un point virgule.""",
            ang = """Sets the values of the vertical diffusion of tracers.
These values may have a significant effect on the evolution of
tracers in time.
Since version 7.1, it has been an array, with one value per tracer,
separated by semicolons.""",
        ),
    ),
)
# -----------------------------------------------------------------------
SEDIMENT_INFO = PROC(nom= "SEDIMENT_INFO",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    SEDIMENT = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Permet de prendre en compte le transport sedimentaire.""",
        ang = """If YES, sediment transport is modelled.""",
    ),
#   -----------------------------------
    DENSITY_OF_THE_SEDIMENT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 2650.,
        fr = """Fixe la valeur de la masse volumique du sediment (kg/m$^3$).""",
        ang = """Value of the sediment density (kg/m$^3$).""",
    ),
#   -----------------------------------
    TIME_STEP_FOR_CONSOLIDATION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1200.,
        fr = """Valeur du pas de temps pour le modele de consolidation qui
peut etre plus grand que le pas de temps hydrodynamique car le
phenomene est tres lent. Ce parametre est utilise si
\telkey{OPTION DU MODELE DE TASSEMENT} = 1 (Modele multicouches
empirique) ou 2 (Modele de Gibson (Lenormant)).""",
        ang = """Time step for the modelling consolidation, which can
be greater than the hydrodynamic time step. This parameter is
used if \telkey{CONSOLIDATION MODEL} = 1 (Empirical multilayer model)
or 2 (Gibson model (Lenormant)).""",
    ),
#   -----------------------------------
    COHESIVE_SEDIMENT = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Permet de dire si le sediment est cohesif ou non.""",
        ang = """Tells if the sediment is cohesive or not.""",
    ),
#   -----------------------------------
    SHIELDS_PARAMETER = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.047,
        fr = """Utilise pour determiner la valeur de la contrainte critique
d''entrainement.""",
        ang = """Used to determine the critical bed shear stress value.""",
    ),
#   -----------------------------------
    MIXED_SEDIMENT = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Si OUI, calcul en sediments mixtes, il y aura un sediment
cohesif et un sediment non cohesif.""",
        ang = """If YES, calculation of mixed sediment transport, there will be
one cohesive sediment and one non cohesive sediment.""",
    ),
#   -----------------------------------
    NUMBER_OF_SEDIMENT_BED_LAYERS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Structure verticale du lit cohesif, le nombre de couches
doit etre inferieur a 20.""",
        ang = """Number of cohesive sediment bed layers, should be less
than 20.""",
    ),
#   -----------------------------------
    INPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DATA = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            READ_CRITICAL_BED_SHEAR_STRESS_PER_LAYER = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Lecture de la contrainte critique d''erosion a
partir du \telkey{FICHIER DE GEOMETRIE}.""",
                ang = """Decides if erosion shear stress at each layer is
read from \telkey{GEOMETRY FILE}.""",
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
            SEDIMENTOLOGICAL_RESULT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les variables
decrivant le fond vaseux (epaisseurs et concentrations\ldots)
avec la periodicite donnee par le mot cle \telkey{PERIODE POUR
LES SORTIES GRAPHIQUES}.""",
                ang = """Name of the file into which the sedimentological computation
results (thickness and concentration of the mud bed\ldots) shall be
written, the periodicity being given by the keyword
\telkey{GRAPHIC PRINTOUT PERIOD}.""",
            ),
#           -----------------------------------
            SEDIMENTOLOGICAL_RESULT_FILE_BINARY = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['STD','IBM','I3E'],
                defaut = 'STD',
                fr = """Type du binaire utilise pour l''ecriture du fichier
des resultats sedimentologiques.
Ce type depend de la machine sur laquelle le fichier a ete genere.
Les valeurs possibles sont :
\begin{itemize}
\item IBM pour un fichier cree sur IBM ;
\item I3E pour un fichier cree sur HP ;
\item STD.
\end{itemize}
Il s''agit alors d''ordres READ et WRITE normaux.""",
                ang = """Binary file type used for writing the results file.
This type depends on the machine on which the file was generated.
The possible values are as follows:
\begin{itemize}
\item IBM, for a file on an IBM (from a CRAY),
\item I3E, for a file on an HP (from a CRAY),
\item STD, binary type of the machine on which the user is working.
\end{itemize}
In that case, normal READ and WRITE commands are used.""",
            ),
        ),
    ),
#   -----------------------------------
    RESTART = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        PREVIOUS_COMPUTATION_SEDIMENTOLOGICAL_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = ('Fichier','All Files (*)'),
            defaut = '',
            fr = """Nom d''un fichier contenant les variables sedimentologiques
decrivant le fond vaseux, resultats d''un calcul precedent realise
sur le meme maillage et dont le dernier pas de temps enregistre
va fournir les conditions initiales pour une suite de de calcul.""",
            ang = """Name of a file containing the sedimentological parameters
(thickness and concentration of the bed\ldots), results of an earlier
computation which was made on the same mesh. The last recorded time
step will provide the initial conditions for the new computation.""",
        ),
    ),
#   -----------------------------------
    PHYSICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        MEAN_DIAMETER_OF_THE_SEDIMENT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = .01,
            fr = """Valeur du diametre D50 pour les sediments non cohesifs.""",
            ang = """Sets the value of the diameter D50 for non cohesive sediments.""",
        ),
#       -----------------------------------
        FRICTION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            RATIO_BETWEEN_SKIN_FRICTION_AND_MEAN_DIAMETER = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 3.0,
                fr = """Ratio pour le calcul du frottement de peau.
rugosite de peau = ratio $\times$ diametre moyen.""",
                ang = """ Ratio for the computation of skin friction.
skin roughness = ratio $\times$ mean diameter.""",
            ),
#           -----------------------------------
            SKIN_FRICTION_CORRECTION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Prise en compte du frottement de peau :
\begin{itemize}
\item 0 : pas de correction (TAUP = TOB) voir aussi
\telkey{RATIO BETWEEN SKIN FRICTION AND MEAN DIAMETER} :
\telfile{S3D\_KSPRATIO} ;
\item 1 : fond plat (KSP = \telfile{S3D\_KSPRATIO} $\times$
\telfile{S3D\_D50}) ;
\item 2 : prise en compte des rides (non programme).
\end{itemize}""",
                ang = """Formula to predict the skin bed roughness:
\begin{itemize}
\item 0: No correction (TAUP = TOB) see also
\telkey{RATIO ENTRE LA RUGOSITE DE PEAU ET LE DIAMETRE MOYEN}
\telfile{S3D\_KSPRATIO},
\item 1: Flat bed (KSP = \telfile{S3D\_KSPRATIO} $\times$
\telfile{S3D\_D50}),
\item 2: Ripple correction factor (not yet implemented).
\end{itemize}""",
            ),
        ),
    ),
#   -----------------------------------
    INITIALIZATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        INITIAL_PERCENTAGE_OF_NON_COHESIVE_SEDIMENT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Pourcentage initial du sediment non cohesif (mixte).""",
            ang = """Initial percentage of non cohesive sediment (mixed sediments).""",
        ),
#       -----------------------------------
        MUD_CONCENTRATIONS_PER_LAYER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Concentration du lit de vase en g/L (definie par couches) en
commencant par la couche du fond.""",
            ang = """Dry density of the mud-bed layers in g/L starting
form the bottom upwards.""",
        ),
#       -----------------------------------
        CRITICAL_EROSION_SHEAR_STRESS_OF_THE_MUD_LAYERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Taux critique d erosion de la vase (N/m$^2$).
Doit etre defini pour chaque couche en commencant par la couche de
fond.""",
            ang = """Critical erosion shear stress of the mud per layer
(N/m$^2$).
Needs to be defined for each layer (N/m$^2$),
starting from the condolidated bottom layer upwards.""",
        ),
#       -----------------------------------
        INITIAL_THICKNESS_OF_SEDIMENT_LAYERS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Epaisseurs initiales des sediments (m).""",
            ang = """Sediment layers thickness (m) for initialisation.""",
        ),
    ),
#   -----------------------------------
    NUMERICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DIFFUSION = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            SOLVER_FOR_DIFFUSION_OF_THE_SEDIMENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["conjugate gradient","conjugate residual","conjugate gradient on a normal equation","minimum error","squared conjugate gradient","cgstab","gmres","direct solver"],
                defaut = "conjugate gradient on a normal equation",
                fr = """Permet de choisir le solveur utilise pour la resolution de
la diffusion du sediment.
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
                ang = """Choice of the solver for the sediment equation.
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
            b_SOLVER_FOR_DIFFUSION_OF_THE_SEDIMENTG = BLOC(condition="SOLVER_FOR_DIFFUSION_OF_THE_SEDIMENT == 'gmres'",
#           -----------------------------------
#               -----------------------------------
                OPTION_OF_SOLVER_FOR_DIFFUSION_OF_THE_SEDIMENT = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'I',
                    defaut = 5,
                    fr = """Dimension de l''espace de Krylov pour la methode GMRES (7).
Ancienne valeur par defaut = 3 jusqu a la version V8P0.""",
                    ang = """Dimension of Krylov space for the GMRES method (7).
Old default value = 3 until version V8P0.""",
                ),
            ),
#           -----------------------------------
            ACCURACY_FOR_DIFFUSION_OF_SEDIMENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.E-8,
                fr = """Fixe la precision demandee pour le calcul de la diffusion
des sediments.
Ancienne valeur par defaut = 1.E-6 jusqu a la version V8P0.""",
                ang = """Sets the accuracy needed for the computation of the
diffusion of sediments.
Old default value = 1.E-6 until version V8P0.""",
            ),
#           -----------------------------------
            MAXIMUM_NUMBER_OF_ITERATIONS_FOR_DIFFUSION_OF_SEDIMENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 60,
                fr = """Limite le nombre d''iterations du solveur a chaque pas
de temps pour le calcul de la diffusion du sediment.""",
                ang = """Limits the number of solver iterations for the diffusion of
sediment.""",
            ),
#           -----------------------------------
            PRECONDITIONING_FOR_DIFFUSION_OF_THE_SEDIMENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["no preconditioning","diagonal","diagonal condensed","diagonal with absolute values","Crout","Gauss-Seidel EBE","Matrix defined by the user","diagonal and Crout","direct solver on the vertical","diagonal condensed and Crout","diagonal and direct solver on the vertical"],
                defaut = "diagonal",
                fr = """Permet de preconditionner le systeme relatif
a la diffusion du sediment. Les choix possibles sont :
\begin{itemize}
\item 0 : aucun ;
\item 2 : diagonal ;
\item 3 : diagonal avec matrice condensee ;
\item 5 : diagonal avec valeurs absolues ;
\item 7 : Crout ;
\item 11 : Gauss-Seidel EBE ;
\item 13 : matrice fournie par l''utilisateur ;
\item 14 : diagonal et Crout ;
\item 17 : solveur direct sur la verticale ;
\item 21 : diagonal condensee et Crout ;
\item 34 : diagonal et solveur direct sur la verticale.
\end{itemize}
Certains preconditionnements sont cumulables
(les diagonaux 2 ou 3 avec les autres).
Pour cette raison on ne retient que les nombres premiers pour
designer les preconditionnements. Si l''on souhaite en cumuler
plusieurs on formera le produit des options correspondantes.""",
                ang = """Choice of the preconditioning in the sediment diffusion
system that the convergence is speeded up when it is being solved.
Possible choices are:
\begin{itemize}
\item 0: no preconditioning,
\item 2: diagonal,
\item 3: diagonal with the condensed matrix,
\item 5: diagonal with absolute values,
\item 7: Crout,
\item 11: Gauss-Seidel EBE,
\item 13: matrix defined by the user,
\item 14: diagonal and Crout,
\item 17: direct solver on the vertical,
\item 21: diagonal condensed and Crout,
\item 34: diagonal and direct solver on the vertical.
\end{itemize}
Some operations (either 2 or 3 diagonal preconditioning) can be
performed concurrently with the others.
Only prime numbers are therefore kept to denote the preconditioning
operations. When several of them are to be performed concurrently,
the product of relevant options shall be done.""",
            ),
        ),
    ),
#   -----------------------------------
    TIDAL_FLATS_INFO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        THRESHOLD_FOR_SEDIMENT_FLUX_CORRECTION_ON_TIDAL_FLATS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.2,
            fr = """Pour les profondeurs inferieures a cette valeur limite, le flux
sedimentaire sera nul. Voir le sous-programme \telfile{FLUSED}.""",
            ang = """Below this limiting depth, all sediment erosion rates are set
to zero. See subroutine \telfile{FLUSED}.""",
        ),
    ),
#   -----------------------------------
    DEPOSITION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        CRITICAL_SHEAR_STRESS_FOR_DEPOSITION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.2,
            fr = """Fixe la valeur de la contrainte de cisaillement au fond
au dessous de laquelle se produit le depot des sediments cohesifs.""",
            ang = """Value of the critical bottom shear stress under which
deposition of cohesive sediments occurs.""",
        ),
#       -----------------------------------
        NON_COHESIVE_BED_POROSITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.4,
            fr = """La concentration massique du lit \telfile{S3D\_CFDEP} est definie par
\telfile{S3D\_CFDEP} = (1-\telfile{S3D\_XKV}) $\times$
\telfile{S3D\_RHOS}.  Ce parametre est utilise pour les sediments
non-cohesifs.""",
            ang = """The bed volume concentration \telfile{S3D\_CFDEP} =
(1-\telfile{S3D\_XKV}) $\times$ \telfile{S3D\_RHOS} is used to calculate
the bed evolution of non-cohesive sand transport.""",
        ),
    ),
#   -----------------------------------
    EROSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        EROSION_COEFFICIENT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 2.E-3,
            fr = """Valeur du coefficient d''erosion utilise dans la formule
de Partheniades en kg/m$^2$/s.""",
            ang = """Value of the erosion coefficient used in Partheniades
formula in kg/m$^2$/s.""",
        ),
    ),
#   -----------------------------------
    SETTLING_VELOCITY = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        CONSTANT_SEDIMENT_SETTLING_VELOCITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.01,
            fr = """Vitesse de chute constante en m/s (> 0 depuis v6.3).
Valeur imposee si
\telkey{INFLUENCE DE LA TURBULENCE SUR LA VITESSE DE CHUTE}
= NON.""",
            ang = """Constant sediment settling velocity in m/s (>0 since v6.3).
Prescribed value if
\telkey{INFLUENCE OF TURBULENCE ON SETTLING VELOCITY} = NO.""",
        ),
#       -----------------------------------
        SETTLING_VELOCITY_OF_SANDS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Vitesse de chute du sediment non cohesif.""",
            ang = """Non cohesive sediment settling velocity.""",
        ),
#       -----------------------------------
        ADVECTION_DIFFUSION_SCHEME_WITH_SETTLING_VELOCITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Choix de schema vertical pour la diffusion et le depot du
sediment :
\begin{itemize}
\item 0 : Diffusion implicite ;
\item 1 : Schema implicite de convection-diffusion
(tridiagonal matrix solver) ;
\item 2 : Convection faible  \telfile{sed\_fall.f}
\end{itemize}""",
            ang = """Choice of the vertical scheme for diffusion and settling of
sediment:
\begin{itemize}
\item 0: Implicit-diffusion scheme,
\item 1: Implicit-convection scheme (Tridiagonal matrix solver),
\item 2: \telfile{set\_fall.f}
\end{itemize}""",
        ),
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
    ),
#   -----------------------------------
    SUSPENSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        REFERENCE_CONCENTRATION_FORMULA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """\begin{itemize}
\item 1 : formule de Zyserman et Fredsoe, formule d''equilibre ;
\item 3 : formule de Van Rijn (1987).
\end{itemize}""",
            ang = """\begin{itemize}
\item 1: Zyserman and Fredsoe, equilibrium formula,
\item 3: Van Rijn formula (1987).
\end{itemize}""",
        ),
    ),
#   -----------------------------------
    FLOCCULATION_INFO = FACT(statut='o',
#   -----------------------------------
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
    DEPRECATED = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        CLEANING_TO_BE_DONE = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            BED_LAYERS_THICKNESS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 5.E-3,
                fr = """Epaisseur de reference pour creer de nouvelles couches
de vase. Ce parametre est utilise seulement dans le cas
\telkey{OPTION DU MODELE DE TASSEMENT} = 2
(modele de Gibson (Lenormant)). Avec ce modele, le sediment
qui se depose sur le fond est tout d''abord stocke dans une couche
tampon appelee couche des depots frais. C''est seulement quand
l''epaisseur de cette couche atteint la valeur donnee par le mot
cle \telkey{EPAISSEUR DES COUCHES DU FOND VASEUX} qu''une nouvelle
couche est cree au niveau du lit de vase.""",
                ang = """Reference thickness considered for the creation of new
bed layers.
This parameter is used if \telkey{CONSOLIDATION MODEL} = 2
(Gibson model (Lenormant)).
With this model, the sediment which settles on
the bottom arrives at first in the fresh deposit layer. When
the thickness of this layer is equal to the
\telkey{BED LAYERS THICKNESS},
a new mud layer is added to the mud bed.""",
            ),
#           -----------------------------------
            MAXIMUM_CONCENTRATION_OF_THE_CONSOLIDATED_MUD = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 500.,
                fr = """Concentration maximale pouvant etre atteinte par une couche
de vase lors du tassement.
Ce parametre est utilise si \telkey{OPTION DU MODELE DE TASSEMENT} = 2
(Modele de Gibson - Lenormant).""",
                ang = """Maximum concentration which may be reached by a mud layer
during consolidation.
This value is used if \telkey{CONSOLIDATION MODEL} = 2
(Gibson model (Lenormant)).""",
            ),
#           -----------------------------------
            RESIDENCE_TIME_FOR_MUD = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=30, max=30,
                fr = """Tableau contenant les temps de sejour en heure et centieme
relatifs a chacune des couches discretisant le fond vaseux
(la premiere valeur correspond a la couche du fond et la derniere
correspond a la couche superficielle).
Valeurs necessaires si \telkey{OPTION DU MODELE DE TASSEMENT} = 1
(Modele multicouches empirique).""",
                ang = """Array which contains the residence times of the mud bed
layers (the first value is related to the bottom layer and the
last one to the top layer).
These values are needed when \telkey{CONSOLIDATION MODEL} = 1
(Empirical multilayer model).""",
            ),
        ),
#       -----------------------------------
        TO_BE_CHECKED = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            CONSOLIDATION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Logique pour la prise en compte du tassement des depots vaseux
a l''aide d''un modele multicouches : les couches discretisant le fond
sont caracterisees par leur temps de sejour, temps au bout duquel la
vase presente dans cette couche bascule dans la couche suivante plus
consolidee.""",
                ang = """If this key word is equal to YES, consolidation is simulated
thanks to a multi-layers model: the bed layers are characterized by
their residence time which is the time after which the quantity of
mud which remains in a layer goes into a more consolidated layer.""",
            ),
#           -----------------------------------
            CONSOLIDATION_MODEL = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Choix du modele de tassement :
\begin{itemize}
\item 1 : Modele multicouches empirique ;
\item 2 : Modele de Gibson (Lenormant).
\end{itemize}""",
                ang = """Choice of the consolidation model:
\begin{itemize}
\item 1: Empirical multilayer model,
\item 2: Gibson model (Lenormant).
\end{itemize}""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
COUPLING = PROC(nom= "COUPLING",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    COUPLING_WITH = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        into = ['','SISYPHE','TOMAWAC','TOMAWACT3D','WAQTEL','DELWAQ','GAIA'],
        defaut = '',
        fr = """Liste des codes avec lesquels on couple \telemac{3d} :
\begin{itemize}
\item \sisyphe : couplage interne avec \sisyphe ;
\item \tomawac : couplage interne avec \tomawac les forces
induites par les vagues sont constantes sur la hauteur ;
\item TOMAWACT3D : couplage interne avec \tomawac les forces
induites par les vagues sont 3D ;
\item \waqtel : couplage interne avec \waqtel ;
\item DELWAQ : sortie de fichiers de resultats pour Delwaq ;
\item \gaia : couplage interne avec \gaia.
\end{itemize}""",
        ang = """List of codes to be coupled with \telemac{3d}:
\begin{itemize}
\item \sisyphe: internal coupling with \sisyphe,
\item \tomawac: internal coupling with \tomawac, forces are constant
along the depth,
\item TOMAWACT3D: internal coupling with \tomawac, forces are 3D,
\item \waqtel: internal coupling with \waqtel,
\item DELWAQ: will yield results file for DELWAQ,
\item \gaia: internal coupling with \gaia.
\end{itemize}""",
    ),
#   -----------------------------------
    b_COUPLING_WITHG = BLOC(condition="COUPLING_WITH == 'SISYPHE'",
#   -----------------------------------
#       -----------------------------------
        SISYPHE_STEERING_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres de \sisyphe en cas de couplage
interne.""",
            ang = """\sisyphe parameter file in case of internal coupling.""",
        ),
#       -----------------------------------
        COUPLING_PERIOD_FOR_SISYPHE = SIMP(statut ='o',
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
    b_COUPLING_WITHH = BLOC(condition="COUPLING_WITH in ['TOMAWAC', 'TOMAWACT3D']",
#   -----------------------------------
#       -----------------------------------
        TOMAWAC_STEERING_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres de \tomawac en cas de couplage
interne.""",
            ang = """\tomawac parameter file in case of internal coupling.""",
        ),
#       -----------------------------------
        COUPLING_PERIOD_FOR_TOMAWAC = SIMP(statut ='o',
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
    b_COUPLING_WITHI = BLOC(condition="COUPLING_WITH == 'WAQTEL'",
#   -----------------------------------
#       -----------------------------------
        WAQTEL_STEERING_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', max='**',
            defaut = '',
            fr = """Fichier des parametres physiques pour les processus de qualite d eau
(internes, pas ceux de DELWAQ).""",
            ang = """File for physical parameters of water quality processes
(local ones of \telemac{3d}-\waqtel not those of DELWAQ).""",
        ),
#       -----------------------------------
        WATER_QUALITY_PROCESS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Donne le numero du processus de qualite d''eau, defini
comme une combinaison multiplicative de nombres premiers
(2, 3, 5, 7, 11, 13, 17 et 19) avec les cas particuliers 0 et 1 :
\begin{itemize}
\item 0 : tous ;
\item 1 : rien ;
\item 2 : O2 ;
\item 3 : BIOMASS ;
\item 5 : EUTRO ;
\item 7 : MICROPOL ;
\item 11 : THERMIC ;
\item 13 : AED2 ;
\item 17 : Loi de degradation ;
\item 19 : Court-circuit temporaire pour les glaces.
\end{itemize}
Exemple : 110 = 2x5x11 activera O2, EUTRO et THERMIC ensemble.
On notera qu''AED2, pour l''instant, doit etre utilise toute seule
sans combinaison possible avec d autre processus de qualite de l eau.""",
            ang = """Gives the water quality process number, defined as
a multiplicative combination of prime numbers (2, 3, 5, 7, 11, 13, 17
and 19) with 0 and 1 having a special role:
\begin{itemize}
\item 0: all,
\item 1: none,
\item 2: O2,
\item 3: BIOMASS,
\item 5: EUTRO,
\item 7: MICROPOL,
\item 11: THERMIC,
\item 13: AED2,
\item 17: Degradation law,
\item 19: Ghost process for ice modelling.
\end{itemize}
Example: 110 = 2x5x11 activates O2, EUTRO and THERMIC together.
It is noted that AED2 should be used on its own, for the time being,
without possible combination with other processes.""",
        ),
    ),
#   -----------------------------------
    b_COUPLING_WITHJ = BLOC(condition="COUPLING_WITH == 'DELWAQ'",
#   -----------------------------------
#       -----------------------------------
        DELWAQ_STEERING_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
            defaut = '',
            fr = """Fichier de commande pour le chainage avec DELWAQ.""",
            ang = """Steering file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        DELWAQ_PRINTOUT_PERIOD = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Periode de sortie des resultats pour DELWAQ.""",
            ang = """Printout period for DELWAQ files.""",
        ),
#       -----------------------------------
        EXCHANGES_BETWEEN_NODES_DELWAQ_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        NODES_DISTANCES_DELWAQ_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        BOTTOM_SURFACES_DELWAQ_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        VOLUMES_DELWAQ_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        EXCHANGE_AREAS_DELWAQ_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        VERTICAL_FLUXES_DELWAQ_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
            defaut = '',
            fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
            ang = """Results file for chaining with DELWAQ.""",
        ),
#       -----------------------------------
        VELOCITY_FOR_DELWAQ = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie de la vitesse pour DELWAQ.""",
            ang = """Triggers the output of velocity for DELWAQ.""",
        ),
#       -----------------------------------
        b_VELOCITY_FOR_DELWAQG = BLOC(condition="VELOCITY_FOR_DELWAQ == True",
#       -----------------------------------
#           -----------------------------------
            VELOCITY_DELWAQ_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
                ang = """Results file for chaining with DELWAQ.""",
            ),
        ),
#       -----------------------------------
        DIFFUSION_FOR_DELWAQ = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie de la diffusion pour DELWAQ.""",
            ang = """Triggers the output of diffusion for DELWAQ.""",
        ),
#       -----------------------------------
        b_DIFFUSION_FOR_DELWAQG = BLOC(condition="DIFFUSION_FOR_DELWAQ == True",
#       -----------------------------------
#           -----------------------------------
            DIFFUSIVITY_DELWAQ_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
                ang = """Results file for chaining with DELWAQ.""",
            ),
        ),
#       -----------------------------------
        TEMPERATURE_FOR_DELWAQ = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie de la temperature pour DELWAQ.""",
            ang = """Triggers the output of temperature for DELWAQ.""",
        ),
#       -----------------------------------
        b_TEMPERATURE_FOR_DELWAQG = BLOC(condition="TEMPERATURE_FOR_DELWAQ == True",
#       -----------------------------------
#           -----------------------------------
            TEMPERATURE_DELWAQ_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
                ang = """Results file for chaining with DELWAQ.""",
            ),
        ),
#       -----------------------------------
        SALINITY_FOR_DELWAQ = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Decide de la sortie de la salinite pour DELWAQ.""",
            ang = """Triggers the output of salinity for DELWAQ.""",
        ),
#       -----------------------------------
        b_SALINITY_FOR_DELWAQG = BLOC(condition="SALINITY_FOR_DELWAQ == True",
#       -----------------------------------
#           -----------------------------------
            SALINITY_DELWAQ_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                defaut = '',
                fr = """Fichier de resultats pour le chainage avec DELWAQ.""",
                ang = """Results file for chaining with DELWAQ.""",
            ),
        ),
    ),
#   -----------------------------------
    SISYPHE = FACT(statut='f',
#   -----------------------------------
    ),
#   -----------------------------------
    TOMAWAC = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        BOTTOM_FRICTION_DUE_TO_WAVES = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Permet de prendre en compte les efforts dus aux vagues sur
le fond dans le cas du couplage 3D. Ceci necessite d avoir
un maillage fin au voisinage du fond pour etre suffisamment precis.""",
            ang = """Allows to take into account the momentum lost by waves due to
bottom friction. You need a fine mesh around the bottom to
be precise.""",
        ),
    ),
#   -----------------------------------
    WAQTEL = FACT(statut='f',
#   -----------------------------------
    ),
#   -----------------------------------
    DELWAQ = FACT(statut='f',
#   -----------------------------------
    ),
#   -----------------------------------
    GAIA = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        GAIA_STEERING_FILE = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Fichier des parametres de \gaia en cas de couplage
interne.""",
            ang = """\gaia parameter file in case of internal coupling.""",
        ),
    ),
)
# -----------------------------------------------------------------------
AUTOMATIC_DIFFERENTIATION = PROC(nom= "AUTOMATIC_DIFFERENTIATION",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    AD_NUMBER_OF_DERIVATIVES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 0,
        fr = """Definit le nombre de derivees utilisateurs, dans le cadre
de la differentiation algorithmique.""",
        ang = """Defines the number of user derivatives, within the framework
of the algorithmic differentiation.""",
    ),
#   -----------------------------------
    AD_NAMES_OF_DERIVATIVES = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM', min= 2, max= 2,
        fr = """Noms des derivees utilisateurs en 32 caracteres,
         16 pour le nom, 16 pour l''unite.""",
        ang = """Name of user derivatives in 32 characters,
         16 for the name, 16 for the unit.""",
    ),
#   -----------------------------------
    AD_NUMBER_OF_DIRECTIONS = SIMP(statut ='f',
#   -----------------------------------
        typ = 'I',
        defaut = 1,
        fr = """Definit le nombre de directions de differentiateurs.""",
        ang = """Defines the number of directions for the differentiators.""",
    ),
#   -----------------------------------
    AD_SYMBOLIC_LINEAR_SOLVER = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = False,
        fr = """Permet le solveur lineaire symbolique pour l AD.""",
        ang = """Enables the symbolic linear solver for AD.""",
    ),
#   -----------------------------------
    AD_LINEAR_SOLVER_RESET_DERIVATIVES = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Remet a zero les derivees pour l AD.""",
        ang = """Resets the derivatives for AD.""",
    ),
#   -----------------------------------
    AD_LINEAR_SOLVER_DERIVATIVE_CONVERGENCE = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Solveur lineaire iteratif : test de convergence des derivees
pour l AD.""",
        ang = """Iterative linear solvers: derivative convergence test for AD.""",
    ),
)
# -----------------------------------------------------------------------
INTERNAL = PROC(nom= "INTERNAL",op = None,
# -----------------------------------------------------------------------
    UIinfo = {"groupes": ("CACHE")},
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
    STEERING_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Nom du fichier contenant les parametres du calcul a realiser.""",
        ang = """Name of the file containing the parameters of the computation.
Written by the user.""",
    ),
#   -----------------------------------
    DICTIONARY = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = 'telemac3d.dico',
        fr = """Dictionnaire des mots cles.""",
        ang = """Key word dictionary.""",
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
GENERAL_PARAMETERS();\
VERTICAL();\
NUMERICAL_PARAMETERS();\
HYDRODYNAMICS();\
TURBULENCE();\
"
Ordre_Des_Commandes = (
'COMPUTATION_ENVIRONMENT',
'GENERAL_PARAMETERS',
'VERTICAL',
'NUMERICAL_PARAMETERS',
'HYDRODYNAMICS',
'TURBULENCE',
'TRACERS',
'SEDIMENT_INFO',
'COUPLING',
'AUTOMATIC_DIFFERENTIATION',
'INTERNAL')
try:
    import TelApy
    source = "eficas"
except Exception as excpt:
    source = "Telemac"
enum = source+'.telemac3d_enum_auto'
dicoCasEn = source+'.telemac3d_dicoCasEnToCata'
dicoCasFr = source+'.telemac3d_dicoCasFrToCata'
