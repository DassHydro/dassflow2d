
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



JdC = JDC_CATA (code = 'TOMAWAC',
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
    OUTPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        LISTING = FACT(statut='f',
#       -----------------------------------
#           -----------------------------------
            PERIOD_FOR_LISTING_PRINTOUTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periodicite,
en nombre de pas de temps
des messages du code sur le fichier listing""",
                ang = """Determines the period,
in number of time step
of the software messages in the listing file.""",
            ),
        ),
#       -----------------------------------
        RESULTS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            VARIABLES_FOR_2D_GRAPHIC_PRINTOUTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM', min=0, max='**',
                into = ["Total variance  (m2)","Spectral significant wave height  (m)","Mean wave direction  (deg)","Mean directional spreading  (deg)","Sea bottom level  (m)","Water depth  (m)","Current along X  (m/s)","Current along Y  (m/s)","Wind along X  (m/s)","Wind along Y  (m/s)","Driving force along X  (m/s2)","Driving force along Y  (m/s2)","Radiation stress along xx  (m3/s2)","Radiation stress along yy  (m3/s2)","Radiation stress along xy  (m3/s2)","Bottom celerity  (m/s)","Wave power (per meter along wave crest)  (kW/m)","Mean frequency FMOY  (Hz)","Mean frequency FM01  (Hz)","Mean frequency FM02  (Hz)","Discrete peak frequency  (Hz)","Peak frequency by Read method of order 5  (Hz)","Peak frequency by Read method of order 8  (Hz)","Surface friction velocity u*  (m/s)","Surface drag coefficient CD  (-)","Surface roughness length Z0  (m)","Surface wave stress  (kg/(m.s2))","Mean period Tmoy  (s)","Mean period Tm01  (s)","Mean period Tm02  (s)","Discrete peak period  (s)","Peak period by Read method of order 5  (s)","Peak period by Read method of order 8  (s)","Private table  (?)","Breaking waves coefficient  (-)"],
                defaut = ["Spectral significant wave height  (m)","Mean wave direction  (deg)"],
                fr = """Noms des variables que l''utilisateur veut ecrire dans
le FICHIER DES RESULTATS 2D. Les variables disponibles sont :
M0   : variance totale
HM0  : Hauteur significative spectrale
DMOY : Direction moyenne de houle
SPD  : Etalement directionnel moyen
ZF   : Cote du fond
WD   : Hauteur d''eau
UX   : Courant suivant X
UY   : Courant suivant Y
VX   : Vent suivant X
VY   : Vent suivant Y
FX   : Force motrice suivant X
FY   : Force motrice suivant Y
SXX  : Contrainte de radiation suivant xx
SYY  : Contrainte de radiation suivant yy
SXY  : Contrainte de radiation suivant xy
UWB  : Vitesse orbitale au fond
POW  : Puissance lineique de houle (par metre de crete)
FMOY : Frequence moyenne FMOY
FM01 : Frequence moyenne FM01
FM02 : Frequence moyenne FM02
FPD  : Frequence de pic discrete
FPR5 : Frequence de pic de Read ordre 5
FPR8 : Frequence de pic de Read ordre 8
US   : Vitesse de frottement en surface u*
CD   : Coefficient de trainee en surface CD
Z0   : Longueur de rugosite en surface Z0
WS   : Contrainte de houle en surface
TMOY : Periode moyenne Tmoy
TM01 : Periode moyenne Tm01
TM02 : Periode moyenne Tm02
TPD  : Periode de pic discrete
TPR5 : Periode de pic de Read ordre 5
TPR8 : Periode de pic de Read ordre 8
PRI  : tableau prive
BETA : coefficient de deferlement
**Mots-cles associes :**
FICHIER DES RESULTATS 2D
NUMERO DE LA PREMIERE ITERATION POUR LES SORTIES GRAPHIQUES
PERIODE POUR LES SORTIES GRAPHIQUES""",
                ang = """Codes of the variables the user wants to write into the 2D
RESULTS FILE. The available variables are as follows
\begin{itemize}
\item M0   : Total variance
\item HM0  : Spectral significant wave height
\item DMOY : Mean wave direction
\item SPD  : Mean directional spreading
\item ZF   : Sea bottom level
\item WD   : Water depth
\item UX   : Current along X
\item UY   : Current along Y
\item VX   : Wind along X
\item VY   : Wind along Y
\item FX   : Driving force along X
\item FY   : Driving force along Y
\item SXX  : Radiation stress along xx
\item SYY  : Radiation stress along yy
\item SXY  : Radiation stress along xy
\item UWB  : Bottom celerity
\item POW  : Wave power (per meter along wave crest)
\item FMOY : Mean frequency FMOY
\item FM01 : Mean frequency FM01
\item FM02 : Mean frequency FM02
\item FPD  : Discrete peak frequency
\item FPR5 : Peak frequency by Read method of order 5
\item FPR8 : Peak frequency by Read method of order 8
\item US   : Surface friction velocity u*
\item CD   : Surface drag coefficient CD
\item Z0   : Surface roughness length Z0
\item WS   : Surface wave stress
\item TMOY : Mean period Tmoy
\item TM01 : Mean period Tm01
\item TM02 : Mean period Tm02
\item TPD  : Discrete peak period
\item TPR5 : Peak period by Read method of order 5
\item TPR8 : Peak period by Read method of order 8
\item PRI  : Private table
\item BETA : Breaking waves coefficient
\end{itemize}
 \begin{CommentBlock}{Related keywords}
2D RESULTS FILE\\
NUMBER OF FIRST ITERATION FOR GRAPHIC PRINTOUTS\\
PERIOD FOR GRAPHIC PRINTOUTS
\end{CommentBlock}""",
            ),
#           -----------------------------------
            PERIOD_FOR_GRAPHIC_PRINTOUTS = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Determine la periode d''impression,
en nombre de pas de temps
des VARIABLES POUR LES SORTIES GRAPHIQUES dans le FICHIER DES
RESULTATS 2D et le FICHIER DES RESULTATS PONCTUELS.
**Mots-cles associes :**
VARIABLES POUR LES SORTIES GRAPHIQUES 2D
ABSCISSES DES POINTS DE SORTIE DU SPECTRE
ORDONNEES DES POINTS DE SORTIE DU SPECTRE
FICHIER DES RESULTATS 2D
FICHIER DES RESULTATS PONCTUELS
NUMERO DE LA PREMIERE ITERATION POUR LES SORTIES GRAPHIQUES""",
                ang = """Determines the printing period,
in number of time step
of the VARIABLES FOR 2D GRAPHIC PRINTOUTS in the 2D
RESULTS FILE and the PUNCTUAL RESULTS FILE.
\\
 \begin{CommentBlock}{Related keywords}
VARIABLES FOR 2D GRAPHIC PRINTOUTS\\
ABSCISSAE OF SPECTRUM PRINTOUT POINTS\\
ORDINATES OF SPECTRUM PRINTOUT POINTS\\
2D RESULTS FILE\\
PUNCTUAL RESULTS FILE\\
NUMBER OF FIRST ITERATION FOR GRAPHIC PRINTOUTS\\
\end{CommentBlock}""",
            ),
#           -----------------------------------
            NUMBER_OF_FIRST_ITERATION_FOR_GRAPHICS_PRINTOUTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Determine le nombre d''iterations sur le pas de temps
ou sur la frequence moyenne a partir duquel debute
l ecriture des resultats dans le FICHIER DES RESULTATS 2D et le
FICHIER DES RESULTATS PONCTUELS
**Mots-cles associes :**
PERIODE POUR LES SORTIES GRAPHIQUES
VARIABLES POUR LES SORTIES GRAPHIQUES 2D
ABSCISSES DES POINTS DE SORTIE DU SPECTRE
ORDONNEES DES POINTS DE SORTIE DU SPECTRE
FICHIER DES RESULTATS 2D
FICHIER DES RESULTATS PONCTUELS""",
                ang = """Determines the number of iterations over mean angular
frequency from which the results are first written into the 2D RESULTS
FILE and the PUNCTUAL RESULTS FILE.
\\
 \begin{CommentBlock}{Related keywords}
PERIOD FOR GRAPHIC PRINTOUTS\\
VARIABLES FOR 2D GRAPHIC PRINTOUTS\\
ABSCISSAE OF SPECTRUM PRINTOUT POINTS\\
ORDINATES OF SPECTRUM PRINTOUT POINTS\\
2D RESULTS FILE\\
PUNCTUAL RESULTS FILE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            ED_RESULTS_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les resultats du calcul
bidimensionnels.
Mots-cles associes
VARIABLES POUR LES SORTIES GRAPHIQUES 2D
PERIODE OUR LES SORTIES GRAPHIQUES
NUMERO DE LA PREMIERE ITERATION POUR LES SORTIES GRAPHIQUES""",
                ang = """Name of the file into which the results of the two-dimensional
computation will be written.
\\
 \begin{CommentBlock}{Related keywords}
VARIABLES FOR 2D GRAPHIC PRINTOUTS\\
PERIOD FOR GRAPHIC PRINTOUTS\\
NUMBER OF FIRST ITERATION FOR GRAPHIC PRINTOUTS
\end{CommentBlock}""",
            ),
#           -----------------------------------
            ED_RESULTS_FILE_FORMAT = SIMP(statut ='o',
#           -----------------------------------
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
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            PUNCTUAL_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les spectres
ponctuels.
**Mots-cles associes :**
FORMAT DU FICHIER DE SPECTRE
ABSCISSES DES POINTS DE SORTIE DU SPECTRE
ORDONNEES DES POINTS DE SORTIE DU SPECTRE
PERIODE POUR LES SORTIES GRAPHIQUES
NUMERO DE LA PREMIERE ITERATION POUR LES SORTIES GRAPHIQUES""",
                ang = """Name of the file into which the punctual spectra will be
written.
\\
 \begin{CommentBlock}{Related keywords}
SPECTRUM FILE FORMAT\\
ABSCISSAE OF SPECTRUM PRINTOUT POINTS\\
ORDINATES OF SPECTRUM PRINTOUT POINTS\\
PERIOD FOR GRAPHIC PRINTOUTS\\
NUMBER OF FIRST ITERATION FOR GRAPHIC PRINTOUTS
\end{CommentBlock}""",
            ),
#           -----------------------------------
            SPECTRUM_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier de spectre
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                ang = """Spectrum results file format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            ABSCISSAE_OF_SPECTRUM_PRINTOUT_POINTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Tableau donnant les abscisses des points de sortie Seraphin du
spectre et de dimension maximale 99. Les points de sortie du spectre
sont les points 2D les plus proches des coordonnees specifiees
**Mots-cles associes :**
ORDONNEES DES POINTS DE SORTIE DU SPECTRE
FICHIER DES RESULTATS PONCTUELS""",
                ang = """Array providing the abscissae of the Seraphin spectrum printout
points with a maximum dimension of 99. The chosen spectrum points are
the closest 2D points to the specified co-ordinates.
\\
 \begin{CommentBlock}{Related keywords}
ORDINATES OF SPECTRUM PRINTOUT POINTS\\
PUNCTUAL RESULTS FILE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            ORDINATES_OF_SPECTRUM_PRINTOUT_POINTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min=0, max='**',
                fr = """Tableau donnant les ordonnees des points de sortie Seraphin du
spectre et de dimension max 99.Les points de sortie du spectre
sont les points 2D les plus proches des coordonnees specifiees
**Mots-cles associes :**
ABSCISSES DES POINTS DE SORTIE DU SPECTRE
FICHIER DES RESULTATS PONCTUELS""",
                ang = """Array providing the ordinates of the Seraphin spectrum printout
points with a maximum dimension of 99. The spectrum printout points are
the closest 2D points to the specified co-ordinates
\\
 \begin{CommentBlock}{Related keywords}
ABSCISSAE OF SPECTRUM PRINTOUT POINTS\\
PUNCTUAL RESULT FILE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            ZD_SPECTRA_RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel seront ecrits les spectres
frequentiels ponctuels (integres selon les directions).
**Mots-cles associes :**
FICHIER DES RESULTATS PONCTUELS
ABSCISSES DES POINTS DE SORTIE DU SPECTRE
ORDONNEES DES POINTS DE SORTIE DU SPECTRE
PERIODE POUR LES SORTIES GRAPHIQUES
NUMERO DE LA PREMIERE ITERATION POUR LES SORTIES GRAPHIQUES""",
                ang = """Name of the file into which the frequential punctual
spectra (integrated according to the directions) will be written.
\\
 \begin{CommentBlock}{Related keywords}
PUNCTUAL RESULTS FILE\\
ABSCISSAE OF SPECTRUM PRINTOUT POINTS\\
ORDINATES OF SPECTRUM PRINTOUT POINTS\\
PERIOD FOR GRAPHIC PRINTOUTS\\
NUMBER OF FIRST ITERATION FOR GRAPHIC PRINTOUTS
\end{CommentBlock}""",
            ),
#           -----------------------------------
            GLOBAL_RESULT_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'),
                defaut = '',
                fr = """Nom du fichier dans lequel sera ecrit le tableau F (spectre
de variance) en fin de calcul dans le but de faire une suite de calcul.
**Mots-cles associes :**
FORMAT DU FICHIER DES RESULTATS GLOBAUX""",
                ang = """Name of the file in which the table F (density spectrum)
is written at the end of the computation in order to realise a next
computation.
\\
 \begin{CommentBlock}{Related keywords}
GLOBAL RESULT FILE FORMAT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            GLOBAL_RESULT_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
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
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            FILE_WITH_COORDINATES_OF_SPECTRA_TO_WRITE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier texte au format DAT de Salome avec les coordonnées
des spectres a ecrire.""",
                ang = """Text file following the DAT format of Salome with the
coordinates of the spectra to output.""",
            ),
        ),
    ),
#   -----------------------------------
    INPUT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NAMES_OF_VARIABLES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM', min= 5, max= 5,
            defaut = 'VELOCITY U      M/S;VELOCITY V      M/S;WIND ALONG X    M/S;WIND ALONG Y    M/S;WATER DEPTH     M',
            fr = """Nom des variables dans les fichiers au format SERAFIN
        1: Vitesse U
        2: Vitesse V
        3: Vitesse du vent suivant X
        4: Vitesse du vent suivant Y
        5: Profondeur""",
            ang = """Names of variables in SERAFIN format files
\begin{itemize}
       \item 1: Velocity U
       \item 2: Velocity V
       \item 3: Wind velocity along X
       \item 4: Wind velocity along Y
       \item 5: Depth
\end{itemize}""",
        ),
#       -----------------------------------
        CURRENT_FILE = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            CURRENTS_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["Selafin, TELEMAC type","User format (couuti.f)"],
                defaut = ["Selafin, TELEMAC type"],
                fr = """Choix du type de format du fichier des courants :
\begin{itemize}
\item 3 = selafin du type TELEMAC
\item 4 = format utilisateur (Modifier alors  la procedure couuti.f)
\end{itemize}
**Mots-cles associes :**
FICHIER DES COURANTS BINAIRE
FICHIER DES COURANTS FORMATE""",
                ang = """Selection of the type of currents file format :
\begin{itemize}
\item 3 = selafin, TELEMAC type
\item 4 = user format (the couuti.f procedure should then be amended)
\end{itemize}
\begin{CommentBlock}{Related keywords}
CURRENTS BINARY FILE\\
CURRENTS FORMATTED FILE\\
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TIME_INCREMENT_NUMBER_IN_TELEMAC_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Indique le numero du pas de temps du fichier de resultats
TELEMAC (fichier des courants) correspondant a l instant desire
pour recuperer la donnee.
**Mots-cles associes :**
RANG DE LA DONNEE TELEMAC A RECUPERER
RECUPERATION DE DONNEES TELEMAC""",
                ang = """Indicates the number of the time increment in the TELEMAC
results file (currents file) corresponding to the desired time for data
recovery.
\\
 \begin{CommentBlock}{Related keywords}
RANK OF THE TELEMAC DATA ITEM TO BE RECOVERED\\
RECOVERY OF TELEMAC DATA ITEM
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TIME_UNIT_IN_CURRENTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Unite donnee en secondes, par exemple 3600. si le temps
est donne en heures""",
                ang = """Unit given in seconds, for example 3600. if time
is given in hours""",
            ),
#           -----------------------------------
            TIME_SHIFT_IN_CURRENTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Sera retranche au temps lu dans le fichier.
L''unite est celle du fichier""",
                ang = """Will be withdrawn from the time read in the file.
 The unit is that of the file""",
            ),
#           -----------------------------------
            BINARY_CURRENTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de donnees de courant (si binaire).
**Mots-cles associes :**
PRISE EN COMPTE D''UN COURANT STATIONNAIRE
PRISE EN COMPTE DE LA MAREE
FICHIER DES COURANTS FORMATE
FORMAT DU FICHIER DES COURANTS""",
                ang = """Name of the current data file (if binary).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF A STATIONARY CURRENT\\
CONSIDERATION OF TIDE\\
FORMATTED CURRENTS FILE\\
CURRENTS FILE FORMAT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            FORMATTED_CURRENTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de donnees de courant (si formate).
**Mots-cles associes :**
PRISE EN COMPTE D''UN COURANT STATIONNAIRE
PRISE EN COMPTE DE LA MAREE
FICHIER DES COURANTS BINAIRE
FORMAT DU FICHIER DES COURANTS""",
                ang = """Name of the current data file (if formatted).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF A STATIONARY CURRENT\\
CONSIDERATION OF TIDE\\
BINARY CURRENTS FILE\\
CURRENTS FILE FORMAT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            BINARY_CURRENTS_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier binaire des courants.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                ang = """Currents binary file format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
        ),
#       -----------------------------------
        WIND_FILE = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            WINDS_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["Selafin, TELEMAC type","User format (venuti.f)"],
                defaut = ["Selafin, TELEMAC type"],
                fr = """Choix du type de format du fichier des vents :
3 = selafin du type TELEMAC
4 = format utilisateur (Modifier alors  la procedure venuti.f)
**Mots-cles associes :**
FICHIER DES VENTS BINAIRE
FICHIER DES VENTS FORMATE""",
                ang = """Selection of winds file format type :
\begin{itemize}
\item 3 = selafin, TELEMAC type
\item 4 = user format (the venuti.f procedure should then be amended)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
WINDS FILE TYPE\\
WINDS FILE\\
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TIME_UNIT_IN_WINDS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Unite donnee en secondes, par exemple 3600. si le temps
est donne en heures""",
                ang = """Unit given in seconds, for example 3600. if time
is given in hours""",
            ),
#           -----------------------------------
            TIME_SHIFT_IN_WINDS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Sera retranche au temps lu dans le fichier.
L''unite est celle du fichier""",
                ang = """Will be withdrawn from the time read in the file.
 The unit is that of the file""",
            ),
#           -----------------------------------
            BINARY_WINDS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de donnees de vent (si binaire).
**Mots-cles associes :**
PRISE EN COMPTE DU VENT
FICHIER DES VENTS FORMATE
FORMAT DU FICHIER DES VENTS""",
                ang = """Name of wind data file (if binary).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF WIND\\
FORMATTED WINDS FILE\\
WINDS FILE FORMAT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            BINARY_WINDS_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier binaire des vents.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                ang = """wind data binary file format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            FORMATTED_WINDS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de donnees de vent (si formate).
**Mots-cles associes :**
PRISE EN COMPTE DU VENT
FICHIER DES VENTS BINAIRE
FORMAT DU FICHIER DES VENTS""",
                ang = """Name of wind data file (if formatted).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF WIND\\
BINARY WINDS FILE\\
WINDS FILE FORMAT
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        TIDAL_FILE = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            TIDAL_WATER_LEVEL_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ["Selafin, TELEMAC type","User format (maruti.f)"],
                defaut = ["Selafin, TELEMAC type"],
                fr = """Choix du type de format du fichier du niveau de la maree :
3 = selafin du type TELEMAC
4 = format utilisateur (Modifier alors  la procedure maruti.f)
**Mots-cles associes :**
PRISE EN COMPTE DE LA MAREE
FICHIER DU NIVEAU DE LA MAREE BINAIRE
FICHIER DU NIVEAU DE LA MAREE FORMATE
PERIODE D ACTUALISATION DE LA MAREE""",
                ang = """Selection of the type of tidal water level file format :
\begin{itemize}
\item 3 = selafin, TELEMAC type
\item 4 = user format (the maruti.f procedure should then be amended)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF TIDE\\
BINARY TIDAL WATER LEVEL FILE\\
FORMATTED TIDAL WATER LEVEL FILE\\
TIDE REFRESHING PERIOD
\end{CommentBlock}""",
            ),
#           -----------------------------------
            RANK_OF_THE_WATER_LEVEL_DATA_IN_THE_TELEMAC_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 4,
                fr = """Rang de la variable donnant le niveau de la maree
dans le fichier TELEMAC
**Mots-cles associes :**
PRISE EN COMPTE DE LA MAREE
FICHIER DU NIVEAU DE LA MAREE BINAIRE
FICHIER DU NIVEAU DE LA MAREE FORMATE
PERIODE D ACTUALISATION DE LA MAREE""",
                ang = """Rank of the water level data in the TELEMAC file
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF TIDE\\
BINARY TIDAL WATER LEVEL FILE\\
FORMATTED TIDAL WATER LEVEL FILE\\
TIDE REFRESHING PERIOD
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TIME_UNIT_IN_TIDAL_WATER_LEVEL_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Unite donnee en secondes, par exemple 3600. si le temps
est donne en heures""",
                ang = """Unit given in seconds, for example 3600. if time
is given in hours""",
            ),
#           -----------------------------------
            TIME_SHIFT_IN_TIDAL_WATER_LEVEL_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Sera retranche au temps lu dans le fichier.
L''unite est celle du fichier""",
                ang = """Will be withdrawn from the time read in the file.
 The unit is that of the file""",
            ),
#           -----------------------------------
            BINARY_TIDAL_WATER_LEVEL_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de donnees du niveau d''eau (si binaire).
**Mots-cles associes :**
PRISE EN COMPTE DE LA MAREE
FICHIER DU NIVEAU DE LA MAREE FORMATE
FORMAT DU FICHIER DU NIVEAU DE LA MAREE
PERIODE D ACTUALISATION DE LA MAREE""",
                ang = """Name of the water level data file (if binary).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF TIDE\\
FORMATTED TIDAL WATER LEVEL FILE\\
TIDAL WATER LEVEL FILE FORMAT\\
TIDE REFRESHING PERIOD
\end{CommentBlock}""",
            ),
#           -----------------------------------
            BINARY_TIDAL_WATER_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier de la maree binaire.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                ang = """binary tidal water file format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            FORMATTED_TIDAL_WATER_LEVEL_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de donnees du niveau d''eau (si formate).
**Mots-cles associes :**
PRISE EN COMPTE D UN COURANT
FICHIER DES COURANTS BINAIRE
FORMAT DU FICHIER DES COURANTS
PERIODE D ACTUALISATION DE LA MAREE""",
                ang = """Name of the tidal data file (if formatted).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF TIDE\\
BINARY TIDAL WATER LEVEL FILE\\
TIDAL WATER LEVEL FILE FORMAT\\
TIDE REFRESHING PERIOD
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        DATA = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            GEOMETRY_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant le maillage du calcul a realiser.
**Mots-cles associes :**
FORMAT DU FICHIER DE GEOMETRIE""",
                ang = """Name of the file containing the mesh of the computation to be
made.
\\
 \begin{CommentBlock}{Related keywords}
GEOMETRY FILE FORMAT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            GEOMETRY_FILE_FORMAT = SIMP(statut ='o',
#           -----------------------------------
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
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            FORTRAN_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'FichierOuRepertoire',
                defaut = '',
                fr = """Nom du fichier FORTRAN a soumettre.""",
                ang = """Name of FORTRAN file to be submitted.""",
            ),
#           -----------------------------------
            BOUNDARY_CONDITIONS_FILE = SIMP(statut ='o',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les types de conditions aux limites.
Ce fichier est rempli de facon automatique par le mailleur au moyen de
couleurs affectees aux noeuds des frontieres du domaine de calcul.""",
                ang = """Name of the file containing the types of boundary conditions.
This file is automatically filled by the grid generator by means of
colours that are assigned to the boundary nodes in the computational
domain.""",
            ),
#           -----------------------------------
            BOTTOM_TOPOGRAPHY_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier eventuel contenant la bathymetrie associee au
maillage au format SINUSX.
Si ce mot-cle est utilise; c''est cette bathymetrie qui sera utilisee
pour le calcul.""",
                ang = """Name of any file containing the bathymetric data associated to
the SINUSX-formatted grid. It this keyword is used, these bathymetric
data shall be used for the computation.""",
            ),
#           -----------------------------------
            BOTTOM_SMOOTHINGS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = 0,
                fr = """Nombre de lissages effectues sur la topographie.
Chaque lissage, effectue a l''aide d''une matrice de masse,
est conservatif. A utiliser lorsque les donnees de bathymetrie
donnent des resultats trop irreguliers apres interpolation.
Voir aussi le sous-programme CORFON""",
                ang = """Number of smoothings made on bottom features. Each smoothing,
being made by means of a mass matrix, is conservative. To be used when
the bathymetric data yield too irregular data after interpolation.
Also refer to the CORFON subroutine.""",
            ),
#           -----------------------------------
            NEXT_COMPUTATION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Indique si on fait une suite de calcul.
**Mots-cles associes :**
FICHIER DU CALCUL PRECEDENT""",
                ang = """Indicates whether a next compution is done.
\\
 \begin{CommentBlock}{Related keywords}
PREVIOUS COMPUTATION FILE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            b_NEXT_COMPUTATIONG = BLOC(condition="NEXT_COMPUTATION == True",
#           -----------------------------------
#               -----------------------------------
                PREVIOUS_COMPUTATION_FILE = SIMP(statut ='o',
#               -----------------------------------
                    typ = ('Fichier','All Files (*)'),
                    defaut = '',
                    fr = """Nom d''un fichier contenant les resultats d''un calcul precedent
realise sur le meme maillage et qui va fournir les conditions
initiales pour une suite de calcul.
**Mots-cles associes :**
BINAIRE DU FICHIER DU CALCUL PRECEDENT""",
                    ang = """Name of the file containing the global results of a previous
computation realised with the same mesh. This file gives the initial
conditions for a next computation.
\\
 \begin{CommentBlock}{Related keywords}
BINARY OF THE PREVIOUS COMPUTATION FILE
\end{CommentBlock}""",
                ),
#               -----------------------------------
                PREVIOUS_COMPUTATION_FILE_FORMAT = SIMP(statut ='o',
#               -----------------------------------
                    typ = 'TXM',
                    into = ['SERAFIN','SERAFIND','MED'],
                    defaut = 'SERAFIN',
                    fr = """Format du fichier de resultats du calcul precedent.
Les valeurs possibles sont seulement:
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                    ang = """Previous computation results file format.
Possible values are only:
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
                ),
            ),
#           -----------------------------------
            REFERENCE_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier de reference en cas de validation.
**Mots-cles associes :**
VALIDATION""",
                ang = """Name of validation data file
\\
 \begin{CommentBlock}{Related keywords}
VALIDATION
\end{CommentBlock}""",
            ),
#           -----------------------------------
            REFERENCE_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
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
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            BINARY_FILE_1 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees code en binaire mis a la disposition
de l''utilisateur.""",
                ang = """Binary-coded data file made available to the user.""",
            ),
#           -----------------------------------
            BINARY_DATA_FILE_1_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['','SERAFIN','SERAFIND','MED'],
                defaut = '',
                fr = """Format du fichier des donnees binaires.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                ang = """binary data file format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}""",
            ),
#           -----------------------------------
            FORMATTED_FILE_1 = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier de donnees formate mis a la disposition
de l''utilisateur.""",
                ang = """Formatted data file made available to the user.""",
            ),
        ),
#       -----------------------------------
        BOUNDARY_CONDITION_FILE = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            FILE_WITH_COORDINATES_OF_SPECTRA_TO_IMPOSE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Fichier texte au format DAT de Salome avec les coordonnées
des spectres qui vont etre impose en conditions limites.""",
                ang = """Text file following the DAT format of Salome with the
coordinates of the spectra that will be imposed on the boundary.
\\
 \begin{CommentBlock}{Related keywords}
IMPOSED SPECTRA FILE
IMPOSED SPECTRA FILE FORMAT
TIME UNIT OF IMPOSED SPECTRA FILE
TIME SHIFT OF IMPOSED SPECTRA FILE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            IMPOSED_SPECTRA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant le maillage des spectres a imposer.""",
                ang = """Name of the file containing the mesh with the imposed spectra.
\\
 \begin{CommentBlock}{Related keywords}
IMPOSED SPECTRA FILE FORMAT
TIME UNIT OF IMPOSED SPECTRA FILE
TIME SHIFT OF IMPOSED SPECTRA FILE
FILE WITH COORDINATES OF SPECTRA TO IMPOSE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            IMPOSED_SPECTRA_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier des spectres imposes.
Les valeurs possibles sont :
- SERAFIN : format standard simple precision pour Telemac;
- SERAFIND: format standard double precision pour Telemac;
- MED     : format MED base sur HDF5""",
                ang = """Imposed spectra file format.
Possible values are:
\begin{itemize}
\item SERAFIN : classical single precision format in Telemac;
\item SERAFIND: classical double precision format in Telemac;
\item MED     : MED format based on HDF5
\end{itemize}

\begin{CommentBlock}{Related keywords}
IMPOSED SPECTRA FILE
TIME UNIT OF IMPOSED SPECTRA FILE
TIME SHIFT OF IMPOSED SPECTRA FILE
FILE WITH COORDINATES OF SPECTRA TO IMPOSE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TIME_UNIT_OF_IMPOSED_SPECTRA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [1.],
                fr = """Unite donnee en secondes, par exemple 3600. si le temps
est donne en heures.""",
                ang = """Unit given in seconds, for example 3600. if time
is given in hours.
\\
 \begin{CommentBlock}{Related keywords}
IMPOSED SPECTRA FILE
IMPOSED SPECTRA FILE FORMAT
TIME SHIFT OF IMPOSED SPECTRA FILE
FILE WITH COORDINATES OF SPECTRA TO IMPOSE
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TIME_SHIFT_OF_IMPOSED_SPECTRA_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Sera retranche au temps lu dans le fichier.
L''unite est celle du fichier.""",
                ang = """Will be withdrawn from the time read in the file.
 The unit is that of the file.
\\
 \begin{CommentBlock}{Related keywords}
IMPOSED SPECTRA FILE
IMPOSED SPECTRA FILE FORMAT
TIME UNIT OF IMPOSED SPECTRA FILE
FILE WITH COORDINATES OF SPECTRA TO IMPOSE
\end{CommentBlock}""",
            ),
        ),
    ),
)
# -----------------------------------------------------------------------
GENERAL_PARAMETERS = PROC(nom= "GENERAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    TIME = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_TIME_STEP = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            fr = """Definit le nombre de pas de temps effectues lors de
l''execution du code.
**Mots-cles associes :**
PAS DE TEMPS""",
            ang = """Define the number of time step.
\\
 \begin{CommentBlock}{Related keywords}
TIME STEP
\end{CommentBlock}""",
        ),
#       -----------------------------------
        TIME_STEP = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            fr = """Definit le pas de temps en secondes.
**Mots-cles associes :**
NOMBRE DE PAS DE TEMPS""",
            ang = """Define the time step in seconds.
\\
 \begin{CommentBlock}{Related keywords}
NUMBER OF TIME STEPS
\end{CommentBlock}""",
        ),
#       -----------------------------------
        DATE_OF_COMPUTATION_BEGINNING = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0,
            fr = """Donne la date du debut du calcul au format aaaammjjhhmm
ce qui permet de se reperer par rapport au fichier des vents,
199310241524 represente le 24 octobre 93 a 15h24.
**Mots-cles associes :**
FICHIER DES VENTS BINAIRE
FICHIER DES VENTS FORMATE
FORMAT DU FICHIER DES VENTS""",
            ang = """Gives the date of the computation beginning. The format
is yyyymmddhhmm,
as an exemple 199310241524 means the 24 october 93
at 15h24. This date gives a reference for reading the
wind file.
\\
 \begin{CommentBlock}{Related keywords}
BINARY WIND FILE\\
FORMATTED WIND FILE\\
WIND FILE FORMAT
\end{CommentBlock}""",
        ),
#       -----------------------------------
        INITIAL_TIME_SET_TO_ZERO = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Remet le temps a zero en cas de suite de calcul""",
            ang = """Initial time set to zero in case of restart""",
        ),
    ),
#   -----------------------------------
    SPECTRAL_DISCRETISATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_DIRECTIONS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 12,
            fr = """Definit le nombre de directions de propagation de la houle.
Les directions de propagation sont regulierement reparties entre 0
et 360 degres.""",
            ang = """Defines the number of wave propagation directions. The
propagation directions are evenly distributed from 0 to 360 degrees.""",
        ),
#       -----------------------------------
        NUMBER_OF_FREQUENCIES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 15,
            fr = """Definit le nombre de frequences de propagation de la houle.
Les frequences sont reparties geometriquement en fonction de la
FREQUENCE MINIMALE et RAISON FREQUENTIELLE.
**Mots-cles associes :**
FREQUENCE MINIMALE
RAISON FREQUENTIELLE
FACTEUR DE QUEUE DU SPECTRE""",
            ang = """Defines the number of wave propagation frequencies. The
propagation frequencies are geometrically distributed as a fonction
of the MINIMAL FREQUENCY OF THE COMPUTATION and the FREQUENTIAL
REASON
\\
 \begin{CommentBlock}{Related keywords}
FREQUENTIAL RATIO\\
SPECTRUM TAIL FACTOR
\end{CommentBlock}""",
        ),
#       -----------------------------------
        MINIMAL_FREQUENCY = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            fr = """Definit la frequence minimale en Hz.
On obtient les autres frequences grace a la RAISON FREQUENTIELLE r et le
NOMBRE DE FREQUENCES NF par $f=f0*r^{k-1}$ k=1,NF.
**Mots-cles associes :**
RAISON FREQUENTIELLE
NOMBRE DE FREQUENCES
FACTEUR DE QUEUE DU SPECTRE""",
            ang = """Define the minimal frequency in Hz. The discretised frequencies
are computed from the FREQUENTIAL RATIO r and the NUMBER OF FREQUENCIES
NF by the relation $f=f_0*r^{k-1}$ k=1,NF.
\\
 \begin{CommentBlock}{Related keywords}
FREQUENTIAL RATIO\\
NUMBER OF FREQUENCIES\\
SPECTRUM TAIL FACTOR
\end{CommentBlock}""",
        ),
#       -----------------------------------
        FREQUENTIAL_RATIO = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.1,
            fr = """Definit le rapport entre 2 frequences de discretisation
successives. On obtient les autres frequences grace a la FREQUENCE
MINIMALE f0 et le NOMBRE DE FREQUENCES NF par $f=f0*r^{k-1}$ k=1,NF.
**Mots-cles associes :**
FREQUENCE MINIMALE
NOMBRE DE FREQUENCES
FACTEUR DE QUEUE DU SPECTRE""",
            ang = """Define the ratio between 2 successive discretised
frequencies
\\
 \begin{CommentBlock}{Related keywords}
MINIMAL FREQUENCY\\
NUMBER OF FREQUENCIES\\
SPECTRUM TAIL FACTOR
\end{CommentBlock}""",
        ),
#       -----------------------------------
        SPECTRUM_TAIL_FACTOR = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 5.,
            fr = """Utilise pour calculer l energie des frequences superieures a la
**Mots-cles associes :**
NOMBRE DE FREQUENCES
FREQUENCE MINIMALE
RAISON FREQUENTIELLE""",
            ang = """decay order of the hight frequencies (Beyond cut off frequency)
\\
 \begin{CommentBlock}{Related keywords}
NUMBER OF FREQUENCIES\\
FREQUENTIAL RATIO
\end{CommentBlock}""",
        ),
#       -----------------------------------
        SPECTRUM_ENERGY_THRESHOLD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-30,
            fr = """ En condition initiale, une energie inferieure a ce seuil
pour un couple frequence-direction donnee est prise nulle. Utile
surtout pour les comparaisons avec WAM cycle 4.""",
            ang = """For initial conditions, the energy on a frequency-direction
component lower to this threshold is taken to 0.
Useful for comparisons with WAM cycle 4.""",
        ),
#       -----------------------------------
        OPTION_FOR_DIAGNOSTIC_TAIL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """Permet de choisir la methode de correction de la
partie diagnostique du spectre.
\begin{itemize}
\item 0 : Pas de queue diagnostique
\item 1 : Une decroissance en $f^{-TAILF}$ est imposee au dessus de
max(4fPM;2.5fmoy)
\end{itemize}
**Mots-cles associes :**
FACTEUR DE QUEUE DU SPECTRE
NOMBRE DE FREQUENCES
RAISON FREQUENTIELLE""",
            ang = """Option to treat the spectrum diagnotic tail.
\begin{itemize}
\item 0 : No diagnostic tail
\item 1 : A decrease in $f^{-TAILF}$ is imposed beyond
max(4fPM;2.5fmoy)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
SPECTRUM TAIL FACTOR\\
NUMBER OF FREQUENCIES\\
FREQUENTIAL RATIO
\end{CommentBlock}""",
        ),
    ),
#   -----------------------------------
    METEO = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TIDE_REFRESHING_PERIOD = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Determine la periode, en nombre d''iterations,
d''actualisation de la profondeur d''eau et des courants de maree.
**Mots-cles associes :**
PRISE EN COMPTE DE LA MAREE
FICHIER DU NIVEAU DE LA MAREE BINAIRE
FICHIER DU NIVEAU DE LA MAREE FORMATE
FORMAT DU FICHIER DU NIVEAU DE LA MAREE""",
            ang = """Determines the period in number of iterations to
update the tidal currents and the water depth.
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF TIDE\\
BINARY TIDAL WATER LEVEL FILE\\
FORMATTED TIDAL WATER LEVEL FILE\\
FORMAT DU FICHIER DU NIVEAU DE LA MAREE\\
\end{CommentBlock}""",
        ),
#       -----------------------------------
        WIND_VELOCITY_ALONG_X = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Vitesse du vent suivant X, constante et homogene (en m/s)
**Mots-cles associes :**
PRISE EN COMPTE DU VENT""",
            ang = """Wind velocity along X axis, constant and homogeneous (m/s)
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF A WIND
\end{CommentBlock}""",
        ),
#       -----------------------------------
        WIND_VELOCITY_ALONG_Y = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.,
            fr = """Vitesse du vent suivant Y, constante et homogene (en m/s)
**Mots-cles associes :**
PRISE EN COMPTE DU VENT""",
            ang = """Wind velocity along Y axis, constant and homogeneous (m/s)
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF A WIND
\end{CommentBlock}""",
        ),
#       -----------------------------------
        CONSIDERATION_OF_A_STATIONARY_CURRENT = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Indique si on prend en compte un courant, dans un fichier ou
dans condiw.f.
**Mots-cles associes :**
FICHIER DES COURANTS""",
            ang = """Indicates whether a stationary current is taken into account,
either in a file or in condiw.f.
\\
 \begin{CommentBlock}{Related keywords}
CURRENTS FILE
\end{CommentBlock}""",
        ),
#       -----------------------------------
        CONSIDERATION_OF_A_WIND = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Indique si on prend en compte un vent, dans un fichier ou dans
venuti.f.
**Mots-cles associes :**
FICHIER DES VENTS""",
            ang = """Indicates whether a wind is taken into account, either in
a file or in venuti.f
\\
 \begin{CommentBlock}{Related keywords}
WINDS FILE
\end{CommentBlock}""",
        ),
#       -----------------------------------
        STATIONARY_WIND = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [True ],
            fr = """Indique si le vent evolue dans le temps et doit etre mis a jour
**Mots-cles associes :**
PRISE EN COMPTE DU VENT""",
            ang = """Indicates whether the wind evolves temporally and requires
to be updated
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF A WIND
\end{CommentBlock}""",
        ),
#       -----------------------------------
        CONSIDERATION_OF_TIDE = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Indique si on prend en compte l''influence de la maree,
c''est-a-dire, prise en compte d''un niveau d''eau et de courants
instationnaires.
**Mots-cles associes :**
FICHIER DU NIVEAU DE LA MAREE FORMATE
FICHIER DU NIVEAU DE LA MAREE BINAIRE
FORMAT DU FICHIER DU NIVEAU DE LA MAREE
PERIODE D ACTUALISATION DE LA MAREE""",
            ang = """Indicates whether a current is taken into account, either in
a file or in cdicow.f.
\\
 \begin{CommentBlock}{Related keywords}
FORMATTED TIDAL WATER LEVEL FILE\\
BINARY TIDAL WATER LEVEL FILE\\
TIDAL WATER LEVEL FILE FORMAT\\
TIDE REFRESHING PERIOD
\end{CommentBlock}""",
        ),
    ),
#   -----------------------------------
    MISCELLANEOUS = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_PRIVATE_ARRAYS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Nombre de tableaux utilises en variables privees par
l''utilisateur""",
            ang = """Number of private arrays used by the user""",
        ),
#       -----------------------------------
        DEBUGGER = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Pour imprimer la sequence des appels, mettre 1""",
            ang = """If 1, calls of subroutines will be printed in the listing""",
        ),
#       -----------------------------------
        PARALLEL_PROCESSORS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
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
        FINITE_ELEMENT_ASSEMBLY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """1 : normal 2 : avec des entiers I8""",
            ang = """1: normal 2: with I8 integers""",
        ),
#       -----------------------------------
        TITLE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = 'SET A TITLE !!!',
            fr = """Titre du cas etudie""",
            ang = """Title of the case being studied.""",
        ),
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
\begin{itemize}
\item 1 : METIS
\item 2 : SCOTCH
\item 3 : PARMETIS
\item 4 : PTSCOTCH
\end{itemize}""",
        ),
#       -----------------------------------
        RECOVERY_OF_TELEMAC_DATA_ITEM = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Indique si on recupere des donnees TELEMAC dans LECDON.f
Si oui il faut veiller a utiliser un FICHIER DES COURANTS au bon format
et donner le rang de la variable en question dans le fichier TELEMAC.
**Mots-cles associes :**
FICHIER DES COURANTS BINAIRE
FICHIER DES COURANTS FORMATE
FORMAT DU FICHIER DES COURANTS
RANG DE LA DONNEE TELEMAC A RECUPERER
NUMERO DU PAS DE TEMPS DU FICHIER TELEMAC""",
            ang = """Indicates whether TELEMAC data are recovered in LECDON. If so,
a proper-formatted CURRENTS FILE should be used and the rank of the
respective variable should be entered into the TELEMAC file.
\\
 \begin{CommentBlock}{Related keywords}
BINARY CURRENTS FILE\\
FORMATTED CURRENTS FILE\\
CURRENTS FILE TYPE\\
RANK OF THE TELEMAC DATA ITEM TO BE RECOVERED\\
TIME INCREMENT NUMBER IN TELEMAC FILE
\end{CommentBlock}""",
        ),
#       -----------------------------------
        CONSIDERATION_OF_PROPAGATION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [True ],
            fr = """Indique si on prend en compte la propagation""",
            ang = """Indicates whether propagation is taken into account.""",
        ),
#       -----------------------------------
        VALIDATION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Logique indiquant si on effectue un calcul de validation
**Mots-cles associes :**
FICHIER DE REFERENCE""",
            ang = """True if the computation is a validation
\\
 \begin{CommentBlock}{Related keywords}
REFERENCE FILE
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_VALIDATIONG = BLOC(condition="VALIDATION == True",
#       -----------------------------------
        ),
#       -----------------------------------
        CHECKING_THE_MESH = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si oui on appelle le sous-programme checkmesh qui verifie
la coherence du maillage, points superposes, etc.""",
            ang = """if this key word is equal to yes, a call to subroutine
checkmesh will look for errors in the mesh, superimposed points, etc.""",
        ),
    ),
#   -----------------------------------
    OTHER_DOMAIN_DEFINITIONS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        ORIGIN_COORDINATES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I', min= 2, max= 2,
            defaut = [0,0],
            fr = """Valeur en metres, utilise pour eviter les trop grands nombres,
transmis dans le format Selafin mais pas d''autre traitement pour
l''instant""",
            ang = """Value in metres, used to avoid large real numbers,
added in Selafin format, but so far no other treatment""",
        ),
#       -----------------------------------
        MINIMUM_WATER_DEPTH = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 0.1,
            fr = """Definit la profondeur d eau minimale en dessous de laquelle
les fonds sont supposes emerges.""",
            ang = """Defines the minimum water depth below which bottom elevations
are regarded as dry.""",
        ),
#       -----------------------------------
        SPHERICAL_COORDINATES = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Indique si on se place ou non en coordonnes spheriques.
ATTENTION, en coordonnees cartesiennes, les coordonnees sont
exprimees em m alors que ce sont des degres en coordonnees
spheriques.""",
            ang = """Indicates whether the coordinates are spherical (unit=
degree) or cartesian (unit = meter).""",
        ),
#       -----------------------------------
        INFINITE_DEPTH = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Indique si on se place dans l hypothese de profondeur infinie.
Cette option inhibe les frottements sur le fond.""",
            ang = """Indicates whether an infinite depth is assumed. If so, bottom
friction is inhibited.""",
        ),
#       -----------------------------------
        TRIGONOMETRICAL_CONVENTION = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Logique indiquant si les directions de propagation
de la houle sont definies dans le sens trigonometrique a
partir de l axe des x positifs, ou definies dans le sens
des aiguilles d une montre a partir du nord geographique""",
            ang = """True if the wave directions are measured
counterclockwise from the positive x-axis, false if
they are measured clockwise fron geographic North""",
        ),
    ),
)
# -----------------------------------------------------------------------
SOURCE_TERMS = PROC(nom= "SOURCE_TERMS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    BAJ_MODELING = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        into = ["classical modelisation","BAJ Modeling"],
        defaut = ["classical modelisation"],
        fr = """Type de calcul de la frequence centrale
Si sa valeur est 0, on prend la formulation classique
Si sa valeur est 1 on pren la modelisation BAJ proposee
par Laugel (2013) pour le calcul de la frequence centrale
et le limiteur de croissance
**Mots-cles associes :**
PRISE EN COMPTE DES TERMES SOURCES""",
        ang = """Choice of the calculus of centrale frequency
if its value is 0, classical choice
if its value is 1, BAJ choice proposed by Laugel (2013).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF SOURCE TERMS
\end{CommentBlock}""",
    ),
#   -----------------------------------
    CONSIDERATION_OF_SOURCE_TERMS = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = True,
        fr = """Indique la prise en compte ou non de l ensemble des termes
sources.
**Mots-cles associes :**
APPORTS DUS AU VENT
DISSIPATION PAR FROTTEMENT SUR LE FOND
DISSIPATION PAR MOUTONNEMENT
DISSIPATION PAR DEFERLEMENT
DISSIPATION PAR WAVE BLOCKING
TRANSFERTS NON LINEAIRES INTER-FREQUENCES
TRANSFERTS ENTRE TRIPLETS DE FREQUENCES""",
        ang = """Indicates whether the source terms are taken into
account or not.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION\\
BOTTOM FRICTION DISSIPATION\\
WHITE CAPPING DISSIPATION\\
DEPTH-INDUCED BREAKING DISSIPATION\\
WAVE BLOCKING DISSIPATION\\
NON-LINEAR TRANSFERS BETWEEN FREQUENCIES\\
TRIAD INTERACTION
\end{CommentBlock}""",
    ),
#   -----------------------------------
    TAKING_INTO_ACCOUNT_SOURCE_TERMS_ON_IMPOSED_BOUNDARIES = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [True ],
        fr = """Indique si on prend en compte les termes sources sur les
frontieres imposees.""",
        ang = """Indicates whether source terms are taken into account on
imposed boundaries.""",
    ),
#   -----------------------------------
    WIND = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WIND_GENERATION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["No wind generation","Wind generation in accordance with WAM cycle 4","Wind generation in accordance with WAM cycle 3","Wind generation in accordance with Yan expression (1987)"],
            defaut = ["No wind generation"],
            fr = """Choix du type de modelisation du terme source d apport
par le vent. Si sa valeur est 0, on ne prend pas en compte les
apports dus au vent, si sa valeur est 1, ils sont integres selon
la formule utilisee dans WAM cycle 4;
si sa valeur est 2, ils sont integres selon la formule utilisee
dans WAM cycle 3; si sa valeur est 3, ils sont integres selon
l expression de Yan (1987).
**Mots-cles associes :**
PRISE EN COMPTE DU VENT
FICHIER DES VENTS
DENSITE DE L''AIR
DENSITE DE L''EAU
COEFFICIENT DE GENERATION PAR LE VENT
CONSTANTE DE VON KARMAN
CONSTANTE DE CHARNOCK
DECALAGE COURBE DE CROISSANCE DUE AU VENT
COTE DE MESURE DES VENTS
COEFFICIENT DE TRAINEE DE VENT
COEFFICIENT A DE GENERATION PAR LE VENT
COEFFICIENT B DE GENERATION PAR LE VENT
COEFFICIENT C DE GENERATION PAR LE VENT
COEFFICIENT D DE GENERATION PAR LE VENT
COEFFICIENT TM DE GENERATION PAR LE VENT""",
            ang = """Selection of the type of modelling of the wind generation
source term. If its value is 0, the wind generation is ignored; if
its value is 1, it is integrated in accordance with the WAM cycle 4
formula; if its value is 2, it is integrated in accordance with the
WAM cycle 3 formula; if its value is 3, it is integrated in accordance
with the Yan (1987) expression.
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF A WIND\\
WINDS FILE\\
AIR DENSITY\\
WATER DENSITY\\
WIND GENERATION COEFFICIENT\\
VON KARMAN CONSTANT\\
CHARNOCK CONSTANT\\
SHIFT GROWING CURVE DUE TO WIND\\
WIND MEASUREMENTS LEVEL\\
WIND DRAG COEFFICIENT\\
WIND GENERATION COEFFICIENT A\\
WIND GENERATION COEFFICIENT B\\
WIND GENERATION COEFFICIENT C\\
WIND GENERATION COEFFICIENT D\\
WIND GENERATION COEFFICIENT TM\\
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_WIND_GENERATIONG = BLOC(condition="WIND_GENERATION == 'Wind generation in accordance with WAM cycle 4'",
#       -----------------------------------
#           -----------------------------------
            WIND_GENERATION_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.2,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT""",
                ang = """Constant used in the wind source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION
\end{CommentBlock}""",
            ),
#           -----------------------------------
            AIR_DENSITY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.225,
                fr = """Le rapport Roair/Roeau est utilise dans le terme source
de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT
DENSITE DE L''EAU""",
                ang = """The ratio ROAIR/ROEAU is used in the wind generation
source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION\\
WATER DENSITY
\end{CommentBlock}""",
            ),
#           -----------------------------------
            WATER_DENSITY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1000.,
                fr = """ Le rapport Roair/Roeau est utilise dans le terme source
de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT
DENSITE DE L''AIR""",
                ang = """The ratio ROAIR/ROEAU is used in the wind generation
source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION\\
AIR DENSITY
\end{CommentBlock}""",
            ),
#           -----------------------------------
            CHARNOCK_CONSTANT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.01,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT""",
                ang = """Constant used in the wind source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION
\end{CommentBlock}""",
            ),
#           -----------------------------------
            WIND_DRAG_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.2875E-3,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT""",
                ang = """Constant used in the wind source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION
\end{CommentBlock}""",
            ),
#           -----------------------------------
            VON_KARMAN_CONSTANT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.41,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT""",
                ang = """Constant used in the wind source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION
\end{CommentBlock}""",
            ),
#           -----------------------------------
            SHIFT_GROWING_CURVE_DUE_TO_WIND = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.011,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT""",
                ang = """Constant used in the wind source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION
\end{CommentBlock}""",
            ),
#           -----------------------------------
            WIND_MEASUREMENTS_LEVEL = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 10.,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent.
**Mots-cles associes :**
APPORTS DUS AU VENT""",
                ang = """Constant used in the wind source term.
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        b_WIND_GENERATIONH = BLOC(condition="WIND_GENERATION == 'Wind generation in accordance with WAM cycle 3'",
#       -----------------------------------
        ),
#       -----------------------------------
        b_WIND_GENERATIONI = BLOC(condition="WIND_GENERATION == 'Wind generation in accordance with Yan expression (1987)'",
#       -----------------------------------
#           -----------------------------------
            YAN_GENERATION_COEFFICIENT_D = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.04,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent de Yan (1987).
**Mots-cles associes :**
APPORTS DUS AU VENT
COEFFICIENT DE GENERATION DE YAN E
COEFFICIENT DE GENERATION DE YAN F
COEFFICIENT DE GENERATION DE YAN H""",
                ang = """Constant used in the wind source term of Yan (1987).
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION\\
YAN GENERATION COEFFICIENT E\\
YAN GENERATION COEFFICIENT F\\
YAN GENERATION COEFFICIENT H
\end{CommentBlock}""",
            ),
#           -----------------------------------
            YAN_GENERATION_COEFFICIENT_E = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.00552,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent de Yan (1987).
**Mots-cles associes :**
APPORTS DUS AU VENT
COEFFICIENT DE GENERATION DE YAN D
COEFFICIENT DE GENERATION DE YAN F
COEFFICIENT DE GENERATION DE YAN H""",
                ang = """Constant used in the wind source term of Yan (1987).
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION\\
YAN GENERATION COEFFICIENT D\\
YAN GENERATION COEFFICIENT F\\
YAN GENERATION COEFFICIENT H
\end{CommentBlock}""",
            ),
#           -----------------------------------
            YAN_GENERATION_COEFFICIENT_F = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.000052,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent de Yan (1987).
**Mots-cles associes :**
APPORTS DUS AU VENT
COEFFICIENT DE GENERATION DE YAN D
COEFFICIENT DE GENERATION DE YAN E
COEFFICIENT DE GENERATION DE YAN H""",
                ang = """Constant used in the wind source term of Yan (1987).
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION\\
YAN GENERATION COEFFICIENT D\\
YAN GENERATION COEFFICIENT E\\
YAN GENERATION COEFFICIENT H
\end{CommentBlock}""",
            ),
#           -----------------------------------
            YAN_GENERATION_COEFFICIENT_H = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = -0.000302,
                fr = """Fait partie de l ensemble des constantes utilisees dans
le terme source de generation par le vent de Yan (1987).
**Mots-cles associes :**
APPORTS DUS AU VENT
COEFFICIENT DE GENERATION DE YAN D
COEFFICIENT DE GENERATION DE YAN E
COEFFICIENT DE GENERATION DE YAN F""",
                ang = """Constant used in the wind source term of Yan (1987).
\\
 \begin{CommentBlock}{Related keywords}
WIND GENERATION\\
YAN GENERATION COEFFICIENT D\\
YAN GENERATION COEFFICIENT E\\
YAN GENERATION COEFFICIENT F
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        LINEAR_WAVE_GROWTH = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["the linear wave growth is ignored","linear wave growth as in Cavaleri and Malanotte-Rizzoli (1981)"],
            defaut = ["the linear wave growth is ignored"],
            fr = """Possibilite d ajouter au terme source d apport
par le vent un terme de croissance lineaire.
Si sa valeur est 0, on ne prend pas en compte le
terme de croissance lineaire, si sa valeur est 1,
il est ajoute au terme source, selon la formule
de Cavaleri et Malanotte-Rizzoli (1981)
**Mots-cles associes :**
PRISE EN COMPTE DU VENT
FICHIER DES VENTS""",
            ang = """Possibility to add a linear wave growth term
to the wind generation source term.
If its value is 0, the linear wave growth is ignored;
if its value is 1, it is added to the source term,
as in the formula of Cavaleri and Malanotte-Rizzoli (1981).
\\
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF A WIND\\
WINDS FILE
\end{CommentBlock}""",
        ),
    ),
#   -----------------------------------
    WHITE_CAPPING = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WHITE_CAPPING_DISSIPATION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["No white capping dissipation","Dissipation in accordance with WAM cycle 4","Dissipation in accordance with Van des Westhuysen(2007)"],
            defaut = ["No white capping dissipation"],
            fr = """Choix du type de modelisation du terme source de dissipation
par moutonnement. Si sa valeur est 0, on ne prend pas en compte la
dissipation par moutonnement, si sa valeur est 1, elle est integree
selon la formule utilisee dans WAM cycle 4; si sa valeur est 2, elle
est integree selon la formule de van der Westhuysen (2007).
**Mots-cles associes :**
COEFFICIENT DE DISSIPATION PAR MOUTONNEMENT
COEFFICIENT DE PONDERATION POUR LE MOUTONNEMENT
COEFFICIENT DE DISSIPATION DE WESTHUYSEN
DISSIPATION PAR MOUTONNEMENT
SEUIL DE SATURATION POUR LA DISSIPATION
DISSIPATION PAR MOUTONNEMENT DE WESTHUYSEN
COEFFICIENT DE PONDERATION DE WESTHUYSEN""",
            ang = """Selection of the modelling type of the white capping source
term. If its value is 0, the white capping dissipation is ignored;
if its value is 1, it is integrated in accordance with a formula that
is similar to that of WAM cycle 4; if its value is 2, it is integrated
in accordance with the formula of van der Westhuysen (2007).
\\
 \begin{CommentBlock}{Related keywords}
WHITE CAPPING DISSIPATION COEFFICIENT\\
WHITE CAPPING WEIGHTING COEFFICIENT\\
WESTHUYSEN DISSIPATION COEFFICIENT\\
SATURATION THRESHOLD FOR THE DISSIPATION\\
WESTHUYSEN WHITE CAPPING DISSIPATION\\
WESTHUYSEN WEIGHTING COEFFICIENT
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_WHITE_CAPPING_DISSIPATIONG = BLOC(condition="WHITE_CAPPING_DISSIPATION == 'Dissipation in accordance with WAM cycle 4'",
#       -----------------------------------
#           -----------------------------------
            WHITE_CAPPING_DISSIPATION_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 4.5,
                fr = """Coefficient de dissipation par moutonnement.
**Mots-cles associes :**
DISSIPATION PAR MOUTONNEMENT
COEFFICIENT DE PONDERATION POUR LE MOUTONNEMENT""",
                ang = """White capping dissipation coefficient .
\\
 \begin{CommentBlock}{Related keywords}
WHITE CAPPING DISSIPATION\\
WHITE CAPPING WEIGHTING COEFFICIENT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            WHITE_CAPPING_WEIGHTING_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.5,
                fr = """ Coefficient de ponderation pour le moutonnement.
**Mots-cles associes :**
DISSIPATION PAR MOUTONNEMENT
COEFFICIENT DE DISSIPATION PAR MOUTONNEMENT""",
                ang = """White capping weighting coefficient.
\\
 \begin{CommentBlock}{Related keywords}
WHITE CAPPING DISSIPATION\\
WHITE CAPPING DISSIPATION COEFFICIENT
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        b_WHITE_CAPPING_DISSIPATIONH = BLOC(condition="WHITE_CAPPING_DISSIPATION == 'Dissipation in accordance with Van des Westhuysen(2007)'",
#       -----------------------------------
#           -----------------------------------
            SATURATION_THRESHOLD_FOR_THE_DISSIPATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.00175,
                fr = """Coefficient de dissipation par moutonnement de l expression
de van der Westhuysen (2007): Br (seuil de saturation).
**Mots-cles associes :**
DISSIPATION PAR MOUTONNEMENT
COEFFICIENT DE DISSIPATION DE WESTHUYSEN
DISSIPATION PAR MOUTONNEMENT DE WESTHUYSEN
COEFFICIENT DE PONDERATION DE WESTHUYSEN""",
                ang = """White capping dissipation coefficient of
van der Westhuysen (2007): Br (saturation threshold).
\\
 \begin{CommentBlock}{Related keywords}
WHITE CAPPING DISSIPATION\\
WESTHUYSEN DISSIPATION COEFFICIENT\\
WESTHUYSEN WHITE CAPPING DISSIPATION\\
WESTHUYSEN WEIGHTING COEFFICIENT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            WESTHUYSEN_WHITE_CAPPING_DISSIPATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 3.29,
                fr = """Coefficient de dissipation par moutonnement de l expression
de van der Westhuysen (2007): Cdis,non-break.
**Mots-cles associes :**
DISSIPATION PAR MOUTONNEMENT
COEFFICIENT DE DISSIPATION DE WESTHUYSEN
SEUIL DE SATURATION POUR LA DISSIPATION
COEFFICIENT DE PONDERATION DE WESTHUYSEN""",
                ang = """White capping dissipation coefficient of
van der Westhuysen (2007): Cdis,non-break.
\\
 \begin{CommentBlock}{Related keywords}
WHITE CAPPING DISSIPATION\\
WESTHUYSEN DISSIPATION COEFFICIENT\\
SATURATION THRESHOLD FOR THE DISSIPATION\\
WESTHUYSEN WEIGHTING COEFFICIENT
\end{CommentBlock}""",
            ),
#           -----------------------------------
            WESTHUYSEN_WEIGHTING_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.0,
                fr = """Coefficient de dissipation par moutonnement de l expression
de van der Westhuysen (2007): delta.
**Mots-cles associes :**
DISSIPATION PAR MOUTONNEMENT
COEFFICIENT DE DISSIPATION DE WESTHUYSEN
SEUIL DE SATURATION POUR LA DISSIPATION
DISSIPATION PAR MOUTONNEMENT DE WESTHUYSEN""",
                ang = """White capping dissipation coefficient of
van der Westhuysen (2007): delta.
\\
 \begin{CommentBlock}{Related keywords}
WHITE CAPPING DISSIPATION\\
WESTHUYSEN DISSIPATION COEFFICIENT\\
SATURATION THRESHOLD FOR THE DISSIPATION\\
WESTHUYSEN WHITE CAPPING DISSIPATION
\end{CommentBlock}""",
            ),
#           -----------------------------------
            WESTHUYSEN_DISSIPATION_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.00005,
                fr = """Coefficient de dissipation par moutonnement de l expression
de van der Westhuysen (2007): Cdis,break.
**Mots-cles associes :**
DISSIPATION PAR MOUTONNEMENT
SEUIL DE SATURATION POUR LA DISSIPATION
DISSIPATION PAR MOUTONNEMENT DE WESTHUYSEN
COEFFICIENT DE PONDERATION DE WESTHUYSEN""",
                ang = """White capping dissipation coefficient of
van der Westhuysen (2007): Cdis,break.
\\
 \begin{CommentBlock}{Related keywords}
WHITE CAPPING DISSIPATION\\
SATURATION THRESHOLD FOR THE DISSIPATION\\
WESTHUYSEN WHITE CAPPING DISSIPATION\\
WESTHUYSEN WEIGHTING COEFFICIENT
\end{CommentBlock}""",
            ),
        ),
    ),
#   -----------------------------------
    BOTTOM_FRICTION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        BOTTOM_FRICTION_DISSIPATION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["No bottom friction dissipation","Dissipation in accordance with WAM cycle 4"],
            defaut = ["No bottom friction dissipation"],
            fr = """Choix du type de modelisation du terme source de dissipation
sur le fond. Si sa valeur est 0, on ne prend pas en compte la
dissipation par frottement, si sa valeur est 1, elle est integree
selon la formule equivalente a celle utilisee dans WAM cycle 4.
**Mots-cles associes :**
PROFONDEUR INFINIE
COEFFICIENT DE FROTTEMENT SUR LE FOND""",
            ang = """Selection of the modelling type of the bottom friction source
term. If its value is 0, the bottom friction dissipation is ignored;
if its value is 1, it is integrated in accordance with a formula that
is similar to that of WAM cycle 4.
\\
 \begin{CommentBlock}{Related keywords}
INFINITE DEPTH\\
BOTTOM FRICTION COEFFICIENT
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_BOTTOM_FRICTION_DISSIPATIONG = BLOC(condition="BOTTOM_FRICTION_DISSIPATION == 'Dissipation in accordance with WAM cycle 4'",
#       -----------------------------------
#           -----------------------------------
            BOTTOM_FRICTION_COEFFICIENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.038,
                fr = """Coefficient de frottement sur le fond.
**Mots-cles associes :**
PROFONDEUR INFINIE
DISSIPATION PAR FROTTEMENT SUR LE FOND""",
                ang = """Bottom friction coefficient.
\\
 \begin{CommentBlock}{Related keywords}
INFINITE DEPTH\\
BOTTOM FRICTION-INDUCED DISSIPATION
\end{CommentBlock}""",
            ),
        ),
    ),
#   -----------------------------------
    QUADRUPLET_INTERACTIONS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        NON_LINEAR_TRANSFERS_BETWEEN_FREQUENCIES = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["No non-linear transfers term","Non-linear transfers term with WAM cycle 4 (DIA Method)","Non-linear transfers term in accordance with MDIA Method","Non-linear transfers term calculated with exact GQM Method"],
            defaut = ["No non-linear transfers term"],
            fr = """Choix du type de modelisation du terme de transfert non lineaire
inter-frequences. Si sa valeur est 0, on ne prend pas en compte les
transferts non lineaires inter-frequences, si sa valeur est 1, ils sont
integres selon la formule utilisee dans WAM cycle 4 (methode DIA), si sa
valeur est 2, la methode MDIA (Multiple DIA) est utilisee pour calculer
le terme de transfert non lineaire, si sa valeur est 3, le terme de
transfert non lineaire est calcule de maniere exacte avec la methode
GQM.
**Mots-cles associes :**
PARAMETRE DE LA CONFIGURATION STANDARD
REGLAGE POUR INTEGRATION SUR OMEGA1
REGLAGE POUR INTEGRATION SUR THETA1
REGLAGE POUR INTEGRATION SUR OMEGA2
SEUIL0 ELIMINATION DE CONFIGURATIONS
SEUIL1 ELIMINATION DE CONFIGURATIONS
SEUIL2 ELIMINATION DE CONFIGURATIONS""",
            ang = """Selection of the modelling type of the non-linear transfert
source term. If its value is 0, the non-linear transfers are ignored;
if its value is 1, they are integrated in accordance with the formula
of WAM cycle 4 (DIA method), if its value is 2, the MDIA
(Multiple DIA) method is used to calculate the non linear transfer
term, if its value is 3, the non linear transfer term is calculated
with the exact GQM method.
\\
 \begin{CommentBlock}{Related keywords}
STANDARD CONFIGURATION PARAMETER\\
SETTING FOR INTEGRATION ON OMEGA1\\
SETTING FOR INTEGRATION ON THETA1\\
SETTING FOR INTEGRATION ON OMEGA2\\
THRESHOLD0 FOR CONFIGURATIONS ELIMINATION\\
THRESHOLD1 FOR CONFIGURATIONS ELIMINATION\\
THRESHOLD2 FOR CONFIGURATIONS ELIMINATION
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_NON_LINEAR_TRANSFERS_BETWEEN_FREQUENCIESG = BLOC(condition="NON_LINEAR_TRANSFERS_BETWEEN_FREQUENCIES == 'Non-linear transfers term with WAM cycle 4 (DIA Method)'",
#       -----------------------------------
#           -----------------------------------
            STANDARD_CONFIGURATION_PARAMETER = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.25,
                fr = """Parametre definissant la configuration d interaction
standard pour les quadruplets dans la methode DIA.
**Mots-cles associes :**
TRANSFERTS NON LINEAIRES INTER-FREQUENCES""",
                ang = """Parameter defining the standard configuration for
the quadruplet interactions in the DIA method.
\\
 \begin{CommentBlock}{Related keywords}
NON-LINEAR TRANSFERS BETWEEN FREQUENCIES
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        b_NON_LINEAR_TRANSFERS_BETWEEN_FREQUENCIESH = BLOC(condition="NON_LINEAR_TRANSFERS_BETWEEN_FREQUENCIES == 'Non-linear transfers term calculated with exact GQM Method'",
#       -----------------------------------
#           -----------------------------------
            SETTING_FOR_INTEGRATION_ON_OMEGA1 = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["medium","fine","rough"],
                defaut = ["rough"],
                fr = """ Choix du reglage donnant le nombre de point d integration sur
omega1 lorsque le terme de transfert non lineaire est calcule de
maniere exacte (methode GQM): grossier 3;moyen 1 ; fin 2""",
                ang = """Choice of setting giving the number of integration points on
omega1 when the non linear transfer term is calculated with the exact
GQM method: rough 3 ; medium 1 ; fine 2""",
            ),
#           -----------------------------------
            SETTING_FOR_INTEGRATION_ON_THETA1 = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["rough","medium","fine"],
                defaut = ["rough"],
                fr = """ Choix du reglage donnant le nombre de point d integration sur
theta1 (nombre de points d integration = 2*NQ\_TE1) lorsque le terme
de transfert non lineaire est calcule de maniere exacte (methode GQM):
grossier 3 ; moyen 4 ; fin 8""",
                ang = """Choice of setting giving the number of integration points on
theta1 (number of integration points= 2*NQ\_TE1) when the non linear
transfer term is calculated with the exact GQM method:
rough 3 ; medium 4 ; fine 8""",
            ),
#           -----------------------------------
            SETTING_FOR_INTEGRATION_ON_OMEGA2 = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["rough","medium","fine"],
                defaut = ["rough"],
                fr = """ Nombre de point d integration sur omega2 lorsque le terme de
transfert non lineaire est calcule de maniere exacte (methode GQM):
grossier 6 ; moyen 8 ; fin 12""",
                ang = """Number of integration points on omega2 when the non linear
transfer term is calculated with the exact GQM method:
rough 6 ; medium 8 ; fine 12""",
            ),
#           -----------------------------------
            THRESHOLD0_FOR_CONFIGURATIONS_ELIMINATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.00,
                fr = """ Choix du seuil pour l elimination de configurations lorsque
le terme de transfert non lineaire est calcule de maniere exacte
(methode GQM)
**Mots-cles associes :**
SEUIL1 ELIMINATION DE CONFIGURATIONS
SEUIL2 ELIMINATION DE CONFIGURATIONS
TRANSFERTS NON LINEAIRES INTER-FREQUENCES""",
                ang = """Choice of threshold for configurations elimination when the
non linear transfer term is calculated with the exact GQM method
\\
 \begin{CommentBlock}{Related keywords}
THRESHOLD1 FOR CONFIGURATIONS ELIMINATION\\
THRESHOLD2 FOR CONFIGURATIONS ELIMINATION\\
NON-LINEAR TRANSFERS BETWEEN FREQUENCIES
\end{CommentBlock}""",
            ),
#           -----------------------------------
            THRESHOLD1_FOR_CONFIGURATIONS_ELIMINATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 10000000000.0,
                fr = """ Choix du seuil1 pour l elimination de configurations lorsque
le terme de transfert non lineaire est calcule de maniere exacte
(methode GQM)
**Mots-cles associes :**
SEUIL0 ELIMINATION DE CONFIGURATIONS
SEUIL2 ELIMINATION DE CONFIGURATIONS
TRANSFERTS NON LINEAIRES INTER-FREQUENCES""",
                ang = """Choice of threshold1 for configurations elimination when the
non linear transfer term is calculated with the exact GQM method
\\
 \begin{CommentBlock}{Related keywords}
THRESHOLD0 FOR CONFIGURATIONS ELIMINATION\\
THRESHOLD2 FOR CONFIGURATIONS ELIMINATION\\
NON-LINEAR TRANSFERS BETWEEN FREQUENCIES
\end{CommentBlock}""",
            ),
#           -----------------------------------
            THRESHOLD2_FOR_CONFIGURATIONS_ELIMINATION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.15,
                fr = """ Choix du seuil2 pour l elimination de configurations lorsque
le terme de transfert non lineaire est calcule de maniere exacte
(methode GQM) : grossier 0.15 ; moyen 0.01 ; fin 0.001
**Mots-cles associes :**
SEUIL0 ELIMINATION DE CONFIGURATIONS
SEUIL1 ELIMINATION DE CONFIGURATIONS
TRANSFERTS NON LINEAIRES INTER-FREQUENCES""",
                ang = """Choice of threshold2 for configurations elimination when the
non linear transfer term is calculated with the exact GQM method:
rough 0.15 ; medium 0.01 ; fine 0.001
\\
 \begin{CommentBlock}{Related keywords}
THRESHOLD0 FOR CONFIGURATIONS ELIMINATION\\
THRESHOLD1 FOR CONFIGURATIONS ELIMINATION\\
NON-LINEAR TRANSFERS BETWEEN FREQUENCIES
\end{CommentBlock}""",
            ),
        ),
    ),
#   -----------------------------------
    NUMERICAL_PARAMETERS = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        NUMBER_OF_ITERATIONS_FOR_THE_SOURCE_TERMS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Nombre de sous-iterations pour le calcul des
contributions des termes sources. Le pas de temps considere
pour les termes sources est le rapport entre le PAS DE TEMPS
et le NOMBRE DE SOUS-ITERATIONS POUR LES TERMES SOURCES.
**Mots-cles associes :**
PAS DE TEMPS""",
            ang = """Number of sub-iterations for the computation of the source terms.
The time step considered in the integration of the source terms
is the ratio between the TIME STEP and the NUMBER OF SUB-ITERATIONS
FOR THE SOURCE TERMS
\\
 \begin{CommentBlock}{Related keywords}
TIME STEP
\end{CommentBlock}""",
        ),
#       -----------------------------------
        IMPLICITATION_COEFFICIENT_FOR_SOURCE_TERMS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.5,
            fr = """Coefficient controlant l implicitation dans le schema
d integration des termes sources, compris entre 0 et 1.
CIMPLI=0.  : explicite
CIMPLI=0.5 : semi-implicite
CIMPLI=1.  : implicite.
**Mots-cles associes :**
PRISE EN COMPTE DES TERMES SOURCES""",
            ang = """Implicitation coefficient for the source terms integration,
included between 0 et 1.
\begin{itemize}
\item CIMPLI=0.  : explicit
\item CIMPLI=0.5 : semi-implicit
\item CIMPLI=1.  : implicit.
\end{itemize}
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF SOURCE TERMS
\end{CommentBlock}""",
        ),
    ),
#   -----------------------------------
    STRONG_CURRENT = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DISSIPATION_BY_STRONG_CURRENT = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["No wave-blocking","Spectrum limitation by a Phillips shape","Dissipation in accordance with Van der Westhuysen(2012)"],
            defaut = ["No wave-blocking"],
            fr = """Lorsque les effets de wave-blocking ou vagues stoppees
par un fort courant sont presents, deux options sont proposees.
Si sa valeur est 1, une limitation est imposee au spectre,
en utilisant une forme d equilibre de Phillips (1977).
Si sa valeur est 2, on utilise le terme de dissipation propose
par Van der Westhuysen (2012).
**Mots-cles associes :**
COEFFICIENT DE DISSIPATION PAR FORT COURANT""",
            ang = """When wave-blocking effects are present (wave stopped by
a strong opposing current), two options are possible.
If its value is 1, an upper limit is imposed to the spectrum,
using a Phillips (1977) shape.
If its value is 2, a dissipative term is added, following
Van der Westhuysen (2012).
\\
 \begin{CommentBlock}{Related keywords}
DISSIPATION COEFFICIENT FOR STRONG CURRENT
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_DISSIPATION_BY_STRONG_CURRENTG = BLOC(condition="DISSIPATION_BY_STRONG_CURRENT == 'Dissipation in accordance with Van der Westhuysen(2012)'",
#       -----------------------------------
#           -----------------------------------
            DISSIPATION_COEFFICIENT_FOR_STRONG_CURRENT = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.65,
                fr = """Coefficient de dissipation pour des vagues stoppees
par un courant fort adverse (effets de wave-blocking).
Expression de van der Westhuysen (2012): Cds,cur.
**Mots-cles associes :**
DISSIPATION PAR FORT COURANT""",
                ang = """Dissipation coefficient for waves stopped
by a strong opposing current (wave blocking effects).
Van der Westhuysen (2012) expression: Cds,cur.
\\
 \begin{CommentBlock}{Related keywords}
DISSIPATION BY STRONG CURRENT
\end{CommentBlock}""",
            ),
        ),
    ),
#   -----------------------------------
    BREAKING = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DEPTH_INDUCED_BREAKING_DISSIPATION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["No breaking","Dissipation in accordance with Battjes et Janssen (1978)","Dissipation in accordance with Thornton et Guza (1983)","Dissipation in accordance with Roelvink (1993)","Dissipation in accordance with Izumiya et Horikawa (1984)"],
            defaut = ["No breaking"],
            fr = """Choix du type de modelisation du terme source de dissipation
par deferlement du a la bathymetrie :
 0 : Pas de prise en compte du deferlement.
 1 : Modele de Battjes et Janssen (1978).
 2 : Modele de Thornton et Guza  (1983).
 3 : Modele de Roelvink (1993).
 4 : Modele de Izumiya et Horikawa (1984).
**Mots-cles associes :**
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 1 (BJ) CONSTANTE ALPHA
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2
DEFERLEMENT 2 (TG) FONCTION DE PONDERATION
DEFERLEMENT 2 (TG) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 2 (TG) CONSTANTE B
DEFERLEMENT 2 (TG) CONSTANTE GAMMA
DEFERLEMENT 3 (RO) DISTRIBUTION DES HAUTEURS DE HOULE
DEFERLEMENT 3 (RO) EXPOSANT FONCTION DE PONDERATION
DEFERLEMENT 3 (RO) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 3 (RO) CONSTANTE ALPHA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA2
DEFERLEMENT 4 (IH) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 4 (IH) CONSTANTE BETA0
DEFERLEMENT 4 (IH) CONSTANTE M2STAR""",
            ang = """Selection of the modelling type of the bathymetric-induced
breaking dissipation source term :
\begin{itemize}
\item 0 : Breaking is ignored.
\item 1 : Battjes and Janssen model (1978).
\item 2 : Thornton and Guza model (1983).
\item 3 : Roelvink model (1993).
\item 4 : Izumiya and Horikawa model (1984).
\end{itemize}
 \begin{CommentBlock}{Related keywords}
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 1 (BJ) QB COMPUTATION METHOD\\
DEPTH-INDUCED BREAKING 1 (BJ) HM COMPUTATION METHOD\\
DEPTH-INDUCED BREAKING 1 (BJ) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA1\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA2\\
DEPTH-INDUCED BREAKING 2 (TG) WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 2 (TG) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT B\\
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT GAMMA\\
DEPTH-INDUCED BREAKING 3 (RO) WAVE HEIGHT DISTRIBUTION\\
DEPTH-INDUCED BREAKING 3 (RO) EXPONENT WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 3 (RO) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA2\\
DEPTH-INDUCED BREAKING 4 (IH) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 4 (IH) COEFFICIENT BETA0\\
DEPTH-INDUCED BREAKING 4 (IH) COEFFICIENT M2STAR
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_DEPTH_INDUCED_BREAKING_DISSIPATIONG = BLOC(condition="DEPTH_INDUCED_BREAKING_DISSIPATION == 'Dissipation in accordance with Battjes et Janssen (1978)'",
#       -----------------------------------
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_1__BJ__QB_COMPUTATION_METHOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["SOLVES BY DICHOTOMY","EXPLICIT INSPIRED FROM CREDIZ VERSION 1","EXPLICIT INSPIRED FROM CREDIZ VERSION 2","EXPLICIT INSPIRED FROM CREDIZ VERSION 3"],
                defaut = ["EXPLICIT INSPIRED FROM CREDIZ VERSION 2"],
                fr = """Choix du mode de resolution de l equation implicite donnant Qb.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 1 (BJ) CONSTANTE ALPHA
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2""",
                ang = """Selection of the method for the resolution of the implicit
equation for QB.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 1 (BJ) HM COMPUTATION METHOD\\
DEPTH-INDUCED BREAKING 1 (BJ) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA1\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_1__BJ__HM_COMPUTATION_METHOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Hm = GAMMA*D","Hm given by the Miche criterium"],
                defaut = ["Hm = GAMMA*D"],
                fr = """Choix du critere de deferlement donnant la hauteur de houle
de deferlement (1 : Hm = GAMMA*D ; 2 : Hm par critere de Miche).
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2""",
                ang = """Selection of the depth-induced breaking criterium
giving the breaking wave height (1 : Hm = GAMMA*D ; 2 : Hm given
the Miche criterium).
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 1 (BJ) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA1\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA2\\
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_1__BJ__CHARACTERISTIC_FREQUENCY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Frequency Fmoy","Frequency F01","Frequency F02","Frequency Fpic","Frequency Fread ordre 5","Frequency Fread ordre 8"],
                defaut = ["Frequency F01"],
                fr = """Choix de la frequence caracteristique du spectre de houle
1 : Frequence Fmoy
2 : Frequence F01 (definie par les moments d ordre 0 et 1 du spectre)
3 : Frequence F02 (definie par les moments d ordre 0 et 2 du spectre)
4 : Frequence Fpic (frequence d echantillonage correspondant au max)
5 : Frequence Fread ordre 5 (frequence de pic methode Read ordre 5)
6 : Frequence Fread ordre 8 (frequence de pic methode Read ordre 8)
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM
DEFERLEMENT 1 (BJ) CONSTANTE ALPHA
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2""",
                ang = """Selection of the characteristic frequency of the wave spectrum
\begin{itemize}
\item 1 : Frequency Fmoy
\item 2 : Frequency F01 (defined by the moments of order 0 and 1 of the
spectrum)
\item 3 : Frequency F02 (defined by the moments of order 0 and 2 of the
spectrum)
\item 4 : Frequency Fpic (sampling frequency corresponding to the max)
\item 5 : Frequency Fread ordre 5 (peak frequency, 5th order Read
method)
\item 6 : Frequency Fread ordre 8 (peak frequency, 8th order Read
method)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 1 (BJ) QB COMPUTATION METHOD\\
DEPTH-INDUCED BREAKING 1 (BJ) HM COMPUTATION METHOD\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA1\\
DEPTH-INDUCED BREAKING 1 (BJ) COEFFICIENT GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_1__BJ__COEFFICIENT_ALPHA = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Constante ALPHA du modele de deferlement de Battjes et Janssen.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2""",
                ang = """ALPHA constant for the Battjes and Janssen model.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB\\
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM\\
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE\\
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1\\
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_1__BJ__COEFFICIENT_GAMMA1 = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.88,
                fr = """Constante GAMMA1 du modele de deferlement de Battjes et Janssen.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 1 (BJ) CONSTANTE ALPHA
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2""",
                ang = """GAMMA1 constant of the Battjes and Janssen model.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB\\
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM\\
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE\\
DEFERLEMENT 1 (BJ) CONSTANTE ALPHA\\
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_1__BJ__COEFFICIENT_GAMMA2 = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.8,
                fr = """Constante GAMMA2 du modele de deferlement de Battjes et Janssen.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 1 (BJ) CONSTANTE ALPHA
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1""",
                ang = """GAMMA1 constant of the Battjes and Janssen model.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE QB\\
DEFERLEMENT 1 (BJ) MODE DE CALCUL DE HM\\
DEFERLEMENT 1 (BJ) CHOIX FREQUENCE CARACTERISTIQUE\\
DEFERLEMENT 1 (BJ) CONSTANTE ALPHA\\
DEFERLEMENT 1 (BJ) CONSTANTE GAMMA1
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        b_DEPTH_INDUCED_BREAKING_DISSIPATIONH = BLOC(condition="DEPTH_INDUCED_BREAKING_DISSIPATION == 'Dissipation in accordance with Thornton et Guza (1983)'",
#       -----------------------------------
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_2__TG__WEIGHTING_FUNCTION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Weight function 1","Weight function 2"],
                defaut = ["Weight function 2"],
                fr = """Choix de l expression de la fonction de ponderation basee
sur une distribution de probabilite des hauteurs de houle.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 2 (TG) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 2 (TG) CONSTANTE B
DEFERLEMENT 2 (TG) CONSTANTE GAMMA""",
                ang = """Selection of the expression for the weighting function
based on a probability distribution of the wave heights.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 2 (TG) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT B\\
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT GAMMA
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_2__TG__CHARACTERISTIC_FREQUENCY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Frequency Fmoy","Frequency F01","Frequency F02","Frequency Fpic","Frequency Fread ordre 5","Frequency Fread ordre 8"],
                defaut = ["Frequency Fread ordre 5"],
                fr = """Choix de la frequence caracteristique du spectre de houle
\begin{itemize}
\item 1 : Frequence Fmoy
\item 2 : Frequence F01 (definie par les moments d ordre 0 et 1 du
spectre)
\item 3 : Frequence F02 (definie par les moments d ordre 0 et 2 du
spectre)
\item 4 : Frequence Fpic (frequence d echantillonage correspondant au
max)
\item 5 : Frequence Fread ordre 5 (frequence de pic methode Read ordre
5)
\item 6 : Frequence Fread ordre 8 (frequence de pic methode Read ordre
8)
\end{itemize}
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 2 (TG) FONCTION DE PONDERATION
DEFERLEMENT 2 (TG) CONSTANTE B
DEFERLEMENT 2 (TG) CONSTANTE GAMMA""",
                ang = """Selection of the characteristic frequency of the wave spectrum
\begin{itemize}
\item 1 : Frequency Fmoy
\item 2 : Frequency F01 (defined by the moments of order 0 and 1 of the
spectrum)
\item 3 : Frequency F02 (defined by the moments of order 0 and 2 of the
spectrum)
\item 4 : Frequency Fpic (sampling frequency corresponding to the max)
\item 5 : Frequency Fread ordre 5 (peak frequency, 5th order Read
method)
\item 6 : Frequency Fread ordre 8 (peak frequency, 8th order Read
method)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION
DEPTH-INDUCED BREAKING 2 (TG) WEIGHTING FUNCTION
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT B
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT GAMMA
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_2__TG__COEFFICIENT_B = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.0,
                fr = """Constante B du modele de deferlement de Thornton et Guza.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 2 (TG) FONCTION DE PONDERATION
DEFERLEMENT 2 (TG) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 2 (TG) CONSTANTE GAMMA""",
                ang = """Coefficient B of the Thornton and Guza model.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 2 (TG) WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 2 (TG) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT GAMMA
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_2__TG__COEFFICIENT_GAMMA = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.42,
                fr = """Constante GAMMA du modele de deferlement de Thornton et Guza.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 2 (TG) FONCTION DE PONDERATION
DEFERLEMENT 2 (TG) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 2 (TG) CONSTANTE B""",
                ang = """Coefficient GAMMA of the Thornton and Guza model.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 2 (TG) WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 2 (TG) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 2 (TG) COEFFICIENT B
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        b_DEPTH_INDUCED_BREAKING_DISSIPATIONI = BLOC(condition="DEPTH_INDUCED_BREAKING_DISSIPATION == 'Dissipation in accordance with Roelvink (1993)'",
#       -----------------------------------
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_3__RO__WAVE_HEIGHT_DISTRIBUTION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Weibull","Rayleigh"],
                defaut = ["Weibull"],
                fr = """Choix de la distribution des hauteurs de houle pour le
modele de deferlement de Roelvink :
   1...Weibull,
   2...Rayleigh.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 3 (RO) EXPOSANT FONCTION DE PONDERATION
DEFERLEMENT 3 (RO) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 3 (RO) CONSTANTE ALPHA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA2""",
                ang = """Selection of the wave height distribution for the
Roelvink breaking model :
   1...Weibull,
   2...Rayleigh.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 3 (RO) EXPONENT WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 3 (RO) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_3__RO__EXPONENT_WEIGHTING_FUNCTION = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 10,
                fr = """Exposant n de la fonction de ponderation utilisee par
le modele de deferlement de Roelvink.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 3 (RO) DISTRIBUTION DES HAUTEURS DE HOULE
DEFERLEMENT 3 (RO) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 3 (RO) CONSTANTE ALPHA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA2""",
                ang = """n exponent of the weighting function used in the Roelvink
breaking model.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 3 (RO) WAVE HEIGHT DISTRIBUTION\\
DEPTH-INDUCED BREAKING 3 (RO) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_3__RO__CHARACTERISTIC_FREQUENCY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Frequency Fmoy","Frequency F01","Frequency F02","Frequency Fpic","Frequency Fread ordre 5","Frequency Fread ordre 8"],
                defaut = ["Frequency Fread ordre 5"],
                fr = """Choix de la frequence caracteristique du spectre de houle
1 : Frequence Fmoy
2 : Frequence F01 (definie par les moments d ordre 0 et 1 du spectre)
3 : Frequence F02 (definie par les moments d ordre 0 et 2 du spectre)
4 : Frequence Fpic (frequence d echantillonage correspondant au max)
5 : Frequence Fread ordre 5 (frequence de pic methode Read ordre 5)
6 : Frequence Fread ordre 8 (frequence de pic methode Read ordre 8)
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 3 (RO) DISTRIBUTION DES HAUTEURS DE HOULE
DEFERLEMENT 3 (RO) EXPOSANT FONCTION DE PONDERATION
DEFERLEMENT 3 (RO) CONSTANTE ALPHA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA2""",
                ang = """Selection of the characteristic frequency of the wave spectrum
\begin{itemize}
\item 1 : Frequency Fmoy
\item 2 : Frequency F01 (defined by the moments of order 0 and 1 of the
spectrum)
\item 3 : Frequency F02 (defined by the moments of order 0 and 2 of the
spectrum)
\item 4 : Frequency Fpic (sampling frequency corresponding to the max)
\item 5 : Frequency Fread ordre 5 (peak frequency, 5th order Read
method)
\item 6 : Frequency Fread ordre 8 (peak frequency, 8th order Read
method)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 3 (RO) WAVE HEIGHT DISTRIBUTION\\
DEPTH-INDUCED BREAKING 3 (RO) EXPONENT WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_3__RO__COEFFICIENT_ALPHA = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.,
                fr = """Constante ALPHA du modele de deferlement de Roelvink (1993).
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 3 (RO) DISTRIBUTION DES HAUTEURS DE HOULE
DEFERLEMENT 3 (RO) EXPOSANT FONCTION DE PONDERATION
DEFERLEMENT 3 (RO) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 3 (RO) CONSTANTE GAMMA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA2""",
                ang = """Coefficient ALPHA of the Roelvink model (1993).
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 3 (RO) WAVE HEIGHT DISTRIBUTION\\
DEPTH-INDUCED BREAKING 3 (RO) EXPONENT WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 3 (RO) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_3__RO__COEFFICIENT_GAMMA = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.54,
                fr = """Constante GAMMA du modele de deferlement de Roelvink (1993).
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 3 (RO) DISTRIBUTION DES HAUTEURS DE HOULE
DEFERLEMENT 3 (RO) EXPOSANT FONCTION DE PONDERATION
DEFERLEMENT 3 (RO) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 3 (RO) CONSTANTE ALPHA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA2""",
                ang = """Coefficient GAMMA of the Roelvink model (1993).
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 3 (RO) WAVE HEIGHT DISTRIBUTION\\
DEPTH-INDUCED BREAKING 3 (RO) EXPONENT WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 3 (RO) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA2
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_3__RO__COEFFICIENT_GAMMA2 = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.65,
                fr = """Constante GAMMA2 du modele de deferlement de Roelvink (1993).
N est utilisee que pour la distribution de Weibull.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 3 (RO) DISTRIBUTION DES HAUTEURS DE HOULE
DEFERLEMENT 3 (RO) EXPOSANT FONCTION DE PONDERATION
DEFERLEMENT 3 (RO) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 3 (RO) CONSTANTE ALPHA
DEFERLEMENT 3 (RO) CONSTANTE GAMMA""",
                ang = """Coefficient GAMMA2 of the Roelvink model (1993).
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 3 (RO) WAVE HEIGHT DISTRIBUTION\\
DEPTH-INDUCED BREAKING 3 (RO) EXPONENT WEIGHTING FUNCTION\\
DEPTH-INDUCED BREAKING 3 (RO) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT ALPHA\\
DEPTH-INDUCED BREAKING 3 (RO) COEFFICIENT GAMMA
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        b_DEPTH_INDUCED_BREAKING_DISSIPATIONJ = BLOC(condition="DEPTH_INDUCED_BREAKING_DISSIPATION == 'Dissipation in accordance with Izumiya et Horikawa (1984)'",
#       -----------------------------------
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_4__IH__CHARACTERISTIC_FREQUENCY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'TXM',
                into = ["Frequency Fmoy","Frequency F01","Frequency F02","Frequency Fpic","Frequency Fread ordre 5","Frequency Fread ordre 8"],
                defaut = ["Frequency Fread ordre 5"],
                fr = """Choix de la frequence caracteristique du spectre de houle
1 : Frequence Fmoy
2 : Frequence F01 (definie par les moments d ordre 0 et 1 du spectre)
3 : Frequence F02 (definie par les moments d ordre 0 et 2 du spectre)
4 : Frequence Fpic (frequence d echantillonage correspondant au max)
5 : Frequence Fread ordre 5 (frequence de pic methode Read ordre 5)
6 : Frequence Fread ordre 8 (frequence de pic methode Read ordre 8)
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
DEFERLEMENT 4 (IH) CONSTANTE BETA0
DEFERLEMENT 4 (IH) CONSTANTE M2STAR""",
                ang = """Selection of the characteristic frequency of the wave spectrum
\begin{itemize}
\item 1 : Frequency Fmoy
\item 2 : Frequency F01 (defined by the moments of order 0 and 1 of the
spectrum)
\item 3 : Frequency F02 (defined by the moments of order 0 and 2 of the
spectrum)
\item 4 : Frequency Fpic (sampling frequency corresponding to the max)
\item 5 : Frequency Fread ordre 5 (peak frequency, 5th order Read
method)
\item 6 : Frequency Fread ordre 8 (peak frequency, 8th order Read
method)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
DEPTH-INDUCED BREAKING 4 (IH) COEFFICIENT BETA0\\
DEPTH-INDUCED BREAKING 4 (IH) COEFFICIENT M2STAR
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_4__IH__COEFFICIENT_BETA0 = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 1.8,
                fr = """Constante BETA0 du modele de deferlement de Izumiya et
Horikawa (1984).
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 4 (IH) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 4 (IH) CONSTANTE M2STAR""",
                ang = """coefficient BETA0 of the Izumiya and Horikawa
model (1984).
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 4 (IH) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 4 (IH) COEFFICIENT M2STAR
\end{CommentBlock}""",
            ),
#           -----------------------------------
            DEPTH_INDUCED_BREAKING_4__IH__COEFFICIENT_M2STAR = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.009,
                fr = """Constante M2STAR du modele de deferlement de Izumiya et
Horikawa (1984).
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT
DEFERLEMENT 4 (IH) CHOIX FREQUENCE CARACTERISTIQUE
DEFERLEMENT 4 (IH) CONSTANTE BETA0""",
                ang = """coefficient M2STAR of the Izumiya and Horikawa
model (1984).
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS\\
DEPTH-INDUCED BREAKING 4 (IH) CHARACTERISTIC FREQUENCY\\
DEPTH-INDUCED BREAKING 4 (IH) COEFFICIENT BETA0
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        NUMBER_OF_BREAKING_TIME_STEPS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Nombre de sous-pas de temps pour la prise en compte de la
dissipation d energie par deferlement. Ces sous-pas de temps sont
en progression geometrique.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
COEFFICIENT POUR LES SOUS-PAS DE TEMPS POUR LE DEFERLEMENT""",
            ang = """Number of time steps for the breaking source term.
These time steps are in a geometric progression
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
COEFFICIENT FOR THE BREAKING TIME STEPS
\end{CommentBlock}""",
        ),
#       -----------------------------------
        MAXIMUM_VALUE_OF_THE_RATIO_HM0_ON_D = SIMP(statut ='o',
#       -----------------------------------
            typ = 'R',
            defaut = 1.,
            fr = """En debut de prise en compte des termes sources, la hauteur
de houle est ecretee de facon a satisfaire le critere specifie.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT""",
            ang = """At the beginning of the integration of the source terms,
the wave height is lopped in order to satisfy the specified
criterium.
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION
\end{CommentBlock}""",
        ),
#       -----------------------------------
        COEFFICIENT_OF_THE_TIME_SUB_INCREMENTS_FOR_BREAKING = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.45,
            fr = """Raison de la suite geometrique des sous-pas de temps pour
le deferlement.
**Mots-cles associes :**
DISSIPATION PAR DEFERLEMENT
NOMBRE DE SOUS-PAS DE TEMPS POUR LE DEFERLEMENT""",
            ang = """Geometrical ratio of the time sub-increments for the
depth-induced breaking
\\
 \begin{CommentBlock}{Related keywords}
DEPTH-INDUCED BREAKING DISSIPATION\\
NUMBER OF BREAKING TIME STEPS
\end{CommentBlock}""",
        ),
    ),
#   -----------------------------------
    LIMITER = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WAVE_GROWTH_LIMITER = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["no wave growth limiter","WAM 4 original limiter","Hersbach et Janssen (1999) limiter","Laugel-BAJ limiter"],
            defaut = ["WAM 4 original limiter"],
            fr = """Choix du type de limiteur de croissance.
Si la valeur est 0, pas de limiteur.
Si la valeur est 1, limiteur type WAM 4 original.
Si la valeur est 2, limiteur de Hersbach et Janssen (1999).
Si la valeur est 3, limiteur de BAJ Laugel.
**Mots-cles associes :**
PRISE EN COMPTE DES TERMES SOURCES""",
            ang = """Choice of the wave growth limiter.
\begin{itemize}
\item If LIMIT=0, no wave growth limiter.
\item If LIMIT=1, WAM 4 original limiter.
\item If LIMIT=2, Hersbach et Janssen (1999) limiter.
\item If LIMIT=3, Laugel BAJ limiter.
\end{itemize}
 \begin{CommentBlock}{Related keywords}
CONSIDERATION OF SOURCE TERMS
\end{CommentBlock}""",
        ),
    ),
#   -----------------------------------
    TRIAD_TRANSFERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        TRIAD_INTERACTIONS = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM',
            into = ["no triad interactions","LTA model (Eldeberky, 1995)","SPB model (Becq, 1998)"],
            defaut = ["no triad interactions"],
            fr = """Choix du type de modelisation du terme de transfert non
lineaire entre triplets de frequences. Si sa valeur est 0, on ne
prend pas en compte les transferts non lineaires entre triplets
de frequences, si sa valeur est 1, ils sont integres selon le
modele LTA, si sa valeur est 2, ils sont integres selon le
modele SPB.
**Mots-cles associes :**
TRIADS 1 (LTA) CONSTANTE ALPHA
TRIADS 1 (LTA) CONSTANTE RFMLTA
TRIADS 2 (SPB) CONSTANTE K
TRIADS 2 (SPB) BORNE DIRECTIONNELLE INFERIEURE
TRIADS 2 (SPB) BORNE DIRECTIONNELLE SUPERIEURE""",
            ang = """Selection of the triad interaction model:
\begin{itemize}
\item 0 : no triad interactions
\item 1 : LTA model (Eldeberky, 1996)
\item 2 : SPB model (Becq, 1998)
\end{itemize}
 \begin{CommentBlock}{Related keywords}
TRIADS 1 (LTA) COEFFICIENT ALPHA\\
TRIADS 1 (LTA) COEFFICIENT RFMLTA\\
TRIADS 2 (SPB) COEFFICIENT K\\
TRIADS 2 (SPB) LOWER DIRECTIONAL BOUNDARY\\
TRIADS 2 (SPB) UPPER DIRECTIONAL BOUNDARY\\
\end{CommentBlock}""",
        ),
#       -----------------------------------
        b_TRIAD_INTERACTIONSG = BLOC(condition="TRIAD_INTERACTIONS == 'LTA model (Eldeberky, 1996)'",
#       -----------------------------------
#           -----------------------------------
            TRIADS_1__LTA__COEFFICIENT_ALPHA = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.5,
                fr = """Constante alpha du modele LTA propose par Eldeberky (1996).
Si alpha = 0, les transferts d energie entre frequences seront nuls.
L intensite de ces transferts augmente avec la valeur de alpha.
**Mots-cles associes :**
TRANSFERTS ENTRE TRIPLETS DE FREQUENCES
TRIADS 1 (LTA) CONSTANTE RFMLTA""",
                ang = """Coefficient alpha of the LTA model proposed by Eldeberky(1996).
If alpha=0, no energy transfers. The energy transfers increase
with alpha.
\\
 \begin{CommentBlock}{Related keywords}
TRIAD INTERACTIONS\\
TRIADS 1 (LTA) COEFFICIENT RFMLTA
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TRIADS_1__LTA__COEFFICIENT_RFMLTA = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 2.5,
                fr = """RFMLTA determine la frequence maximale vers laquelle les
transferts d energie peuvent avoir lieu. La frequence maximale est
calculee comme le produit entre la constante RFMLTA et la frequence
de pic du spectre.
**Mots-cles associes :**
TRANSFERTS ENTRE TRIPLETS DE FREQUENCES
TRIADS 1 (LTA) CONSTANTE ALPHA""",
                ang = """RFMLTA determines the upper frequency on which the energy
transfers may occur. The maximal frequency is calculated as the
product of the constant RFMLTA by the peak frequency of the spectrum.
\\
 \begin{CommentBlock}{Related keywords}
TRIAD INTERACTIONS\\
TRIADS 1 (LTA) COEFFICIENT ALPHA
\end{CommentBlock}""",
            ),
        ),
#       -----------------------------------
        b_TRIAD_INTERACTIONSH = BLOC(condition="TRIAD_INTERACTIONS == 'SPB model (Becq, 1998)'",
#       -----------------------------------
#           -----------------------------------
            TRIADS_2__SPB__COEFFICIENT_K = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.34,
                fr = """Parametre d ajustement du modele SPB
**Mots-cles associes :**
TRANSFERTS ENTRE TRIPLETS DE FREQUENCES
TRIADS 2 (SPB) BORNE DIRECTIONNELLE INFERIEURE
TRIADS 2 (SPB) BORNE DIRECTIONNELLE SUPERIEURE""",
                ang = """coefficient K of the SPB model
\\
 \begin{CommentBlock}{Related keywords}
TRIAD INTERACTIONS\\
TRIADS 2 (SPB) LOWER DIRECTIONAL BOUNDARY\\
TRIADS 2 (SPB) UPPER DIRECTIONAL BOUNDARY
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TRIADS_2__SPB__LOWER_DIRECTIONAL_BOUNDARY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Borne directionnelle inferieure du modele SPB
**Mots-cles associes :**
TRANSFERTS ENTRE TRIPLETS DE FREQUENCES
TRIADS 2 (SPB) CONSTANTE K
TRIADS 2 (SPB) BORNE DIRECTIONNELLE SUPERIEURE""",
                ang = """Lower directional boundary of the SPB model
\\
 \begin{CommentBlock}{Related keywords}
TRIAD INTERACTIONS\\
TRIADS 2 (SPB) COEFFICIENT K\\
TRIADS 2 (SPB) UPPER DIRECTIONAL BOUNDARY
\end{CommentBlock}""",
            ),
#           -----------------------------------
            TRIADS_2__SPB__UPPER_DIRECTIONAL_BOUNDARY = SIMP(statut ='o',
#           -----------------------------------
                typ = 'R',
                defaut = 360.,
                fr = """Borne directionnelle superieure du modele SPB
**Mots-cles associes :**
TRANSFERTS ENTRE TRIPLETS DE FREQUENCES
TRIADS 2 (SPB) CONSTANTE K
TRIADS 2 (SPB) BORNE DIRECTIONNELLE INFERIEURE""",
                ang = """Upper directional boundary of the SPB model
\\
 \begin{CommentBlock}{Related keywords}
TRIAD INTERACTIONS\\
TRIADS 2 (SPB) COEFFICIENT K\\
TRIADS 2 (SPB) LOWER DIRECTIONAL BOUNDARY
\end{CommentBlock}""",
            ),
        ),
    ),
#   -----------------------------------
    VEGETATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        VEGETATION_TAKEN_INTO_ACCOUNT = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si oui, on appelle le sous-programme QVEG, dans lequel des
donnees sur la vegetation doivent etre renseignees""",
            ang = """If YES, subroutine QVEG will be called, it contains data
on vegetation that are case-specific and must thus be modified""",
        ),
    ),
#   -----------------------------------
    POROUS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        POROUS_MEDIA = SIMP(statut ='o',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si oui, on appelle le sous-programme QPOROS, dans lequel des
donnees sur le milieu poreux doivent etre renseignees""",
            ang = """If YES, subroutine QPOROS will be called, it contains data
on POROUS MEDIA that are case-specific and must thus be modified""",
        ),
    ),
)
# -----------------------------------------------------------------------
INITIAL_CONDITIONS = PROC(nom= "INITIAL_CONDITIONS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    TYPE_OF_INITIAL_DIRECTIONAL_SPECTRUM = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        into = ["non-existent spectrum","JONSWAP spectrum 1","JONSWAP spectrum 2","JONSWAP spectrum 3","JONSWAP spectrum 4","JONSWAP spectrum 5","JONSWAP spectrum 6","TMA spectrum"],
        defaut = ["non-existent spectrum"],
        fr = """Si ce mot-cle est pris egal a 0, on specifie un spectre
initial nul. Si il est pris egal entre 1 et 6
un spectre de type JONSWAP est specifie du champ de vent initial
et ou des valeurs des autres mots cles. Si il est pris egal a 7
un spectre de type TMA parametre est specifie.
**Mots-cles associes :**
HAUTEUR SIGNIFICATIVE INITIALE
FREQUENCE DE PIC INITIALE
FACTEUR DE PIC INITIAL
VALEUR INITIALE DE SIGMA-A POUR SPECTRE
VALEUR INITIALE DE SIGMA-B POUR SPECTRE
CONSTANTE DE PHILLIPS INITIALE
VALEUR MOYENNE DU FETCH INITIAL
FREQUENCE DE PIC MAXIMALE
DIRECTION PRINCIPALE 1 INITIALE
ETALEMENT DIRECTIONNEL 1 INITIAL
DIRECTION PRINCIPALE 2 INITIALE
ETALEMENT DIRECTIONNEL 2 INITIAL
FACTEUR DE PONDERATION POUR FRA INITIALE""",
        ang = """If this keyword is set to 0, a non-existent spectrum is speci
fied at the initial time step. If it ranges from 1 to 7, a
JONSWAP (or TMA)-typed spectrum is specified as a
function of the initial wind field and/or of the values of the
following keywords
\\
 \begin{CommentBlock}{Related keywords}
INITIAL SIGNIFICANT WAVE HEIGHT\\
INITIAL PEAK FREQUENCY\\
INITIAL PEAK FACTOR\\
INITIAL VALUE OF SIGMA-A FOR SPECTRUM\\
INITIAL VALUE OF SIGMA-B FOR SPECTRUM\\
INITIAL PHILLIPS CONSTANT\\
INITIAL MEAN FETCH VALUE\\
INITIAL MAXIMUM PEAK FREQUENCY\\
INITIAL MAIN DIRECTION 1\\
INITIAL DIRECTIONAL SPREAD 1\\
INITIAL MAIN DIRECTION 2\\
INITIAL DIRECTIONAL SPREAD 2\\
INITIAL WEIGHTING FACTOR FOR ADF
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_ANGULAR_DISTRIBUTION_FUNCTION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["$cos{^2s}(T-T0)$, T in [T0-pi/2,T0+pi/2]","$exp(-0.5((T-T0)/s)^2)$, T in [T0-pi/2,T0+pi/2]","$cos{^2s}((T-T0)/2)$ (de type Mitsuyasu)","Goda Formula allowing a frequence dependancy"],
        defaut = ["$cos{^2s}(T-T0)$, T in [T0-pi/2,T0+pi/2]"],
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel. Permet le calcul de
la fonction de repartition angulaire bimodale pour une serie de
directions
1 : modele en $cos{^2s}(T-T0)$ ; T dans [T0-pi/2;T0+pi/2]
2 : modele en $exp(-0.5((T-T0)/s)^2)$ ; T dans [T0-pi/2;T0+pi/2]
3 : modele en $cos{^2s}((T-T0)/2)$ (de type Mitsuyasu)
4 : Formule de Goda : Mitsuyasu avec $s=(F/fp)^5$ or $s=(F/fp)^{-2.5}$
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL INITIAL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum. Allow the computation of the
angular distribution function
\begin{itemize}
\item 1 : $cos^{2s}(T-T0)$ ; with T in [T0-pi/2;T0+pi/2]
\item 2 : $exp(-0.5((T-T0)/s)^2)$ ; with T in [T0-pi/2;T0+pi/2]
\item 3 : $cos^{2s}((T-T0)/2)$ (of type Mitsuyasu)
\item 4 : Mitsuyasu with $s=(F/fp)^5$ or $s=(F/fp)^{-2.5}$
\end{itemize}
where $s$ is the boundary directionnal spread (\telkey{SPRED1}
or \telkey{SPRED2})
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_SIGNIFICANT_WAVE_HEIGHT = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 1.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL INITIAL""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_PEAK_FREQUENCY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.067,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_PEAK_FACTOR = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 3.3,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_VALUE_OF_SIGMA_A_FOR_SPECTRUM = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.07,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_VALUE_OF_SIGMA_B_FOR_SPECTRUM = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.09,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_PHILLIPS_CONSTANT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.0081,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initiale directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_MEAN_FETCH_VALUE = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 30000.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_MAXIMUM_PEAK_FREQUENCY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.2,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_MAIN_DIRECTION_1 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial. Exprimee en degres
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_DIRECTIONAL_SPREAD_1 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 2.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_MAIN_DIRECTION_2 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial. Exprimee en degres.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_DIRECTIONAL_SPREAD_2 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 2.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_WEIGHTING_FACTOR_FOR_ADF = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 1.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial.
**Mots-cles associes :**
INITIALISATION DU SPECTRE DIRECTIONNEL""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF INITIAL DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    GODA_COEFFICIENT_FOR_ANGULAR_SPEADING = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 1.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
l initialisation du spectre directionnel en fonction du champ
de vent initial. Coefficient intervenant dans la loi de repartition
angulaire. Selon Goda 10(Wind waves), 25 swell short decay, 75 long deca
mais des valeurs autres sont acceptees.
**Mots-cles associes :**
FONCTION DE REPARTITION ANGULAIRE INITIALE""",
        ang = """Is part of the set of constants used for computing the
initial directional spectrum as a function of the wind field.
This coefficient is in the formula of Goda formulation.
It should be 10Wind waves, 25 swell short decay, 75 long decay
but other values are accepted.
\\
 \begin{CommentBlock}{Related keywords}
INITIAL ANGULAR DISTRIBUTION FUNCTION
\end{CommentBlock}""",
    ),
#   -----------------------------------
    INITIAL_STILL_WATER_LEVEL = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Parametre permettant de calculer la profondeur
initiale du plan d eau (DEPTH) a partir de la cote du
fond (ZF) : DEPTH=ZREPOS-ZF.""",
        ang = """Parameter used in the computation of the initial
water DEPTH : DEPTH=ZREPOS-ZF.""",
    ),
)
# -----------------------------------------------------------------------
INTERNAL = PROC(nom= "INTERNAL",op = None,
# -----------------------------------------------------------------------
    UIinfo = {"groupes": ("CACHE")},
#   -----------------------------------
    VECTOR_LENGTH = SIMP(statut ='o',
#   -----------------------------------
        typ = 'I',
        defaut = [1],
        fr = """Indique la longueur du vecteur de la machine vectorielle
utilisee.""",
        ang = """Indicates the vector length of the vectorial machine
being used.""",
    ),
#   -----------------------------------
    STEERING_FILE = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Nom du fichier contenant les parametres du calcul a realiser.""",
        ang = """Name of the file containing the parameters of the computation
to be made.""",
    ),
#   -----------------------------------
    DICTIONARY = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = 'tomawac.dico',
        fr = """Dictionnaire des mots cles.""",
        ang = """Key word dictionary.""",
    ),
#   -----------------------------------
    BOUNDARY_CONDITION_BEFORE_TIME_STEP = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Indique si on tient compte de la condition à la limite
avant ou apres le pas de temps""",
        ang = """Indicates whether boundary condition is taken into account
before or after time step""",
    ),
#   -----------------------------------
    ECRET_FOR_SMALL_HEIGHT = SIMP(statut ='f',
#   -----------------------------------
        typ = bool,
        defaut = [True ],
        fr = """Indique si on ecrete la hauteur de vague pour les
hauteurs d eau quasi nulle""",
        ang = """Indicates whether wave height is null for very
small heigth of water""",
    ),
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
)
# -----------------------------------------------------------------------
BOUNDARY_CONDITIONS = PROC(nom= "BOUNDARY_CONDITIONS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    TYPE_OF_BOUNDARY_DIRECTIONAL_SPECTRUM = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        into = ["non-existent spectrum","JONSWAP spectrum 1","JONSWAP spectrum 2","JONSWAP spectrum 3","JONSWAP spectrum 4","JONSWAP spectrum 5","JONSWAP spectrum 6","TMA spectrum"],
        defaut = ["non-existent spectrum"],
        fr = """Si ce mot-cle est pris egal a 0, on specifie un spectre
nul aux limites. Si il est pris egal entre 1 et 6
un spectre de type JONSWAP est specifie en tout point des limites
d entree en fonction du champ de vent aux limites et ou des
valeurs des autres mots cles. Si il est pris egal a 7 un spectre de
type TMA parametre est specifie aux limites.
**Mots-cles associes :**
HAUTEUR SIGNIFICATIVE AUX LIMITES
FREQUENCE DE PIC AUX LIMITES
FACTEUR DE PIC AUX LIMITES
VALEUR AUX LIMITES DE SIGMA-A POUR SPECTRE
VALEUR AUX LIMITES DE SIGMA-B POUR SPECTRE
CONSTANTE DE PHILLIPS AUX LIMITES
VALEUR MOYENNE DU FETCH AUX LIMITES
FREQUENCE DE PIC MAXIMALE
DIRECTION PRINCIPALE 1 AUX LIMITES
ETALEMENT DIRECTIONNEL 1 AUX LIMITES
DIRECTION PRINCIPALE 2 AUX LIMITES
ETALEMENT DIRECTIONNEL 2 AUX LIMITES
FACTEUR DE PONDERATION POUR FRA AUX LIMITES""",
        ang = """If this keyword is set to 0, a non-existent spectrum is speci
fied at the inlet boundaries of the domain. If it ranges from 1 to 7, a
JONSWAP (or TMA) -typed spectrum is specified at these very points as a
function of the initial wind field and/or of the values of the following
keywords
\\
 \begin{CommentBlock}{Related keywords}
BOUNDARY SIGNIFICANT WAVE HEIGHT\\
BOUNDARY PEAK FREQUENCY\\
BOUNDARY PEAK FACTOR\\
BOUNDARY VALUE OF SIGMA-A FOR SPECTRUM\\
BOUNDARY VALUE OF SIGMA-B FOR SPECTRUM\\
BOUNDARY PHILLIPS CONSTANT\\
BOUNDARY MEAN FETCH VALUE\\
BOUNDARY MAXIMUM PEAK FREQUENCY\\
BOUNDARY MAIN DIRECTION 1\\
BOUNDARY DIRECTIONAL SPREAD 1\\
BOUNDARY MAIN DIRECTION 2\\
BOUNDARY DIRECTIONAL SPREAD 2\\
BOUNDARY WEIGHTING FACTOR FOR ADF\\
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_ANGULAR_DISTRIBUTION_FUNCTION = SIMP(statut ='f',
#   -----------------------------------
        typ = 'TXM',
        into = ["$cos{^2s}(T-T0)$, T in [T0-pi/2,T0+pi/2]","$exp(-0.5((T-T0)/s)^2)$, T in [T0-pi/2,T0+pi/2]","$cos{^2s}((T-T0)/2)$ (de type Mitsuyasu)","Goda Formula allowing a frequence dependancy"],
        defaut = ["$cos{^2s}(T-T0)$, T in [T0-pi/2,T0+pi/2]"],
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites. Permet le calcul de
la fonction de repartition angulaire bimodale pour une serie de
directions
1 : modele en $cos{^2s}(T-T0)$ ; T dans [T0-pi/2;T0+pi/2]
2 : modele en $exp(-0.5((T-T0)/s)^2)$ ; T dans [T0-pi/2;T0+pi/2]
3 : modele en $cos{^2s}((T-T0)/2)$ (de type Mitsuyasu)
4 : Formule de Goda : Mitsuyasu avec $s=(F/fp)^5$ or $s=(F/fp)^{-2.5}$
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum. Allow the computation of the
angular distribution function
\begin{itemize}
\item 1 : $cos^{2s}(T-T0)$ ; with T in [T0-pi/2;T0+pi/2]
\item 2 : $exp(-0.5((T-T0)/s)^2)$ ; with T in [T0-pi/2;T0+pi/2]
\item 3 : $cos^{2s}((T-T0)/2)$ (of type Mitsuyasu)
\item 4 : Mitsuyasu with $s=(F/fp)^5$ or $s=(F/fp)^{-2.5}$
\end{itemize}
where $s$ is the boundary directionnal spread (\telkey{SPRE1L}
or \telkey{SPRE2L})
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_SIGNIFICANT_WAVE_HEIGHT = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 1.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent.
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_PEAK_FREQUENCY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.067,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_SPECTRUM_VALUE_OF_SIGMA_A = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.07,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_SPECTRUM_VALUE_OF_SIGMA_B = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.09,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_PHILLIPS_CONSTANT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.0081,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_MEAN_FETCH_VALUE = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 30000.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_MAXIMUM_PEAK_FREQUENCY = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.2,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_MAIN_DIRECTION_1 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_DIRECTIONAL_SPREAD_1 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 2.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_MAIN_DIRECTION_2 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 0.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_DIRECTIONAL_SPREAD_2 = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 2.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_WEIGHTING_FACTOR_FOR_ADF = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 1.,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    BOUNDARY_PEAK_FACTOR = SIMP(statut ='o',
#   -----------------------------------
        typ = 'R',
        defaut = 3.3,
        fr = """Fait partie de l ensemble des constantes utilisees dans
le calcul du spectre directionnel aux limites en fonction du champ
de vent .
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Is part of the set of constants used for computing the
boundary directional spectrum as a function of the wind field.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
#   -----------------------------------
    LIMIT_SPECTRUM_MODIFIED_BY_USER = SIMP(statut ='o',
#   -----------------------------------
        typ = bool,
        defaut = [False],
        fr = """Indique si l''utilisateur desire modifier le spectre
aux limites. Il doit alors rapatrier chez lui le sous programme
limwac.f dans le cas d une discretisation frequentielle du
spectre et le sous programme spelim.f sinon.
**Mots-cles associes :**
TYPE DE SPECTRE DIRECTIONNEL AUX LIMITES""",
        ang = """Indicates whether the user wants to modify the boundary
spectrum. He should then retrieve the limwac.f subroutine, if the
spectrum is frequency discretized, or the spelim.f subroutine,
otherwise.
\\
 \begin{CommentBlock}{Related keywords}
TYPE OF BOUNDARY DIRECTIONAL SPECTRUM
\end{CommentBlock}""",
    ),
)
# -----------------------------------------------------------------------
TRANSPORT = PROC(nom= "TRANSPORT",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    DIFFRACTION_PARAMETERS = FACT(statut='f',
#   -----------------------------------
#       -----------------------------------
        DIFFRACTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["Diffraction is not taken into account","Mild Slope Equation model (Berkhoff - 1972)","Revised Mild Slope Equation model (Porter - 2003)"],
            defaut = ["Diffraction is not taken into account"],
            fr = """Choix du model pour representer la diffraction :
 0 : Pas de prise en compte de la diffraction
 1 : Mild Slope Equation de Berkhoff (1972)
 2 : Revised Mild Slope Equation de Porter (2003)
La diffraction est modelisee suivant l approche proposee
par Holthuijsen (2003)
**Mots-cles associes :**
PAS DE TEMPS DEBUT DIFFRACTION
SEUIL DE VARIANCE CONSIDEREE POUR DIFFRACTION
FILTRE POUR DIFFRACTION""",
            ang = """Caution : We do not guarantee the modele of diffraction.
Selection of the model used to represent the diffraction :
\begin{itemize}
\item 0 : Diffraction is not taken into account
\item  1 : Mild Slope Equation model (Berkhoff - 1972)
\item  2 : Revised Mild Slope Equation model (Porter - 2003)
\end{itemize}
The phase-decoupled approach proposed by Holthuijsen (2003) is
used to simulate diffraction in TOMAWAC
\\
 \begin{CommentBlock}{Related keywords}
STARTING TIME STEP FOR DIFFRACTION\\
VARIANCE THRESHOLD FOR DIFFRACTION\\
DIFFRACTION FILTER
\end{CommentBlock}""",
        ),
#       -----------------------------------
        STARTING_TIME_STEP_FOR_DIFFRACTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 1,
            fr = """Numero du pas de temps a partir duquel la diffraction est
prise en compte dans la simulation.
**Mots-cles associes :**
DIFFRACTION
SEUIL DE VARIANCE CONSIDEREE POUR DIFFRACTION
FILTRE POUR DIFFRACTION""",
            ang = """Number of the time step from which the diffraction
is taken into account until the end of the simulation.
\\
 \begin{CommentBlock}{Related keywords}
DIFFRACTION\\
VARIANCE THRESHOLD FOR DIFFRACTION\\
DIFFRACTION FILTER
\end{CommentBlock}""",
        ),
#       -----------------------------------
        VARIANCE_THRESHOLD_FOR_DIFFRACTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-12,
            fr = """Seuil minimum de variance spectrale pris en compte
dans la diffraction
**Mots-cles associes :**
DIFFRACTION
PAS DE TEMPS DEBUT DIFFRACTION
FILTRE POUR DIFFRACTION""",
            ang = """Minimum spectral variance threshold taken into account
when diffraction is considered
\\
 \begin{CommentBlock}{Related keywords}
DIFFRACTION\\
STARTING TIME STEP FOR DIFFRACTION\\
DIFFRACTION FILTER
\end{CommentBlock}""",
        ),
#       -----------------------------------
        OPTION_FOR_SECOND_DERIVATIVES = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = [1],
            fr = """1 : methode Freemesh 2 : deux derivees simples""",
            ang = """1: Freemesh method 2: two simple derivatives""",
        ),
#       -----------------------------------
        DIFFRACTION_FILTER = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = [False],
            fr = """Si la diffraction est prise en compte, le logique
indique si les amplitudes calculees sont filtrees
pour le calcul du parametre de diffraction et des vitesses
de transfert.
**Mots-cles associes :**
DIFFRACTION
SEUIL DE VARIANCE CONSIDEREE POUR DIFFRACTION
PAS DE TEMPS DEBUT DIFFRACTION""",
            ang = """If diffraction is considered, the keyword indicates
whether the local amplitudes of the directional spectra are
filtered to compute the diffraction parameter and the transfer
rates.
\\
 \begin{CommentBlock}{Related keywords}
DIFFRACTION\\
VARIANCE THRESHOLD FOR DIFFRACTION\\
STARTING TIME STEP FOR DIFFFRACTION
\end{CommentBlock}""",
        ),
    ),
)
TEXTE_NEW_JDC = "\
COMPUTATION_ENVIRONMENT();\
GENERAL_PARAMETERS();\
SOURCE_TERMS();\
INITIAL_CONDITIONS();\
BOUNDARY_CONDITIONS();\
"
Ordre_Des_Commandes = (
'COMPUTATION_ENVIRONMENT',
'GENERAL_PARAMETERS',
'SOURCE_TERMS',
'INITIAL_CONDITIONS',
'INTERNAL',
'BOUNDARY_CONDITIONS',
'TRANSPORT')
try:
    import TelApy
    source = "eficas"
except Exception as excpt:
    source = "Telemac"
enum = source+'.tomawac_enum_auto'
dicoCasEn = source+'.tomawac_dicoCasEnToCata'
dicoCasFr = source+'.tomawac_dicoCasFrToCata'
