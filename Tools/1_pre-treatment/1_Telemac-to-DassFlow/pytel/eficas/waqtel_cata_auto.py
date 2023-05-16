
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



JdC = JDC_CATA (code = 'WAQTEL',
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
        WAQ_CASE_TITLE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            defaut = '',
            fr = """Titre du cas etudie. Ce titre sera inscrit dans les sorties.""",
            ang = """Title of the case being considered.
This title will be marked on the printouts.""",
        ),
#       -----------------------------------
        PARALLEL_PROCESSORS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'I',
            defaut = 0,
            fr = """Nombre de processeurs en calcul parallele :
\begin{itemize}
\item 0 : 1 machine, compilation sans bibliotheque de parallelisme ;
\item 1 : 1 machine, compilation avec bibliotheque de parallelisme ;
\item 2 : 2 processeurs ou machines en parallele etc.
\end{itemize}""",
            ang = """Number of processors for parallel processing:
\begin{itemize}
\item 0: 1 machine, compiling without parallel library,
\item 1: 1 machine, compiling with a parallel library,
\item 2: 2 processors or machines in parallel etc.
\end{itemize}""",
        ),
    ),
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
                defaut = '',
                fr = """Fichier de geometrie, pareil que celui de TELEMAC.""",
                ang = """Geometry file same as the TELEMAC one.""",
            ),
#           -----------------------------------
            GEOMETRY_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier de geometrie.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Geometry file format.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            BOUNDARY_CONDITIONS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                fr = """Nom du fichier contenant les types de conditions aux limites.
Ce fichier est rempli de facon automatique par le mailleur au moyen de
couleurs affectees aux noeuds des frontieres du domaine de calcul.""",
                ang = """Name of the file containing the types of boundary conditions.
This file is filled automatically by the mesh generator through
colours that are assigned to the boundary nodes.""",
            ),
#           -----------------------------------
            FORTRAN_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Nom du fichier FORTRAN a soumettre.""",
                ang = """Name of FORTRAN file to be submitted.""",
            ),
#           -----------------------------------
            HYDRODYNAMIC_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'), max='**',
                defaut = '',
                fr = """Fichier des donnees hydrodynamiques provenant de TELEMAC.
Lu mais pas utilise pour le moment.""",
                ang = """Hydrodynamic data file coming from TELEMAC.
Read but not used at the moment.""",
            ),
#           -----------------------------------
            HYDRODYNAMIC_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier hydrodynamique.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Hydrodynamic file format.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            VALIDATION = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Option utilisee principalement pour le dossier de validation.
Si ce mot-cle vaut OUI, les resultats du calcul vont alors etre
compares aux valeurs du fichier de reference.
Lu mais pas utilise pour le moment.""",
                ang = """This option is primarily used for the validation
documents. If this keyword is equal to YES, the REFERENCE FILE
is then considered as a reference which the computation is
going to be compared with.
Read but not used at the moment.""",
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
                    fr = """Format du fichier hydrodynamique.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                    ang = """Hydrodynamic file format.
Possible choices are:
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
                    fr = """Nom du fichier servant a valider le calcul.
Si \telkey{VALIDATION} = OUI, les resultats du calcul vont etre
comparees aux valeurs contenues dans ce fichier.
La comparaison est effectuee par le sous-programme VALIDA.
a implementer.""",
                    ang = """Name of the file used to validate the computation.
If \telkey{VALIDATION} = YES, the results of the computation will
be compared with the values of this file. The comparison is
made by the subroutine BIEF\_VALIDA. (not implemented yet).""",
                ),
            ),
        ),
#       -----------------------------------
        AED2 = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            AED2_STEERING_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les parametres AED2 du calcul
QE a realiser.""",
                ang = """Name of the file containing AED2 parameters of the WAQ
computation.""",
            ),
#           -----------------------------------
            AED2_PHYTOPLANKTON_STEERING_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les parametres phytoplancton AED2 du
calcul QE a realiser.""",
                ang = """Name of the file containing AED2 phytoplankton parameters of
the WAQ computation.""",
            ),
#           -----------------------------------
            AED2_ZOOPLANKTON_STEERING_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les parametres zooplancton AED2 du
calcul QE a realiser.""",
                ang = """Name of the file containing AED2 zooplankton parameters of the
WAQ computation.""",
            ),
#           -----------------------------------
            AED2_PATHOGEN_STEERING_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les parametres pathogenes AED2 du
calcul QE a realiser.""",
                ang = """Name of the file containing AED2 pathogen parameters of the
WAQ computation.""",
            ),
#           -----------------------------------
            AED2_BIVALVE_STEERING_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)'),
                defaut = '',
                fr = """Nom du fichier contenant les parametres bivalves AED2 du
calcul QE a realiser.""",
                ang = """Name of the file containing AED2 bivalve parameters of the
WAQ computation.""",
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
            RESULTS_FILE = SIMP(statut ='f',
#           -----------------------------------
                typ = ('Fichier','All Files (*)','Sauvegarde'), max='**',
                fr = """Nom du fichier dans lequel seront ecrits les resultats avec
une periodicite donnee par le mot cle
\telkey{PERIODE POUR LES SORTIES QUALITE D EAU}.""",
                ang = """Name of the file into which the computation results will be
written, the periodicity being given by the keyword
\telkey{WATER QUALITY PRINTOUT PERIOD}.""",
            ),
#           -----------------------------------
            RESULTS_FILE_FORMAT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM',
                into = ['SERAFIN','SERAFIND','MED'],
                defaut = 'SERAFIN',
                fr = """Format du fichier des resultats.
Les valeurs possibles sont :
\begin{itemize}
\item SERAFIN : format standard simple precision pour \tel ;
\item SERAFIND: format standard double precision pour \tel ;
\item MED     : format MED double precision base sur HDF5.
\end{itemize}""",
                ang = """Results file format.
Possible choices are:
\begin{itemize}
\item SERAFIN : classical single precision format in \tel,
\item SERAFIND: classical double precision format in \tel,
\item MED     : MED double precision format based on HDF5.
\end{itemize}""",
            ),
#           -----------------------------------
            WATER_QUALITY_PRINTOUT_PERIOD = SIMP(statut ='o',
#           -----------------------------------
                typ = 'I',
                defaut = 1,
                fr = """Periode pour les sorties graphiques QE.
Lu mais pas utilise pour le moment.""",
                ang = """Graphic outputs period for WAQ.
Read but not used at the moment.""",
            ),
#           -----------------------------------
            VARIABLES_FOR_WAQ_PRINTOUTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', max='**',
                into = ['','to edit !!!'],
                defaut = '',
                fr = """Noms des variables que l''utilisateur veut ecrire dans
le fichier des resultats QE.
Chaque variable est representee par une lettre.
Ne semble pas utilise pour le moment.""",
                ang = """Names of variables the user wants to write
into the graphic results file.
Does not seem to be used at the moment.""",
            ),
        ),
#       -----------------------------------
        LISTING = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            MASS_BALANCE = SIMP(statut ='f',
#           -----------------------------------
                typ = bool,
                defaut = False,
                fr = """Determine si oui ou non le bilan de masse est realise.
Lu mais pas utilise pour le moment.""",
                ang = """Determines whether a check of the mass-balance over the domain
is made or not.
Read but not used at the moment.""",
            ),
#           -----------------------------------
            WAQ_VARIABLES_TO_BE_PRINTED = SIMP(statut ='f',
#           -----------------------------------
                typ = 'TXM', max='**',
                defaut = '',
                fr = """Nom des variables que l''utilisateur desire ecrire sur
le listing. Memes possibilites que pour les sorties graphiques.
Ne semble pas utilise pour le moment.""",
                ang = """Names of variables the user wants to write on the listing.
Each variable is represented by a letter in the same manner as
it is done in the graphic results file.
Does not seem to be used at the moment.""",
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
        defaut = [0],
        fr = """Pour imprimer la sequence des appels, mettre 1.""",
        ang = """If 1, calls of subroutines will be printed in the listing.""",
    ),
)
# -----------------------------------------------------------------------
HYDRODYNAMICS = PROC(nom= "HYDRODYNAMICS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    PHYSICAL_PARAMETERS = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WATER_DENSITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 999.972,
            fr = """Fixe la valeur de la masse volumique de l''eau.""",
            ang = """Sets the value of water density.""",
        ),
#       -----------------------------------
        KINEMATIC_WATER_VISCOSITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 1.E-6,
            fr = """Definit la viscosite cinematique de l''eau.
En m/s$^2$. Lu mais pas utilise pour le moment.""",
            ang = """Specifies the water kinematic viscosity.
En m/s$^2$. Read but not used at the moment.""",
        ),
    ),
#   -----------------------------------
    SUSPENSION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        DISPERSION_ALONG_THE_FLOW = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1.E-2],
            fr = """Lu mais pas utilise pour le moment.""",
            ang = """Read but not used at the moment.""",
        ),
#       -----------------------------------
        DISPERSION_ACROSS_THE_FLOW = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1.E-2],
            fr = """Lu mais pas utilise pour le moment.""",
            ang = """Read but not used at the moment.""",
        ),
    ),
)
# -----------------------------------------------------------------------
PHYSICAL_PARAMETERS = PROC(nom= "PHYSICAL_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    EVAPORATION_RATE = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = [0.],
        fr = """Taux d evaporation - meme unite que la pluie en m$^3$/s/m$^2$.""",
        ang = """Rate of evaporation - same unit as rainfall in m$^3$/s/m$^2$.""",
    ),
)
# -----------------------------------------------------------------------
WAQ_PARAMETERS = PROC(nom= "WAQ_PARAMETERS",op = None,
# -----------------------------------------------------------------------
#   -----------------------------------
    METHOD_OF_COMPUTATION_OF_RAY_EXTINCTION_COEFFICIENT = SIMP(statut ='o',
#   -----------------------------------
        typ = 'TXM',
        into = ["ATKINS FORMULA","MOSS FORMULA","GIVEN CONSTANT"],
        defaut = "ATKINS FORMULA",
        fr = """Choix de la methode de calcul du coefficient d extinction
du rayonnement solaire dans l eau $k_e$ en m$^{-1}$. Les choix sont:
\begin{itemize}
\item 1 : formule d Atkins (1.7/Secchi) ;
\item 2 : formule de Moss si la profondeur de Secchi est inconnue ;
\item 3 : constante fournie par l utilisateur.
\end{itemize}
Pour EUTRO+BIOMASS et THERMIC en 3D.""",
        ang = """Choice of the method of calculation of the extinction of
sun ray in water $k_e$ in m$^{-1}$. The choices are :
\begin{itemize}
\item 1: Atkins formula (1.7/Secchi),
\item 2: Moss formula if Secchi depth is unknown,
\item 3: constant provided by the user with the keyword
\telkey{LIGHT EXTINCTION COEFFICIENT}.
\end{itemize}
For EUTRO+BIOMASS and THERMIC in 3D.""",
    ),
#   -----------------------------------
    SECCHI_DEPTH = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.9,
        fr = """En m. Utilise pour EUTRO+BIOMAS et THERMIC en 3D.""",
        ang = """In m. Used for EUTRO+BIOMAS, and THERMIC in 3D.""",
    ),
#   -----------------------------------
    LIGHT_EXTINCTION_COEFFICIENT = SIMP(statut ='f',
#   -----------------------------------
        typ = 'R',
        defaut = 0.2,
        fr = """En m$^{-1}$. Utilise pour EUTRO+BIOMAS et THERMIC en 3D.""",
        ang = """In m$^{-1}$. Used for EUTRO+BIOMAS, and THERMIC in 3D.""",
    ),
#   -----------------------------------
    EUTROPHICATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WATER_TEMPERATURE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [7.],
            fr = """En $^{\circ}$C, temperature moyenne de l eau
necessaire pour calculer les valeurs de $C_s$.""",
            ang = """In $^{\circ}$C, mean temperature necessary for
computing different values of $C_s$.""",
        ),
#       -----------------------------------
        O2 = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            VEGETAL_RESPIRATION_R = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.06],
                fr = """En mgO$_2$/j/l. Variable $R$ dans la documentation.
Pour O2 seulement.""",
                ang = """In mgO$_2$/d/l. Variable $R$ in the documentation.
For O2 only.""",
            ),
#           -----------------------------------
            CONSTANT_OF_DEGRADATION_OF_ORGANIC_LOAD_K1 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.25],
                fr = """En j$^{-1}$. Pour O2 seulement.""",
                ang = """In d$^{-1}$. For O2 only.""",
            ),
#           -----------------------------------
            CONSTANT_OF_NITRIFICATION_KINETIC_K4 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.35],
                fr = """En j$^{-1}$. Pour O2 seulement.""",
                ang = """In d$^{-1}$. For O2 only.""",
            ),
        ),
#       -----------------------------------
        EUTRO_AND_O2 = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            BENTHIC_DEMAND = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.1],
                fr = """En gO$_2$/m$^2$/j. Variable $BEN$ dans la documentation.
Pour EUTRO+O2.""",
                ang = """In gO$_2$/m$^2$/d. Variable $BEN$ in the documentation.
For EUTRO+O2.""",
            ),
#           -----------------------------------
            PHOTOSYNTHESIS_P = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [1.],
                fr = """En mgO$_2$/j/l. Variable $P$ dans la documentation.
Son ordre de grandeur est compris entre 0.3 et 9~mgO$_2$/j/l
selon la nature du cours d eau.
Pour EUTRO+O2.""",
                ang = """In mgO$_2$/d/l. Variable $P$ in the documentation.
Between 0.3 and 9~mgO$_2$/d/l depending on the type of river.
For EUTRO+O2.""",
            ),
#           -----------------------------------
            O2_SATURATION_DENSITY_OF_WATER__CS_ = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [11.],
                fr = """En mgO$_2$/l. Variable $C_s$ = 9~mgO$_2$/l at 20$^{\circ}$C.
Pour EUTRO+O2.""",
                ang = """In mgO$_2$/l. Variable $C_s$ = 9~mgO$_2$/l at 20$^{\circ}$C.
For EUTRO+O2.""",
            ),
#           -----------------------------------
            FORMULA_FOR_COMPUTING_K2 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = [1],
                fr = """Donne le choix de calcul de la formule de $k_2$ de la
reaeration naturelle, les options sont les suivantes :
\begin{itemize}
\item 0 : $k_2$ constant, valeur de $k_2$ = 0.9 ;
\item 1 : formule de Tenessee Valley Authority ;
\item 2 : formule de Owens et al. ;
\item 3 : formule de Churchill et al. ;
\item 4 : formule de O Connor \& Dobbins ;
\item 5 : formule combinant les formules 2, 3 et 4.
\end{itemize}
Pour EUTRO+O2.""",
                ang = """Gives how to compute the reaeration coefficient $k_2$
options are:
\begin{itemize}
\item 0: $k_2$ constant, in this case $k_2$ = 0.9,
\item 1: formula of The Tenessee Valley Authority,
\item 2: formula of Owens et al.,
\item 3: formula of Churchill et al.,
\item 4: formula of O Connor \& Dobbins,
\item 5: formula combining the formulae 2, 3 et 4.
\end{itemize}
For EUTRO+O2.""",
            ),
#           -----------------------------------
            K2_REAERATION_COEFFICIENT = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.9],
                fr = """Ou coefficient d echange gazeux eau-atmosphere a 20$^{\circ}$C
En j$^{-1}$. Valeur comprise entre 0.1 et > 1.15~j$^{-1}$.
Variable $k_2$ dans la documentation. Pour EUTRO+O2.""",
                ang = """In d$^{-1}$. Value between 0.1 and > 1.15~d$^{-1}$.
Variable $k_2$ in the documentation. For EUTRO+O2.""",
            ),
#           -----------------------------------
            FORMULA_FOR_COMPUTING_CS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = [0],
                fr = """En j$^{-1}$, les options sont les suivantes :
\begin{itemize}
\item 0 : constante ;
\item 1 : formule de Elmore \& Hayes ;
\item 2 : formule de Montgomery.
\end{itemize}
Pour EUTRO+O2.""",
                ang = """In d$^{-1}$, here are the available options:
\begin{itemize}
\item 0: constant,
\item 1: Elmore \& Hayes formula,
\item 2: Montgomery formula.
\end{itemize}
For EUTRO+O2.""",
            ),
        ),
#       -----------------------------------
        EUTRO_AND_BIOMASS = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            MAXIMUM_ALGAL_GROWTH_RATE_AT_20C = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 2.,
                fr = """Taux de croissance algale maximale a 20$^{\circ}$C.
Pour EUTRO+BIOMASS.""",
                ang = """Maximum algal growth rate at 20$^{\circ}$C.
For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            ALGAL_TOXICITY_COEFFICIENTS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min= 2, max= 2,
                defaut = [1.,0.],
                fr = """$\alpha_1$ et $\alpha_2$ dans la documentation.
$\alpha_i = 1$ signifie absence de toxicite. Pour EUTRO+BIOMASS.""",
                ang = """$\alpha_1$ and $\alpha_2$ in the documentation.
$\alpha_i = 1$ means no toxicity. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            VEGETAL_TURBIDITY_COEFFICIENT_WITHOUT_PHYTO = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """En m$^{-1}$. Pour EUTRO+BIOMASS.""",
                ang = """In m$^{-1}$. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            PARAMETER_OF_CALIBRATION_OF_SMITH_FORMULA = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [120.],
                fr = """En W/m$^2$. Ordre de grandeur = 100~W/m$^2$.
Pour EUTRO+BIOMASS.""",
                ang = """In W/m$^2$. Around 100~W/m$^2$. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            SUNSHINE_FLUX_DENSITY_ON_WATER_SURFACE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = 0.,
                fr = """Densite de flux du rayonnement solaire a la surface de l eau
en W/m$^2$. Pour EUTRO+BIOMASS.""",
                ang = """Density of sunshine flux on the water surface
in W/m$^2$. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            CONSTANT_OF_HALF_SATURATION_WITH_PHOSPHATE = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.005],
                fr = """En mgP/l. Environ 0.005~mgP/l. Pour EUTRO+BIOMASS.""",
                ang = """In mgP/l. Around 0.005~mgP/l. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            CONSTANT_OF_HALF_SATURATION_WITH_NITROGEN = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.03],
                fr = """En mgN/l. Environ 0.03~mgN/l. Pour EUTRO+BIOMASS.""",
                ang = """In mgN/l. Around 0.03~mgN/l. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            RESPIRATION_RATE_OF_ALGAL_BIOMASS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.05],
                fr = """En j$^{-1}$, a 20$^{\circ}$C. Variable $RP$ dans la
documentation. Pour EUTRO+BIOMASS.""",
                ang = """In d$^{-1}$, at 20$^{\circ}$C. Variable $RP$ in
the documentation. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            COEFFICIENTS_OF_ALGAL_MORTALITY_AT_20C = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min= 2, max= 2,
                defaut = [0.1,0.003],
                fr = """Variables $M_1$ et $M_2$ dans la documentation.
Pour EUTRO+BIOMASS.""",
                ang = """Variables $M_1$ and $M_2$ in the documentation.
For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            PROPORTION_OF_PHOSPHORUS_WITHIN_PHYTO_CELLS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.0025],
                fr = """En mgP/$\mu$gChlA. Variable $fp$ dans la documentation.
Pour EUTRO+BIOMASS.""",
                ang = """In mgP/$\mu$gChlA. Variable $fp$ in the documentation.
For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            PERCENTAGE_OF_PHOSPHORUS_ASSIMILABLE_IN_DEAD_PHYTO = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.5],
                fr = """En $\%$. Pour EUTRO+BIOMASS.""",
                ang = """In $\%$. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            RATE_OF_TRANSFORMATION_OF_POR_TO_PO4 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.03],
                fr = """Taux de transformation du phosphore degradable non assimilable
par le phytoplancton (POR) en phosphore mineral dissous assimilable par
le phyto (PO4) par le biais de la mineralisation bacterienne
en j$^{-1}$.
Variables $k_1$ pour BIOMAS et $k_{320}$ pour EUTRO (a 20$^{\circ}$C)
dans la documentation. Pour EUTRO+BIOMASS.""",
                ang = """Rate of transformation of degradable and non assimilable
phosphorus (POR) to mineral dissolved assimilable phosphorus (PO4)
by means of bacteria mineralization in d$^{-1}$.
Variables $k_1$ for BIOMAS and $k_{320}$ for EUTRO (at 20$^{\circ}$C)
in the documentation. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            RATE_OF_TRANSFORMATION_OF_NOR_TO_NO3 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Taux de transformation de l azote degradable
non assimilable par le phytoplancton (NOR) en azote mineral dissous
assimilable par le phytoplancton NO3 par le biais de la
mineralisation bacterienne en j$^{-1}$.
Variables $k_2$ pour BIOMAS et $k_{620}$ pour EUTRO (a 20$^{\circ}$C)
dans la documentation. Pour EUTRO+BIOMASS.""",
                ang = """Rate of transformation of degradable and non assimilable
nitrogen (NOR) to mineral dissolved assimilable nitrogen (NO3)
by bacteria mineralization in d$^{-1}$.
Variables $k_2$ for BIOMAS and $k_{620}$ for EUTRO (at 20$^{\circ}$C)
in the documentation. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            PROPORTION_OF_NITROGEN_WITHIN_PHYTO_CELLS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.0035],
                fr = """En mgN/$\mu$gChlA. Variable $fn$ dans la documentation.
Pour EUTRO+BIOMASS.""",
                ang = """In mgN/$\mu$gChlA. Variable $fn$ in the documentation.
For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            PERCENTAGE_OF_NITROGEN_ASSIMILABLE_IN_DEAD_PHYTO = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.5],
                fr = """En $\%$. Variable $dtn$ dans la documentation.
Pour EUTRO+BIOMASS.""",
                ang = """In $\%$. Variable $dtn$ in the documentation.
For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            SEDIMENTATION_VELOCITY_OF_ORGANIC_PHOSPHORUS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Vitesse de sedimentation du phosphore organique non algal
en m/s. Pour EUTRO+BIOMASS.""",
                ang = """Sedimentation velocity of non algal organic phosphorus
in m/s. For EUTRO+BIOMASS.""",
            ),
#           -----------------------------------
            SEDIMENTATION_VELOCITY_OF_NON_ALGAL_NITROGEN = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """Vitesse de sedimentation de l azote organique non algal
en m/s. Pour EUTRO+BIOMASS.""",
                ang = """Sedimentation velocity of non algal organic nitrogen
in m/s. For EUTRO+BIOMASS.""",
            ),
        ),
#       -----------------------------------
        EUTRO = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            CONSUMED_OXYGEN_BY_NITRIFICATION = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [5.2],
                fr = """En mgO$_2$/mgNH$_4$. Variable $n$ in the documentation.
Pour EUTRO seulement.""",
                ang = """In mgO$_2$/mgNH$_4$. Variable $n$ in the documentation.
For EUTRO only.""",
            ),
#           -----------------------------------
            CONSTANT_FOR_THE_NITRIFICATION_KINETIC_K520 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.35],
                fr = """Constante de la cinetique de nitrification a 20$^{\circ}$C.
En j$^{-1}$. Pour EUTRO seulement.""",
                ang = """Constant for the nitrification kinetic at 20$^{\circ}$C.
In d$^{-1}$. For EUTRO only.""",
            ),
#           -----------------------------------
            CONSTANT_OF_DEGRADATION_OF_ORGANIC_LOAD_K120 = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.35],
                fr = """Constante de cinetique de degradation de la charge organique
a 20$^{\circ}$C. En j$^{-1}$. Pour EUTRO seulement.""",
                ang = """Constant of degradation kinetic of organic load
at 20$^{\circ}$C. In d$^{-1}$. For EUTRO only.""",
            ),
#           -----------------------------------
            SEDIMENTATION_VELOCITY_OF_ORGANIC_LOAD = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.],
                fr = """En m/s. Pour EUTRO seulement.""",
                ang = """In m/s. For EUTRO only.""",
            ),
#           -----------------------------------
            OXYGEN_PRODUCED_BY_PHOTOSYNTHESIS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [0.15],
                fr = """En mgO$_2$/$\mu$gChlA. Variable $f$ dans la documentation.
Pour EUTRO seulement.""",
                ang = """In mgO$_2$/$\mu$gChlA. Variable $f$ in the documentation.
For EUTRO only.""",
            ),
        ),
#       -----------------------------------
        SOURCES = FACT(statut='o',
#       -----------------------------------
#           -----------------------------------
            FORMULA_FOR_COMPUTING_RS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'I',
                defaut = [0],
                fr = """Donne le choix de calcul de la formule de RS de la
reaeration naturelle au niveau des seuils,
les options sont les suivantes :
\begin{itemize}
\item 0 : RS constant, valeur de RS = 1.0 ;
\item 1 : formule de Gameson 1 ;
\item 2 : formule de Gameson 2 ;
\item 3 : formule de WRL1 ;
\item 4 : formule de WRL2.
\end{itemize}
Lu mais pas utilise pour le moment.""",
                ang = """Gives how to cumpute the weir reaeration coefficient RS
options are:
\begin{itemize}
\item 0: RS constant, in this case RS = 1.0,
\item 1: formula of Gameson 1,
\item 2: formula of Gameson 2,
\item 3: formula of WRL1,
\item 4: formula of WRL2.
\end{itemize}
Read but not used at the moment.""",
            ),
#           -----------------------------------
            WEIR_REAERATION_COEFFICIENT_RS = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R',
                defaut = [1.0],
                fr = """En j$^{-1}$. Lu mais pas utilise pour le moment.""",
                ang = """In d$^{-1}$. Read but not used at the moment.""",
            ),
#           -----------------------------------
            COEFFICIENTS_A_AND_B_FOR_RS_FORMULA = SIMP(statut ='f',
#           -----------------------------------
                typ = 'R', min= 2, max= 2,
                defaut = [1.2,0.7],
                fr = """Coefficients intervenant dans le calcul de RS :
$a$ est entre 0.65 (eau tres polluee) et 1.8 (eau tres claire)
et $b$ varie beaucoup (voir tableau dans la documentation).
Lu mais pas utilise pour le moment.""",
                ang = """Coefficients needed for the calculation of RS:
$a$ is between 0.65 (very polluted water and 1.8 (very clear water)
and $b$ varies a lot (see array in the documentation).
Read but not used at the moment.""",
            ),
        ),
    ),
#   -----------------------------------
    MICROPOL = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        SEDIMENTATION_CRITICAL_STRESS = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [5.],
            fr = """Tension de cisaillement critique de sedimentation
en Pa. Pour MICROPOL seulement.""",
            ang = """Sedimentation critical shear stress
in Pa. For MICROPOL only.""",
        ),
#       -----------------------------------
        SEDIMENT_SETTLING_VELOCITY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [6.E-6],
            fr = """Vitesse de sedimentation en m/s. Variable $w$ dans la
documentation. Pour MICROPOL seulement.""",
            ang = """Sediment velocity in m/s. Variable $w$ in the documentation.
For MICROPOL only.""",
        ),
#       -----------------------------------
        CRITICAL_STRESS_OF_RESUSPENSION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1000.],
            fr = """Tension de cisaillement critique de remise en
suspension en Pa. Variable $\tau_r$ dans la documentation.
Pour MICROPOL seulement.""",
            ang = """Sedimentation critical shear stress
in Pa.  Variable $\tau_r$ in the documentation. For MICROPOL only.""",
        ),
#       -----------------------------------
        EROSION_RATE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [0.],
            fr = """Taux d erosion caracteristique des MES deposees
ou encore appele constante de Partheniades.
Variable $e$ dans la documentation.
Pour MICROPOL seulement.""",
            ang = """ Characteristic erosion rate  of deposited SPM
or also called Partheniades s constant.
Variable $e$ in the documentation.
For MICROPOL only.""",
        ),
#       -----------------------------------
        CONSTANT_OF_DESORPTION_KINETIC = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [2.5E-7],
            fr = """En s$^{-1}$. Variable $k_{-1}$ dans la documentation.
Pour MICROPOL seulement.""",
            ang = """In s$^{-1}$. Variable $k_{-1}$ in the documentation.
For MICROPOL only.""",
        ),
#       -----------------------------------
        COEFFICIENT_OF_DISTRIBUTION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1775.],
            fr = """En m$^3$/kg ou l/g. Variable $K_d$ dans la documentation.
Pour MICROPOL seulement.""",
            ang = """In m$^3$/kg or l/g. Variable $K_d$ in the documentation.
For MICROPOL only.""",
        ),
#       -----------------------------------
        EXPONENTIAL_DESINTEGRATION_CONSTANT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1.13E-7],
            fr = """En s$^{-1}$, loi de decroissance exponentielle comme celle de
la radioactivite. Variable $L$ dans la documentation.
Pour MICROPOL seulement.""",
            ang = """In s$^{-1}$, exponential decrease law like the one of
radioactivity. Variable $L$ in the documentation. For MICROPOL only.""",
        ),
    ),
#   -----------------------------------
    THERMIC = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        WATER_SPECIFIC_HEAT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [4180.],
            fr = """En J/kg/$^{\circ}$C. Pour THERMIC seulement.""",
            ang = """In J/kg/$^{\circ}$C. For THERMIC only.""",
        ),
#       -----------------------------------
        COEFFICIENTS_FOR_CALIBRATING_ATMOSPHERIC_RADIATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.97,
            fr = """Variable $e_{\rm{air}}$ dans la documentation.
Pour THERMIC seulement.""",
            ang = """Variable $e_{\rm{air}}$ in the documentation.
For THERMIC only.""",
        ),
#       -----------------------------------
        COEFFICIENT_OF_CLOUDING_RATE = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.17,
            fr = """Coefficient qui depend du type de nuages :
\begin{itemize}
\item Cirrus = 0.04 ;
\item Cirro stratus = 0.08 ;
\item Alto cumulus = 0.17 ;
\item Alto stratus = 0.2 ;
\item Stratus = 0.24
\end{itemize}
Alto Cumulus (valeur moyenne) est utilise habituellement
(T.V.A. 1972), valeur par defaut en 3D.
En 3D, uniquement utilise avec la formule de Swinbank (1963)
pour le calcul du rayonnement atmospherique.
Variable $k$ dans la documentation. Pour THERMIC seulement.
Ancienne valeur par defaut = 0.2 jusqu a la version V8P1.""",
            ang = """Coefficient depending on the type of clouds:
\begin{itemize}
\item Cirrus = 0.04,
\item Cirro stratus = 0.08,
\item Alto cumulus = 0.17,
\item Alto stratus = 0.2,
\item Stratus = 0.24.
\end{itemize}
Alto Cumulus (mean value) is usually used (T.V.A. 1972),
default value in 3D.
In 3D, only used with Swinbank formula (1963)
for the computation of the atmospheric radiation.
Variable $k$ in the documentation.
For THERMIC only.
Old default value = 0.2 until release V8P1.""",
        ),
#       -----------------------------------
        COEFFICIENTS_FOR_CALIBRATING_SURFACE_WATER_RADIATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = 0.97,
            fr = """Il depend du site et des obstacles entourant le plan d eau.
Pour une riviere etroite bordee d arbres, ca serait voisin de 0.97
et pour un terrain largement decouvert, ca serait voisin de 0.92.
Variable $e_{\rm{eau}}$ dans la documentation.
Pour THERMIC seulement.""",
            ang = """It depends on the location and the obstacles around the
water. For a narrow river bordered with trees, it would be around
0.97 and for a widely oopen field, it would be around 0.92.
Variable $e_{\rm{eau}}$ in the documentation.
For THERMIC only.""",
        ),
#       -----------------------------------
        AIR_SPECIFIC_HEAT = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [1005.],
            fr = """En J/kg/$^{\circ}$C. Pour THERMIC seulement.""",
            ang = """In J/kg/$^{\circ}$C. For THERMIC only.""",
        ),
#       -----------------------------------
        COEFFICIENTS_OF_AERATION_FORMULA = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min= 2, max= 2,
            defaut = [0.002,0.0012],
            fr = """Couple de coefficients de calage pour la fonction de vent
dans les modeles d echanges eau-atmosphere.
Leurs valeurs tres proches valent environ 0.0025.""",
            ang = """Couple of calibration coefficients for the wind function
of the atmosphere-water exchange models.
Their close values are around 0.0025.""",
        ),
#       -----------------------------------
        ATMOSPHERE_WATER_EXCHANGE_MODEL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["NO MODEL","LINEARISED FORMULA AT THE FREE SURFACE","MODEL WITH COMPLETE BALANCE"],
            defaut = ["NO MODEL"],
            fr = """Choix du modele d echanges entre l eau et l atmosphere.
\begin{itemize}
\item 0 : pas de mode d echanges eau-atmosphere (defaut) ;
\item formule linearisee a la surface (3D seulement) ;
\item modele a bilan complet.
\end{itemize}
En 2D, si autre processus que THERMIC, laisser obligatoirement
a 0.""",
            ang = """Choice of the atmosphere-water exchange model.
\begin{itemize}
\item 0: no model (default),
\item 1: linearised formula at the free surface,
\item 2: model with complete balance.
\end{itemize}
In 2D, if another processus than THERMIC, mandatory to let to 0.""",
        ),
#       -----------------------------------
        COEFFICIENT_TO_CALIBRATE_THE_ATMOSPHERE_WATER_EXCHANGE_MODEL = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R',
            defaut = [0.0025],
            fr = """Valeur du coefficient de calage pour la fonction de vent
dans les modeles d echanges eau-atmosphere
(formule linearisee a la surface ou bilan complet).
Une valeur comprise entre 0.0017 et 0.0035 est conseillee.
Seulement pour le THERMIC 3D.""",
            ang = """Value of the calibration coefficient for the wind function
of the atmosphere-water exchange models
(linearised formula at the free surface or complete balance).
A value between 0.0017 and 0.0035 is advised.
Only for THERMIC 3D.""",
        ),
#       -----------------------------------
        FORMULA_OF_ATMOSPHERIC_RADIATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["IDSO AND JACKSON (1969)","SWINBANK (1963)","BRUTSAERT (1975)","YAJIMA TONO DAM (2014)"],
            defaut = "SWINBANK (1963)",
            fr = """Formule au choix pour le calcul du rayonnement atmospherique.
Voir GLM.
\begin{itemize}
\item 1 : Idso et Jackson (1969) ;
\item 2 : Swibank (1963), defaut ;
\item 3 : Brutsaert (1975) ;
\item 4 : Yajima Tono Dam (2014).
\end{itemize}
Seulement pour le THERMIC 3D.""",
            ang = """Formula to be chosen to compute the atmospheric radiation.
See GLM.
\begin{itemize}
\item 1: Idso and Jackson (1969),
\item 2: Swibank (1963), default,
\item 3: Brutsaert (1975),
\item 4: Yajima Tono Dam (2014).
\end{itemize}
Only for THERMIC 3D.""",
        ),
#       -----------------------------------
        LIGHTNESS_OF_THE_SKY = SIMP(statut ='f',
#       -----------------------------------
            typ = 'TXM',
            into = ["VERY BRIGHT, PURE SKY","MODERATELY BRIGHT SKY","FOGGY LIKE THE SKY OF INDUSTRIAL AREA "],
            defaut = "MODERATELY BRIGHT SKY",
            fr = """Degre de clarte (purete) du ciel.
Les choix possibles sont :
\begin{itemize}
\item ciel tres clair, tres pur ;
\item ciel moyennement clair ;
\item ciel d une zone industrielle, opaque.
\end{itemize}
Pour THERMIC en 3D seulement.""",
            ang = """How the sky is bright (pure).
Possible choices are:
\begin{itemize}
\item very bright, pure sky,
\item moderately bright sky,
\item foggy like the sky of industrial area.
\end{itemize}
For THERMIC in 3D only.""",
        ),
#       -----------------------------------
        SOLAR_RADIATION_READ_IN_METEO_FILE = SIMP(statut ='f',
#       -----------------------------------
            typ = bool,
            defaut = False,
            fr = """Si la donnee meteo de rayonnement solaire est disponible,
elle peut etre lue dans le
\telkey{FICHIER ASCII DE DONNEES ATMOSPHERIQUES} de \telemac{2d} ou
\telemac{3d} au lieu d etre calculee par le module d echanges
eau-atmosphere en activant ce mot cle a OUI.
Pour THERMIC en 3D seulement. This is mandatory in 2D.""",
            ang = """If solar radiation data is available, it can be read in the
\telkey{ASCII ATMOSPHERIC DATA FILE} of \telemac{2d} or \telemac{3d}
instead of been computed by the heat exchange with atmosphere module
by activating this keyword to YES.
For THERMIC in 3D only. This is mandatory in 2D.""",
        ),
    ),
#   -----------------------------------
    DEGRADATION = FACT(statut='o',
#   -----------------------------------
#       -----------------------------------
        LAW_OF_TRACERS_DEGRADATION = SIMP(statut ='o',
#       -----------------------------------
            typ = 'TXM', min=0, max='**',
            into = ["NO DEGRADATION","F(T90) LAW","FIRST ORDER KINETIC LAW, CONSTANT OF DEGRADATION IN HOURS","FIRST ORDER KINETIC LAW, CONSTANT OF DEGRADATION IN DAYS","LAW IMPLEMENTED BY USER"],
            defaut = ["NO DEGRADATION","NO DEGRADATION"],
            fr = """Prise en compte d''une loi de decroissance des traceurs. Valeur 0 si
pas de prise en compte, 1 si loi 1... Pour chaque valeur entree, un
nom de traceur correspondant doit etre present dans
\telkey{NOMS DES TRACEURS}, auquel sera appliquee la loi.
Les choix possibles sont :
\begin{itemize}
\item 0 : pas de degradation ;
\item 1 : loi de degradation bacteriologique avec coefficient $T_{90}$ ;
\item 2 : loi de degradation d ordre 1, constante de degradation
en h$^{-1}$ ;
\item 3 : loi de degradation d ordre 1, constante de degradation
en jours$^{-1}$ ;
\item 4 : loi programmee par l utilisateur.
\end{itemize}""",
            ang = """Take in account a law for tracers decrease. Value 0 if not taken into
account, 1 if law 1... For each value entered, a corresponding name
should be present in the keyword \telkey{NAMES OF TRACERS}, so that the
decrease law is applied to the correct tracer(s).
Possible choices are:
\begin{itemize}
\item 0: no degradation,
\item 1: law for bacterial degradation with $T_{90}$ coefficient,
\item 2: degradation law of first order, constant of tracer kinetic
degradation in h$^{-1}$,
\item 3: degradation law of first order, constant of tracer kinetic
degradation in d$^{-1}$,
\item 4: law implemented by user.
\end{itemize}""",
        ),
#       -----------------------------------
        COEFFICIENT_1_FOR_LAW_OF_TRACERS_DEGRADATION = SIMP(statut ='f',
#       -----------------------------------
            typ = 'R', min=0, max='**',
            fr = """Coefficient 1 de la loi de decroissance des traceurs. Voir aussi la
correspondance avec \telkey{NOMS DES TRACEURS} et l''aide du mot cle
\telkey{LOI DE DEGRADATION DES TRACEURS}.""",
            ang = """Coefficient 1 of law for tracers decrease. Check also the relation
between the keywords \telkey{NAMES OF TRACERS} and
\telkey{LAW OF TRACERS DEGRADATION}.""",
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
        typ = ('Fichier','All Files (*)'),
        defaut = '',
        fr = """Nom du fichier contenant les parametres du calcul
QE a realiser.""",
        ang = """Name of the file containing parameters of the WAQ
computation Written by the user.""",
    ),
#   -----------------------------------
    DICTIONARY = SIMP(statut ='f',
#   -----------------------------------
        typ = ('Fichier','All Files (*)'),
        defaut = 'waqtel.dico',
        fr = """Dictionnaire des mots cles.""",
        ang = """Key word dictionary.""",
    ),
)
TEXTE_NEW_JDC = "\
COMPUTATION_ENVIRONMENT();\
GENERAL_PARAMETERS();\
HYDRODYNAMICS();\
PHYSICAL_PARAMETERS();\
WAQ_PARAMETERS();\
"
Ordre_Des_Commandes = (
'COMPUTATION_ENVIRONMENT',
'GENERAL_PARAMETERS',
'HYDRODYNAMICS',
'PHYSICAL_PARAMETERS',
'WAQ_PARAMETERS',
'INTERNAL')
try:
    import TelApy
    source = "eficas"
except Exception as excpt:
    source = "Telemac"
enum = source+'.waqtel_enum_auto'
dicoCasEn = source+'.waqtel_dicoCasEnToCata'
dicoCasFr = source+'.waqtel_dicoCasFrToCata'
