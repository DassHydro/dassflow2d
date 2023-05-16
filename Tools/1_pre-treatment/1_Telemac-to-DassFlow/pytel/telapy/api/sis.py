# -*- coding: utf-8 -*-
"""
    Python wrapper to the Fortran APIs of Sisyphe

    Author(s): Fabrice Zaoui, Yoann Audouin, Cedric Goeury, Renaud Barate

    Copyright EDF 2016
"""

from __future__ import print_function
import os
from telapy.api.api_module import ApiModule
from utils.exceptions import TelemacException


class Sisyphe(ApiModule):
    """The Sisyphe Python class for APIs"""
    _instanciated = False

    def __new__(cls, *args, **kwargs):
        if cls._instanciated:
            raise TelemacException("a Sisyphe instance already exists")
        instance = ApiModule.__new__(cls)
        cls._instanciated = True
        return instance

    def __init__(self, casfile,
                 user_fortran=None,
                 dicofile=None,
                 lang=2, stdout=6,
                 comm=None,
                 log_lvl='INFO',
                 recompile=True):
        """
        Constructor for Sisyphe
        @param casfile (string) Name of the steering file
        @param user_fortran (string) Name of the user Fortran
        @param dicofile (string) Path to the dictionary
        @param lang (int) Language for ouput (1: French, 2:English)
        @param stdout (int) Where to put the listing
        @param comm (MPI.Comm) MPI communicator
        @param recompile (boolean) If true recompiling the API
        @param log_lvl (string) Logger level
        """
        if dicofile is None:
            hometel = os.getenv("HOMETEL")
            if hometel is not None:
                default_dicofile = os.path.join(os.getenv("HOMETEL"),
                                                "sources",
                                                "sisyphe",
                                                "sisyphe.dico")
            else:
                default_dicofile = 'sisyphe.dico'
            dicofile = default_dicofile
        super(Sisyphe, self).__init__("sis", casfile, user_fortran,
                                      dicofile, lang, stdout, comm,
                                      recompile, log_lvl=log_lvl)

    def __del__(self):
        """
        Destructor
        """
        Sisyphe._instanciated = False
