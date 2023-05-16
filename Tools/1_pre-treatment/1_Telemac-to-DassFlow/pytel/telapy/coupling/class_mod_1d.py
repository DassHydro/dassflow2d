# coding: utf-8
"""
ClassMod1D class
=================
"""
import json
from telapy.tools.study_masc_driven import MascaretStudy

from telapy.coupling.class_comm_1d import ClassComm1D
from telapy.coupling.class_cpl_1d import ClassCpl1D


class ClassMod1D:
    """
    Relevant variables and operations to define a coupled 1D model
    """

    def __init__(self):
        """
        Constructor

        1. Initialize the coupling parallel context under
        """

        # Initialize fields (as a mnemonic content list)
        self.freq_res = 0  # Output frequence (in seconds)
        self.maxiter = 0  # Max. nr of iterates for cpl convergence
        self.cplsteps = 0  # Nr of coupling time steps
        self.cplstp = 0  # Coupling interval  counter
        self.iter = 0  # Coupling iteration counter
        self.d_t = 0.0
        self.t_i = 0.0
        self.t_f = 0.0
        self.initst = 0.0
        self.endtst = 0.0
        self.ti_ini = 0.0
        self.inirun = False
        self.nit_2d = 0
        self.nb_1d_models = 0
        self.cplsect = 0
        self.cploutnam = ''
        self.coupling_law = 0
        self.dt_2d = 0.0
        self.tstincpl = 0
        self.tbcincpl = 0
        self.cplfreq = 0.0
        self.config = ''
        self.cplmethod = 'multiplicativeschwarz'

        self.converged = False
        self.cploutbc = list()
        self.times_bc = None
        self.times2d = None
        self.nbbc = 0
        self.tbcinstp = 0
        self.tstinstp = 0
        self.instance = 0
        self.current_bc = None
        self.bc2 = None
        self.tstart = 0
        self.state_id = 0

        # Initialize the coupling parallel context
        self.comm = ClassComm1D(self)
        self.runcfg = {}
        try:
            with open('CurrentRunDef.json') as cfgfile:
                self.runcfg = json.load(cfgfile)
        except (ValueError, IOError):
            raise ValueError('Current Run Def isn\'t a correct json file.\n'
                             'Validate it with https://jsonlint.com\nSTOP')

        self.get_cpl_config()
        self.jsonbc = 'bc1D_restart_' + self.config + '_' + \
            str(int(self.ti_ini)) + '.json'

        # Initialize the model driver as an instance of MascaretStudy
        # from a JSON setting file

        dico_mas = self.runcfg['config_{}'.format(self.config)]

        self.study = MascaretStudy(dico_mas,
                                   lig=self.restartfile,
                                   log_lvl='CRITICAL',
                                   iprint=0,
                                   working_directory='study_par_'
                                   + self.config)
        # Complete the coupling definition once the model is loaded
        self.complete_cpl_config()
        # Initialize the coupling service utilities
        self.cpl = ClassCpl1D(self)

        # Initialize the communications
        self.comm.init_comms()
        self.cpl.init_conv()

        # Initialize the result output file (fixed variable list in a .opt)
        self.opt_init()

    def get_cpl_config(self):
        """
        Read in the coupling definition from json file.
        """

        # Number of criteria
        self.nbcriteria = 6

        self.jsoncfg = self.runcfg['coupling_def']
        self.cplsteps = self.runcfg["Run"]["CplTSteps"]
        self.ti_ini = self.runcfg["Run"]["StartTime"]
        self.inirun = self.runcfg["Run"]["InitialRun"]

        # Load instance independent parameters
        cpldef = self.jsoncfg["Coupling"]
        mods1d = self.jsoncfg["1D"]
        mod2d = self.jsoncfg["2D"][list(self.jsoncfg["2D"].keys())[0]]
        itface = self.jsoncfg["Interfaces"]
        self.cplfreq = cpldef["TimeStep"]
        self.maxiter = cpldef["MaxIter"]
        self.dt_2d = mod2d["TimeStep"]
        self.nit_2d = int(self.cplfreq / self.dt_2d) + 1
        self.nb_1d_models = len(mods1d)
        self.crit_arret_h = [i["ConvCriteria"]["Height"] for i in itface]
        self.crit_arret_v = [i["ConvCriteria"]["Velocity"] for i in itface]
        if "CplStepRestart1D" in cpldef:
            self.cpl1dstart = cpldef["CplStepRestart1D"].lower()
        else:
            self.cpl1dstart = 'persist2d'

        # Load current instance
        ordered_cfg = list(mods1d.keys())
        ordered_cfg.sort()
        self.config = ordered_cfg[self.instance]
        mod1d = mods1d[self.config]
        self.d_t = mod1d["TimeStep"]
        self.freq_res = max(int(mod1d["OutputFreq"] / self.cplfreq), 1)
        self.restartfile = self.runcfg["1D"][self.config]["WaterLineFile"]
        self.cplmethod = cpldef["Method"].lower()

    def complete_cpl_config(self):
        """
        Complete the coupling definition with the model dependent data.
        NOTICE THAT IT SHOULD BE INTERFACE RELATED
         NOT ANYMORE INDEXED ON MODELS
        """
        masc = self.study.masc
        itface = self.jsoncfg["Interfaces"]
        myitfc = [i for i in itface if i["Id1D"] == self.config]
        # From now on, we rely on the 1 interface per model hypothesis
        nbextr, _, _ = masc.get_var_size('Model.Boundary.Name')
        for extr in range(nbextr):
            if myitfc[0]["IdExtr1D"].lower() == \
                    masc.get('Model.Boundary.Name', extr).lower():
                self.cplsect = \
                    masc.get('Model.Connect.NodeNumFreeOutflow', extr) - 1
                self.coupling_law = masc.get('Model.Boundary.GraphNum', extr)

    def opt_init(self):
        """
        Initialize the .opt output file

        1. Create the header with the fixed choice of fields
        2. Initializes the list of the sections and the section
            - bief association
        """

        self.opt = open(self.study.settings['files']['res'], 'w', buffering=1)

        self.opt.write('[variables]\n')
        self.opt.write('"Cote du fond";"ZREF";"m";4\n')
        self.opt.write('"Cote de l eau";"Z";"m";3\n')
        self.opt.write('"Debit mineur";"QMIN";"m3/s";3\n')
        self.opt.write('"Debit majeur";"QMAJ";"m3/s";3\n')
        self.opt.write(
            '"Coefficient de frottement mineur";"KMIN";"m1/3/s";0\n')
        self.opt.write(
            '"Coefficient de frottement majeur";"KMAJ";"m1/3/s";0\n')
        self.opt.write('"Nombre de Froude";"FR";"";5\n')
        self.opt.write('"Debit total";"Q";"m3/s";3\n')
        self.opt.write('[resultats]\n')

        opt_stock = self.study.masc.get('Model.RecOption')
        if opt_stock == 2:
            nb_res, _, _ = self.study.masc.get_var_size('Model.NodeRes')
            if nb_res > 0:
                self.noderes = \
                    [self.study.masc.get('Model.NodeRes', node)
                     for node in range(nb_res)]
            else:
                self.noderes = []
        elif opt_stock == 1:
            nb_res, _, _ = self.study.masc.get_var_size('Model.X')
            self.noderes = [node + 1 for node in range(nb_res)]
        else:
            self.noderes = []

        nbbf = self.study.masc.get_var_size('Model.Connect.FirstNdNum')[0]
        oribf = [self.study.masc.get('Model.Connect.FirstNdNum', i) - 1
                 for i in range(nbbf)]
        endbf = [self.study.masc.get('Model.Connect.LastNdNum', i)
                 for i in range(nbbf)]
        self.num_bief_s = [ib + 1 for ib in range(nbbf)
                           for i in range(oribf[ib], endbf[ib])]
        del oribf
        del endbf

    def opt_final(self):
        """
        Close .opt output file
        """
        self.opt.close()

    def opt_store(self, time):
        """
        Add an entry to the .opt output file

        1. output of T,X,Zbot,Z,Q1,Q2,CF1,CF2,Froude,Q

        @param time (float) : time for output
        """
        for node in self.noderes:
            self.opt.write(
                '%12.1f;"%2i";"%5i"' % (time, self.num_bief_s[node - 1], node))
            self.opt.write(
                ';%11.2f' % self.study.masc.get('Model.X', node - 1))
            self.opt.write(
                ';%13.4f' % self.study.masc.get('Model.Zbot', node - 1))
            self.opt.write(
                ';%12.3f' % self.study.masc.get('State.Z', node - 1))
            self.opt.write(
                ';%12.3f' % self.study.masc.get('State.Q1', node - 1))
            self.opt.write(
                ';%12.3f' % self.study.masc.get('State.Q2', node - 1))
            self.opt.write(
                ';%9.0f' % self.study.masc.get(
                    'Model.FricCoefMainCh', node - 1))
            self.opt.write(
                ';%9.0f' % self.study.masc.get('Model.FricCoefFP', node - 1))
            self.opt.write(
                ';%14.5f' % self.study.masc.get('State.Froude', node - 1))
            self.opt.write(
                ';%12.3f' % self.study.masc.get('State.Q', node - 1))
            self.opt.write('\n')

    def handle_restart(self):
        """
        Manage the restart.

        1. At the first iterate store the initial model state
            and clean old ones
        2. At successive iterates, restore it
        """

        if self.iter == 0:
            self.study.masc.free_all_saved_states()
            self.state_id = self.study.masc.save_state()
        else:
            self.study.masc.set_state(self.state_id, 0)

    def __call__(self):
        """
        Iterator.

        1. Store or reload the restart
        2. Store the initial output boundary conditions
        3. Integrate the model on the coupling step window
            and store the output BC's
        4. Compute the convergence criteria at the coupling sections
        """

        # Store or reload the restart
        self.handle_restart()

        # Store initial 1d coupling bc's
        self.cploutbc[0] = self.study.masc.get(self.cploutnam, self.cplsect)

        for tst in range(self.tstinstp):
            self.tstart = self.t_i + tst * self.d_t
            self.study.masc.compute_bc(self.tstart, self.tstart + self.d_t,
                                       self.d_t,
                                       self.times_bc, self.tbcinstp, self.nbbc,
                                       self.current_bc, self.bc2)

            # Store current 1d coupling bc's
            self.cploutbc[tst + 1] = \
                self.study.masc.get(self.cploutnam, self.cplsect)
        # Compute convergence criteria
        self.cpl.compute_criteria()

    def selective_output(self):
        """
        Selective output at prescribed output frequency.

        1. Add a line to the opthyca file
        2. Store a restart waterline
        """

        if self.endtst % self.freq_res == 0:
            # Add a record to the .opt
            self.opt_store(self.t_f)

    def finalize(self):
        """
        Finalize a model integration

        1. Send last bc to 2D if the method is additive
        2. Store the restart information for the 1D bc
        3. Close the output file
        """

        if self.cplmethod == 'additiveschwarz':
            self.cpl.transmit_bc()
        # Save the restart
        filelig = \
            'WaterLine_' + self.config + '_' + str(int(self.t_f)) + '.lig'
        self.study.save_lig_restart(filelig)
        self.cpl.store_bc()
        self.opt_final()
