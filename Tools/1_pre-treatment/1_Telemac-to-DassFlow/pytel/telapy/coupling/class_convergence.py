# coding: utf-8
"""
ClassConvergence class
================
"""


class ClassConvergence:
    """
    Check coupling convergence
    """

    def __init__(self, nb_model1d, do_stop, output=False):
        """
        Constructor
        @param nb_model1d (int) :  number of 1D model
        @param do_stop (int) : Max iteration of schwarz loop
        @param output (bool) : boolean indicating the creation
         of convergence file
        """
        self.nb_model1d = nb_model1d
        self.output = output
        self.do_stop = do_stop

        self.crit_h = None
        self.crit_v = None

        self.h_1d = list()
        self.h_2d = list()
        self.v_1d = list()
        self.v_2d = list()

        # variables initialisation
        self.id_cv = False
        val = [9999] * self.nb_model1d
        self.list_var = ['v', 'h', 'q', 'i', 'j', 'sm']
        self.list_out = ['v', 'sm', 'q', 'h', 'i', 'j']
        self.list_wdt = [7, 10, 9, 7, 7, 7]
        for var in self.list_var:
            for idx in ['_1d', '_2d']:
                setattr(self, var + idx, val)
            setattr(self, 'crit_' + var, val)
        if self.output:
            self.file = open('Convergence_criteria.out', 'w', buffering=1)
            self.filex = open('Convergence_criteria.csv', 'w', buffering=1)
            hdr = 'STEP, ITER, MODEL, '
            for ivar in range(len(self.list_out)):
                var = self.list_out[ivar]
                if ivar == len(self.list_out) - 1:
                    hdr += var.upper() + '\n'
                else:
                    hdr += var.upper() + ', '
            self.filex.write(hdr)
        else:
            self.file = None
            self.filex = None

    def set_criteria(self, height, velocity):
        """
        function allows to change h and v variables
        @param height (list) : convergence criterion on the height
        @param velocity (list): convergence criterion on the velocity
        """
        self.crit_h = height
        self.crit_v = velocity

    def get_value(self, var_1d, var_2d):
        """
        Get computed convergence variables
        @param var_1d (np.array): 1D model variables
        @param var_2d (np.array) :  2D model variables
        """

        for i, var in enumerate(self.list_var):
            setattr(self, var + '_1d', var_1d[:, i])
            setattr(self, var + '_2d', var_2d[:, i])

    def check_conv(self, ib_do):
        """
        Check convergence function
        @param ib_do (int): index of Schwarz iterration
        """

        self.id_cv = True
        for i in range(self.nb_model1d):
            self.id_cv = self.id_cv and \
                         abs(self.h_1d[i] - self.h_2d[i]) <= self.crit_h[i] \
                         and abs(self.v_1d[i] - self.v_2d[i]) <= self.crit_v[i]
        self.id_cv = self.id_cv or ib_do >= self.do_stop - 1

    def write_outfile(self, ib_tps, ib_do):
        """
        Write output file
        param ib_tps (int)  time step index
        param ib_do (int) Schwarz step index

        """
        if self.id_cv and self.file is not None:
            for i in range(self.nb_model1d):
                msg = 'id_modele : {},  ib_tps : {}, ib_do : {}\n' \
                    .format(i + 1, ib_tps + 1, ib_do + 1)
                self.file.write(msg)
                val_1d = ''
                val_2d = ''
                for ivar in range(len(self.list_out)):
                    var = self.list_out[ivar]
                    wdt = self.list_wdt[ivar]
                    val_1d += '{}_1D : {:{wdt}.3f}' \
                              .format(var.upper(),
                                      eval('self.{}_1d[{}]'.format(var, i)),
                                      wdt=wdt)
                    val_2d += '{}_2D : {:{wdt}.3f}' \
                              .format(var.upper(),
                                      eval('self.{}_2d[{}]'.format(var, i)),
                                      wdt=wdt)
                    if ivar == len(self.list_out) - 1:
                        val_1d += '\n'
                        val_2d += '\n'
                    else:
                        val_1d += ', '
                        val_2d += ', '

                self.file.write(val_1d + val_2d + '\n')
                line1d = '{}, {}, 1D.{}, '.format(ib_tps + 1, ib_do + 1, i + 1)
                line2d = '{}, {}, 2D.{}, '.format(ib_tps + 1, ib_do + 1, i + 1)
                for ivar in range(len(self.list_out)):
                    var = self.list_out[ivar]
                    wdt = self.list_wdt[ivar]
                    line1d += '{:{wdt}.3f}' \
                              .format(eval('self.{}_1d[{}]'.format(var, i)),
                                      wdt=wdt)
                    line2d += '{:{wdt}.3f}' \
                              .format(eval('self.{}_2d[{}]'.format(var, i)),
                                      wdt=wdt)
                    if ivar == len(self.list_out) - 1:
                        line1d += '\n'
                        line2d += '\n'
                    else:
                        line1d += ', '
                        line2d += ', '
                self.filex.write(line1d + line2d)

    def finalize(self):
        """
        Close file function
        """
        if self.output:
            self.file.close()
            self.filex.close()

    def main(self, ib_do, var_1d, var_2d, ib_tps=None):
        """
        Main function to check convergence:
            - Get values
            - Check convergence
            - Write or not output file
            - Return if convergence or not
        @param ib_do (int) : Schwarz step index
        @param var_1d (np.array): computed convergence variables for 1D model
        @param var_2d (np.array) :  computed convergence variables for 2D model
        @param ib_tps (int) : Time step index
        @return (bool) : if convergence or not
        """
        self.get_value(var_1d, var_2d)
        self.check_conv(ib_do)
        if self.output:
            self.write_outfile(ib_tps, ib_do)
        return self.id_cv
