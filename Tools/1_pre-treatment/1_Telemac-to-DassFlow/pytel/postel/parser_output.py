r"""@author TELEMAC-MASCARET Consortium

    @brief
            Functions to read Telemac result file
"""
import re
from os import path, walk
from fnmatch import fnmatch
from utils.files import get_file_content
from utils.exceptions import TelemacException
import numpy as np


def get_latest_output_files(input_file):
    """
    Inspired from matchSafe in utils.py
    Follows first the simple template casFile+'_*??h??min??s.sortie'
    and then look for the parallel equivalent
    """
    # ~~> list all entries
    dir_path, _, filenames = next(walk(path.dirname(input_file)))
    # ~~> match expression
    exnames = []  # [ path.basename(input_file) ]
    for fout in filenames:
        if fnmatch(fout, path.basename(input_file)+'_*??h??min??s.sortie'):
            exnames.append(fout)
    if exnames == []:
        return []
    casbase = sorted(exnames).pop()
    exnames = [path.join(dir_path, casbase)]
    casbase = path.splitext(casbase)[0]
    for fout in filenames:
        if fnmatch(fout, casbase+'_*.sortie'):
            exnames.append(path.join(dir_path, fout))
    return exnames


class OutputFileData():
    """
    @brief : Read and store data result file
    :return:
    """
    def __init__(self, file_name=''):
        self.output = []
        if file_name != '':
            self._get_file_content(file_name)

    def _get_file_content(self, file_name):
        """
        @brief Get the content of the file
                 if file does not exist then system exit
        """
        if not path.exists(file_name):
            raise TelemacException(
                '... could not find your .sortie file: '
                '{}'.format(file_name))
        self.output = get_file_content(file_name)

    def get_time_profile(self):
        """
        @brief : Returns the time profile in iteration and in seconds,
                         read from the TELEMAC-* output file
                         Also sets the xLabel to either 'Time (s)'
                         or 'Iteration #'
        @returns (str, list), (str, list): Iteration name, list and same of the
        times
        """
        form = re.compile(r'\s*ITERATION\s+(?P<iteration>\d+)'
                          + r'\s+(TIME)[\s:]*(?P<others>.*?)'
                          + r'(?P<number>\b((?:(\d+)\b)|(?:(\d+(|\.)'
                          + r'\d*[dDeE](\+|\-)?\d+'
                          + r'|\d+\.\d+))))\s+S\s*(|\))\s*\Z', re.I)
        itr = []
        time = []
        for line in self.output:
            proc = re.match(form, line)
            if proc:
                itr.append(int(proc.group('iteration')))
                time.append(float(proc.group('number')))
        return ('Iteration #', np.array(itr)), ('Time (s)', np.array(time, dtype=np.float))

    def get_exec_time(self):
        """
        Returns the execution time (in seconds)

        @returns (float) The time
        """
        form = re.compile(r'\A\s*('\
                          r'((?P<days>[0-9]+)\s*DAYS)'\
                          r'|((?P<hours>[0-9]+)\s*HOURS)'\
                          r'|((?P<minutes>[0-9]+)\s*MINUTES)'\
                          r'|((?P<seconds>[0-9]+)\s*SECONDS)'\
                          r')\s*\Z'\
                          )

        exec_time = 0.0
        for line in self.output:
            proc = re.match(form, line)
            if proc:
                if proc.group('days') is not None:
                    exec_time += 24*60*60*int(proc.group('days'))
                if proc.group('hours') is not None:
                    exec_time += 60*60*int(proc.group('hours'))
                if proc.group('minutes') is not None:
                    exec_time += 60*int(proc.group('minutes'))
                if proc.group('seconds') is not None:
                    exec_time += int(proc.group('seconds'))

        return exec_time

    # ~~ Name of the study ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Returns the name of the study, read from the TELEMAC output file
    #   +> If no name, returns NO NAME
    #   +> If not found, returns NAME NO FOUND
    def get_name_of_study(self):
        """
        the name of the study, read from the TELEMAC output file

        @returns (string) returns NO NAME OF STUDY if no name and NAME OF STUDY
        NO FOUND if not found
        """
        form = re.compile(r'\s*(?P<before>.*?)(NAME OF THE STUDY'
                          + r')[\s:]*(?P<after>.'
                          + r'*?)\s*\Z', re.I)
        for line in range(len(self.output)):
            proc = re.match(form, self.output[line])
            if proc:
                if self.output[line+1].strip() == '':
                    return 'NO NAME OF STUDY'
                return self.output[line+1].strip()
        return 'NAME OF STUDY NOT FOUND'

    def get_volume_profile(self):
        """
        Returns the time series of Values, read from
        the TELEMAC-2D output file
        volumes, fluxes, errors, etc.
        Assumptions:
        - if VOLUME ... is not found, will not try
        to read FLUX and ERROR
        - for every VOLUME instance, it will advance
         to find FLUX and ERROR
        Also sets the yLabel to either 'Volume (m3/s)'
        or 'Fluxes (-)' or 'Error(-)'
        """
        form_liqnumbers = re.compile(
            r'\s*(THERE IS)\s+(?P<number>\d+)'
            + r'\s+(LIQUID BOUNDARIES:)\s*\Z', re.I)
        form_liqnumberp = re.compile(
            r'\s*(NUMBER OF LIQUID BOUNDARIES)'
            + r'\s+(?P<number>\d+)\s*\Z', re.I)
        form_volinitial = re.compile(
            r'\s*(INITIAL VOLUME )'
            + r'[\s:]*\s+(?P<value>\b([-+]|)((?:(\d+)\b)|'
            + r'(?:(\d+(|\.)'
            + r'\d*[dDeE](\+|\-)?\d+|\d+\.\d+))))\s+'
            + r'(?P<after>.*?)\s*\Z', re.I)
        form_voltotal = re.compile(
            r'\s*(VOLUME IN THE DOMAIN)[\s:]*'
            + r'\s+(?P<value>\b([-+]|)((?:(\d+)\b)|'
            + r'(?:(\d+(|\.)'
            + r'\d*[dDeE](\+|\-)?\d+|\d+\.\d+))))\s+'
            + r'(?P<after>.*?)\s*\Z', re.I)
        form_volfluxes = re.compile(
            r'\s*(FLUX BOUNDARY)\s+'
            + r'(?P<number>\d+)\s*:\s*'
            + r'(?P<value>[+-]*\b((?:(\d+)\b)|(?:(\d+(|\.)'
            + r'\d*[dDeE](\+|\-)?\d+|\d+\.\d+))))(.\s|\s)+'
            + r'(?P<after>.*?)\s*\Z', re.I)
        form_volerror = re.compile(
            r'\s*(RELATIVE ERROR IN VOLUME AT T '
            + r'=)\s+'
            + r'(?P<at>[+-]*\b((?:(\d+)\b)|(?:(\d+(|\.)'
            + r'\d*[dDeE](\+|\-)?\d+|\d+\.\d+))))\s+S :\s+'
            + r'(?P<value>[+-]*\b((?:(\d+)\b)|(?:(\d+(|\.)'
            + r'\d*[dDeE](\+|\-)?\d+|\d+\.\d+))))'
            + r'\s*\Z', re.I)
        iline = 0
        # ~~ Searches for number of liquid boundaries ~~~~~~~~~~~~~~~~~~~
        fluxes_prof = []
        bound_names = []
        liqnumber = 0
        while iline < len(self.output):
            proc = re.match(form_liqnumbers, self.output[iline])
            if not proc:
                proc = re.match(form_liqnumberp, self.output[iline])
            if proc:
                liqnumber = int(proc.group('number'))
                for i in range(liqnumber):
                    fluxes_prof.append([])
                    bound_names.append('Boundary ' + str(i+1))
                break
            iline = iline + 1
        # ~~ Initiates profiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        volumes_prof = []
        volumes_name = 'Volumes (m3/s)'
        errors_prof = []
        errors_name = 'Error (-)'
        fluxes_name = 'Fluxes (-)'
        # ~~ Reads the rest of time profiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        while iline < len(self.output):
            if re.match(form_volinitial, self.output[iline]):
                break
            proc = re.match(form_voltotal, self.output[iline])
            if proc:
                volumes_prof.append(float(proc.group('value')))
                for i in range(liqnumber):
                    iline = iline + 1
                    proc = re.match(form_volfluxes, self.output[iline])
                    while not proc:
                        iline = iline + 1
                        if iline >= len(self.output):
                            raise TelemacException(
                                '... Could not parse FLUXES'
                                'FOR BOUNDARY {}'.format(str(i+1)))
                        proc = re.match(form_volfluxes, self.output[iline])
                    fluxes_prof[i].append(float(proc.group('value')))
                iline = iline + 1
                proc = re.match(form_volerror, self.output[iline])
                while not proc:
                    iline = iline + 1
                    if iline >= len(self.output):
                        raise TelemacException(
                            '... Could not parse RELATIVE ERROR IN VOLUME ')
                errors_prof.append(float(proc.group('value')))

            iline = iline + 1

        # ~~ Adds initial volume ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        while iline < len(self.output):
            proc = re.match(form_volinitial, self.output[iline])
            if proc:
                volumes_prof.insert(0, float(proc.group('value')))
                break
            iline = iline + 1
        # ~~ Adds initial error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        errors_prof.insert(0, 0.0)  # assumed
        # ~~ Adds initial fluxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for i in range(liqnumber):      # 0.0 may not be the correct value
            fluxes_prof[i].insert(0, 0.0)

        return (volumes_name, volumes_prof), \
               (fluxes_name, bound_names, fluxes_prof),\
               (errors_name, errors_prof)

    def get_value_history_output(self, vrs):
        """
        @brief Read values from the TELEMAC-2D output file
        @param vrs (string) voltotal: extract total volume
                            volfluxes: extract boundary fluxes
                            volerror: extract error in volume
        @return (list, list) x,y arrays for plotting
        """
        # ~~ Extract data
        _, time = self.get_time_profile()
        voltot, volflu, volerr = self.get_volume_profile()
        volinfo = []
        for varname in vrs.split(';'):
            # ~~ y-axis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if varname == "voltotal":
                volinfo.append(voltot)
            elif varname == "volfluxes":
                volinfo.append(volflu)
            elif varname == "volerror":
                volinfo.append(volerr)
            else:
                raise TelemacException(
                    '... do not know how to extract: '
                    '{} '.format(varname))
        return time, volinfo

    def get_sediment_mass_profile(self):
        """
        Return the total and lost mass of sediment for each class for each
        listing iteration

        @returns (dict) Keys are the id of the sediment and each contains a
        dictionary with eky total and lost
        """
        # Regular expression for sediment number
        sed_class = re.compile(r'\s*(SEDIMENT CLASS NUMBER)\s+(=)\s*(?P<sed_class>\d+)', re.I)
        all_class = re.compile(r'\s*( GAIA MASS-BALANCE OF SEDIMENTS OVER ALL CLASSES:)\s*', re.I)
        # Regular expression for total mass
        total_mass = re.compile(r'\s*(TOTAL MASS)\s+(=)\s*(?P<mass>[-0-9de.+]+)', re.I)
        # Regular expression for lost mass
        lost_mass = re.compile(r'\s*(LOST MASS)\s+(=)\s*(?P<mass>[-0-9de.+]+)', re.I)
        # Regular expression for lost mass
        cumu_lost_mass = re.compile(r'\s*(CUMULATED LOST MASS)\s+(=)\s*(?P<mass>[-0-9de.+]+)', re.I)
        # Regular expression for relative error to initial act layer
        re_ial_mass = re.compile(r'\s*(RELATIVE ERROR TO INITIAL ACT LAYER MASS)\s+(=)\s*(?P<mass>[-0-9de.+]+)', re.I)
        # extration data

        mass = {}
        iclass = -1
        for line in self.output:
            # Identify for which class is the information
            proc = re.match(sed_class, line)
            if proc:
                iclass = int(proc.group('sed_class'))
                if iclass not in mass:
                    mass[iclass] = {'total':[], 'lost':[], 'cumu_lost':[], 're_ial':[]}
            # When more than one class mass over all classes is given
            proc = re.match(all_class, line)
            if proc:
                iclass = 0
                if iclass not in mass:
                    mass[iclass] = {'total':[], 'lost':[], 'cumu_lost':[], 're_ial':[]}

            # Getting the total mass
            proc = re.match(total_mass, line)
            if proc:
                mass[iclass]['total'].append(float(proc.group('mass')))

            # Getting the lost_mass
            proc = re.match(lost_mass, line)
            if proc:
                mass[iclass]['lost'].append(float(proc.group('mass')))

            # Getting the cumu_lost_mass
            proc = re.match(cumu_lost_mass, line)
            if proc:
                mass[iclass]['cumu_lost'].append(float(proc.group('mass')))

            # Getting the re_ial_mass
            proc = re.match(re_ial_mass, line)
            if proc:
                mass[iclass]['re_ial'].append(float(proc.group('mass')))

        for iclass in mass:
            for ttype in ['total', 'lost', 'cumu_lost', 're_ial']:
                if ttype in ['total', 're_ial']:
                    # Removing last value as it is the one from the final bilan
                    array = np.array(mass[iclass][ttype][:-1], dtype=np.float)
                else:
                    array = np.array(mass[iclass][ttype], dtype=np.float)
                mass[iclass][ttype] = array

        return mass

    def get_user_defined_output(self, user_form, line_num=False):
        """
        Read user defined values from the TELEMAC-* output file

        @param user_form (string) user gives a regular expression
                                  to find in the file.
        @param line_num (boolean) return or not the line number

        @return (list, list) or (list) file line (and line number if
        line_num==True)
        """
        iline = 0
        user_value = []
        num_line = []
        while iline < len(self.output):
            proc = re.match(user_form, self.output[iline])
            if proc:
                user_value.append(self.output[iline].strip())
                num_line.append(iline)
            iline += 1
        if line_num:
            return user_value, num_line
        else:
            return user_value

    def get_tracer_mass_profile(self):
        """
        Return the total and lost mass of suspended tracer for each class for
        each listing iteration

        @returns (dict) Keys are the id of the tracer and each contains a
        dictionary with eky total and lost
        """
        # Regular expression for water mass header
        water_mass_header = re.compile(r'\s*(WATER)', re.I)
        # Regular expression for tracer number mass header
        tracer_class = re.compile(r'\s*(TRACER)\s*(?P<tracer_class>\d+):', re.I)
        # Regular expression for total mass
        total_mass = re.compile(r'\s*(MASS AT THE PRESENT TIME STEP)\s+(:)\s*(?P<mass>[-0-9de.+]+)', re.I)
        # Regular expression for lost mass
        lost_mass = re.compile(r'\s*(ERROR ON THE MASS DURING THIS TIME STEP)\s+(:)\s*(?P<mass>[-0-9de.+]+)', re.I)

        # extract data
        mass = {}
        iclass = -1
        tracer_listing = True
        for line in self.output:

            # Check for the water mass section (so we ignore it)
            proc = re.match(water_mass_header, line)
            if proc:
                tracer_listing = False

            # Identify for which class is the information
            proc = re.match(tracer_class, line)
            if proc:
                tracer_listing = True
                iclass = int(proc.group('tracer_class'))
                if iclass not in mass:
                    mass[iclass] = {'total':[], 'lost':[]}

            # Get the total mass (if not in water mass section)
            if tracer_listing:
                proc = re.match(total_mass, line)
                if proc:
                    mass[iclass]['total'].append(float(proc.group('mass')))

            # Getting the lost_mass
                proc = re.match(lost_mass, line)
                if proc:
                    mass[iclass]['lost'].append(float(proc.group('mass')))


        for iclass in mass:
            for ttype in ['total', 'lost']:
                array = np.array(mass[iclass][ttype], dtype=np.float)
                mass[iclass][ttype] = array

        return mass
