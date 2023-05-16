"""
Verifications on steering file content
"""
from datetime import datetime, timedelta
import re
import numpy as np
from data_manip.extraction.telemac_file import TelemacFile
from utils.exceptions import TelemacException

def check_cas(module, cas):
    """
    Does some check on the steering file

    @param cas (TelemacCas) steering file
    """
    if cas.lang == 'fr':
        default = 'DEFAUT'
    else:
        default = 'DEFAUT1'

    print("\n  ~>  Checking default values:\n")
    for key in cas.values:
        data = cas.dico.data[key]
        if default in data:
            if cas.values[key] == data[default]:
                print('{} set to default value {}'.format(key, data[default]))

    check_time(module, cas)

def get_simulation_date(module, cas):
    """
    Get the beginning and ending time of the simulation

    @param cas (TelemacCas) Steerinf file
    """

    start_date = datetime(1900, 1, 1)
    end_date = datetime(1900, 1, 1)

    if module in ['telemac2d', 'telemac3d']:

        ntimestep = cas.get("NUMBER OF TIME STEPS")
        time_step = cas.get("TIME STEP")
        duration = cas.get("DURATION", 0.0)

        ntimestep = max(ntimestep, int(duration/time_step + 0.5))

        previous_comp = cas.get('COMPUTATION CONTINUED')
        tel_date = get_cas_date(module, cas)

        if previous_comp:
            # Get time from previous computation file
            previous_file = cas.get("PREVIOUS COMPUTATION FILE")
            record_prev = cas.get('RECORD NUMBER FOR RESTART')
            # Record are starting from 0 in python 1 n the steering file
            # 0 in steering file means last time step
            record_prev -= 1
            date, _, time = get_file_date(previous_file, last=record_prev)
            start_date = date + timedelta(seconds=time)
        else:
            start_date = tel_date

        end_date = start_date + timedelta(seconds=ntimestep*time_step)

    return start_date, end_date


def get_file_date(file_name, last=-1):
    """
    Get date, first time last time from file

    @param file_name (str) Name of the file

    @returns (datetime, float, float)
    """
    res = TelemacFile(file_name)
    date = res.datetime
    # Default date
    if all(date[0:3] == 0):
        date[0:3] = [1900, 1, 1]
    date = datetime(date[0], date[1], date[2],
                    date[3], date[4], date[5])
    begin = res.times[0]
    end = res.times[last]
    res.close()

    return date, begin, end

def get_cas_date(module, cas):
    """
    Get date from steering file

    @param module (str) Name of the telemac module
    @param cas (TelemacCas) Steerfin file

    @returns (datetime)
    """
    tel_date = None
    if module in ['telemac2d', 'telemac3d']:
        # If no date using 1900-01-01
        tel_date = cas.get('ORIGINAL DATE OF TIME')
        if tel_date == [0, 0, 0]:
            tel_date = [1900, 1, 1]
        tel_time = cas.get('ORIGINAL HOUR OF TIME')
        tel_date = datetime(tel_date[0], tel_date[1], tel_date[2],
                            tel_time[0], tel_time[1], tel_time[2])

    return tel_date

def get_lid_date(file_name):
    """
    Get start time, end time and date from a liq boundary file like file

    @param file_name (str) name of the file

    @returns
    """
    data = np.loadtxt(file_name, dtype=str)

    # Identify column containing the time (not always the first one)
    data_col = -1
    for row in data:
        for i, val in enumerate(row):
            if val.upper() == 'T':
                data_col = i
                break
        if data_col != -1:
            break

    if data_col == -1:
        data_col = 0

    is_number = re.compile(r"-?\d+(.\d+)?")
    data_row = 0
    for i, val in enumerate(data[:, data_col]):
        if is_number.match(val):
            data_row = i
            break

    time_start = float(data[data_row, data_col])
    time_end = float(data[-1, data_col])

    # check in the head if there is a reference date
    date = datetime(1900, 1, 1)
    with open(file_name, 'r') as f:
        for line in f.readlines():
            if line[0] != '#':
                break
            if line.startswith('#REFDATE'):
                date = datetime.strptime(line[9:].strip(),
                                         "%Y-%m-%d %H:%M:%S")

    return time_start, time_end, date

def check_previous_comp(module, cas, comp_cont="COMPUTATION CONTINUED",
                        prev_file="PREVIOUS COMPUTATION FILE"):
    """
    Check for compuration continued
    """
    # Checking that same date in steering file and previous computation file
    previous_comp = cas.get(comp_cont)


    if previous_comp:
        print("  ~> Checking {} coherence".format(comp_cont.lower()))

        tel_date = get_cas_date(module, cas)
        reset_time = cas.get('INITIAL TIME SET TO ZERO')
        record_prev = cas.get('RECORD NUMBER FOR RESTART')
        # Record are starting from 0 in python 1 n the steering file
        # 0 in steering file means last time step
        record_prev -= 1

        previous_file = cas.get(prev_file)
        prev_date, _, end = get_file_date(previous_file, last=record_prev)

        prev_date = prev_date + timedelta(seconds=end)


        if not reset_time:
            tel_date = tel_date + timedelta(seconds=end)

        if tel_date != prev_date:
            raise TelemacException(
                "Warning:\nMissmatch between previous computation date ({}) "
                "and steering file date ({})\n"
                "The good values should be:\n"
                "if reset time to zero:\n"
                "steering date = previous_file date + last_timestep"
                "else\n"
                "steering date = previous_file date"
                .format(prev_date, tel_date))



def check_time(module, cas):
    """
    Checking time coherence between input files and steering file info

    @param cas (TelemacCas) Steering file
    """

    passed = True

    if module in ['telemac2d', 'telemac3d']:
        check_previous_comp(module, cas)

    if module == 'telemac3d':
        # Checkgin for 2d continuation
        check_previous_comp(module, cas, comp_cont="2D CONTINUATION",
                            prev_file="FILE FOR 2D CONTINUATION")

    print("  ~> Displaying simulation date")
    start_date, end_date = get_simulation_date(module, cas)
    print("Starting date: {}".format(start_date))
    print("Ending date:   {}".format(end_date))

    atmospheric_file = cas.get('BINARY ATMOSPHERIC DATA FILE')
    if atmospheric_file != '' and module in ['telemac2d', 'telemac3d']:
        print('  ~> Checking atmo binary file time coherence')
        meteo_date, time_start, time_end = get_file_date(atmospheric_file)
        meteo_start = meteo_date + timedelta(seconds=time_start)
        meteo_end = meteo_date + timedelta(seconds=time_end)
        print("Time range of the atmospheric data  {} to {}"
              .format(meteo_start, meteo_end))

        if start_date < meteo_start:
            print("The simulation start at {} wheras the atmo binary "
                  "file start at {}"
                  .format(start_date, meteo_start))
            passed = False
        if end_date > meteo_end:
            print("The simulation ends at {} wheras the atmo binary "
                  "file ends at {}"
                  .format(end_date, meteo_end))
            passed = False

        if passed:
            print("  ~> OK")
        else:
            print("  ~> Failed")

    ascii_meteo_file = cas.get('ASCII ATMOSPHERIC DATA FILE')
    if ascii_meteo_file != '' and module in ['telemac2d', 'telemac3d']:
        print('  ~> Checking atmo ascii file time coherence')
        time_start, time_end, date = get_lid_date(ascii_meteo_file)

        meteo_start = date + timedelta(seconds=time_start)
        meteo_end = date + timedelta(seconds=time_end)
        print("Time range of the atmospheric data  {} to {}"
              .format(meteo_start, meteo_end))

        if start_date < meteo_start:
            print("The simulation start at {} wheras the atmo ascii "
                  "file start at {}"
                  .format(start_date, meteo_start))
            passed = False
        if end_date > meteo_end:
            print("The simulation ends at {} wheras the atmo ascii "
                  "file ends at {}"
                  .format(end_date, meteo_end))
            passed = False

        if passed:
            print("  ~> OK")
        else:
            print("  ~> Failed")

    liq_bnd_file = cas.get('LIQUID BOUNDARIES FILE')
    if liq_bnd_file != '' and module in ['telemac2d', 'telemac3d']:
        print('  ~> Checking liquid boundaries file time coherence')
        time_start, time_end, date = get_lid_date(liq_bnd_file)

        liq_start = date + timedelta(seconds=time_start)
        liq_end = date + timedelta(seconds=time_end)
        print("Time range of the liquid boundary data  {} to {}"
              .format(liq_start, liq_end))

        if start_date < liq_start:
            print("The simulation start at {} wheras the liquid boundaries "
                  "file start at {}"
                  .format(start_date, liq_start))
            passed = False
        if end_date > liq_end:
            print("The simulation ends at {} wheras the liquid boundaries "
                  "file ends at {}"
                  .format(end_date, liq_end))
            passed = False

        if passed:
            print("  ~> OK")
        else:
            print("  ~> Failed")


    src_file = cas.get('SOURCES FILE')
    if src_file != '' and module in ['telemac2d', 'telemac3d']:
        print('  ~> Checking liquid boundaries file time coherence')
        time_start, time_end, date = get_lid_date(src_file)

        src_start = date + timedelta(seconds=time_start)
        src_end = date + timedelta(seconds=time_end)
        print("Time range of the source data  {} to {}"
              .format(src_start, src_end))

        if start_date < src_start:
            print("The simulation start at {} wheras the sources "
                  "file start at {}"
                  .format(start_date, src_start))
            passed = False
        if end_date > src_end:
            print("The simulation ends at {} wheras the sources "
                  "file ends at {}"
                  .format(end_date, src_end))
            passed = False

        if passed:
            print("  ~> OK")
        else:
            print("  ~> Failed")

    bnd_bin_file = cas.get('BINARY BOUNDARY DATA FILE', '')
    if bnd_bin_file != '' and module in ['telemac3d']:
        print('  ~> Checking atmo binary file time coherence')
        bnd_bin_date, time_start, time_end = get_file_date(bnd_bin_file)
        bnd_bin_start = bnd_bin_date + timedelta(seconds=time_start)
        bnd_bin_end = bnd_bin_date + timedelta(seconds=time_end)
        print("Time range of the atmospheric data  {} to {}"
              .format(bnd_bin_start, bnd_bin_end))

        if start_date < bnd_bin_start:
            print("The simulation start at {} wheras the binary boundary"
                  "file start at {}"
                  .format(start_date, bnd_bin_start))
            passed = False
        if end_date > bnd_bin_end:
            print("The simulation ends at {} wheras the binary boundary "
                  "file ends at {}"
                  .format(end_date, bnd_bin_end))
            passed = False

        if passed:
            print("  ~> OK")
        else:
            print("  ~> Failed")
