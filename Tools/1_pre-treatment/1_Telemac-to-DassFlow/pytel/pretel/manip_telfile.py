"""@author TELEMAC-MASCARET Consortium

   Manipulation of Telemac files (mesh, results files)
"""
from os import path, remove
from data_manip.extraction.telemac_file import TelemacFile
from utils.exceptions import TelemacException


def scan(tel_file, bnd_file, data):
    """
    Ascii dump of content of TelemacFile

    @param tel_file (str) Name of the telemac file
    @param bnd_file (str) Name of the boundary file
    @param data (boolean) If True display data information as well
    """

    res = TelemacFile(tel_file, bnd_file=bnd_file)

    res.print_info(full=data)

def diff(tel_file1, tel_file2, diff_file):
    """
    Creates a file that is the difference between two files (file1 - file2) of
    the same shape (same mesh same number of variable and same number of
    records)

    @param tel_file1 (str) Name of the first file
    @param tel_file1 (str) Name of the second file
    @param diff_file (str) Name of the file containing the difference
    """
    res1 = TelemacFile(tel_file1)
    res2 = TelemacFile(tel_file2)
    if path.exists(diff_file):
        print(" ~> Deleting eixisting diff_file: {}".format(diff_file))
        remove(diff_file)
    dif = TelemacFile(diff_file, access='w')

    res1.read()
    res2.read()

    if res1._values.shape != res2._values.shape:
        raise TelemacException(
            "Different shape between files (ntimestep, nvar, npoin):\n"
            "{}: {}\n{}: {}"
            .format(tel_file1, res1._values.shape,
                    tel_file2, res2._values.shape))

    dif.read(res1)

    for time in range(res1.ntimestep):
        for var in range(res1.nvar):
            dif._values[time, var, :] -= res2._values[time, var, :]

    dif.write()
