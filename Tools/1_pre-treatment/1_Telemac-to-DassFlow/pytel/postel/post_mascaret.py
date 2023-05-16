r"""@author TELEMAC-MASCARET Consortium

   @brief
         Functions to read rubens and opthyca mascaret result file
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from struct import unpack
# ~~> dependencies towards other pytel/modules
from utils.exceptions import TelemacException
import numpy as np


def lit_header_opt(resultfile):
    """
    @brief : read the header of an opthyca result file
    :param resultfile: opthyca result file
    :return: list of variables available (name, short name, unit, indicator)
    """
    var_liste = []

    try:
        line = resultfile.readline()
    except OSError:
        print("\n\nlit_header_opt: Wrong file format\n\n")
        return

    if '[variables]' not in line:
        print("\n\nlit_header_opt: Wrong first line\n\n")
        return

    line = resultfile.readline()
    while '[resultats]' not in line:
        var_liste.append(line.split(';', -1))
        line = resultfile.readline()

    return var_liste


def lit_res_opt(fle):
    """
    @brief : read the opthyca result file
    :param fle: opthyca result file
    :return: list of time, reach, section and PK available and values
    """
    time_records = []
    reach_idx = []
    section_idx = []
    section_pk = []

    res_file = open(fle, 'r')
    _ = lit_header_opt(res_file)
    lines = res_file.readlines()
    for line in lines:
        time_record, reach, sec1, sec2, _ = line.split(';', 4)
        if time_record not in time_records:
            time_records.append(time_record)
        if reach not in reach_idx:
            reach_idx.append(reach)
        if sec1 not in section_idx:
            section_idx.append(sec1)
        if sec2 not in section_pk:
            section_pk.append(sec2)
    res_file.close()
    return time_records, reach_idx, section_idx, section_pk


def get_endian_from_char(input_file, nchar):
    """
    @brief : get the endianess of the file
    :param input_file: rubens result file
    :param nchar:
    :return: endianess
    """
    pointer = input_file.tell()
    endian = ">"  # "<" means little-endian, ">" means big-endian
    ll, _, chk = unpack(endian+'i'+str(nchar)+'si', input_file.read(4+nchar+4))
    if chk != nchar:
        endian = "<"
        input_file.seek(pointer)
        ll, _, chk = \
            unpack(endian+'i'+str(nchar)+'si', input_file.read(4+nchar+4))
    if ll != chk:
        raise TelemacException(
           '... Cannot read characters from your binary file'
           '     +> Maybe it is the wrong file format ? {}'.format(str(nchar)))
    input_file.seek(pointer)
    return endian


def lit_res_rub(fle):
    """
    @brief : read the rubens result file
    :param fle: rubens result file
    :return: list of time, reach, section and PK available and values
    """

    res_file = open(fle, 'rb')
    ftype = get_endian_from_char(res_file, 72)

    _, _, _ = unpack(ftype + 'i72si', res_file.read(4 + 72 + 4))
    _, _, _ = unpack(ftype + 'i72si', res_file.read(4 + 72 + 4))
    _, _, _ = unpack(ftype + 'i72si', res_file.read(4 + 72 + 4))
    _, _, _ = unpack(ftype + 'i4si', res_file.read(4 + 4 + 4))
    _, _, _ = unpack(ftype + 'i4si', res_file.read(4 + 4 + 4))
    _, cfin, _ = unpack(ftype + 'i4si', res_file.read(4 + 4 + 4))
    _, nbief1, _, _ = unpack(ftype + 'iiii', res_file.read(4 + 4 + 4 + 4))

    orig_bief = []
    for _ in range(nbief1):
        _, ilu, _ = unpack(ftype + 'iii', res_file.read(4 + 4 + 4))
        orig_bief.append(ilu)
    fin_bief = []
    for _ in range(nbief1):
        _, ilu, _ = unpack(ftype + 'iii', res_file.read(4 + 4 + 4))
        fin_bief.append(ilu)

    nomvar_indep = []
    _, clu, _ = unpack(ftype + 'i4si', res_file.read(4 + 4 + 4))
    while clu != cfin:
        nomvar_indep.append(clu)
        _, clu, _ = unpack(ftype + 'i4si', res_file.read(4 + 4 + 4))
    _, nsto1, _, _ = unpack(ftype + 'iiii', res_file.read(4 + 4 + 4 + 4))

    valvar_indep = []
    for _ in enumerate(nomvar_indep):
        res_file.seek(4, 1)
        varlu = np.array(unpack(ftype+str(nsto1)+'f',
                         res_file.read(nsto1 * 4)))
        res_file.seek(4, 1)
        valvar_indep.append(varlu)

    nomvar_dep = []
    _, clu, _ = unpack(ftype + 'i4si', res_file.read(4 + 4 + 4))
    while clu != cfin:
        nomvar_dep.append(clu)
        _, clu, _ = unpack(ftype + 'i4si', res_file.read(4 + 4 + 4))

    itemps = []
    temps = []
    valvar_dep = []
    while True:
        try:
            _, ilu, ilu, _ = \
                unpack(ftype + 'iiii', res_file.read(4 + 4 + 4 + 4))
            itemps.append(ilu)
            res_file.seek(4, 1)
            flu, flu = unpack(ftype+'ff', res_file.read(4 + 4))
            temps.append(flu)
            res_file.seek(4, 1)
            _, nsto1, _, _ = unpack(ftype + 'iiii',
                                    res_file.read(4 + 4 + 4 + 4))
            for _ in enumerate(nomvar_dep):
                res_file.seek(4, 1)
                varlu = np.array(unpack(ftype+str(nsto1)+'f',
                                        res_file.read(nsto1 * 4)))
                res_file.seek(4, 1)
                valvar_dep.append(varlu)
        except OSError:
            break

    res_file.close()

    return


def main():
    """
    @brief : Main function of post_mascaret
    """
#    time_res, reach_idx, section_idx, section_pk =\
#    lit_res_opt('mascaret_imp_ecr.opt')
    lit_res_rub('mascaret.rub')
#    print(time_res)
#    print(reach_idx)
#    print(section_idx)
#    print(section_pk)
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')


if __name__ == "__main__":
    main()
