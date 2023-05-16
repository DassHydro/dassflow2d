""" script to compare steering files """

from collections import OrderedDict
from utils.exceptions import TelemacException
from execution.telemac_cas import TelemacCas, get_dico

def diff_cas(cas_names, csv_name, dico_name, delimiter=',', allkey=False):
    """ Write a .csv file of keyword Damocles differences between .cas files

    @param cas_names (list/tuple) list of .cas Telemac steering files
    @param csv_name (string) name of the .csv to write to show the differences
    @param module (str) name of the common Telemac module for all the .cas
    @param all (bool) If True all keyword added in csv otherwise only the ones
    where there is a difference
    @param delimiter (str) separator for the csv file
    """
    import pandas as pd
    # Test the presence of .cas files
    if (not isinstance(cas_names, list)) & (not isinstance(cas_names, tuple)):
        raise TelemacException("The parameter 'cas_names' is not a list/tuple")

    if len(cas_names) < 2:
        raise TelemacException('At least two .cas files are expected')

    for i, cas_name in enumerate(cas_names):
        if not isinstance(cas_name, str):
            raise TelemacException('cas_names[%d] is not string' % i)

    # Test csv_name
    if not isinstance(csv_name, str):
        raise TelemacException('csv_name is not string')

    # Save info on the current directory
    # curdir = getcwd()

    # Read all the .cas files
    cases = list()
    for file_name in cas_names:
        cases.append(TelemacCas(file_name, dico_name, check_files=False).values)

    # Merge all the steering files keywords
    merged_keyword = OrderedDict()
    listall = set(cases[0])
    for k in cases[1:]:
        listall = listall | set(k)
    listall = sorted(listall)
    for k in listall:
        merged_keyword[k] = list()
    for cas in cases:
        for k in listall:
            merged_keyword[k].append(cas.get(k, None))

    # Identify the rows of equal values
    merged_list = list(merged_keyword.items())
    line = list()

    # Identify lines where values are all the same
    if not allkey:
        for k, valist in enumerate(merged_list):
            if all(x == valist[1][0] for x in valist[1]):
                line.append(k)

    # Write the .csv file using a DataFrame
    diffarray = pd.DataFrame.from_dict(merged_keyword, orient='index',
                                       columns=cas_names)

    if not allkey:
        # Remove rows with no difference
        diffarray = diffarray.drop(diffarray.index[line])

    # diffarray.drop_duplicates(inplace=True, keep=False)
    diffarray.to_csv(csv_name)
