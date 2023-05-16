#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author TELEMAC-MASCARET Consortium
@brief Convert a listing courlis binary file into ASCII Opthyca file (with
sediment budget file)
"""

import argparse
from data_manip.formats.mascaret_file import MascaretFile

def convert_listing_courlis(name, output_name, write_opt, budget,
                            txt_format, csv_format, xlsx_format,
                            time_check):
    """
    Converting binary listing Courlis file into ASCII Opthyca file

    @param name (str) Name of the listing
    @param output_name (str) Name of the output Opthyca file
    @param write_opt (bool) If true writting opt file
    @param budget (bool) If true writting budget file
    @param txt_format (bool) If true writting sediment budget in text format
    @param csv_format (bool) If true writting sediment budget in text format
    @param xlsx_format (bool) If true writting sediment budget in Excel format
    @param time_check (bool) If true  display time when reading/writting
    """
    short_name = ".".join(name.split(".")[:-1])

    if output_name.split(".")[-1:] == ["opt"]:
        output_name = short_name + output_name

    listing_courlis = MascaretFile(name, fformat='listingcourlis')

    if write_opt:
        print('\nConversion from '+name+' to '+output_name+'\n')
        listing_courlis.write_optfile(output_name, timecheck=time_check)

    if budget:
        if txt_format:
            print('\nGenerating ' + short_name + '_sediment_budget.txt\n')
            listing_courlis.export_sediment_budget_to_txt(\
                    short_name + '_sediment_budget.txt')

        if csv_format:
            print('\nGenerating ' + short_name + '_sediment_budget.csv\n')
            listing_courlis.export_sediment_budget_to_csv(\
                    short_name + '_sediment_budget.csv')

        if xlsx_format:
            print('\nGenerating ' + short_name + '_sediment_budget.xslx\n')

            try:
                import pandas as pd
            except ImportError as xcpt:
                print("Pandas is mandatory to genreate an Excel file")
                raise xcpt
            try:
                import openpyxl
            except ImportError as xcpt:
                print("openpyxl is mandatory to genreate an Excel file")
                raise xcpt
            mud, sand, total = listing_courlis.sediment_budget(pandas=True)
            writer = pd.ExcelWriter(short_name+'_sediment_budget.xlsx')

            mud.to_excel(writer, sheet_name='mud')
            sand.to_excel(writer, sheet_name='sand')
            total.to_excel(writer, sheet_name='total')

            writer.save()

def convert_courlis_parser(subparser):
    """
    Building parser

    @param subparser (ArgumentParser) The parser

    @returns (ArgumentaParser) The udpated parser
    """
    parser = subparser.add_parser('convert_listing_courlis',
                                  help=\
        'convert binary listing '
        'Courlis file into ASCII Opthyca file with sediment budget file')
    parser.add_argument("args", metavar='listing Courlis file')
    parser.add_argument("-o", "--output",
                        dest="outputFileName", default='_listingCourlis.opt',
                        help="Option to give a name to the output Opthyca file")
    parser.add_argument("-no-opt", "--no-opthyca",
                        dest="write_opt", default=True, action="store_false",
                        help="Option to disable writing of Opthyca file")
    parser.add_argument("-no-b", "--no-budget",
                        dest="budget", default=True, action="store_false",
                        help="Option to disable writing of budget file")
    parser.add_argument("-no-csv", "--no-budget_file-csv",
                        dest="csv_format", default=True, action="store_false",
                        help="Option to export sediment budget as csv file "
                             "(sediment_budget.csv)")
    parser.add_argument("-txt", "--budget_file-txt",
                        dest="txt_format", default=False, action="store_true",
                        help="Option to export sediment budget as txt file "
                             "(sediment_budget.txt)")
    parser.add_argument("-xlsx", "--budget_file-xlsx",
                        dest="xlsx_format", default=False, action="store_true",
                        help="Option to export sediment budget as xlsx Exclel "
                             "file (sediment_budget.xlsx)")
    parser.add_argument("-time", "--time-check",
                        dest="time_check", default=False, action="store_true",
                        help="Option to display reading and writing time step")
    return subparser

def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        description=
        'convert_listing_courlis is a Python script converting binary listing '
        'Courlis file into ASCII Opthyca file with sediment budget file')

    parser.add_argument("args", metavar='listing Courlis file')
    parser.add_argument("-o", "--output",
                        dest="outputFileName", default='_listingCourlis.opt',
                        help="Option to give a name to the output Opthyca file")
    parser.add_argument("-no-opt", "--no-opthyca",
                        dest="write_opt", default=True, action="store_false",
                        help="Option to disable writing of Opthyca file")
    parser.add_argument("-no-b", "--no-budget",
                        dest="budget", default=True, action="store_false",
                        help="Option to disable writing of budget file")
    parser.add_argument("-no-csv", "--no-budget_file-csv",
                        dest="csv_format", default=True, action="store_false",
                        help="Option to export sediment budget as csv file "
                             "(sediment_budget.csv)")
    parser.add_argument("-txt", "--budget_file-txt",
                        dest="txt_format", default=False, action="store_true",
                        help="Option to export sediment budget as txt file "
                             "(sediment_budget.txt)")
    parser.add_argument("-xlsx", "--budget_file-xlsx",
                        dest="xlsx_format", default=False, action="store_true",
                        help="Option to export sediment budget as xlsx Exclel "
                             "file (sediment_budget.xlsx)")
    parser.add_argument("-time", "--time-check",
                        dest="time_check", default=False, action="store_true",
                        help="Option to display reading and writing time step")

    options = parser.parse_args()

    convert_listing_courlis(options.args, options.outputFileName,
                            options.write_opt, options.budget,
                            options.txt_format, options.csv_format,
                            options.xlsx_format, options.time_check)

if __name__ == "__main__":
    main()
