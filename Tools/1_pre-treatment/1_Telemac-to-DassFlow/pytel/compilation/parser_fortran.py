r"""
@author Sebastien E. Bourban and Noemie Durand

@history 28/04/2011 -- Sebastien E. Bourban
    Now supports SYSTELCFG as a directory (old Perl version, to which
    systel.cfg is added) or as a file.

@history 30/04/2011 -- Sebastien E. Bourban
    Upgrade made to config parsing to include the option to reset the
    version and the root from the command line option:
    -v <version>, reset the version read in the config file with this
    -r <root>, reset the root path read in the config file with this

@history 04/12/2012 -- Juliette Parisi and Sebastien E. Bourban
Simplifying call to parseConfigFile, which now takes two arguments
    options.configFile, and options.configName and return one or more
    valid configurations in an array. Testing for validity is now done
    within config.py

@history 10/01/2013 -- Yoann Audouin
ScanSources goes through subdirectories as well now ignore
hidden directories
Adding scan of .F and .F90 files as well

@history 06/02/2013 -- Sebastien E. Bourban
Adding the functionality of displaying changes (html/diff) made
    to a PRINCI file by comparing individual subroutines to their
    original version.
Further capability to compare changes made between 2 PRINCI files.

@history 01/07/2013 -- Sebastien E. Bourban and Yoann Audoin
Upgrade to the new structure

@history 13/07/2013 -- Sebastien E. Bourban
Now deals with DECLARATIONS first before identifying unkonwn externals

@history 23/09/2014 -- Sebastien E. Bourban and Yoann Audoin
The content of the log files from GRETEL and PARTEL are now reported
in the error report.

@history 11/06/2018 -- Sebastien Bourban and Jean-Romain Delisle
Move of the compilation to new python3 structure

@brief
"""
from __future__ import print_function

# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import re
import sys
from optparse import Values
from copy import deepcopy
from os import path, sep

# ~~> dependencies towards other pytel/modules
from utils.files import is_newer, add_to_list, \
    get_file_content, put_file_content, diff_text_files
from utils.exceptions import TelemacException
from utils.progressbar import ProgressBar
#from config import add_config_argument, update_config, CFGS

DEBUG = False

# _____               ______________________________________________
# ____/ Instructions /_____________________________________________/
#

INSTRUCTIONS = ['ALLOCATE', 'ASSIGN',
                'BACKSPACE', 'BLOCK DATA',
                'CALL', 'CASE', 'CLOSE', 'COMMON', 'CYCLE', 'CONTINUE',
                'DATA', 'DEALLOCATE', 'DEFAULT', 'DO',
                'ELSE', 'ELSEIF', 'ENDIF', 'ENDDO', 'ENDSELECT', 'END',
                'ENDFILE',
                'EQUIVALENCE', 'EXIT',
                'FORMAT',
                'GO', 'TO', 'GOTO',
                'IF', 'IMPLICIT NONE', 'INCLUDE', 'INQUIRE', 'INTERFACE',
                'MULTIPLE',
                'NAMELIST', 'NULLIFY',
                'OPEN',
                'PRINT',
                'READ', 'REWIND', 'RETURN',
                'SELECT', 'STOP', 'SAVE',
                'THEN', 'USE',
                'WHILE', 'WHERE', 'WRITE']

INTRINSICS = [
    'ABS', 'ACCESS', 'ACHAR', 'ACOS', 'ACOSH', 'ADJUSTL', 'ADJUSTR', 'AIMAG',
    'AINT', 'ALARM', 'ALL',
    'ALLOCATED', 'AND', 'ANINT', 'ANY', 'ASIN', 'ASINH', 'ASSOCIATED', 'ATAN',
    'ATAN2', 'ATANH',
    'BESJ0', 'BESJ1', 'BESJN', 'BESY0', 'BESY1', 'BESYN', 'BIT_SIZE', 'BTEST',
    'CEILING', 'CHAR', 'CHDIR', 'CHMOD', 'CMPLX', 'COMMAND_ARGUMENT_COUNT',
    'CONJG', 'COS', 'COSH', 'COUNT', 'CPU_TIME', 'CSHIFT', 'CTIME',
    'DATE_AND_TIME', 'DBLE', 'DCMPLX', 'DFLOAT', 'DIGITS', 'DIM',
    'DOT_PRODUCT', 'DPROD', 'DREAL', 'DTIME',
    'DMAX1', 'DMIN1', 'DMOD', 'DSQRT', 'DSIN', 'DCOS', 'DTAN', 'DABS',
    'DATAN', 'DATAN2', 'DEXP', 'DLOG', 'DSINH', 'DCOSH', 'DTANH',
    'EOSHIFT', 'EPSILON', 'ERF', 'ERFC', 'ETIME', 'EXIT', 'EXP', 'EXPONENT',
    'FDATE', 'FGET', 'FGETC', 'FLOAT', 'FLOOR', 'FLUSH', 'FNUM', 'FPUT',
    'FPUTC', 'FRACTION', 'FREE', 'FSEEK', 'FSTAT', 'FTELL',
    'GERROR', 'GETARG', 'GET_COMMAND', 'GET_COMMAND_ARGUMENT', 'GETCWD',
    'GETENV', 'GET_ENVIRONMENT_VARIABLE', 'GETGID', 'GETLOG', 'GETPID',
    'GETUID', 'GMTIME',
    'HOSTNM', 'HUGE',
    'IACHAR', 'IAND', 'IARGC', 'IBCLR', 'IBITS', 'IBSET', 'ICHAR', 'IDATE',
    'IEOR', 'IERRNO',
    'INDEX', 'IDINT', 'INT', 'INT2', 'INT8', 'IOR', 'IRAND', 'ISATTY',
    'ISHFT', 'ISHFTC', 'ITIME',
    'KILL', 'KIND',
    'LBOUND', 'LEN', 'LEN_TRIM', 'LGE', 'LGT', 'LINK', 'LLE', 'LLT', 'LNBLNK',
    'LOC', 'LOG', 'LOG10', 'LOGICAL', 'LONG', 'LSHIFT', 'LSTAT', 'LTIME',
    'MALLOC', 'MATMUL', 'MAX', 'MAX0', 'MAXEXPONENT', 'MAXLOC', 'MAXVAL',
    'MCLOCK', 'MCLOCK8', 'MERGE', 'MIN', 'MIN0', 'MINEXPONENT', 'MINLOC',
    'MINVAL', 'MOD', 'MODULO', 'MOVE_ALLOC', 'MVBITS',
    'NEAREST', 'NEW_LINE', 'NINT', 'NOT', 'NULL',
    'OR',
    'PACK', 'PERROR', 'PRECISION', 'PRESENT', 'PRODUCT',
    'RADIX', 'RANDOM_NUMBER', 'RANDOM_SEED', 'RAND', 'RANGE', 'RAN', 'REAL',
    'RENAME', 'REPEAT', 'RESHAPE', 'RRSPACING', 'RSHIFT',
    'SCALE', 'SCAN', 'SECNDS', 'SECOND', 'SELECTED_INT_KIND',
    'SELECTED_REAL_KIND', 'SET_EXPONENT', 'SHAPE', 'SIGN', 'SIGNAL', 'SIN',
    'SINH', 'SIZE',
    'SLEEP', 'SNGL', 'SPACING', 'SPREAD', 'SQRT', 'SRAND', 'STAT', 'SUM',
    'SYMLNK', 'SYSTEM', 'SYSTEM_CLOCK',
    'tN', 'TANH', 'TIME', 'TIME8', 'TINY', 'TRANSFER', 'TRANSPOSE',
    'TRIM', 'TTYNAM',
    'UBOUND', 'UMASK', 'UNLINK', 'UNPACK', 'VERIFY', 'XOR']

# _____                             ________________________________
# ____/ Global Regular Expressions /_______________________________/
#
"""
  Softer regular expression to ensure scanning of any sources
"""

EMPTY_LINE = re.compile(r'\Z')

F77_COMMENT = re.compile(r'[C!#*]')
F77_CONTINU2 = re.compile(r'(\s{5}\S)(?P<line>.*)')
F90_COMMENT = re.compile(r'(?P<line>([^"]*"[^"]*"[^"!]*|[^\']*\'[^\']*\'[^\'!]*|[^!]*))!{1}(?P<rest>[^\b]*)')
F90_CONTINU1 = re.compile(r'(?P<line>.*)&\Z')
F90_CONTINU2 = re.compile(r'(&\s?)(?P<line>.*)&\Z')
F90_CONTINU3 = re.compile(r'(&\s?)(?P<line>.*)')

VAR_DQUOTS = re.compile(r'(?P<dquot>".*?")')
VAR_SQUOTS = re.compile(r"(?P<squot>'.*?')")
VAR_BRACKS = re.compile(r'(?P<brack>\([\w,*\s+-/:%]*?\))')

# no (\+|\-)? to capture the sign if there ... different from the utils version
VAR_DOUBLEP = re.compile(\
    r'(?P<before>.*?)(?P<number>\b(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')
VAR_INTEGER = re.compile(
    r'(?P<before>.*?)(?P<number>\b(|[^a-zA-Z(,])(?:(\d+)(_\w+|\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')
F90_LOGIC = '(FALSE|TRUE|EQV|EQ|NE|GT|LT|GE|LE|OR|AND|XOR)'
VAR_LOGICAL = re.compile(r'.*?(?P<logic>\.\s?%s\s?\.)' % (F90_LOGIC))
VAR_POINTER = re.compile(r'(?P<before>.*?)(?P<this>%\w+)(?P<after>.*?)\Z')
VAR_WORD = re.compile(r'(?P<before>.*?)(?P<word>\b\w+?\b)(?P<after>.*?)\Z')

VAR_ONLY = re.compile(r'\s*ONLY\s*:(?P<after>.*?)\Z')
VAR_ASSOCS = re.compile(r'(?P<before>.*?)(?P<name>\s?\b\w+?\b\s?=\s?)(?P<after>.*?)\Z')
VAR_EQUALS = re.compile(r'(?P<before>.*?)(?P<name>\s?\b\w+?\b)(?P<value>\s?=[^,]*)(?P<after>.*?)\Z')
VAR_OPERAND = re.compile(r'(?P<before>.*?)(?=[^\s])\W+(?P<after>.*?)\Z')
ARGNAMES = re.compile(r'(?P<name>\b\w+\b)\s*(|\((?P<args>.*)\))\Z')

ITF_TITLE = re.compile(r'\bINTERFACE\b(|\s(?P<name>\w+).*?)\Z')
ITF_ULESS = re.compile(r'\bMODULE\sPROCEDURE\s(?P<name>\w+).*?\Z')
ITF_CLOSE = re.compile(r'\bEND\sINTERFACE\b(|\s(?P<name>\w+).*?)\Z')
USE_TITLE = re.compile(r'\bUSE\b\s(?P<name>\b\w+\b)\s?(|,\s?(?P<after>.*?))\Z')
XTN_TITLE = re.compile(r'.*?\bEXTERNAL\b(.*?::)?\s?(?P<vars>.*?)\Z')
ITZ_TITLE = re.compile(r'\bINTRINSIC\b([^:]?::)?\s?(?P<vars>.*?)\Z')
EQV_TITLE = re.compile(r'\bEQUIVALENCE\b([^:]?::)?\s?(?P<vars>.*?)\Z')
IMPLICIT_NONE = re.compile(r'\bIMPLICIT\sNONE\b\Z')
INC_TITLE = re.compile(r'\bINCLUDE\b\s?(?P<file>.*?)\Z')
CMN_TITLE = re.compile(r'\bCOMMON\b\s?/\s?(?P<name>\w*)\s?/\s?(?P<after>.*?)\Z')
CTN_TITLE = re.compile(r'\bCONTAINS\b\Z')
DEF_TITLE = re.compile(r'\bTYPE\b([^:]?::)?\s(?P<name>\b\w+\b)\Z')
DEF_CLOSE = re.compile(r'\bEND TYPE\b(|\s(?P<name>\b\w+\b))\Z')

ALS_CORE = re.compile(r'(?P<before>.*?)(?P<alias>\b\w([\w%]|\([\w(,:)%]*\))*)\s?=>\s?(?P<link>\b\w+?\b)(?P<after>.*?)\Z')
VAR_FORMAT = re.compile(r'(?P<before>.*?)\d+?\s\bFORMAT\b(?P<after>.*?)\Z')

CLS_TITLE = re.compile(r".*?\bCALL\b\s(?P<name>\b\w+\b)\s?(|\((?P<args>[\w\s\*\+\-\/=,%('.:)]*)\))\Z")  # , re.I)
FCT_TITLE = re.compile(r"(?P<before>.*?)(?P<name>\b\w+\b)\s?\((?P<args>[\w\s\*\+\-\/=,%('.:)]*)\)(?P<after>.*?)\Z")
F90_TYPES = r'(CHARACTER|LOGICAL|INTEGER|REAL|COMPLEX|DOUBLE\s?(PRECISION\s?(COMPLEX|)|COMPLEX))\s?(\**\s?\d+|\**\(.*?\))?|TYPE\s?\([\w\s,=(*)]*?\)'
F90_XPORT = '(PUBLIC|PRIVATE|SAVE|PARAMETER|DATA|SEQUENCE)'
F95_NAME = re.compile(r"(?P<name>\b\w+\b)\s?:\s?(?P<after>.*?)\Z")
TYP_ARGS = re.compile(r'(.*?::)?\s?(?P<vars>.*?)\Z')
TYP_NAME = re.compile(r'(?P<type>(%s))\Z' % (F90_TYPES))
TYP_TITLE = re.compile(r'(?P<type>(%s))\s?(?P<after>.*?)\Z' % (F90_TYPES))
TYP_XPORT = re.compile(r'(?P<type>(%s))\s?(?P<after>.*?)\Z' % (F90_XPORT))

PCL_TITLE = re.compile(r'((?P<type>[\w\s(=*+-/)]*?)|)\b(?P<object>(PROGRAM|FUNCTION|PURE FUNCTION|RECURSIVE FUNCTION|SUBROUTINE|RECURSIVE SUBROUTINE|MODULE|BLOCK DATA))\b\s+(?P<after>.*?)\s*(\bRESULT\b[\s\(]*(?P<result>\w+[\w\s]*)\)?|)\Z')
PCL_CLOSE = re.compile(r'\bEND\b(|\s(?P<object>(PROGRAM|FUNCTION|RECURSIVE FUNCTION|RECURSIVE SUBROUTINE|SUBROUTINE|MODULE))(|\s(?P<name>\w+?)))\Z')

"""
  Ruled regular expression to ensure proper coding for the TELEMAC system
"""

# _____                         ____________________________________
# ____/ FORTRAN Parser Toolbox /___________________________________/
#

def clean_spaces(istr):
    """
    Remove double-spaces, trailing spaces and spaces at the begining of a line
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    return istr.strip().replace('  ', ' ').replace('  ', ' ')


def clean_assoc(istr):
    """
    Deal with '=' associations, and clean spaces
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    while 1:
        ostr = istr
        proc = re.match(VAR_ASSOCS, ostr)
        if proc:
            istr = proc.group('before') + proc.group('after')
        if ostr == istr:
            break
    return clean_spaces(ostr)


def clean_bracks(istr):
    """
    Remove the string between '(' and ')', and clean spaces
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    while 1:
        ostr = istr
        for brack in re.findall(VAR_BRACKS, ostr):
            istr = istr.replace(brack, '')
        if ostr == istr:
            return clean_spaces(ostr)


def clean_equals(istr):
    """
    Deal with '=' associations, and clean spaces
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    while 1:
        ostr = istr
        proc = re.match(VAR_EQUALS, ostr)
        if proc:
            istr = proc.group('before') + proc.group('name') + \
                   proc.group('after')
        if ostr == istr:
            break
    return clean_spaces(ostr)


def clean_formatted(istr):
    """
    Remove the FORMAT statements.
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    proc = re.match(VAR_FORMAT, istr)
    if proc:
        istr = ''
    return istr.strip()


def clean_instruction(istr, list_tag):
    """
    Remove FORTRAN instructions.

    @param istr (string) line to be processed
    @param list_tag (list) List of tag not to clean

    @return a modified (string) line processed
    """
    ostr = ''
    while 1:
        proc = re.match(VAR_WORD, istr)
        if proc:
            ostr = ostr + ' ' + proc.group('before')
            istr = proc.group('after')
            instr = proc.group('word')
            if instr not in list_tag:
                ostr = ostr + ' ' + instr
        else:
            break
    return clean_spaces(ostr)


def clean_logicals(istr):
    """
    Remove LOGICAL instructions. (FALSE|TRUE|EQV|EQ|NE|GT|LT|GE|LE|OR|AND|XOR)
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    while 1:
        ostr = istr
        proc = re.match(VAR_LOGICAL, ostr)
        if proc:
            istr = istr.replace(proc.group('logic'), ' ')
        if ostr == istr:
            break
    return clean_spaces(ostr)


def clean_numbers(istr):
    """
    Remove numbers.
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    while 1:
        ostr = istr
        proc = re.match(VAR_DOUBLEP, ostr)
        if proc:
            istr = proc.group('before') + proc.group('after')
        if ostr == istr:
            break
    while 1:
        ostr = istr
        proc = re.match(VAR_INTEGER, ostr)
        if proc:
            istr = proc.group('before') + proc.group('after')
        if ostr == istr:
            break
    return clean_spaces(ostr)


def clean_operators(istr):
    """
    Remove operators.
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    while 1:
        ostr = istr
        proc = re.match(VAR_OPERAND, ostr)
        if proc:
            istr = proc.group('before') + ' ' + proc.group('after')
        if ostr == istr:
            break
    return clean_spaces(ostr)


def clean_pointers(istr):
    """
    Deal with pointers.
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    while 1:
        ostr = istr
        proc = re.match(VAR_POINTER, ostr)
        if proc:
            istr = proc.group('before') + proc.group('after')
        if ostr == istr:
            break
    return clean_spaces(ostr)


def clean_quotes(istr):
    """
    Clean strings within double or single quotes, and also remove
    internal apostrophees.
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    istr = istr.replace("'''", "'").replace('"""', "'")
    # ~~> Deal with single quotes
    while 1:
        ostr = istr
        for quote in re.findall(VAR_SQUOTS, ostr):
            istr = istr.replace(quote, "''")
        if ostr == istr:
            break
    # ~~> Deal with double quotes (replace by sigle quotes)
    while 1:
        ostr = istr
        for quote in re.findall(VAR_DQUOTS, ostr):
            istr = istr.replace(quote, "''")
        if ostr == istr:
            break
    # ~~> Remove the internal apostrophies
    ostr = ostr.replace("''''", "''").replace("''''", "''")
    return ostr


def clean_empty_quotes(istr):
    """
    Remove empty quotes
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    return istr.replace("''", "").replace('"', "")


def clean_empty_spaces(istr):
    """
    Remove spaces between commas
    @param istr (string) line to be processed
    @return a modified (string) line processed
    """
    return istr.replace(',,', ',')


def refactor_name(istr, nami, namo):
    """
    Rename name (nami) into another name (namo)
    @param istr (string) line to be processed
    @param nami (string) name replaced
    @param namo (string) replaced name
    @return a modified (string) line processed
    """
    n = re.subn(r"(\b)(%s)(\b)" % (nami), namo, istr)[1]
    while n > 0:
        istr, n = re.subn(r"(\b)(%s)(\b)" % (nami), namo, istr)
    return istr


def parse_alias(lines):
    """
    Deal with aliases (over several lines)
    @param lines (list of strings) Fortran lines
    @return modified lines and list of aliases found
    """
    list_alias = []
    count = 0
    core = []
    core.extend(lines)
    for line in lines:
        line = clean_quotes(line)
        proc = re.match(ALS_CORE, line)
        if proc:
            alias = proc.group('alias').strip()
            if not re.match(VAR_POINTER, alias):
                if alias not in list_alias:
                    list_alias.append(alias)
            core[count] = proc.group('before') + ' ' + proc.group('link') + \
                          proc.group('after')
        count = count + 1
    return core, list_alias


def parse_args(ilist):
    """
    Return a cleaned list of variables (arguments) for a call
    @param ilist (string) to be processed
    @return the (list of strings) of individual arguments
    """
    return clean_pointers(\
            clean_numbers(\
             clean_logicals(\
              clean_equals(\
               clean_bracks(\
                clean_quotes(ilist)))))).replace(' ', '').split(',')


def clean_implicit_none(lines):
    """
    Clean IMPLICIT NONE statements
    @param lines (list of strings) of Fortran code to be processed
    @return the (list of strings) processed
    """
    core = []
    core.extend(lines)
    if lines == []:
        return lines, False
    line = lines[0]
    proc = re.match(IMPLICIT_NONE, line)
    if proc:
        core.pop(0)
        return core, True
    return core, False


def parse_declarations(lines):
    """
    Parse and sort out variable declarations whether these are part of a
    COMMON, INTRINSIC, EXTERNAL, ALIASES, predefined TYPEs and local
    variable declarations
    @param lines (list of strings) of Fortran code to be processed
    @return (dict) defined by the keys 'cmn', 'dec', 'itz', 'xtn',
    'als' and 'typ'.
    """
    list_common = []
    list_declar = []
    list_intrinsic = []
    list_external = []
    list_type = {}
    core = []
    core.extend(lines)
    ignoreblock = False
    for line in lines:
        headline = False

        # ~~ TYPE Definition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if not ignoreblock:
            proc = re.match(DEF_TITLE, line)
            if proc:
                ignoreblock = True
                name = proc.group('name')
                list_type.update({name: []})
                core.pop(0)
                headline = True
                continue
        else:
            # /!\ assuming here you have no COMMON, INTRINSIC, INCLUDE,
            #      EXTERNAL, etc in the TYPE block
            proc = re.match(TYP_TITLE, core[0])
            if proc:
                if not re.match(XTN_TITLE, proc.group('after')):
                    proc = re.match(TYP_ARGS, proc.group('after'))
                    if proc:
                        if proc.group('vars') is not None:
                            for k in parse_args(
                                    clean_empty_quotes(proc.group('vars'))):
                                if k != '':
                                    list_type[name].append(k)
            proc = re.match(DEF_CLOSE, line)
            if proc:
                ignoreblock = False
                if proc.group('name') and DEBUG:
                    print('No mane at TYPE END of TYPE ' + name)
            core.pop(0)
            headline = True
            continue

        # ~~ Common Declaration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        proc = re.match(CMN_TITLE, line)
        if proc:
            list_common.append([proc.group('name'),
                                parse_args(proc.group('after'))])
            core.pop(0)
            headline = True
            continue

        # ~~ Private/Public Declarations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        proc = re.match(TYP_XPORT, line)
        if proc:
            # list_common.append([proc.group('name'), \
            #   parse_args(proc.group('after'))])
            core.pop(0)
            headline = True
            continue

        # ~~ INCLUDE Declarations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # you should parse the content from MainWrap
        proc = re.match(INC_TITLE, line)
        if proc:
            core.pop(0)
            headline = True
            continue

        # ~~ Type Declaration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        proc = re.match(TYP_TITLE, line)
        if proc:
            if not re.match(XTN_TITLE, proc.group('after')):
                proc = re.match(TYP_ARGS, proc.group('after'))
                if proc:
                    if proc.group('vars') is not None:
                        for k in parse_args(
                                clean_empty_quotes(proc.group('vars'))):
                            if k != '':
                                list_declar.append(k)
                    core.pop(0)
                    headline = True
                    continue

        # ~~ Intrinsic Declaration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        proc = re.match(ITZ_TITLE, line)
        if proc:
            if proc.group('vars') is not None:
                list_intrinsic.extend(parse_args(proc.group('vars')))
            core.pop(0)
            headline = True
            continue

        # ~~ Equivalence Declaration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        proc = re.match(EQV_TITLE, line)
        if proc:
            core.pop(0)
            headline = True
            continue

        # ~~ External Declaration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        proc = re.match(XTN_TITLE, line)
        if proc:
            if proc.group('vars') is not None:
                list_external.extend(parse_args(proc.group('vars')))
            core.pop(0)
            headline = True
            continue

        # ~~ Reached main core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if not headline:
            break

    return core, {'cmn': list_common, 'dec': list_declar, 'itz': list_intrinsic,
                  'xtn': list_external, 'als': [], 'typ': list_type}


def parse_uses(lines):
    """
    Sort out USEs statements
    @param lines (list of strings) of Fortran code to be processed
    @return (dict) use dependencies
    """
    list_uses = {}
    core = []
    core.extend(lines)
    for line in lines:
        proc = re.match(USE_TITLE, line)
        if proc:
            name = proc.group('name')  # You should not find None here
            args = ''
            # Do you need to distinguish the ONLYs ?
            if proc.group('after') is not None:
                p_after = proc.group('after')
                if p_after != '':
                    proc = re.match(VAR_ONLY, p_after)
                    if proc:
                        args = proc.group('after').strip()
                    else:
                        args = p_after.strip()
            core.pop(0)
            add_to_list(list_uses, name, args)
        else:
            break
    return core, list_uses


def parser_vars(ilist):
    """
    Parse list or variables
    @param ilist (string) to be processed
    @return (dict) vars and their dependencies
    """
    return clean_pointers(\
            clean_numbers(\
             clean_logicals(\
              clean_assoc(\
               clean_bracks(\
                clean_quotes(ilist)))))).replace(' ', ',').split(',')


def parse_calls(lines):
    """
    Parse calls to subroutines
    @param lines (string) to be processed
    @return (dict) subroutiness and their dependencies
    """
    list_calls = {}
    core = []
    core.extend(lines)
    for line in lines:
        # you might not want to clean Logical at this stage /!\
        proc = re.match(CLS_TITLE, clean_pointers(clean_quotes(line)))
        if proc:
            name = proc.group('name')
            args = []
            if proc.group('args') is not None:
                args = parser_vars(proc.group('args'))
            add_to_list(list_calls, name, args)
    return core, list_calls


def parse_functions(lines):
    """
    Parse call to functions
    @param lines (string) to be processed
    @return (dict) functions and their dependencies
    """
    list_fcts = set([])
    list_tags = []
    for line in lines:
        line = clean_quotes(line)
        line = clean_pointers(line)
        proc = re.match(CLS_TITLE, clean_logicals(line))
        if proc:
            line = ''
            if proc.group('args') is not None:
                line = ' '.join(parser_vars(proc.group('args')))
        line = clean_formatted(line)
        proc = re.match(FCT_TITLE, line)
        if proc:
            line = proc.group('before') + proc.group('name') + '(' + \
                   clean_assoc(proc.group('args')) + ')' + proc.group('after')
        line = clean_logicals(line)
        line = clean_numbers(line)
        proc = re.match(F95_NAME, line)
        if proc:
            line = proc.group('after')
            list_tags.append(proc.group('name'))
        line = clean_instruction(line, list_tags)
        line = clean_operators(line)
        line = clean_instruction(line, INSTRUCTIONS)
        line = clean_instruction(line, INTRINSICS)
        if line != '':
            list_fcts = list_fcts | set(line.split())
    return list_fcts


def parse_principal_wrap(lines):
    """
    Parse the principal / high level structure of a Fotran file
    @param lines (string) of Fortran code to be processed
    @return (list of string) as a split of the top level of the Fortran
        file (excluding decalrations), the main name, its arguments and return,
        the interfaces if present na dthe main care
    """
    # you could be parsing the INTERFACE / END INTERFACE as well
    # and in the case of a module ... declarations / definitions
    core = []
    core.extend(lines)
    face = []
    proc = re.match(PCL_TITLE, lines[0])
    if proc:
        name = proc.group('after')
        resu = proc.group('result')
        objt = proc.group('object').split()[-1]
        typ = proc.group('type').strip()
        if typ != '':
            if not re.match(TYP_NAME, typ):
                raise TelemacException(\
                        'Invalid header type {} {} {}\n|    {}\n| ...'
                        ''.format(typ, objt, name, lines[:5]))
        proc = re.match(ARGNAMES, name)
        if proc:
            name = proc.group('name')
            args = []
            if proc.group('args') is not None:
                args = parse_args(proc.group('args'))
            # ~~> header completed
            count = 0
            block = 0
            ctain = 0
            ltain = False
            lface = False
            for line in lines[1:]:  # /!\ does not work with empty lines
                count = count + 1
                # ~~> interface
                if lface:
                    proc = re.match(ITF_CLOSE, line)
                    if proc:
                        lface = False
                    else:
                        proc = re.match(ITF_ULESS, line)
                        if not proc:
                            face.append(line)
                    core.pop(count)
                    count = count - 1
                    continue  # THIS IS TO IGNORE THE LOCAL VARIABLES
                else:
                    proc = re.match(ITF_TITLE, line)
                    if proc:
                        lface = True
                        core.pop(count)
                        count = count - 1
                # ~~> contains
                if ltain:
                    ctain = ctain - 1
                if not ltain:
                    if re.match(CTN_TITLE, line):
                        ltain = True
                        ctain = ctain - 1
                proc = re.match(PCL_CLOSE, line)
                if proc:
                    block = block - 1
                else:
                    proc = re.match(PCL_TITLE, line)
                    if proc:
                        p_type = proc.group('type').strip()
                        if p_type != '':
                            if not re.match(TYP_NAME, p_type):
                                continue
                        block = block + 1
                if block < 0:
                    if ltain:
                        ctain = ctain + 1
                    if proc.group('name') != name:
                        if proc.group('name') is not None:
                            print('Different name at END ' + objt + ' ' + name)
                    if proc.group('object') != objt:
                        if proc.group('object') is not None and \
                                proc.group('object') not in objt:
                            print('Different type at END ' + objt)
                    return core[1:count + ctain], \
                           [objt[0:1], name, args, resu], \
                           face, core[count + ctain + 1:count], core[count + 1:]
            # ~~> Acount for empty MODULE (including only INTERFACE and
            # CONTAINS)
            if ltain and block == 0 and count + 1 >= len(core) - 1:
                if proc.group('name') != name:
                    if DEBUG:
                        print('Different name at END ' + objt + ' ' + name)
                if proc.group('object') != objt:
                    if DEBUG and proc.group('object') not in objt:
                        print('Different type at END ' + objt)
                return core[1:count + ctain], [objt[0:1], name, args, resu], \
                       face, core[count + ctain + 1:count], core[count + 2:]
        else:
            raise TelemacException('Invalid header type for first line '
                                   '{}'.format(lines[0]))

    raise TelemacException(\
            'Invalid header type for first line {}'.format(lines[0]))
    #return  # /!\ this return is just for python parsing


def parse_principal_main(lines, who, typ, name, args, resu):
    """
    Parse the main core of a Fortran file, sorting out its high level structure
    @param lines (string) of Fortran code to be processed
    @param who (string) of Fortran code to be processed
    @param typ (string) of Fortran code to be processed
    @param name (string) of Fortran code to be processed
    @param args (string) of Fortran code to be processed
    @param resu (string) of Fortran code to be processed
    @return (list of string) as a split between the top part of the Fortran
        file (including decalrations), the main name, its arguments and return,
        the interfaces if present na dthe main care
    """
    core = []
    core.extend(lines)
    whi = deepcopy(who)
    whi['uses'] = {}
    whi['vars'] = {}
    whi['calls'] = {}
    whi['called'] = []
    whi['type'] = typ
    whi['name'] = name
    whi['args'] = args
    whi['resu'] = resu

    # ~~ Lists aliases in File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core, alias = parse_alias(core)

    # ~~ Lists uses in File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core, uses = parse_uses(core)
    for k in uses:
        whi['uses'].update({k: []})
        for var in uses[k]:
            add_to_list(whi['uses'], k, var)

    # ~~ Imposes IMPLICIT NONE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core, is_li = clean_implicit_none(core)
    if not is_li and whi['type'] != 'M' and DEBUG:
        print('No IMPLICIT NONE in ' + whi['name'] + ' ' + whi['file'])

    # ~~ Lists declarations in File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core, decs = parse_declarations(core)
    for dec in decs['xtn']:
        if dec in decs['dec']:
            decs['dec'].remove(dec)
    # ~~> remove common variables from decs (keps them in common)
    for dec in decs['cmn']:
        for k in dec[1]:
            if k in decs['dec']:
                decs['dec'].remove(k)
    # ~~> remove arguments from decs (keps them as arguments)
    for dec in args:
        if dec in decs['dec']:
            decs['dec'].remove(dec)
    for k in decs:
        whi['vars'][k] = []
        if isinstance(decs[k], dict):
            add_to_list(whi['vars'], k, '')
            whi['vars'][k] = decs[k]
        else:
            for var in decs[k]:
                add_to_list(whi['vars'], k, var)
    whi['vars']['als'] = alias

    # ~~ Lists calls in File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core, calls = parse_calls(core)
    for k in calls:
        if k in INTRINSICS:
            continue
        whi['calls'][k] = []
        for var in calls[k]:
            add_to_list(whi['calls'], k, var)  # still includes xtn calls

    # ~~ Lists functions in File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fcts = []
    for fct in parse_functions(core):
        if fct not in args:
            if fct not in whi['vars']['dec']:
                if fct not in whi['vars']['als']:
                    found = False
                    for cmn in whi['vars']['cmn']:
                        if fct in cmn[1]:
                            found = True
                    if not found and fct != name:
                        fcts.append(fct)
    whi['functions'] = fcts  # still includes xtn

    return name, whi, core


def del_continueds_f77(lines):
    """
        In order for multiple line of code to be interpreted, remove the
        continuation symbols so every line is self contained --
        This takes into account the possibility for f77
        continuation -- assumes that in between comments have been
        removed already
        Return the set of lines without continuation

        @param lines (string) lines

        @return cmds (list) the set of lines without continuation
    """
    # ~~ Assumes code without comments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cmds = []
    cmd = ''
    for line in lines:
        proc1 = re.match(F77_CONTINU2, line)
        if proc1:
            # /!\ Looking to save on upper space
            cmd = (cmd.rstrip() + proc1.group('line')).strip() \
                .replace('  ', ' ') \
                .replace('  ', ' ')
        else:
            if cmd != '':
                cmds.append(cmd)
            # /!\ Looking to save on upper space
            cmd = line.strip().replace('  ', ' ').replace('  ', ' ')
    cmds.append(cmd)
    return cmds


def del_continueds_f90(lines):
    """
    In order for multiple line of code to be interpreted, remove the
    continuation symbols so every line is self contained --
    This takes into account the possibility for f90
    continuation -- assumes that in between comments have been
    removed already
    Return the set of lines without continuation

    @param lines (string) lines

    @return cmds (list) the set of lines without continuation
    """
    # ~~ Assumes code without comments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cmds = []
    cmd = ''
    cnt = False
    for line in lines:
        proc2 = re.match(F90_CONTINU1, line)
        proc3 = re.match(F90_CONTINU2, line)
        proc4 = re.match(F90_CONTINU3, line)
        if proc3:
            # /!\ Looking to save on upper space
            cmd = clean_spaces(cmd + proc3.group('line'))
            cnt = True
        elif proc2:
            if not cnt:
                if cmd != '':
                    cmds.append(cmd)
                cmd = ''
            # /!\ Looking to save on upper space
            cmd = clean_spaces(cmd + proc2.group('line'))
            cnt = True
        elif proc4:
            # /!\ Looking to save on upper space
            cmd = clean_spaces(cmd + proc4.group('line'))
            cnt = False
        else:
            if cnt:
                # /!\ Looking to save on upper space
                cmd = clean_spaces(cmd + line)
            else:
                if cmd != '':
                    cmds.append(cmd)
                # /!\ Looking to save on upper space
                cmd = clean_spaces(line)
            cnt = False
    if cmd != '':
        cmds.append(cmd)
    return cmds


def del_comments(lines):
    """
    In order for multiple lines of code to be interpreted, remove all
    comments form these, whether they include f77 or f90 comments --
    Return the command lines (without empty lines)

    @param lines (string) lines

    @return cmds (list) the command lines (without empty lines)
    """
    # ~~ Also removes end lines and sets to UPPERCASE ~~~~~~~~~~~~~~~
    cmds = []
    for line in lines:
        line = clean_quotes(line).rstrip().upper()
        proc1 = re.match(F77_CONTINU2, line)  # no strip here
        proc = re.match(F77_COMMENT, line)  # no strip here
        if proc and not proc1:
            cmd = ''
        else:
            while 1:
                cmd = line
                proc = re.match(F90_COMMENT, cmd)
                if proc:
                    line = proc.group('line').rstrip()
                if cmd == line:
                    break
        if cmd != '':
            proc = re.match(EMPTY_LINE, cmd)
            if not proc:
                cmds.append(cmd.replace('\n', '').replace('\r', ''))
    return cmds


def put_scan_content(cfg, file_name, wcw):
    """
    Print to file the the scanning processes as an XML tree.

    @param cfg (dict) configuration dictionary
    @param file_name (string) file name
    @param wcw (dictionary) ?????
    """

    # ~~ High level root ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xml_indents = 0
    xml_sources = ['<' + path.basename(cfg['root']) + '>']

    for mod in wcw:
        # ~~ For each module ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xml_indents += 3
        xml_sources.append(' ' * xml_indents + '<' + mod +
                           ' path="' + wcw[mod]['path'] + '" >')

        for name in wcw[mod]:
            if name == 'path':
                continue
            whoi = wcw[mod][name]
            # ~~ For each file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            xml_indents += 3
            xml_sources.append(' ' * xml_indents + '<' + name + ' type="' +
                               whoi['type'] +
                               '" file="' + whoi['file'] + '">')

            # ~~> uses
            if whoi['uses'] != {}:
                xml_indents += 3
                xml = ' ' * xml_indents + '<uses'
                for use in whoi['uses']:
                    xml += ' ' + use + '="' + ','.join(whoi['uses'][use]) + '"'
                xml += '/>'
                xml_sources.append(xml)
                xml_indents -= 3
            # ~~> declarations
            if whoi['vars']['dec'] != [] or whoi['vars']['xtn'] != [] or \
                    whoi['args'] != []:
                xml_indents += 3
                xml = ' ' * xml_indents + '<vars'
                if whoi['args'] != []:
                    xml += ' args="' + ','.join(whoi['args']) + '"'
                if whoi['vars']['dec'] != []:
                    xml += ' dec="' + ','.join(whoi['vars']['dec']) + '"'
                if whoi['vars']['xtn'] != []:
                    xml += ' xtn="' + ','.join(whoi['vars']['xtn']) + '"'
                xml += '/>'
                xml_sources.append(xml)
                xml_indents -= 3
            # ~~> calls
            if whoi['calls'] != {}:
                xml_indents += 3
                xml = ' ' * xml_indents + '<calls'
                for call in whoi['calls']:
                    xml += ' ' + call + '="' + str(len(whoi['calls'][call])) + \
                           'x"'
                xml += '/>'
                xml_sources.append(xml)
                xml_indents -= 3
            # ~~> calls
            if whoi['called'] != {}:
                xml_indents += 3
                xml = ' ' * xml_indents + '<called ' + \
                      ','.join(whoi['called']) + ' />'
                xml_sources.append(xml)
                xml_indents -= 3

            # ~~> function / unknown
            if whoi['functions'] != []:
                xml_indents += 3
                xml = ' ' * xml_indents + '<unknowns'
                if whoi['functions'] != []:
                    xml += ' functions="' + ','.join(whoi['functions']) + '"'
                xml += '/>'
                xml_sources.append(xml)
                xml_indents -= 3

            # ~~> contains
            if whoi['contains'] != []:
                xml_indents += 3
                xml_sources.append(' ' * xml_indents + '<contains>')

                for contain in whoi['contains']:
                    xml_indents += 3
                    xml_sources.append(' ' * xml_indents + '<'+contain+' />')
                    xml_indents -= 3

                xml_sources.append(' ' * xml_indents + '</contains>')
                xml_indents -= 3

            xml_sources.append(' ' * xml_indents + '</' + name + '>')
            xml_indents -= 3

        xml_sources.append(' ' * xml_indents + '</' + mod + '>')
        xml_indents -= 3

    xml_sources.append('</' + path.basename(cfg['root']) + '>')
    put_file_content(file_name, xml_sources)

    return


def refactor_sources(subset, cmdf, bypass):
    """
    Refactor Sources

    @param subset (list) ???? todo
    @param cmdf (dict) modules dictionary
    @param bypass (boolean) continue with raise exception
    """
    genmod = cmdf['general']['module']
    name = cmdf['general']['name'].split('.')[0]
    # ~~> scanning main names
    refnames = []
    numfiles = 0
    for mod in cmdf:
        if mod == 'general':
            continue
        numfiles += len(cmdf[mod]['files'])
        for fle in cmdf[mod]['files']:
            genpath = path.join(cmdf[mod]['path'], '.' + genmod)
            refnames.append(
                get_principal_wrap_names(path.join(genpath, fle))[0][1])
    # ~~> simplifying subset
    sub = []
    for ssub in subset:
        sub.append(path.basename(sep.join(ssub[0:2])))
    # ~~> replacing main names
    ibar = 0
    pbar = ProgressBar(maxval=numfiles).start()
    for mod in cmdf:
        if mod == 'general':
            continue
        for fle in cmdf[mod]['files']:
            ibar = ibar + 1
            pbar.update(ibar)
            if fle not in sub:
                continue
            genpath = path.join(cmdf[mod]['path'], '.' + genmod)
            # read the file
            f = open(path.join(genpath, fle), 'r')
            srci = ''.join(list(f))
            f.close()
            # refactor the file
            for refname in refnames:
                srci = refactor_name(srci, refname,
                                     refname + '_' + name.upper())
            # write the file
            f = open(path.join(genpath, fle), 'wb')
            f.write(srci)
            f.close()
    pbar.finish()


def scan_sources(cfgdir, cfg, bypass, verbose):
    """
    Scan sources

    @param cfgdir (string) configuration path
    @param cfg (dict)  configuration dictionary
    @param bypass (boolean) continue with raise exception
    @param verbose (boolean) If true display what file is scanned

    @return fic (dict) file dictionary
    @return mdl (dict) ???todo
    @return sbt (dict) ???todo
    @return fct (dict) function dictionary
    @return prg (dict) program dictionary
    @return top (dict) ???todo
    @return odd (dict) ??? todo
    @return wcw (dict) ???todo
    """
    fic = {}
    mdl = {}
    sbt = {}
    fct = {}
    prg = {}
    top = {}
    odd = {}
    wcw = {}

    # ~~ Looking at each file individually ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # TODO: do this in parallel
    for mod in cfg['MODULES'].keys():

        # Dealing with only one set of files
        #  by skipping the duplicate tagged files (active and passive)
        if mod in cfg['ADDONES'] or mod in cfg['ODDONES']:
            continue
        if 'files' not in cfg['MODULES'][mod].keys():
            continue

        # Skipping masacaret module
        if mod == 'mascaret':
            continue

        wcw.update({mod: {'path': cfg['MODULES'][mod]['path']}})
        fic.update({mod: {}})
        # ~~ Scans the sources that are relevant to the model ~~~~~~~~
        src_dir = cfg['MODULES'][mod]['path']  # not anymore

        # In case of subdirectories loop on the subdirectories as well
        file_list = cfg['MODULES'][mod]['files']
        if len(file_list) == 0:
            raise TelemacException('... found an empty module: '
                                   '{}'.format(mod))
        o_dir = path.join(cfg['MODULES'][mod]['path'], cfgdir)

        print('... now scanning ' + path.basename(cfg['MODULES'][mod]['path']))

        ibar = 0
        pbar = ProgressBar(maxval=len(file_list)).start()
        # /!\ if you need to print within this loop, you now need to use
        #     pbar.write(text, ibar) so the progress bar moves along
        for fle in file_list:

            ibar = ibar + 1
            pbar.update(ibar)
            fle = path.join(src_dir, fle)

            if not mod in fic:
                fic.update({mod: {}})
            fic[mod].update({fle: []})
            who = {'path': src_dir,
                   'file': fle.replace(src_dir + sep, '').replace(sep, '|'),
                   'libname': mod,
                   'type': '',
                   'name': '',
                   'args': [],
                   'resu': '',
                   'contains': [],
                   'uses': {},
                   'vars': {},
                   'calls': {},
                   'called': [],
                   'functions': [],
                   'rank': 1,
                   'time': 0}
            if path.isdir(o_dir):
                if cfg['COMPILER']['REBUILD'] == 2:
                    who['time'] = 0
                else:
                    who['time'] = is_newer(
                        fle,
                        path.join(o_dir,
                                  path.splitext(path.basename(fle))[0] +
                                  cfg['SYSTEM']['sfx_obj'].lower()))

            if verbose:
                print("Parsing: ", fle)
            src_file = open(fle, 'r', encoding='utf-8')
            if path.splitext(who['file'])[1].lower() in ['.f90', '.f95']:
                # Strips the F90+ commented lines
                flines = del_continueds_f90(del_comments(src_file))
            else:
                # Strips the F77 commented lines
                flines = del_continueds_f77(del_comments(src_file))
            src_file.close()  # and deletes the continuation characters

            core = flines
            found = False
            # ignores what might be in the file after the main program
            while core != [] and not found:

                # ~~ Parse Main Structure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                code, what, face, ctns, core = parse_principal_wrap(core)
                name, whoi, _ = parse_principal_main(code, who,
                                                     what[0], what[1],
                                                     what[2], what[3])
                if what[0] == 'P':
                    top = add_to_list(top, name, whoi['libname'])
                # others
                if what[0] == 'P':
                    prg = add_to_list(prg, name, whoi['libname'])  # main program
                if what[0] == 'P':
                    found = bypass
                if what[0] == 'S':
                    sbt = add_to_list(sbt, name, whoi['libname'])  # subroutine
                if what[0] == 'M':
                    mdl = add_to_list(mdl, name, whoi['libname'])  # module
                if what[0] == 'F':
                    fct = add_to_list(fct, name, whoi['libname'])  # function
                fic[mod][fle].append(name)
                # ~~> you may not need this afterall or do you ?
                while face != []:
                    fcode, f_what, _, _, face = parse_principal_wrap(face)
                    _, whof, _ = parse_principal_main(fcode, who, f_what[0],
                                                      f_what[1], f_what[2],
                                                      f_what[3])
                    for k in whof['uses']:
                        for var in whof['uses'][k]:
                            add_to_list(whoi['uses'], k, var)
                while ctns != []:  # contains fcts & subs
                    ccode, c_what, _, _, ctns = parse_principal_wrap(ctns)
                    cname, whoc, _ = parse_principal_main(ccode, who, c_what[0],
                                                          c_what[1], c_what[2],
                                                          c_what[3])
                    whoi['contains'].append(cname)
                    # ~~> not unkown anymore
                    if cname in whoi['functions']:
                        whoi['functions'].remove(cname)
                    if cname in whoi['calls']:
                        whoi['calls'].pop(cname)
                    if c_what[0] == 'S':
                        sbt = add_to_list(sbt, cname, whoi['libname']) #subroutine
                    if c_what[0] == 'F':
                        fct = add_to_list(fct, cname, whoi['libname'])  # function
                    # ~~> transfer global context to identify unknowns
                    for k in whoc['uses']:
                        for var in whoc['uses'][k]:
                            add_to_list(whoi['uses'], k, var)
                    for k in deepcopy(whoc['calls']):
                        if k in whoi['vars']['xtn']:
                            whoc['calls'].pop(k)
                        elif k in whoi['contains'] or k in whoi['calls']:
                            whoc['calls'].pop(k)
                        else:
                            whoi['calls'].update({k: whoc['calls'][k]})
                    for k in whoc['calls']:
                        for var in whoc['calls'][k]:
                            add_to_list(whoi['calls'], k, var)
                    for k in deepcopy(whoc['functions']):
                        for var in whoi['vars']['cmn']:
                            if k in var[1]:
                                whoc['functions'].remove(k)
                        if k in whoc['functions']:
                            if k in whoi['args'] or \
                                    k in whoi['vars']['dec'] or \
                                    k in whoi['vars']['xtn']:

                                whoc['functions'].remove(k)
                            elif k in whoi['contains'] or \
                                    k in whoi['functions']:
                                whoc['functions'].remove(k)
                            else:
                                whoi['functions'].append(k)

                whoi['vars'].update({'use': {}})
                wcw[mod].update({name: whoi})  # ~~ All ~~~~~~~~~~

                # separating the top of the tree from the odd ones
                for m in cfg['ODDONES']:
                    if path.basename(fle) in cfg['ODDONES'][m]:
                        odd = add_to_list(odd, name, m + '.' + name.lower())
                        # the combination m+name is unique, and therefore
                        # allows you to specify the path
                        if m not in cfg['ADDONES']:
                            top = add_to_list(top, name, m + '.' + name.lower())
                            wcw.update(
                                {m + '.' + name.lower():
                                     {'path': cfg['MODULES'][mod]['path'],
                                      name: whoi}})
                # /!\ tree tops now also include some of the ODD ONES
                # (if there not
                # also amongst the ADD ONES)
                for m in cfg['ADDONES']:
                    if path.basename(fle) in cfg['ADDONES'][m]:
                        top = add_to_list(top, name, m + '.' + name.lower())
                        # the combination m+name is unique, and therefore
                        # allows you to specify the path
                        wcw.update(
                            {m + '.' + name.lower():
                                 {'path': cfg['MODULES'][mod]['path'],
                                  name: whoi}})

        pbar.finish()

    # ~~: wcw.keys()
    #    holds the list of all primary and secondary libraries
    #    where source files were found. For each library, the list of files
    #    scanned is provided by cfg['MODULES'][mod]['files']
    # ~~: wcw[mod].keys()
    #    holds the 'path' and well as the list of names of all subroutines
    #    / functions / programs / modules contained in that mod
    #    (uppercase names)
    # ~~: wcw[mod][name].keys()
    #    holds a split of the analysis of the content of 'name' with the
    #    following fields:
    #  ['functions', 'libname', 'name', 'vars', 'contains', 'args', 'rank',
    #    'resu', 'uses', 'calls', 'file', 'time', 'path', 'type', 'called']

    # ~~ Cross-referencing CALLS together ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # For those CALLs stored in 'calls' but not part of the system:
    #    move them from 'calls' to 'function' (outsiders remain)
    #
    for mod in wcw.keys():
        iter_key_wcw = list(wcw[mod].keys())
        for name in iter_key_wcw:
            # for name in copy(wcw[mod].keys()):
            if name == 'path':
                continue
            who = wcw[mod][name]
            iter_key_who = list(who['calls'].keys())
            for sub in iter_key_who:
                # for s in copy(who['calls'].keys()):
                # separate items in the top tree by declaring then as external
                if sub in top.keys() and sub not in odd.keys():
                    del wcw[mod][name]['calls'][sub]
                # calls bot in the list of subroutines msut be functions (?)
                elif sub not in sbt:
                    del wcw[mod][name]['calls'][sub]
                    wcw[mod][name]['functions'].append(sub)

    # ~~ Cross-referencing FUNCTIONS together ~~~~~~~~~~~~~~~~~~~~~~~
    for mod in wcw.keys():
        for name in wcw[mod].keys():
            if name == 'path':
                continue
            who = wcw[mod][name]
            f, use = sort_functions(who['functions'], who['vars']['use'], wcw,
                                    mdl, who['uses'])
            # /!\ because this is a dico-list, updating who updates wcw
            who['functions'] = f
            who['vars']['use'].update(use)

    for mod in wcw:
        for name in wcw[mod]:
            if name == 'path':
                continue
            who = wcw[mod][name]
            for f in who['vars']['xtn']:
                if f not in who['functions']:
                    if DEBUG:
                        print(f + ' declared but not used in ' + who['name'])
                    who['functions'].append(f)
            for f in fct:
                while f in who['functions']:
                    who['functions'].remove(f)
                    if f not in top.keys() or f in odd.keys():
                        who['calls'].update({f: [['']]})

    # ~~ Sort out referencing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Fill-in the 'called' category
    for mod in wcw.keys():
        for name in wcw[mod].keys():
            if name == 'path':
                continue
            for call in wcw[mod][name]['calls']:
                if call in fct:
                    for use in fct[call]:
                        if call in wcw[use]:
                            wcw[use][call]['called'].append(name)

    if DEBUG:
        put_scan_content(cfg, cfg['root'] + sep + 'scanTELEMAC.xml', wcw)

    return fic, mdl, sbt, fct, prg, top, odd, wcw


def sort_functions(ifcts, iuses, llist, mods, xuses):
    """
    Change order function list

    @param ifcts (dict) todo
    @param iuses (dict) todo
    @param llist (dict) ??? todo
    @param mods (list) ???todo
    @param xuses (list) ???todo

    @return ifcts (dictionary)
    @return iuses (dictionary)
    """
    ofcts = []
    ofcts.extend(ifcts)
    for fct in ifcts:
        for use in xuses:
            if use not in mods:
                continue
            if fct in llist[mods[use][0]][use]['vars']['dec']:
                ofcts.remove(fct)
                add_to_list(iuses, use, fct)
                break
            if fct in llist[mods[use][0]][use]['vars']['als']:
                ofcts.remove(fct)
                add_to_list(iuses, use, fct)
                break
    ifcts = ofcts
    for use in xuses:
        if use in mods and ifcts != []:
            ifcts, iuses = sort_functions(ifcts, iuses, llist, mods,
                                          llist[mods[use][0]][use]['uses'])

    return ifcts, iuses


# _____                     _____________________________________________
# ____/ Coding Rules  /____________________________________________/
#
def scan_enunciation(name):
    """
    Display differenciating between name and name.cmdf

    @param name (string) name file

    @return name (string) name file
    """
    # ~~ get file content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lines = get_file_content(name)

    # ~~ scan indents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # ~~ scan comments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # ~~ proposed file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dif_file = name + '.mdf'
    put_file_content(dif_file, lines)

    # ~~ differentiation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~> default options
    options = Values()
    options.unified = False
    options.ndiff = False
    options.html = True
    options.ablines = 3
    options.context = False
    # ~~> differentiate
    print('\n\nDifferenciating:\n     +> ' + name + '\nand +> ' + dif_file + \
          '\n' + '~' * 72 + '\n')
    # ~~> use of writelines because diff is a generator
    diff = diff_text_files(name, dif_file, options)
    # ~~> write differentiator
    if options.html:
        html_file = path.splitext(dif_file)[0] + '.html'
        f = open(html_file, 'wb')
        f.writelines(diff)
        f.close()
        print('\n        ~> produced: ' + html_file, '\n')

    return name

# _____                      ____________________________________________
# ____/ DOXYGEN parse  /___________________________________________/
#

DOXY_COMMENT = re.compile(r'!>(?P<doxy>.*)\s*\Z')
DOXY_TAG = re.compile(r'!>\s*?@(?P<tag>[\w\[,\]]+)(?P<after>.*)\s*\Z')
DOXY_TITLE = re.compile(r'\s+(?P<title>.*)\s*\Z')


def get_next_doxy_tag(core):
    """
        Parse a list of entry lines, removing the lines that are
        understood to be part of a tag -- the lines are removed
        from the 'core' iteratively -- Return one complete tag

        @param core (list) list of entry lines

        @return tag (string) tag
        @return name (string) name
        @return ltag (list) tag list
        @return core (list) lines list
    """
    ltag = []
    tag = ''
    name = ''
    found = False
    while core != []:
        proc = re.match(DOXY_TAG, core[0])
        if proc:
            if not found:
                tag = proc.group('tag')
                proc = re.match(DOXY_TITLE, proc.group('after'))
                if proc:
                    name = proc.group('title')
                found = True
            else:
                break
        if found:
            ltag.append(core[0])
        core.pop(0)

    return tag, name, ltag, core


def parse_doxy_tags(core):
    """
    Parse a list of entry lines, removing the doxygen blocks
    listing them out by tag names -- the doxygen blocks are removed
    from the 'core' iteratively -- Return the list of tags by blocks

    @param core (list) list of entry lines

    @return tags (list) list of tags by blocks
    """
    tags = []
    while 1:
        tname, title, field, core = get_next_doxy_tag(core)
        if tname == '':
            break
        tags.append([tname, title, field])
    return tags


# _____                      ____________________________________________
# ____/ FORTRAN parse  /___________________________________________/
#

VAR_COMMENT = re.compile(r'[C!#*]\s*?[|!]\s*?(?P<vars>[\s\w,()]*)(?P<inout>(|[|!->=<\s]*[|!]))(?P<after>(|[^|!]*))\s*[|!]?\s*\Z')


def parse_fort_header(core):
    """
    Parse fortran header

    @param core (list) list of entry lines

    @return docs (list) docs information
    @return vrs (dict) ??? todo
    @return core (list)lines list
    """
    docs = []
    vrs = {}
    while 1:
        line = del_comments([core[0]])
        if line == []:
            proc = re.match(VAR_COMMENT, core[0])
            if proc:
                # print( core[0].rstrip() )
                if proc.group('vars').strip() != '':
                    var = proc.group('vars').strip()
                    val = proc.group('after').strip()
                    ino = proc.group('inout').strip()
                    if ino == '!!' or ino == '||':
                        ino = '<>'
                    vrs.update({var: [ino, [val]]})
                elif proc.group('after').strip() != '':
                    vrs[var][1].append(proc.group('after').strip())
            core.pop(0)
            continue
        line = line[0].rstrip()
        proc = re.match(PCL_TITLE, line)
        if proc:
            typ = proc.group('type').strip()
            if typ != '':
                if not re.match(TYP_NAME, typ):
                    break
            docs.append(core[0].rstrip())
            core.pop(0)
            continue
        proc = re.match(USE_TITLE, line)
        proc = proc or re.match(IMPLICIT_NONE, line)
        proc = proc or re.match(ITF_TITLE, line)
        proc = proc or re.match(DEF_TITLE, line)
        proc = proc or re.match(CMN_TITLE, line)
        proc = proc or re.match(TYP_XPORT, line)
        proc = proc or re.match(INC_TITLE, line)
        proc = proc or re.match(TYP_TITLE, line)
        proc = proc or re.match(ITZ_TITLE, line)
        proc = proc or re.match(XTN_TITLE, line)
        if proc:
            break
        docs.append(core[0].rstrip())
        core.pop(0)

    return docs, vrs, core


def get_principal_wrap_names(dif_file):
    """
    Get name of subroutine/functions defined in file

    @param dif_file (string) name file

    @return pFiles (list) list of [file_name , function_name]
    """
    # Filter most unuseful
    src_file = open(dif_file, 'r', encoding='utf-8')
    if path.splitext(path.basename(dif_file))[1].lower() in ['.f90', '.f95']:
        # Strips the F90+ commented lines
        flines = del_continueds_f90(del_comments(src_file))
    else:
        # Strips the F77 commented lines
        flines = del_continueds_f77(del_comments(src_file))
    src_file.close()  # and deletes the continuation characters
    # Identify main items
    p_files = []
    while flines != []:
        _, what, _, _, flines = parse_principal_wrap(flines)
        p_files.append([what[0], what[1]])
    return p_files


def filter_principal_wrap_names(u_names, s_files):
    """
    Idenfy the name of function associated with files

    @param u_names (list) list to check ???? todo
    @param s_files (list) file names list

    @return o_files (dictionary)   ???? todo
    """
    o_files = {}
    for s_file in s_files:
        src_file = open(s_file, 'r', encoding='utf-8')
        if path.splitext(path.basename(s_file))[1].lower() in ['.f90', '.f95']:
            flines = del_continueds_f90(del_comments(src_file))
        else:
            flines = del_continueds_f77(del_comments(src_file))
        src_file.close()
        _, what, _, _, flines = parse_principal_wrap(flines)
        if what[1] in u_names:
            o_files.update({what[1]: s_file})
    return o_files


# Note that the spaces are kept in the 'after' text for possible formatting
DOXY_TAGS = '(brief|note|warning|history|bug|code)'
DOXY_COMMENT = re.compile(r'\s*!(?P<name>%s\b)(?P<after>.*?)\s*\Z'%(DOXY_TAGS))
DOXY_COMMENT_ADD = re.compile(r"\s*!\+(?P<after>.*?)\s*\Z")


def parse_doxy_header(body):
    """
    Format text for Doxygen header

    @param body (string) text line

    @return name (string) name
    @return tags (list) ??? todo
    """
    # It is now assumed that DOXY_TAGS could be part of the main
    # core of the Fortran.
    tags = []
    tfound = False
    tcount = -1
    core = []
    core.extend(body)
    bcount = -1
    nfound = False
    # ~~ Parse Doxygen Tags ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    while core != []:
        line = core[0].rstrip()
        bcount = bcount + 1

        # ~~> parsing the tags
        proc = re.match(DOXY_COMMENT, line)
        if proc:
            tcount = tcount + 1
            tags.extend([[]])
            tags[tcount].append(proc.group('name'))
            tags[tcount].append([proc.group('after')])
            tfound = True
            core.pop(0)
            body.pop(bcount)
            bcount = bcount - 1
            continue
        proc = re.match(DOXY_COMMENT_ADD, line)
        if proc and tfound:
            tags[tcount][1].append(proc.group('after'))
            core.pop(0)
            body.pop(bcount)
            bcount = bcount - 1
            continue
        tfound = False

        # ~~> parsing the name of the program for later reference
        if not nfound:
            proc = re.match(PCL_TITLE, line)
            if proc:
                proc = re.match(VAR_WORD, proc.group('after'))
                if proc:
                    name = proc.group('word')
                nfound = True
        core.pop(0)

    return name, tags


def parse_doxy_wrap(lines):
    """
        Split a set of lines (content of a file) into a Doxygen header
        the definition of the FORTRAN entity and the text between the
        the definition oan the core of the FORTRAN entity -- The wrap
        iteratively goes through all included entities and sub-entities
        - lines contains the content of the file
        - icount is the number of entities included (functions, subroutines,
        etc.)

        @param lines (string) lines content of a file

        @return doxy(list) Doxygen header
    """
    core = []
    core.extend(lines)
    wrap = []
    wrap.extend([[]])
    blevel = 0
    bcount = 0

    # ~~ Split the content of the file by blocks ~~~~~~~~~~~~~~~~~~~~
    while core != []:
        line = core[0].rstrip()
        wrap[bcount].append(line)
        core.pop(0)

        proc = re.match(PCL_CLOSE, line)
        if proc:
            blevel = blevel - 1
            if blevel == 0:
                bcount = bcount + 1
                wrap.extend([[]])
            continue
        proc = re.match(PCL_TITLE, line)
        if proc:
            blevel = blevel + 1
            continue

    # ~~ Clean the wrap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # remove the the last item of the wrap,
    # which is empty since it gets set just after the PCL_CLOSE
    # >wrap.pop(len(wrap)-1)

    # ~~ Parse the DOXYGEN tags ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    doxy = []
    name = ''
    for what in wrap:
        if what != []:
            iname, tag = parse_doxy_header(what)
            if name == '':
                name = iname  # /!\ here takes the first name of the file
            doxy.extend([(name, iname, tag, what)])

    return doxy


def trim_tree(name, lname, lst, rebuild, maksystel):
    """
    Scan the trre of FORTRAN calls to populate who is "up" and who is "dw"
        (calle and called). Create the order of dependency of the libraries.

    @param name (string) Name of the subroutine or function or module as read
        in the FORTRAN file (in upper case)
    @param lname (string) Name of the library withinwhich the name is
    @param lst (dict) Dictionary containing the result of the parsing of the
        FORTRAN. The principal keys include the name of the dependencies as well
        as "path" (which is the path to the file)
    @param maksystel (dict) Has three keys: "add", "tag" and "deps", where
        "add" include the list of file to be compiled
        "tag" include the list of files to be ignored (already compiled)
        "deps" include the list of files not part of the TELEMAC system
    @param rebuild (integer) Level of reconstruction ?

    @return lst (dict) the tree of dependencies up to that point.
    """
    liborder = []
    #
    # ~~: lrank
    #  provides 'dw' and 'up' lists of calls / called respectively for each
    #  lib in in the tree of dependencies of name. For instance,
    #  >  lrank['special']['dw'] is probably empty []
    #  >  lrank[lname]['up'] should also be empty []
    lrank = {}
    _ = get_tree(name, lname, lst, [], lrank, rebuild, maksystel)
    #
    # ~~: libdws
    #  the list of libs touced by name, incuding lname, but not necessarily
    #  in order of dependencies
    libdws = [lrank.keys()]
    #
    # gradually emptying lrank until none are left, sorting out ups and downs
    while lrank != {}:
        for libdw in libdws:
            for lib in libdw:
                if lrank[lib]['up'] == [] and lib not in liborder:
                    liborder.insert(0, lib)
        if liborder == []:
            raise TelemacException(' ... found recursive loop with '
                                   '{}:{}'.format(name, lrank))
        libdws = []
        for lib in liborder:
            if lib in lrank:
                libdws.append(lrank[lib]['dw'])
                del lrank[lib]
            for ldw in lrank:
                if lib in lrank[ldw]['up']:
                    lrank[ldw]['up'].remove(lib)
    return liborder


# _____                         __________________________________________
# ____/ TOP OF THE TREE  /_________________________________________/
#

def get_tree(name, lname, lst, level, lrank, rebuild, maksystel):
    """
    Recursively scan through the list of calls (subroutines and functions) and
        uses to document who calls who.

    @param name (string) Name of the subroutine or function or module as read
        in the FORTRAN file (in upper case)
    @param lname (string) Name of the library withinwhich the name is
    @param lst (dict) Dictionary containing the result of the parsing of the
        FORTRAN. The principal keys include the name of the dependencies as well
        as "path" (which is the path to the file)
    @param level (list) List of pairs of string, showing the call path from
        the main program, down to the current subroutine/function/module.
    @param lrank (dict) Includes the "up" and "dw" keys, which will be populated
        within trim_tree. Empty for now.
    @param maksystel (dict) Has three keys: "add", "tag" and "deps", where
        "add" include the list of file to be compiled
        "tag" include the list of files to be ignored (already compiled)
        "deps" include the list of files not part of the TELEMAC system
    @param rebuild (int) TODO

    @return lst (dict) the tree of dependencies up to that point.
    """

    # debug = False
    # ~~ Recursive tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if level != []:
        listname = [v[0] for v in level]
        if name in listname:
            print('found recursive loop with ' + name + ':' + \
                  ' => '.join(zip(*level)[0]))
            return lst[lname][name]['time']
    # ~~ New leaf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    level.append([name, lname])
    # ~~ Ranking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if lname not in lrank:
        lrank.update({lname: {'up': [], 'dw': []}})
    for _, lib in level:
        if lib not in lrank[lname]['up'] and lib != lname:
            lrank[lname]['up'].append(lib)
    # ~~ prints the tree to screen:
    time = lst[lname][name]['time']
    for use in sorted(lst[lname][name]['uses']):
        libname = lname
        if use not in lst[lname]:
            for lib in lst:
                if use in lst[lib] and lib != lname:
                    libname = lib
        if use in lst[libname]:
            if libname not in lrank[lname]['dw'] and libname != lname:
                lrank[lname]['dw'].append(libname)
            t_tree = get_tree(use, libname, lst, level, lrank, rebuild,
                              maksystel)
            level.pop()
            if rebuild < 3:
                time = time * t_tree
    for call in sorted(lst[lname][name]['calls']):
        libname = lname
        if call not in lst[lname]:
            for lib in lst:
                if call in lst[lib] and lib != lname:
                    libname = lib
        if call in lst[libname]:
            if libname not in lrank[lname]['dw'] and libname != lname:
                lrank[lname]['dw'].append(libname)
            t_tree = get_tree(call.strip(), libname, lst, level, lrank,
                              rebuild, maksystel)
            level.pop()
            if rebuild < 3:
                time = time * t_tree
    for function in sorted(lst[lname][name]['functions']):
        libname = lname
        if function not in lst[lname]:
            for lib in lst:
                if function in lst[lib] and lib != lname:
                    libname = lib
        if function in lst[libname]:
            if libname not in lrank[lname]['dw'] and libname != lname:
                lrank[lname]['dw'].append(libname)
            t_tree = get_tree(function.strip(), libname, lst, level, lrank,
                              rebuild, maksystel)
            level.pop()
            if rebuild < 3:
                time = time * t_tree
    lst[lname][name]['time'] = time
    if time == 0:
        if [lst[lname][name]['file'], lname] not in maksystel['list']:
            maksystel['list'].append([lst[lname][name]['file'], lname])
    else:
        if [name, lname] not in maksystel['tag']:
            maksystel['tag'].append([name, lname])

    return lst[lname][name]['time']


# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Sebastien E. Bourban; Noemie Durand"
__date__ = "$19-Jul-2010 08:51:29$"


