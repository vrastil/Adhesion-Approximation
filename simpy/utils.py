"""
'utils.py' module contains various help functions
without additional dependencies on simpy modules
"""

###################################
# modules
###################################
from __future__ import print_function

import os
import sys
import fnmatch

import traceback
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import TerminalFormatter


###################################
# colorful exception info
###################################

def print_exception(file=sys.stdout):
    """ print catched exception with colors """
    tbtext = traceback.format_exc()
    lexer = get_lexer_by_name("pytb", stripall=True)
    formatter = TerminalFormatter()

    print("\n")
    print("=" * 110)
    file.write(highlight(tbtext, lexer, formatter))
    print("=" * 110)

###################################
# FIND, SORT, SLICE
###################################

def create_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

def del_duplicate(seq):
    # type: (Seq) -> List
    """remove duplicate elements in sequence while preserving order, return list"""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def get_files_in_traverse_dir(a_dir, patterns):
    # type: (str, str) -> List[str]
    """ return list of all files in directory which matches 'patterns'
    support Unix filename pattern matching ('*', '?', [seq], [!seq])
    and multiple option in 'patterns' (space delimetered) """

    return list(set([ # throw away duplicate files
        os.path.join(root, name) # full file name
        for root, _, files in os.walk(a_dir) # go through all subdirectores
        for pattern in patterns.split() # if multiple patterns given
        for name in fnmatch.filter(files, pattern) # pattern matching
        ]))

def sort_lists(*lists):
    return zip(*sorted(zip(*lists), reverse=True))

def sort_get_z(files, app, skip_init=False):
    # type: (List[str], struct.SimInfo) -> List[str], List[TypeVar(str, float)]
    zs = []
    for a_file in files:
        if app + '_z' in a_file:
            zs.append(float(a_file[a_file.index(app + '_z') + len(app+'_z'):-4]))
        elif app + '_init' in a_file and not skip_init:
            zs.append('init')
        else:
            print("WARNING! Skipping file '%s', unknown format." % a_file)
    return sort_lists(zs, files)

###################################
# COLORED TEXT
###################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

len_info = 0

def print_info(info):
    global len_info
    len_info = len(info)
    print(bcolors.BOLD + '*'*len_info + '\n' + info + '\n' + '*'*len_info + bcolors.ENDC)

def print_info_end():
    global len_info
    print(bcolors.BOLD + '*'*len_info + '\nAll runs analyzed!' + bcolors.ENDC)
    len_info = 0

def print_done():
    print(bcolors.OKGREEN + "[Done]" + bcolors.ENDC)

def print_skipped(done=False):
    if done:
        print(bcolors.OKBLUE + "[Skipped]  (already done)" + bcolors.ENDC)
    else:
        print(bcolors.OKBLUE + "[Skipped]" + bcolors.ENDC)

def print_skipped_miss():
    print(bcolors.WARNING + "[Skipped]  (missing data)" + bcolors.ENDC)