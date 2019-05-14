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
from bson.son import SON

import traceback
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import TerminalFormatter
from IPython.display import display, Math

###################################
# sorted dictionary
###################################
def get_sorted_bson(a_dict):
    res = SON()
    for k, v in sorted(a_dict.items()):
        if isinstance(v, dict):
            res[k] = get_sorted_bson(v)
        else:
            res[k] = v
    return res

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

def get_z_from_file(a_file, app):
    if app + '_z' in a_file:
        return float(a_file[a_file.index(app + '_z') + len(app + '_z'):-4])
    elif app + '_init' in a_file:
        return 'init'
    else:
        return None

def sort_get_z(files, app):
    # type: (List[str], struct.SimInfo) -> List[str], List[TypeVar(str, float)]
    zs = []
    for a_file in files:
        z = get_z_from_file(a_file, app)
        if z is not None:
            zs.append(z)
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

def print_info(info, math_mode=False, app=None):
    global len_info
    len_info = len(info)
    if math_mode:
        len_info += len(math_mode)
        if app == 'CHI':
            len_info -= 25
        elif app:
            len_info -= 10
        print(bcolors.BOLD + '*'*len_info)
        display(Math(r'\textrm{' + info + r'}' + math_mode))
        print(bcolors.BOLD + '*'*len_info + bcolors.ENDC)
    else:
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

def print_warning(warning):
    print(bcolors.FAIL + warning + bcolors.ENDC)

###################################
# ACCESS DICT ITEMS VIA DOT
###################################

class Map(dict):
    """
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.iteritems():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.iteritems():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]