from datetime import datetime
import sys
import gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json

from . import *
from . import plot
from . import data


def compare(sim_info_1, sim_info_2):
    attributes = [a for a in dir(sim_info_1) if not a.startswith('__') and not callable(getattr(sim_info_1, a))
                  and a != "seed" and a != "results" and a != "file" and a != "dir" and a != "res_dir" and a != "cosmo"]
    for att in attributes:
        att1 = getattr(sim_info_1, att)
        att2 = getattr(sim_info_2, att)

        if isinstance(att1, dict):
            for key in att1:
                if att1[key] != att2[key]:
                    return False
        else:
            if att1 != att2:
                return False

    return True


def insert(a_sim_info, sep_files):
    for sep_sim_infos in sep_files:
        if compare(a_sim_info, sep_sim_infos[0]):
            sep_sim_infos.append(a_sim_info)
            return True
    return False


def stack_all(out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/'):
    # get all runs
    files = get_files_in_traverse_dir(out_dir, 'sim_param.json')
    sim_infos = []
    for args in files:
        sim_infos.append(data.SimInfo(*args))

    # separate files according to run parameters
    sep_files = []
    for a_sim_info in sim_infos:
        if not insert(a_sim_info, sep_files):
            sep_files.append([a_sim_info])

    # print info about separated files
    num_all_runs = 0
    num_all_sep_runs = 0
    num_sep_runs = 0
    for sep_sim_infos in sep_files:
        num_all_runs += len(sep_sim_infos)
        if len(sep_sim_infos) > 1:
            num_sep_runs += 1
            num_all_sep_runs += len(sep_sim_infos)
    print "There are in total %i different runs, from which %i share the same parameters, constituting %i group(s) eligible for smoothing:" % (
        num_all_runs, num_all_sep_runs, num_sep_runs)
    for sep_sim_infos in sep_files:
        print "\n" + sep_sim_infos[0].info_tr()
        for a_sim_info in sep_sim_infos:
            print "\t" + a_sim_info.file
