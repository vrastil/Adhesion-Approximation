"""
'database.py' module serves for creating, updating and managinh databse with
all data from simulations
"""

# system packages
from __future__ import print_function
import os
import sys
import pymongo
from bson.son import SON
from getpass import getpass
import json
import datetime
import numpy
try:
    from urllib.parse import urlparse
except ImportError:
     from urlparse import urlparse

# simpy packages
from .struct import RESULTS_DIRS
from .utils import *

# version independent raw_input
if sys.version[0]=="3": raw_input=input

def create_database(host='localhost', port=27017, user='admin'):
    """create database and user admin in it,
    databse should be started without authentication this first time"""

    print("Database server should be running without authentication on host '%s' and port '%i'" % (
        host, port))
    raw_input("Press Enter to continue...")
    client = pymongo.MongoClient(host, port)
    db = client['admin']
    try:
        user = raw_input("Username (default '%s'):" % user) or user
        db.command("createUser", user, pwd=getpass(prompt='Password:'), roles=[{'role':'userAdminAnyDatabase','db':'admin'}, "readWriteAnyDatabase"])
    except pymongo.errors.DuplicateKeyError as e:
        print(e)
    else:
        print("User '%s' successfully created. Restart the server with authentication enabled." % user)

def connect_db(host='localhost', port=27017, user='admin'):
    for _ in range(3):
        try:
            client = pymongo.MongoClient(host, port, username=user, password=getpass(prompt='Password:'))
            client['admin'].list_collection_names()
        except pymongo.errors.OperationFailure as e:
            print(e)
        else:
            print("Successfully conected to the database.")
            return client

def create_indices(db, collection='data'):
    db[collection].create_index([
        ("app", 1),
        ("box_opt.mesh_num", 1),
        ("box_opt.par_num", 1),
        ("box_opt.mesh_num_pwr", 1),
        ("box_opt.box_size", 1)
    ])

def get_id_keys(db, app, collection='data'):
    doc = db[collection].find_one({'app' : app})
    opt_keys = [x for x in doc.keys() if
        x != 'app' and
        x != 'data' and
        x != '_id' and
        x != 'run_opt' and
        x != 'out_opt' and
        x != 'app_opt'
    ]

    all_keys = {}
    for key in opt_keys:
        all_keys[key] = doc[key].keys()

    if app != 'TZA':
        all_keys["cosmo"].remove("smoothing_k")
    
    return all_keys


def get_pipeline(db, app):
    keys = get_id_keys(db, app)
    group = {"_id" : {}}
    for opt_key, opt_key_list in keys.items():
        for key in opt_key_list:
            group["_id"]["%s_%s" % (opt_key, key)] = "$%s.%s" % (opt_key, key)

    group["count"] = {"$sum": 1}
    pipeline = [
        {"$match" : {'app' : app}},
        {"$group" : group },
        {"$sort": SON([
            ("_id.box_opt_mesh_num", 1),
            ("_id.box_opt_mesh_num_pwr", 1),
            ("_id.box_opt_par_num", 1)
            ])}
    ]

    if app == 'CHI':
        pipeline[2]["$sort"]["_id.chi_opt_phi"] = 1
        pipeline[2]["$sort"]["_id.chi_opt_n"] = 1

    return pipeline

def print_db_info(db, collection='data'):
    apps = db.data.distinct('app', {})
    print("There are in total %i number of simulations." % db.data.count_documents({}))
    for app in apps:
        print("%s (%i):" % (app, db[collection].count_documents({'app' : app})))

        pipeline = get_pipeline(db, app)
        for doc in db.data.aggregate(pipeline):
            _id = doc["_id"]
            msg = "\tNm = %i, NM = %i, Np = %i, L = %i" % (
                _id["box_opt_mesh_num"],
                _id["box_opt_mesh_num_pwr"],
                _id["box_opt_par_num"],
                _id["box_opt_box_size"]
                )
            if app == 'CHI':
                msg += ", phi = %.1E, n = %.1f" % (
                    _id["chi_opt_phi"],
                    _id["chi_opt_n"]
                )
                if _id['chi_opt_linear']:
                    msg += ' (linear)'
                else:
                    msg += ' (non-linear)'
            msg += ", num = %i" % doc['count']
            print(msg)

def add_one_sim_data(a_file, db, collection='data', override=False):
    """load simulation parameters with available data and save them into database,
    return true if the record was not in the database before"""

    # open sim_param.json
    with open(a_file) as json_file:
        data = json.loads(json_file.read())

    # check if we already loaded this simulation
    if data["results"] is not None and "database" in data["results"]:
        if not override:
            return False
        elif datetime.datetime.strptime(data["results"]["database"], '%Y-%m-%d %H:%M:%S.%f') > override:
            return False

    # get rid of unwanted data
    data.pop("results")
    data["out_opt"].pop("out_dir")
    data["run_opt"].pop("num_thread")
    
    # load all the files and add data to the database
    data['data'] = {}
    dirs = set(RESULTS_DIRS.values())
    root_dir = os.path.dirname(a_file) + '/'
    # go through all possible directories
    for data_dir in dirs:
        data['data'][data_dir] = {}
        # go through all files in the directory
        for data_file in get_files_in_traverse_dir(root_dir + data_dir, '*'):
            # load data and save them as plain lists
            data['data'][data_dir][data_file] = numpy.loadtxt(data_file).tolist()

    # save eveyrthing into the database
    if db[collection].find_one_and_update(data,  {"$set" : data}, upsert=True) is None:
        new = True
    else:
        new = False

    # update sim_param.json
    with open(a_file, 'r+') as json_file:
        data_orig = json.loads(json_file.read())
        if data_orig["results"] is None:
            data_orig["results"] = {}
        data_orig["results"]["database"] = str(datetime.datetime.utcnow())
        json_file.seek(0)
        json.dump(data_orig, json_file, indent=2)
    return new


def add_many_sim_data(a_dir, db, collection='data', a_file='sim_param.json', verbose=True, override=False):
    """add all simulations in subdirectories of 'a_dir'"""

    if verbose:
        print("Starting to populate the database with simulations...")

    # find all simulations parameters and load data
    new_docs = 0
    counter = 0
    sim_files = get_files_in_traverse_dir(a_dir, a_file)
    length = len(sim_files)
    for sim_file in sim_files:
        counter += 1
        print("Adding simulation %i from %i\r" % (counter, length), end='')
        new_docs += add_one_sim_data(sim_file, db, collection=collection, override=override)

    # print some useful info
    if verbose:
        print("There were added %i new simulations to the database from %i in the directory."
             % (new_docs, length))
        print_db_info(db, collection=collection)

def init_database(a_dir, host='localhost', port=27017, user='admin', init=True, override=False):
    """initialize database, create user with admin rights, add data"""
    if init:
        # start database
        create_database(host=host, port=port, user=user)

        # wait for restart of the database
        raw_input("Press Enter to continue...")

    # connect to database and fill it with data
    db = connect_db(host='localhost', port=27017, user='admin')["fastsim"]
    create_indices(db)
    add_many_sim_data(a_dir, db, override=override)

    return db
