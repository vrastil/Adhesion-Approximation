"""
'database.py' module serves for creating, updating and managinh databse with
all data from simulations
"""

# system packages
from __future__ import print_function
import os
import pymongo
from getpass import getpass
import json
try:
    from urllib.parse import urlparse
except ImportError:
     from urlparse import urlparse

import numpy

# simpy packages
from .struct import RESULTS_DIRS
from .utils import get_files_in_traverse_dir

# for better search in database use integers instead of strings
APPS_IDX = ['ZA', 'TZA', 'FF', 'FPA', 'CHI']
ZA = APPS_IDX.index('ZA')
TZA = APPS_IDX.index('TZA')
FF = APPS_IDX.index('FF')
FPA = APPS_IDX.index('FPA')
CHI = APPS_IDX.index('CHI')

def create_database(host='localhost', port=27017, user='admin'):
    """create database and user admin in it,
    databse should be started without authentication this first time"""

    client = pymongo.MongoClient(host, port)
    db = client['admin']
    try:
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
            return client

def add_sim_info(a_file, db, collection='data', upsert=True):
    """load simulation parameters with available data and save them into database"""

    # open sim_param.json
    with open(a_file) as json_file:
        data = json.loads(json_file.read())

    # change 'app' field from string to integer
    data['app'] = APPS_IDX.index(data['app'])
    
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

    # save eveyrthing into the database -- check for already inserted data
    db[collection].find_one_and_update(data,  {"$set" : data}, upsert=upsert)
