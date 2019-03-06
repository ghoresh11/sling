import os
import utils
import subprocess
import sys
from shutil import copyfile

class Error (Exception): pass

def constuct_req_file(args):
    ''' read the input requirements from user and create a new txt
    file in the data_env that includes the requirements given by the user.
    if a value wasn't provided by the user, use the defaults'''
    print("Summarising the structural requirements...")
    req_dict = {}
    ## 1. Load default params:
    with open(os.path.join(args.sling_dir,"default.txt")) as f: ## get all the values from the req file
        for line in f:
            key, val = line.strip().split()
            if key != "order":
                val = int(val)
            req_dict[key] = val

    orders = ["upstream", "downstream", "either", "both"]

    args_dict = vars(args)
    for key in req_dict:
        if args_dict[key] is not None:
            if key == "order":
                args_dict[key] = args_dict[key].lower()
            req_dict[key] = args_dict[key]

    utils.check_reqs(req_dict)

    with open(os.path.join(args.sling_dir, args.name + ".txt"), "w") as out:
        for key in req_dict:
            out.write(key + "\t" + str(req_dict[key]) + "\n")
    return

def constuct_hmm_dbs(args):
    ''' copy the HMM file to the SLING directory and run hmmpress'''
    print("Constructing the HMM profile collection...")
    if not os.path.isfile(args.hmm_db):
        sys.exit("Error: could not find HMM file: [" + args.hmm_db + "]")
    copyfile(args.hmm_db, os.path.join(args.sling_dir, args.name))
    args.hmm_db = os.path.join(args.sling_dir, args.name)
    subprocess.call([args.hmmpress, args.hmm_db])
    return

def add_db(args):
    ''' add the new database to the list of available DBs'''
    print("Finalising")
    with open(os.path.join(args.sling_dir,"DATABASES"),"a") as out:
        out.write("\n" + args.name)
    print("Successfully complete!")
    return

def run(args):
    args.name = args.name.lower()
    args.hmm_db = os.path.abspath(args.hmm_db)

    ## get SLING directory
    if args.sling_dir is not None:
        args.sling_dir = os.path.join(os.path.abspath(sling_dir),"sling","data")
    else:
        print("#### Warning: path to git repository not provided. Adding collection to local installation of SLING. ####")
        d = os.path.abspath(os.path.dirname(__file__))
        args.sling_dir = os.path.join(d, 'data/')

    ## check the name is OK
    reqs = utils.databases
    if args.name in reqs:
        sys.exit("Error: name given [" + args.name + "] already exists. Please choose different name")

    constuct_req_file(args)
    constuct_hmm_dbs(args)
    add_db(args)
    return
