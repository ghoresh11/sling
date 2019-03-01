import multiprocessing
import os
import utils
import subprocess
import sys
import warnings
from shutil import copyfile
import copy

class Error (Exception):
    pass


def run_hmmpress(args):
    ''' run HMM press on new HMM db provided '''
    args.hmm_db = os.path.abspath(args.hmm_db)
    if not os.path.isfile(args.hmm_db):
        sys.exit(
            "Error: could not find HMM file: [" + args.hmm_db + "]")
    res = subprocess.call([args.hmmpress, args.hmm_db])
    return

def copy_data(args, scan_dir):
    ''' if database exists in SLING, copy its content to the
    current working directory (prevents runtime errors)'''
    d = os.path.abspath(os.path.dirname(__file__))
    data_env = os.path.join(d, 'data/')
    curr_env = os.path.join(scan_dir, "data")
    utils.assure_path_exists(curr_env)

    for f in os.listdir(data_env):
        if args.hmm_db in f or "domains" in f:
            copyfile(os.path.join(data_env, f), os.path.join(curr_env,f))

    args.hmm_db =  os.path.join(scan_dir, 'data', args.hmm_db)
    return

def get_hmmer_version(command):
    ''' get the version of the HMMER command '''
    p = subprocess.Popen([command, "-h"], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate(
        b"input data that is passed to subprocess' stdin")
    rc = p.returncode
    output = output.split()
    for i in range(0, len(output)):
        if output[i] == "HMMER":
            return (command + ": HMMER-" + output[i + 1] + "\n")
    return (command + ": HMMER version not found\n")

def create_jobs_list(args, prep_dir, scan_dir):
    ''' create a list of scanning jobs for HMMsearch
    while generating a LOG file in the SCAN directory'''
    # get version of hmmscan for log file
    log_other = get_hmmer_version(args.hmmsearch)
    log_other = log_other + get_hmmer_version(args.hmmpress)
    # keeping a text file of all the genomes used
    log_other = log_other + "###   INPUT   ### \ncnt\tgenome\tprep_fasta_file\n"

    jobs = []
    for f in os.listdir(prep_dir):
        if f.endswith(".fasta"):
            basename = os.path.basename(f)
            basename = basename.replace(".fasta", "")
            prep_file = os.path.join(prep_dir, f)

            scan_genome = {"basename": basename, "fasta_file": prep_file,
                           "out_dir": scan_dir, "hmm_db": args.hmm_db, "hmmsearch": args.hmmsearch}
            jobs.append(scan_genome)

            log_other = log_other + \
                    str(len(jobs)) + "\t" + basename + "\t" + \
                    prep_file + "\n"

    utils.write_log(os.path.join(scan_dir, "LOG"), "STEP 2 : GENOME SCANNING", vars(args), log_other)
    return jobs

def run_hmmer(args):
    ''' Run one hmmsearch job on a FASTA file '''
    command = map(str, [args["hmmsearch"], "--cpu", "1", "--max", "--noali", "--domtblout",
                        os.path.join(args["out_dir"], args["basename"] + ".result"), args["hmm_db"], args["fasta_file"]])
    res = 1
    attempt = 0
    while attempt < utils.MAX_ATTEMPTS and res != 0:
        res = subprocess.call(command)
        attempt += 1
    if res != 0:
        warnings.warn("Warning: Failed to complete hmmsearch for [" + args["fasta_file"] + "].\
         Please check log files! Continued to next file.")
        return res
    return res


def run(args):
    args = copy.deepcopy(args)
    # create output directory
    out_dir = os.path.abspath(args.out_dir)

    if "prep_id" not in vars(args):
        prep_id = args.id
    elif args.prep_id is None:
        prep_id = args.id
    else:
        prep_id = args.prep_id

    prep_dir = os.path.join(out_dir, prep_id + "_PREPARE")

    scan_dir = os.path.join(out_dir, args.id + "_SCAN")

    utils.assure_path_exists(scan_dir) ## create the output directory

    if args.hmm_db not in utils.databases: ## check if the database exists
        run_hmmpress(args)
    else:  # before starting, copy the data directory to the out direcotry
        copy_data(args, scan_dir)

    ## create list of jobs for HMMSEARCH
    jobs = create_jobs_list(args, prep_dir, scan_dir)
    ## run the pool
    pool = multiprocessing.Pool(args.cpu)
    results = pool.map_async(run_hmmer, tuple(jobs))
    pool.close()
    pool.join()
    return
