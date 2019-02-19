import os
from sling import __version__
import datetime
import multiprocessing
import signal
import sys


# ### get the built-in databases
def get_dbs():
	d = os.path.abspath(os.path.dirname(__file__))
	data_env = os.path.join(d, 'data/')
	return open(os.path.join(data_env,"DATABASES")).read().split()

## every time there's a new database, add here
MAX_ATTEMPTS = 2 # try every process a maximum MAX attempts times
databases = get_dbs()

def assure_path_exists(path):
    if not os.path.exists(path):
    	os.makedirs(path)

## get the order from the requirement file
def get_order(req,data_env):

	if req not in databases:
		## read the requirement file given by the user, if any field is missing -> default value is taken
		req_file = os.path.join(data_env,"default.txt")
	else:
		req_file = os.path.join(data_env,req + ".txt")

	with open(req_file) as f:
		for line in f:
			if line == "": # empty line
				continue
			key, val = line.strip().split()
			if key == "order":
				return val

def write_log(log_file_path, title, params, other):
	log_file = open(log_file_path,"w")
	log_file.write("########## SLING  ########## \n")
	log_file.write("##########"  + title +  "########## \n")
	log_file.write(" SLING version: " + str(__version__))
	log_file.write("\nTime : {:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now()))
	log_file.write("\n### PARAMS ###\n")
	for p in params:
		log_file.write(p + " : " + str(params[p]) + "\n")
	log_file.write(other)
	log_file.close()


def run_pool(jobs, args, fun):
    ''' run the pool of prep workers'''
    # create input handler for using CTRL+C
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = multiprocessing.Pool(processes=args["cpu"])
    signal.signal(signal.SIGINT, original_sigint_handler)

    try:
        results = pool.map_async(fun, tuple(jobs))
        results.get(120000)
    except KeyboardInterrupt as e:
        pool.terminate()
        sys.exit("Terminated by user")
    except ValueError as e:
        pool.terminate()
        sys.exit(
            "Error: please check LOG files")
    else:
        pool.close()
    pool.join()
    return
