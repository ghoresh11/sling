import os
import sys
from sling import __version__ 
import datetime

## every time there's a new database, add here
MAX_ATTEMPTS = 2 # try every process a maximum MAX attempts times
databases = ["toxins", "RND_pump"]


def assure_path_exists(path):
    if not os.path.exists(path):
    	os.makedirs(path)

## get the order from the requirement file
def get_order(req,data_env):

	if req not in databases and req != "flexible":
		## read the requirement file given by the user, if any field is missing -> default value is taken
		req_file = os.path.abspath(req)
	else:
		req_file = os.path.join(data_env,req + ".txt")

	with open(req_file) as f:
		for line in f:
			if line == "": # empty line
				continue
			key, val = line.strip().split()
			if key == "order":
				return val
	sys.exit("Must provide 'order' argument in the requirement file!")

def load_config_file():
	config_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),"CONFIG")
	configs = {}
	with open(config_file) as f:
		for line in f:
			line = line.strip().split()
			configs[line[0]] = line[1]
	return configs


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