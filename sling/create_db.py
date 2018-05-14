import multiprocessing
import os
import utils
import subprocess
import sys
from shutil import copyfile

class Error (Exception): pass

class CreateDB:

	def __init__(self, 
		name,
		hmm_file,
		sling_dir,
		order = None,
		max_diff_avg_length = None,
		min_hit_length = None,
		max_hit_length = None,
		min_upstream_length = None,
		max_upstream_length = None,
		min_downstream_length = None,
		max_downstream_length = None,
		max_distance = None,
		max_overlap = None):

		self.name = name.lower()
		self.hmm_file = os.path.abspath(hmm_file)
		self.sling_dir = sling_dir
		if sling_dir != None:
			self.sling_dir = os.path.join(os.path.abspath(sling_dir),"sling","data")
		else:
			print("####  \n Warning: path to git repository not provided. Adding collection only to local installation of SLING.\n ####")
		
		self.order = order
		self.max_diff_avg_length = max_diff_avg_length
		self.min_hit_length = min_hit_length
		self.max_hit_length = max_hit_length
		self.min_upstream_length = min_upstream_length
		self.max_upstream_length = max_upstream_length
		self.min_downstream_length = min_downstream_length
		self.max_downstream_length = max_downstream_length
		self.max_distance = max_distance
		self.max_overlap = max_overlap
		self.configs = utils.load_config_file() 


	def _add_db(self, data_env):
		with open(os.path.join(data_env,"DATABASES"),"a") as out:
			out.write("\n" + self.name)


	def _constuct_hmm_files(self,data_env):
		if not os.path.isfile(self.hmm_file):
			sys.exit("Error: could not find HMM file: [" + self.hmm_file + "]") 

		copyfile(self.hmm_file, os.path.join(data_env, self.name))
		
		self.hmm_file = os.path.join(data_env, self.name)
		subprocess.call([self.configs["hmmpress"],self.hmm_file])

	def _constuct_req_file(self,data_env):
		req_dict = {}
		## 1. Load default params:
		with open(os.path.join(data_env,"default.txt")) as f: ## get all the values from the req file
			for line in f:
				key, val = line.strip().split()
				if key != "order":
					val = int(val)
				req_dict[key] = val
		## replace default if given by user, and assert the value is ok.
		if self.order != None:
			self.order = self.order.lower()
			orders = ["upstream", "downstream", "either", "both"]
			if self.order not in orders:
				sys.exit("Error: order must be " + str(orders) + ".")
			req_dict["order"] = self.order
		
		if self.max_diff_avg_length != None:
			req_dict["max_diff_avg_length"] = self.max_diff_avg_length

		if self.min_hit_length != None:
			req_dict["min_hit_length"] = self.min_hit_length
		
		if self.max_hit_length != None:
			req_dict["max_hit_length"] = self.max_hit_length
		
		if self.min_upstream_length != None:
			req_dict["min_upstream_length"] = self.min_upstream_length		
		
		if self.max_upstream_length != None:
			req_dict["max_upstream_length"] = self.max_upstream_length

		if self.min_downstream_length != None:
			req_dict["min_downstream_length"] = self.min_downstream_length		
		
		if self.max_downstream_length != None:
			req_dict["max_downstream_length"] = self.max_downstream_length

		if self.max_distance != None:
			req_dict["max_distance"] = self.max_distance		
		
		if self.max_overlap != None:
			req_dict["max_overlap"] = self.max_overlap

		with open(os.path.join(data_env, self.name + ".txt"), "w") as out:
			for key in req_dict:
				out.write(key + "\t" + str(req_dict[key]) + "\n")


	def run(self):
		reqs = utils.databases

		if self.name in reqs:
			sys.exit("Error: name given [" + self.name + "] already exists. Please choose different name")

		d = os.path.abspath(os.path.dirname(__file__))
		data_env = os.path.join(d, 'data/')

		## create the requirements file
		print("Summarising the structural requirements...")
		self._constuct_req_file(data_env)
		if self.sling_dir != None:
			self._constuct_req_file(self.sling_dir)

		## create the hmm files and put them in the data environment
		print("Constructing the HMM profile collection...")
		self._constuct_hmm_files(data_env)
		if self.sling_dir != None:
			self._constuct_hmm_files(self.sling_dir)
		## once everything has succeeded, add to the DATABASES file
		print("Finalising")
		self._add_db(data_env)
		if self.sling_dir != None:
			self._add_db(self.sling_dir)
		
		print("Successfully complete!")

			
			
		


