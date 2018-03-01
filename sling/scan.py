import multiprocessing
import os
import utils
import subprocess
import sys
import math
import signal


class Error (Exception): pass

class Scan:

	def __init__(self,
		prep_id,
		scan_id,
		hmm_db, 
		out_dir = ".",
		cpu = 2):

		self.args = {"prep_dir" : os.path.join(os.path.abspath(out_dir),prep_id + "_PREPARE") ,
		"scan_dir" : os.path.join(os.path.abspath(out_dir),scan_id + "_SCAN") ,
		"cpu" : cpu,
		"hmm_db" : hmm_db,
		"configs" : utils.load_config_file() ,
		"prep_id": prep_id,
		"scan_id": scan_id}


	def _run_hmmpress(self):
		self.args["hmm_db"] = os.path.abspath(self.args["hmm_db"])
		res = subprocess.call([self.args["configs"]["hmmpress"],self.args["hmm_db"]])



	def _get_hmmer_version(self,command):
		p = subprocess.Popen([self.args["configs"][command], "-h"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output, err = p.communicate(b"input data that is passed to subprocess' stdin")
		rc = p.returncode
		output = output.split()
		for i in range(0,len(output)):
			if output[i] == "HMMER":
				return (command + ": HMMER-" + output[i+1] + "\n")
		return (command + ": HMMER version not found\n")


	def run(self):
		
		## create output directory
		utils.assure_path_exists(self.args["scan_dir"])

		if self.args["hmm_db"] not in utils.databases:
			self._run_hmmpress()
		else:## before starting, copy the data directory to the out direcotry
			d = os.path.abspath(os.path.dirname(__file__))
			data_env = os.path.join(d, 'data/')
			os.system("mkdir -p " + os.path.join(self.args["scan_dir"],"data"))
			## Do both??? -> this seems to choose how to behave arbitrarily
			os.system("cp -r " + data_env + " " + os.path.join(self.args["scan_dir"]))
			os.system("cp -r " + data_env + " " + os.path.join(self.args["scan_dir"],"data"))
			self.args["hmm_db"] = os.path.join(self.args["scan_dir"],'data',self.args["hmm_db"])

		

		
		## get version of hmmscan for log file
		log_other =  self._get_hmmer_version("hmmscan")
		log_other = log_other + self._get_hmmer_version("hmmpress")
		log_other = log_other + "###   INPUT   ### \ncnt\tgenome\tfasta_file\tgff_file\n" # keeping a text file of all the genomes used

		jobs = []
		cnt = 1
		for file in os.listdir(self.args["prep_dir"]):
			if file.endswith(".sixframe.fasta"):
				basename = os.path.basename(file)
				basename = basename.replace(".sixframe.fasta","")

				sixframe_file = os.path.join(self.args["prep_dir"], basename + ".sixframe.fasta")
				annotated_file = os.path.join(self.args["prep_dir"], basename + ".annotated.fasta")

				scan_genome = {"basename" :basename, "source": "sixframe", "fasta_file": sixframe_file, "out_dir": self.args["scan_dir"], "hmm_db": self.args["hmm_db"], "hmmscan": self.args["configs"]["hmmscan"]}
				jobs.append(scan_genome)


				if os.path.isfile(annotated_file):
					scan_genome = {"basename" :basename, "source": "annotated", "fasta_file": annotated_file, "out_dir": self.args["scan_dir"], "hmm_db": self.args["hmm_db"], "hmmscan": self.args["configs"]["hmmscan"]}
					jobs.append(scan_genome)
					log_other = log_other + str(cnt) +"\t" + basename +"\t"+ sixframe_file +"\t"+ annotated_file+"\n"
				else:
					log_other = log_other + str(cnt) +"\t" + basename +"\t"+ sixframe_file +"\tnot found\n"
				cnt += 1


		utils.write_log(os.path.join(self.args["scan_dir"], "LOG"), "STEP 2 : GENOME SCANNING", self.args, log_other)

		## Allocated CPUs in an efficient way
		CPUs_per_scan = 5
		
		
		## if there are more CPUs than jobs, or we have less than 5 CPUs, don't use multiprocessing, just give 
		## the maximum number of CPUs to each hmmscan call
		if self.args["cpu"] > len(jobs) or self.args["cpu"] < CPUs_per_scan:
			pool_size = 1
			CPUs_per_scan = self.args["cpu"]
		## otherwise, create pool of size CPUS/5, and each scan will run on 5 CPUs (what to do with CPUs I lose?)
		else:
			pool_size =  self.args["cpu"] / CPUs_per_scan 

		for j in jobs:
			j["scan_cpu"] = CPUs_per_scan
		


		pool = multiprocessing.Pool(processes = pool_size)
		results = pool.map_async(run_scan,tuple(jobs))
		pool.close()
		pool.join()


def run_scan(args):
	command = map(str,[args["hmmscan"],"--cpu",args["scan_cpu"],"--max","--noali","--domtblout",
		os.path.join(args["out_dir"] , args["basename"] + "." + args["source"] + ".result"), args["hmm_db"], args["fasta_file"]])

	res = 1 
	attempt = 0
	while attempt < utils.MAX_ATTEMPTS and res != 0:
		res = subprocess.call(command)
	if res != 0:
		sys.exit("Error: Failed to complete hmmscan for [" + self.basename + "_" + self.source +"]. Please check log files.")

	return res

