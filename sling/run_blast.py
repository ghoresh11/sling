import utils
import glob, os
import subprocess
import sys
import warnings

class Error (Exception): pass

def run_blast(out_dir,file_type,evalue,cpu):	
	configs = utils.load_config_file()
	
	command = map(str,[configs["makeblastdb"],"-in",out_dir + "/" + file_type + ".fasta","-dbtype","prot"] )
	command2 = map(str,[configs["blastp"], "-db",out_dir + "/" + file_type + ".fasta", "-query", out_dir + "/" + file_type + ".fasta","-out",out_dir + "/" + file_type + "_blast_results","-outfmt","6","-evalue",evalue, "-num_threads", cpu])
	attempt = 0
	res = 1 
	while attempt < utils.MAX_ATTEMPTS and res != 0:
		res = subprocess.call(command)
		if res != 0:
			attempt += 1
			continue
		res = subprocess.call(command2)
		attempt += 1
	if res != 0:
		sys.exit("Error: Failed to run BLASTP for [" + file_type +"]. Please check log files.")


class RunBlast:

	def __init__(self, order, 
		filter_id ,
		group_id, 
		out_dir = ".",
		min_blast_evalue = 0.01,
		sep = ",",
		report_unfit = False,
		cpu = 2):
		
		self.results_dir = os.path.join(out_dir,filter_id + "_FILTER")
		self.out_dir = os.path.join(out_dir,group_id + "_GROUP","blast_files")
		
		
		self.min_blast_evalue = min_blast_evalue
		self.order = order
		self.sep = sep
		self.report_unfit = report_unfit
		self.cpu = 2


	def _hits_to_fasta(self): # combine all the hits and write them into a single fasta file to run blast
		utils.assure_path_exists(self.out_dir)
		hits = open(os.path.join(self.out_dir,"hits.fasta"),"w")
		
		if self.order == "both":
			downstream = open(os.path.join(self.out_dir,"downstream.fasta"), "w")
			upstream = open(os.path.join(self.out_dir ,"upstream.fasta"),"w")
		else:
			partners = open(os.path.join(self.out_dir ,"partners.fasta"),"w")


		for file in os.listdir(self.results_dir):
			if not file.endswith(".csv"):
				continue
			strain = os.path.basename(file)
			strain = strain.replace(".csv","")

			with open(os.path.join(self.results_dir,file)) as f:
				line_num = 1
				for line in f:
					if line.startswith("Domain"):
						continue
					toks = line.strip().split(self.sep)
					identifier = ">" + strain + "|" + str(line_num)

					if self.order == "either" or self.order == "both":
						hits.write(identifier + "*hit" + "\n" + toks[15] + "\n")
						if self.order == "either":
							if toks[16] != "":
								partners.write(identifier +"*upstream"+ "\n" + toks[16] + "\n")
							if toks[17] != "":
								partners.write(identifier + "*downstream" + "\n" + toks[17] + "\n")
						else:
							upstream.write(identifier +"*upstream"+ "\n" + toks[16] + "\n")
							downstream.write(identifier + "*downstream" + "\n" + toks[17] + "\n")
					else:
						hits.write(identifier + "*hit" + "\n" + toks[11] + "\n")
						partners.write(identifier + "*" + self.order + "\n" + toks[12] + "\n")
						
					line_num += 1


		## if in the previous step, the unfit were also reported, add them to the "hits" in the network analysis
		if self.report_unfit and not os.path.exists(os.path.join(self.results_dir,"UNFIT")):
			warnings.warn("Could not find UNFIT files from SUMMARISE step. To report unfit, turn on --report_unfit / -u flag in FILTER and run again.")

		if os.path.exists(os.path.join(self.results_dir,"UNFIT")):
			for file in os.listdir(os.path.join(self.results_dir,"UNFIT")):
				if not file.endswith(".csv"):
					continue
				with open(os.path.join(self.results_dir,"UNFIT",file)) as f:
					line_num = 1
					for line in f:
						if line.startswith("Strain"):
							continue
						toks = line.strip().split(",")
						strain = toks[0]
						identifier = ">" + strain + "|" + str(line_num) + "*unfit"
						hits.write(identifier + "\n" + toks[8] + "\n")
						line_num += 1
		
		hits.close()
		
		if self.order == "both":
			upstream.close()
			downstream.close()
		else:
			partners.close()

		
	def run(self):

		## 1. preprocess the files
		self._hits_to_fasta()

		### run blast in multiple threads
		run_blast(self.out_dir,"hits",self.min_blast_evalue,self.cpu)
		
		if self.order == "both":
			run_blast(self.out_dir,"upstream",self.min_blast_evalue, self.cpu)
			run_blast(self.out_dir,"downstream",self.min_blast_evalue, self.cpu)
		else:
			run_blast(self.out_dir,"partners",self.min_blast_evalue, self.cpu)

