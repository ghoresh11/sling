import utils
import run_blast
import group_operons
import sys
import os
import subprocess

class Group:

	def __init__(self, 
		filter_id,
		group_id,
		req,
		out_dir=".",
		order = None,
		min_identity = 30,
		min_blast_evalue = 0.01,
		save_to_ITOL = False,
		sep = ",",
		report_unfit = False,
		cpu = 2):

		self.req = req
		self.order = order
		self.filter_id = filter_id
		self.group_id = group_id
		self.out_dir = os.path.abspath(out_dir)
		self.min_identity = min_identity
		self.min_blast_evalue = min_blast_evalue
		self.save_to_ITOL = save_to_ITOL
		self.sep = sep
		self.report_unfit = report_unfit
		self.cpu = cpu

		if self.order == None: ## user didn't override the order argument
			d = os.path.abspath(os.path.dirname(__file__))
			data_env = os.path.join(d, 'data/')
			self.order = utils.get_order(self.req,data_env)

		## get blast versions
		self.blastp = self._get_blast_version("blastp")
		self.makeblastdb = self._get_blast_version("makeblastdb")

		utils.assure_path_exists(os.path.join(self.out_dir,group_id+ "_GROUP"))
		utils.write_log(os.path.join(self.out_dir,group_id + "_GROUP", "LOG"), "STEP 4 : GROUP ORFs", vars(self), "")

	
	## quick function to get the versions of BLAST used for log files
	def _get_blast_version(self,command):
		configs = utils.load_config_file()
		p = subprocess.Popen([configs[command], "-h"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output, err = p.communicate(b"input data that is passed to subprocess' stdin")
		rc = p.returncode
		output = output.split()
		for i in range(0,len(output)):
			if "version" == output[i] and command=="makeblastdb":
				return ("BLAST-" + output[i+1])
			elif "version" in output[i] and command == "blastp":	
				return ("BLAST-" + output[i+4])
		return ("BLAST version not found\n")

	def run(self):
		


		blast = run_blast.RunBlast(self.order,self.filter_id, self.group_id, out_dir = self.out_dir, 
			min_blast_evalue = self.min_blast_evalue, sep = self.sep, report_unfit = self.report_unfit, cpu = self.cpu)
		blast.run()
		group = group_operons.GroupHits(self.order,self.filter_id, self.group_id,
			out_dir = self.out_dir,
			min_identity = self.min_identity,
			save_to_ITOL=self.save_to_ITOL,sep = self.sep,report_unfit = self.report_unfit)
		group.run()
