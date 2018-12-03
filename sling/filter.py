import utils
import os
import pandas
import sys

class Error (Exception): pass

class Summarise:

	def __init__(self, 
		prep_id,
		scan_id,
		filter_id,
		req,
		out_dir=".",
		sep = ",",
		min_hmmscan_score = 20,
		order = None,
		max_diff_avg_length = None,
		min_hit_length = None,
		max_hit_length = None,
		min_upstream_length = None,
		max_upstream_length = None,
		min_downstream_length = None,
		max_downstream_length = None,
		max_distance = None,
		max_overlap = None, 
		domains_file = "",  # the user can provide his own file of domains with the average length expected for each domain
		ignore_file = "",
		report_unfit = False):

		self.args={"prep_id" : prep_id, "scan_id": scan_id, "filter_id": filter_id}

		self.args["req"] = req
		self.args["out_dir"] = os.path.abspath(out_dir)
		self.args["results_dir"] = os.path.join(self.args["out_dir"],filter_id+"_FILTER")
		
		self.req_dict = {}

		d = os.path.abspath(os.path.dirname(__file__))
		data_env = os.path.join(d, 'data/')
		self._get_requirements(data_env)

		self.args["min_hmmscan_score"] = min_hmmscan_score
	
		### override value in requirement file if want to
		if max_diff_avg_length != None:
			self.req_dict["max_diff_avg_length"] = max_diff_avg_length
		if min_hit_length != None:	
			self.req_dict["min_hit_length"] = min_hit_length
		if max_hit_length != None:
			self.req_dict["max_hit_length"] = max_hit_length

		if min_downstream_length != None:			
			self.req_dict["min_downstream_length"] = min_downstream_length
		if max_downstream_length != None:
			self.req_dict["max_downstream_length"] = max_downstream_length
		if min_upstream_length != None:
			self.req_dict["min_upstream_length"] = min_upstream_length
		if max_upstream_length != None:
			self.req_dict["max_upstream_length"] = max_upstream_length

		if max_distance != None:
			self.req_dict["max_distance"] = max_distance
		if max_overlap != None:
			self.req_dict["max_overlap"] = max_overlap
	
		if order != None:
			order = order.lower()
			if order not in ["upstream","downstream","either","both"]:
				sys.exit("Order can be: 'upstream', 'downstream', 'either' or 'both'")
			self.req_dict["order"] = order

		## if user hasn't provided files, take defaults
		if domains_file == "" and req in utils.databases:
			domains_file = os.path.join(data_env,"domains.txt")
		if ignore_file == "" and req in utils.databases:
			ignore_file = os.path.join(data_env,"domains_to_ignore.txt")

		self.args["sep"] = sep
		self.args["profiles_lengths_file"] = domains_file
		self.args["ignore_file"] = ignore_file
		self.args["report_unfit"] = report_unfit
		self.args["req_dict"] = self.req_dict
		self.args["profile_lengths"] = {}
		self.args["profiles_to_ignore"] = {}


	
	def _get_requirements(self,data_env):
		## path to the requirement document of each built in DB
		reqs = utils.databases

		if self.args["req"] in reqs: ## built in collections
			req_file = os.path.join(data_env , self.args["req"]+ ".txt")
		else: # not a built in collection, take the default values
			req_file = os.path.join(data_env , "default.txt")

		with open(req_file) as f: ## get all the values from the req file
			for line in f:
				key, val = line.strip().split()
				if key != "order":
					val = int(val)
				self.req_dict[key] = val

	def run(self):

		utils.assure_path_exists(self.args["results_dir"])
		## get domains' lengths and domains to ignore
		self._parse_domains()
		self._parse_domains_to_ignore()

		
		log_other = "###   INPUT  ### \ngenome\tfasta_hmmer_result\tgff_hmmer_result\n" # keeping a text file of all the genomes used

		jobs = []
		for file in os.listdir(os.path.join(self.args["out_dir"],self.args["scan_id"] + "_SCAN")):
			if file.endswith(".sixframe.result"):
				basename = os.path.basename(file)
				basename = basename.replace(".sixframe.result","")
				
				scan_summary = dict(self.args)
				scan_summary["basename"] = basename
				scan_summary["out_file"] = os.path.join(self.args["results_dir"], basename + ".csv")
				scan_summary["orf_locs"] = None
				scan_summary["orfs"] = {}
				scan_summary["unfit_orfs"] = {}
				if scan_summary["report_unfit"]:
					scan_summary["report_unfit"] = os.path.join(self.args["results_dir"],"UNFIT",basename + ".csv")

				jobs.append(scan_summary)
		
		for j in jobs:
			sixframe_file, annotation_file = run_summarise(j)
			log_other = log_other + j["basename"] + "\t" + sixframe_file + "\t" + annotation_file + "\n"

		utils.write_log(os.path.join(self.args["results_dir"], "LOG"), "STEP 3 : FILTER HMMER HITS", self.args, log_other)

	def _parse_domains(self):
		if self.args["profiles_lengths_file"] == "":
			return
		self.args["profiles_lengths_file"] = os.path.abspath(self.args["profiles_lengths_file"])
		with open(self.args["profiles_lengths_file"], 'rb') as f:
			for line in f:
				line = line.strip().split()
				self.args["profile_lengths"][line[0]] = float(line[1])

	def _parse_domains_to_ignore(self):
		if self.args["ignore_file"] == "":
			return
		self.args["ignore_file"] = os.path.abspath(self.args["ignore_file"])
		self.args["profiles_to_ignore"] = [line.rstrip() for line in open(self.args["ignore_file"])]



def run_summarise(args):
	
	parse_orf_locs(args) # get location and sequence of all available ORFs

	# parse all the hmmer results, keep only rows that match conditions
	sixframe_file = parse_hmmer_results(args, "sixframe") # on sixframe translation
	annotation_file = parse_hmmer_results(args,"annotation")	# on annotation

	## find antitoxins
	find_candidate_antitoxins(args)
	merge_hits(args) 

	write_results(args)
	report_unfit(args)

	return (sixframe_file,annotation_file)




# if an "unfit" ORF is actually found in the original collection, remove it
def merge_unfit(args):

	## merge with other identical ORFs which were found to be fit
	for o1_k in args["unfit_orfs"]:
		for o2_k in args["orfs"]:
			o1 = args["unfit_orfs"][o1_k]
			o2 = args["orfs"][o2_k]
			if o2.keep and o1 == o2:
				o1.unfit = False
	
	## merge duplicates of unfit ORFs
	for o1_k in args["unfit_orfs"]:
		for o2_k in args["unfit_orfs"]:
			o1 = args["unfit_orfs"][o1_k]
			o2 = args["unfit_orfs"][o2_k]
			if o1.name == o2.name:
				continue
			if o1.unfit == False or o2.unfit == False:
				continue
			if o1 == o2:
				o1.unfit=False
				o2.unfit=True

def report_unfit(args):
	if not args["report_unfit"]:
		return
	merge_unfit(args)
	utils.assure_path_exists(os.path.dirname(args["report_unfit"]))
	report = open(args["report_unfit"],"w")
	report.write(args["sep"].join(["Strain","Domain","HMMER_score","Contig","Strand","Hit_length","Hit_start","Hit_stop","Sequence","Reason1","Reason2"])+ "\n")

	for o in args["unfit_orfs"]:
		if args["unfit_orfs"][o].unfit:
			report.write(args["unfit_orfs"][o]._to_unfit(args["basename"],args["sep"]) + "\n")
	report.close()


def parse_orf_locs(args):
	## files to parse the orf_locs df
	
	sixframe_orf_locs_file = os.path.join(args["out_dir"],args["prep_id"]+"_PREPARE",args["basename"] + ".sixframe.bed")
	
	args["orf_locs"] = pandas.read_table(sixframe_orf_locs_file,header=None)

	annotated_orf_locs_file = os.path.join(args["out_dir"],args["prep_id"]+"_PREPARE",args["basename"] + ".annotated.bed")
	
	## if there's an annotation file, add it
	if os.path.isfile(annotated_orf_locs_file):
		try:
			args["orf_locs"] = pandas.concat([args["orf_locs"],pandas.read_table(annotated_orf_locs_file ,header=None)],
                     ignore_index=True)
		except:
			return



def parse_hmmer_results(args, source):

	results_file = os.path.join(args["out_dir"],args["scan_id"]+"_SCAN", args["basename"] + "."+ source + ".result")
	
	if source == "annotation":
		results_file = os.path.join(args["out_dir"], args["scan_id"]+"_SCAN" ,args["basename"] + ".annotated.result")
		if not os.path.isfile(results_file): ## there is no annotation results file
			return "Not found"

	with open(results_file) as f:
		prev_name = "" # save previous name, to not repeat actions
		for line in f:
			if (line.startswith("#")):
				continue

			curr_line = line.split()

			if len(curr_line) < 8: ## run hasn't completed
				continue
			

			name = curr_line[0] # orf name
			domain_name = curr_line[3]
			score = float(curr_line[7]) # score of current hit

			## if using hmmsearch, the name is toks[0] and the domain is in curr_line[3]
			## check this by checking if name starts with source
			if not name.startswith(source):
				name = curr_line[3]
				domain_name = curr_line[0]

			if score < args["min_hmmscan_score"]: # ignore hits with low scores
				continue

			if domain_name in args["profiles_to_ignore"]:
				continue

			if name != prev_name: # prevents looking up the name when it appears in multiple rows
				curr_orf = args["orf_locs"][ args["orf_locs"][3] == name]
				prev_name = name
			
			sequence = curr_orf[6].tolist()[0]
			length=len(sequence)

			min_hit_length = args["req_dict"]["min_hit_length"]
			max_hit_length = args["req_dict"]["max_hit_length"]

			if domain_name in args["profile_lengths"]:
				min_hit_length = max(args["profile_lengths"][domain_name] - args["req_dict"]["max_diff_avg_length"], args["req_dict"]["min_hit_length"])
				max_hit_length = args["profile_lengths"][domain_name] + args["req_dict"]["max_diff_avg_length"]


			if length >= min_hit_length and length <= max_hit_length:
				if name in args["orfs"]:
					if score > args["orfs"][name].score:
						args["orfs"][name].score = score
						args["orfs"][name].domain = domain_name
				else:
					contig, strand, start, stop = curr_orf[0].tolist()[0], curr_orf[5].tolist()[0] ,curr_orf[1].tolist()[0] ,curr_orf[2].tolist()[0]
					args["orfs"][name] = ORF(name, domain_name, score, sequence, contig, strand, start, stop, length,source)
			else:
				contig, strand, start, stop = curr_orf[0].tolist()[0], curr_orf[5].tolist()[0] ,curr_orf[1].tolist()[0] ,curr_orf[2].tolist()[0]
				args["unfit_orfs"][name] = ORF(name, domain_name, score, sequence, contig, strand, start, stop, length,source)
				args["unfit_orfs"][name].reason = "hit length"
				args["unfit_orfs"][name].unfit = True
	return results_file


def choose_delta(deltas):
	deltas = list(deltas)
	chosen = deltas[0]
	for d in deltas:
		if abs(d) < abs(chosen):
			chosen = d
	return chosen 

### cases in which we have AT - T (upstream for +, downstream for -)
def find_upstream_antitoxin(args, orf, orf_locs, flag = False):

	## this is the downstream gene for the reverse strand
	if (orf.strand == "-" and args["req_dict"]["order"] == "upstream") or (orf.strand == "+" and args["req_dict"]["order"] == "downstream"): 
		return False

	## AT stop >= T start - MAX_distance 
	filtered = orf_locs[orf_locs[2] >= orf.start - args["req_dict"]["max_distance"]] 
	
	## AT stop <= T start + max overlap
	filtered = filtered[filtered[2] <= orf.start + args["req_dict"]["max_overlap"]]


	## couldn't find adjacent genes
	if filtered.shape[0] == 0:
		if orf.name not in args["unfit_orfs"]:
			args["unfit_orfs"][orf.name] = orf
		if orf.reason != "":
			orf.reason += ","
		if orf.strand == "+":
			orf.reason = orf.reason + "No adjacent upstream ORF"
			orf.unfit = True
		else:
			orf.reason = orf.reason + "No adjacent downstream ORF"
			orf.unfit = True
		return False

	if orf.strand == "+": ## this is the upstream gene
		## min_length(aa) *3 <= AT stop - AT start <= max_length(aa)*3 
		filtered = filtered[filtered[2] - filtered[1] >=  args["req_dict"]["min_upstream_length"]*3]
		filtered = filtered[filtered[2] - filtered[1] <=  args["req_dict"]["max_upstream_length"]*3]
	else: ## this is the downstream gene
		filtered = filtered[filtered[2] - filtered[1] >=  args["req_dict"]["min_downstream_length"]*3]
		filtered = filtered[filtered[2] - filtered[1] <=  args["req_dict"]["max_downstream_length"]*3]

	## didn't meet length requirements
	if filtered.shape[0] == 0:
		if orf.name not in args["unfit_orfs"]:
			args["unfit_orfs"][orf.name] = orf
		if orf.reason != "":
			orf.reason += ","
		if orf.strand == "-":
			orf.reason = orf.reason + "Downstream length"
			orf.unfit = True
		else:
			orf.reason = orf.reason + "Upstream length"
			orf.unfit = True
		return False

	filtered = filtered[orf.start - filtered[2]  == choose_delta(orf.start - filtered[2])]
	
	if orf.strand == "+":
		orf.upstream_at_start, orf.upstream_at_stop, orf.upstream_at_sequence = filtered[1].tolist()[0], filtered[2].tolist()[0], filtered[6].tolist()[0]
		orf.delta_up = orf.start - orf.upstream_at_stop
	else:
		orf.downstream_at_start, orf.downstream_at_stop, orf.downstream_at_sequence = filtered[1].tolist()[0], filtered[2].tolist()[0], filtered[6].tolist()[0]
		orf.delta_down = orf.start - orf.downstream_at_stop
	return True


### cases in which we have T - AT (downstream for +, upstream for -)
def find_downstream_antitoxin(args, orf, orf_locs, flag = False):

	## this is the upstream gene for the reverse strand
	if (orf.strand == "+" and args["req_dict"]["order"] == "upstream") or (orf.strand == "-" and args["req_dict"]["order"] == "downstream"): 
		return False
	
	filtered = orf_locs[orf_locs[1] <= orf.stop + args["req_dict"]["max_distance"]] # antitoxin start location
	filtered = filtered[filtered[1] >= orf.stop - args["req_dict"]["max_overlap"]] # antitoxin start location

	if filtered.shape[0] == 0:
		if orf.name not in args["unfit_orfs"]:
			args["unfit_orfs"][orf.name] = orf
		if orf.reason != "":
			orf.reason += ","
		if orf.strand == "-":
			orf.reason = orf.reason + "No adjacent upstream ORF"
			orf.unfit = True
		else:
			orf.reason = orf.reason + "No adjacent downstream ORF"
			orf.unfit = True

		return False

	if orf.strand == "+":
		filtered = filtered[filtered[2] - filtered[1] >=  args["req_dict"]["min_downstream_length"]*3]
		filtered = filtered[filtered[2] - filtered[1] <=  args["req_dict"]["max_downstream_length"]*3]
	else:
		filtered = filtered[filtered[2] - filtered[1] >=  args["req_dict"]["min_upstream_length"]*3]
		filtered = filtered[filtered[2] - filtered[1] <=  args["req_dict"]["max_upstream_length"]*3]
	
	if filtered.shape[0] == 0:
		if orf.name not in args["unfit_orfs"]:
			args["unfit_orfs"][orf.name] = orf
		if orf.reason != "":
			orf.reason += ","
		if orf.strand == "+":
			orf.reason = orf.reason + "Downstream length"
			orf.unfit = True
		else:
			orf.reason = orf.reason + "Upstream length"
			orf.unfit = True
		return False

	filtered = filtered[filtered[1] - orf.stop == choose_delta(filtered[1] - orf.stop)]
	

	if orf.strand == "-":
		orf.upstream_at_start, orf.upstream_at_stop, orf.upstream_at_sequence = filtered[1].tolist()[0], filtered[2].tolist()[0], filtered[6].tolist()[0]
		orf.delta_up = orf.upstream_at_start - orf.stop
	else:
		orf.downstream_at_start, orf.downstream_at_stop, orf.downstream_at_sequence = filtered[1].tolist()[0], filtered[2].tolist()[0], filtered[6].tolist()[0]
		orf.delta_down = orf.downstream_at_start - orf.stop
	return True



## find the antitoxin to all the toxins that matched the search
def find_candidate_antitoxins(args):


	## remove ORFs that have a large number of "XXX" as potentail antitoxins
	args["orf_locs"] = args["orf_locs"][~args["orf_locs"][6].str.contains("XXXXXXXX")] 
	args["orf_locs"]= args["orf_locs"][~args["orf_locs"][6].str.contains("NNNNNNNN")] 
	for o in args["orfs"]:

		curr_orf = args["orfs"][o]
		
		orf_locs = args["orf_locs"][args["orf_locs"][0] == curr_orf.contig] # filter by contig
		orf_locs = orf_locs[orf_locs[5] == curr_orf.strand] # filter by strand
		
		up = find_upstream_antitoxin(args,curr_orf, orf_locs)
		down = find_downstream_antitoxin(args,curr_orf, orf_locs)


		## For 3 component operons, must find up and down
		if args["req_dict"]["order"] == "both" and not (up and down):
			curr_orf.keep = False
		
		# couldn't find a candidate antitoxin
		if not up and not down:
			curr_orf.keep = False


## hierarchy to choose how to merge ORFs that are the same
def merge_hits(args):
	for o1_k in args["orfs"]:
		for o2_k in args["orfs"]:
			o1 = args["orfs"][o1_k]
			o2 = args["orfs"][o2_k]
			if o1.name == o2.name:
				continue # ignore an ORF against itself
			
			if o1.keep == False or o2.keep == False: # have already been checked
				continue

			if o1 == o2:
				## if losing an ORF, prefer to keep both
				if o1.downstream_at_sequence == "-" and o2.downstream_at_sequence != "-" and o2.upstream_at_sequence!="-": # don't lose downstream ORF
					o1.keep = False
					o2.keep = True
					if o1.source != o2.source:
						o2.source = "both"
				elif o2.downstream_at_sequence == "-" and o1.downstream_at_sequence != "-" and o1.upstream_at_sequence != "-":
					o2.keep = False
					o1.keep = True
					if o1.source != o2.source:
						o1.source = "both"
				elif o1.upstream_at_sequence == "-" and o2.upstream_at_sequence != "-" and o2.downstream_at_sequence!="-":
					o1.keep = False
					o2.keep = True
					if o1.source != o2.source:
						o2.source = "both"
				elif o2.upstream_at_sequence == "-" and o1.upstream_at_sequence != "-" and o1.downstream_at_sequence != "-":
					o2.keep = False
					o1.keep = True
					if o1.source != o2.source:
						o1.source = "both"
				elif o1.source != o2.source and o2.source == "annotation":
					o2.keep = True
					o1.keep = False
					o2.source = "both"
				elif o1.source != o2.source:
					o1.keep = True
					o2.keep = False
					o1.source = "both"
				elif o1.score >= o2.score:
					o1.keep = True
					o2.keep = False
				else:
					o2.keep = True
					o1.keep = False


def write_results(args):
	# create output directory if doesn't exist
	utils.assure_path_exists(os.path.dirname(args["out_file"]))

	## write files
	out = open(args["out_file"], "w")
	if args["req_dict"]["order"] == "upstream":
		out.write(args["sep"].join(["Domain","HMMER_score","Contig","Strand","Hit_length","Hit_start","Hit_stop","Upstream_length",
			"Upstream_delta","Upstream_start","Upstream_stop",
			"Hit","Upstream","Source"]) + "\n")
	elif args["req_dict"]["order"] == "downstream":
		out.write(args["sep"].join(["Domain","HMMER_score","Contig","Strand","Hit_length","Hit_start","Hit_stop",
			"Downstream_length","Downstream_delta","Downstream_start","Downstream_stop",
			"Hit","Downstream","Source"]) + "\n")
	else:
		out.write(args["sep"].join(["Domain","HMMER_score","Contig","Strand","Hit_length","Hit_start","Hit_stop","Upstream_length",
			"Upstream_delta","Upstream_start","Upstream_stop","Downstream_length","Downstream_delta","Downstream_start","Downstream_stop",
			"Hit","Upstream","Downstream","Source"]) + "\n")
	for o in args["orfs"]:
		if args["orfs"][o].keep:
			out.write( args["orfs"][o]._to_string(args["sep"],args["req_dict"]["order"]) + "\n")
	out.close()




def ORF(name, domain, score, sequence, contig, strand, start, stop, length, source):
	orf = {}
	orf["name"] = name
	orf["domain"] = domain
	orf["score"]  = score
	orf["sequence"]  = sequence
	orf["contig"]  = contig
	orf["strand"]  = strand 
	orf["start"]  = start
	orf["stop"]  = stop
	orf["length"]  = length
	orf["source"] = source
	orf["reason"]  = ""
	orf["keep"]  = True # boolean to use when merging the hits
	orf["unfit"] = False # boolean to save unfit orfs
	
	### upstream antitoxin
	orf["upstream_at_start"] =""
	orf["upstream_at_stop"] =""
	orf["upstream_at_sequence"] =""
	orf["delta_up"]  = ""

	## downstream antitoxin
	orf["downstream_at_start"] =""
	orf["downstream_at_stop"] =""
	orf["downstream_at_sequence"] =""
	orf["delta_down"]  = ""
	

		
		
class ORF: ## class to keep results of one putative toxin ORF

	def __init__(self, name, domain, score, sequence, contig, strand, start, stop, length,source):
		self.name = name
		self.domain = domain
		self.score = score
		self.sequence = sequence
		self.contig = contig
		self.strand = strand 
		self.start = start
		self.stop = stop
		self.length = length
		self.source = source
		self.reason = ""
		self.keep = True # boolean to use when merging the hits
		self.unfit = False # boolean to save unfit orfs
		
		### upstream antitoxin
		self.upstream_at_start=""
		self.upstream_at_stop=""
		self.upstream_at_sequence=""
		self.delta_up = ""

		## downstream antitoxin
		self.downstream_at_start = ""
		self.downstream_at_stop = ""
		self.downstream_at_sequence = ""
		self.delta_down = ""


	def _to_string(self,sep,order):
		if order == "upstream":
			return (sep.join(map(str,[self.domain,self.score,self.contig,self.strand,self.length,self.start, self.stop,
				len(self.upstream_at_sequence), self.delta_up, self.upstream_at_start, self.upstream_at_stop,
				self.sequence, self.upstream_at_sequence, self.source])))
		if order == "downstream":
			return (sep.join(map(str,[self.domain,self.score,self.contig,self.strand,self.length,self.start, self.stop,
				len(self.downstream_at_sequence), self.delta_down, self.downstream_at_start, self.downstream_at_stop,
				self.sequence, self.downstream_at_sequence, self.source])))
		return (sep.join(map(str,[self.domain,self.score,self.contig,self.strand,self.length,self.start, self.stop,
			len(self.upstream_at_sequence), self.delta_up, self.upstream_at_start, self.upstream_at_stop, len(self.downstream_at_sequence), self.delta_down, self.downstream_at_start, self.downstream_at_stop,
			self.sequence, self.upstream_at_sequence, self.downstream_at_sequence, self.source])))

	def _to_unfit(self,basename,sep):
		return(sep.join(map(str,[basename, self.domain,self.score,self.contig,self.strand,self.length,self.start,self.stop,self.sequence,self.reason])))

	def __eq__(self,other): ## consider the ORFs to be the same if their start and stop ~100bp apart.
		if (other==None):
			return False
		if (self.name == other.name):
			return True
		if (self.strand == other.strand and self.contig==other.contig and self.start-100 <= other.start and other.start <= self.start+100 ):
			return True
		if (self.strand == other.strand and self.contig==other.contig and self.stop-100 <= other.stop and other.stop <= self.stop + 100 ):
			return True
		return False

	def __ne__(self,other):
		return not self.__eq__(other)
		






