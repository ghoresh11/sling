from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
from Bio.Seq import translate
import numpy as np
import sys
import os
from utils import *
import re
import multiprocessing
import warnings
import sys
import signal


class Error (Exception): pass

class Prepare:

	def __init__(self,
		fasta_dir,
		ID,
		out_dir=".",
		fasta_suffix = ".fasta",
		gff_suffix = ".gff", 
		gff_dir = "", 
		min_orf_length = 20,
		start_codons = "atg,gtg,ttg",
		codon_table = "Standard",
		cpu = 1,
	):


		if gff_dir != "":
			gff_dir = os.path.abspath(gff_dir) 
		else:
			gff_dir = os.path.abspath(fasta_dir)

		start_codons = start_codons.split(",")
		if len(start_codons) < 1:
			sys.exit("Must provide at least one start codon! Start codons need to be comma delemited.")
		
		for i in range(0,len(start_codons)):
			s = start_codons[i] 
			start_codons[i] = s.lower()
			if len(s) != 3 or not re.match("^[acgt]*$",s):
				sys.exit("Start codons must be of length 3 and contain only the letters a,c,g or t: " + s)
		
		self.args =  {"fasta_dir" : os.path.abspath(fasta_dir),
		"out_dir" :  os.path.join(os.path.abspath(out_dir),ID + "_PREPARE"),
		"gff_dir" : gff_dir,
		"fasta_suffix" : fasta_suffix,
		"gff_suffix" : gff_suffix,
		"min_orf_length" : min_orf_length,
		"start_codons" : start_codons,
		"codon_table" : codon_table,
		"cpu" : cpu,
		"prep_id": ID}



	def run(self):
		
		original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
		pool = multiprocessing.Pool(processes = self.args["cpu"])
		signal.signal(signal.SIGINT, original_sigint_handler)
		
		assure_path_exists(self.args["out_dir"])

		gff_files = {}
		gff_filenames = os.listdir(self.args["gff_dir"])
		for f in gff_filenames:
			if f.endswith(self.args["gff_suffix"]):
				basename = os.path.basename(f)
				basename = basename.replace(self.args["gff_suffix"],"")
				gff_files[basename] = f


		log_other = "###   INPUT ### \ncnt\tgenome\tfasta_file\tgff_file\n" # keeping a text file of all the genomes used

		jobs = []
		file_cnt = 1
		for file in os.listdir(self.args["fasta_dir"]):
			if file.endswith(self.args["fasta_suffix"]) and not file.startswith("."): # ignore hidden files


				basename = os.path.basename(file)
				basename = basename.replace(self.args["fasta_suffix"],"")

				

				gff_file = ""
				if basename in gff_files:
					gff_file = os.path.join(self.args["gff_dir"],gff_files[basename])

				genome_prep = dict(self.args)
				genome_prep["strain"] = basename
				genome_prep["fasta_file"] = os.path.join(self.args["fasta_dir"],file)
				genome_prep["gff_file"] = gff_file
				
				if gff_file == "":
					gff_file = "Not found"
				log_other = log_other + str(file_cnt) + "\t" + basename + "\t" + genome_prep["fasta_file"] + "\t" + gff_file + "\n"
				file_cnt += 1

				jobs.append(genome_prep)
		
		write_log(os.path.join(self.args["out_dir"], "LOG"), "STEP 1 : GENOME PREP", self.args, log_other)
		
		try:
			results = pool.map_async(run_prepare,tuple(jobs))
			results.get(12000) # give up to 20 mins for each process, if the function takes longer than 10 mins abort	
		except KeyboardInterrupt as e:
			pool.terminate()
		else:
			pool.close()
		pool.join()
		

def run_prepare(args):
	## does it work to initate the class here?
	args["orf_id"] = 0 
	args["contigs"] = {}
	if len(set('[~!@$%^&*()+{}":;\']+$').intersection(args["strain"])) > 0:
			sys.exit("Avoid using special characters in file name! \nFile: " + args["fasta_file"] + "\nRemove chars: " + "\t".join(list(set('[~!@$%^&*()+{}":;\']+$').intersection(args["strain"]))))

	sixframe_translation(args) 
	annotated_orf_locs(args)


def sixframe_translation(args):
	outfile_sixframe = os.path.join(args["out_dir"],args["strain"] + ".sixframe.fasta") # input for hmmscan
	outfile_orf_locs = os.path.join(args["out_dir"],args["strain"] + ".sixframe.bed") # locations of ORFs for later	

	### open output files for writing
	sixframe_out = open(outfile_sixframe,"w")
	orf_locs_out = open(outfile_orf_locs,"w")

	with open(args["fasta_file"]) as handle:
		for values in SimpleFastaParser(handle):

			contig = values[0].split()[0]
			genome = values[1].lower()
			rev_genome = reverse_complement(genome) 

			args["contigs"][contig] = genome # save the contigs for the annotated ORF locs

			## doesn't assume the contigs are circular
			for i in range(0,3):
				# standard
				get_frame_orfs(args, translate(genome[i:], table=args["codon_table"]),genome[i:], i, contig, sixframe_out, orf_locs_out)
				# complement
				get_frame_orfs(args, translate(rev_genome[i:], table = args["codon_table"]),rev_genome[i:], i+3, contig, sixframe_out, orf_locs_out)
	
	sixframe_out.close()
	orf_locs_out.close()

''' get all open reading frames in current frame '''
def get_frame_orfs(args,protein_frame,nuc_frame,i,contig,sixframe_out,orf_locs_out):
	frame_length = len(protein_frame) # get length of the protein sequence

	ORFs = protein_frame.split("*") # split at stop codons
	
	stops = map(len,ORFs) # get the stop codon positions	
	stops = [x+1 for x in stops]
	stops= list(np.cumsum(stops))
	starts = [0] + stops[:-1] # get the start codons positions

	for index in range(0,len(ORFs)): # iterate over all possible ORFs
		nuc_seq = nuc_frame[starts[index] * 3 : stops[index] * 3] # get nuc sequence of that ORF
		check_orf_conditions(args,ORFs[index],nuc_seq,i,contig,starts[index],stops[index],sixframe_out,orf_locs_out,frame_length)
	return 1

''' given an ORF, find the relevant start codon and write to file if meets requirements '''
def check_orf_conditions(args,protein_seq,nuc_seq,strand,contig,start,stop,sixframe_out,orf_locs_out,frame_length):	

	for i in range(0,len(args["start_codons"])):
		res = find_start_codon(args,args["start_codons"][i],protein_seq,nuc_seq,start,stop,strand,frame_length)
		if res != None: # found an ORF with the higher priority codon
			break

	if res == None:
		return 

	orf = res[0] # the actual ORF
	strand_symbol = res[1] # strand
	start = res[2] # start
	stop = res[3] # stop
	## give the ORF a name that can later be useful
	ORF_name = "sixframe|Strain:" + args["strain"] + "|ORF:" +str(args["orf_id"]) + "|Contig:"+ str(contig) + "|Strand:"+strand_symbol+"|Start:" + str(start) + "|Stop:" + str(stop)
	
	## write the results to the fasta file
	sixframe_out.write(">" + ORF_name +"\n" + orf + "\n")

	## write the results to the ORF locations file
	orf_locs_out.write("\t".join([contig,str(start) ,str(stop), ORF_name,"0",strand_symbol,orf]) + "\n")
	
	# increase ORF id by one for next ORF
	args["orf_id"] += 1

''' given an ORF, look for a start codon in the candidate ORF '''
def find_start_codon(args,codon,protein_seq,nuc_seq,start,stop,strand,frame_length):

	aa = translate(codon, args["codon_table"])
	aa_index = protein_seq.find(aa) # find the relevant aa of the start codon
	
	if aa_index>0 and len(protein_seq[aa_index:])>=args["min_orf_length"] and nuc_seq[aa_index*3:aa_index*3 + 3] == codon: # check that all conditions apply
		## use M
		protein_seq = protein_seq[aa_index:]
		start = start + aa_index 
		strand_symbol = "+"
		nuc_seq = nuc_seq[aa_index*3:]
		
		## handle minus strand	
		if strand>2:
			tmp = start
			start = frame_length - stop 
			stop = frame_length - tmp 
			strand_symbol = "-"
		
		## fix coordinates in nucleotide sequence:
		if strand_symbol == "+":
			start = start*3 + 1 + strand%3
			stop = stop*3 - 2 + strand%3
		else:
			### "bandade" to get the cooridantes correct
			add_stop = 0
			add_start = 3
			if strand == 3:
				add_stop = 1
				add_start = 4
			elif strand == 5:
				add_stop = 2
				add_start = 5
			start = start * 3 + add_start
			stop = stop * 3 + add_stop

		protein_seq = "M" + protein_seq[1:] # always put an M at start of protein sequence
		return [protein_seq,strand_symbol,start,stop] # return all values
	
	# couldn't find ORF
	return None


''' get all the open reading frames from the GFF file, save as BED and as FASTA '''
def annotated_orf_locs(args):

	if args["gff_file"]=="":
		return

	args["orf_id"] = 0

	flag = False

	annotated_fasta = os.path.join(args["out_dir"], args["strain"]+".annotated.fasta")
	annotated_orf_locs = os.path.join(args["out_dir"],args["strain"]+".annotated.bed") # locations of ORFs for later

	## create dummy files even if not using annotation
	fasta_out = open(annotated_fasta, "w")
	orf_locs_out = open(annotated_orf_locs,"w")


	with open(args["gff_file"]) as f:
		for line in f:
			if line.startswith("##FASTA"): # end of file
				break

			if line.startswith("#"): # comments
				continue

			line = line.strip().split("\t")

			if len(line)<7: # not full line, ignore
				continue

			if (line[2]=="CDS"): # coding sequence
				strand = line[6]
				start = int(line[3])
				stop= int(line[4])
				contig = line[0]

				if contig not in args["contigs"]:
					sys.exit("Name of contigs in FASTA file must match contigs in GFF file")

				orf_name = "annotation|Strain:" + args["strain"] + "|ORF:" +str(args["orf_id"]) + "|Contig:"+ str(contig) + "|Strand:"+strand+"|Start:" + str(start) + "|Stop:" + str(stop)

				if strand == "-":
					sequence = args["contigs"][contig][start:stop]
					sequence = reverse_complement(sequence)
				else:
					sequence = args["contigs"][contig][start-1:stop-1]

				sequence = translate(sequence, table = args["codon_table"])
				
				## write the results to the fasta file
				fasta_out.write(">" + orf_name +"\n" + sequence + "\n")

				## write the results to the ORF locations file
				orf_locs_out.write("\t".join([contig,str(start) ,str(stop), orf_name,"0",strand,sequence]) + "\n")

				args["orf_id"] += 1
	
	fasta_out.close()
	orf_locs_out.close()