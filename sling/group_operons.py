import utils
import glob, os
import networkx as nx
import csv
import numpy as np
import pandas
import math

class Error (Exception): pass

class GroupHits():

	def __init__(self, order,
		filter_id,
		group_id, 
		out_dir = ".",
		min_identity = 30,
		save_to_ITOL=False, 
		sep = ",",
		report_unfit = True):
		
		self.group_dir = os.path.join(out_dir,group_id + "_GROUP")
		self.blast_dir = os.path.join(self.group_dir,"blast_files")
		self.out_dir = os.path.abspath(out_dir)
		self.min_identity = min_identity
		self.save_to_ITOL = save_to_ITOL
		self.sep = sep

		self.filter_dir = os.path.join(out_dir,filter_id+ "_FILTER")

		self.hits = {}
		self.unfits = {}
		self.strains = []

		self.order = order
		self.report_unfit = report_unfit


	def _parse_hits_files(self,dir_name):
		## parse all the hits files
		results_dir = self.filter_dir
		if dir_name == "UNFIT":
			results_dir = os.path.join(results_dir,"UNFIT")
		for file in os.listdir(results_dir):
			if not file.endswith(".csv"):
				continue
			
			strain = os.path.basename(file)
			strain = strain.replace(".csv","")

			results_table =  pandas.read_table(os.path.join(results_dir,file),sep=self.sep,header=0,na_values="",keep_default_na=True)
			results_table = results_table.replace(np.nan, '', regex=True)

			with open(os.path.join(results_dir ,file)) as f:
				line_num = 1
				for line in f:
					if line.startswith("Domain"):
						continue
					ID = strain + "|" + str(line_num)
					if dir_name == "":
						self.hits[ID] = Hit(ID,results_table[line_num-1:line_num]) ## save a list of all the values of this hit for each hit
					else:
						self.unfits[ID] =  Hit(ID,results_table[line_num-1:line_num]) 
					line_num += 1

	def _write_toxin_output(self,label,curr_nodes, num_copies,domains,scores,toxins_lengths):
		# calc general properties
		utils.assure_path_exists(os.path.join(self.group_dir,"hits_clusters"))
		outfile = open(os.path.join(self.group_dir,"hits_clusters",label +".txt"),"w")
		outfile.write("##  ID: " + label + "\n")
		outfile.write("##  Num_Strains: " + str(num_copies) + "\n")
		outfile.write("##  Domains: " + ",".join(domains)+ "\n")	
		outfile.write("##  Average Hit Length: " +str(np.mean(toxins_lengths))+ "\n")
		outfile.write("##  Average HMMER score: " + str(np.mean(scores)) + "\n")

		outfile.write(self.sep.join(['Strain','Upstream_Cluster','Downstream_Cluster','Domain', 'HMMER_score', 'Contig', 'Strand',
		 'Hit_length', 'Hit_start', 'Hit_stop', 'Upstream_length', 'Upstream_delta', 'Upstream_start', 
			'Upstream_stop', 'Downstream_length', 'Downstream_delta', 'Downstream_start', 
			'Downstream_stop', 'Hit', 'Upstream', 'Downstream', 'source']) + "\n")
		for n in curr_nodes:
			hit_toks = n.split("*")
			hit_ID = hit_toks[0]
			curr_hit = self.hits[hit_ID]
			toks = curr_hit.toks
			strain = hit_ID.split("|")[0]
			outfile.write(self.sep.join(map(str,[strain,curr_hit.clusters["upstream"],curr_hit.clusters["downstream"],toks['Domain'].tolist()[0], toks['HMMER_score'].tolist()[0],
				toks['Contig'].tolist()[0], toks['Strand'].tolist()[0], toks['Hit_length'].tolist()[0], toks['Hit_start'].tolist()[0],toks['Hit_stop'].tolist()[0],
				toks['Upstream_length'].tolist()[0],toks['Upstream_delta'].tolist()[0],
				toks['Upstream_start'].tolist()[0],toks['Upstream_stop'].tolist()[0],toks['Downstream_length'].tolist()[0],
				toks['Downstream_delta'].tolist()[0],toks['Downstream_start'].tolist()[0],toks['Downstream_stop'].tolist()[0],
				toks['Hit'].tolist()[0],toks['Upstream'].tolist()[0],toks['Downstream'].tolist()[0],toks['Source'].tolist()[0]])) + "\n")
		outfile.close()


	def _write_unfit_output(self, label, curr_nodes, num_copies, domains, scores, lengths, reasons):

		utils.assure_path_exists(os.path.join(self.group_dir,"unfit_clusters"))
		outfile = open(os.path.join(self.group_dir,"unfit_clusters" ,label +".txt"),"w")
		outfile.write("##  ID: " + label + "\n")
		outfile.write("##  Num Copies: " + str(num_copies) + "\n")
		outfile.write("##  Domains: " + ",".join(domains)+ "\n")	
		outfile.write("##  Average Hit Length: " +str(np.mean(lengths))+ "\n")
		outfile.write("##  Average HMMER score: " + str(np.mean(scores)) + "\n")
		outfile.write("## Reasons: " + ",".join(map(str,list(reasons))) + "\n")
		## add information on each of the hits seperately
		outfile.write(self.sep.join(["Strain","Unfit_Cluster","Hit_Clusters","Domain","HMMER_score","Contig","Strand","Hit_Length","Hit_Start","Hit_Stop","Sequence","Reason1","Reason2"])+"\n")
		for n in curr_nodes:
			hit_toks = n.split("*")
			hit_ID = hit_toks[0]
			curr_hit = self.unfits[hit_ID]
			toks = curr_hit.toks
			strain = hit_ID.split("|")[0]
			unfit_cluster = curr_hit.clusters["unfit"]
			hit_clusters = curr_hit.clusters["hit"]
			outfile.write(self.sep.join(map(str,[strain, unfit_cluster, hit_clusters ,toks['Domain'].tolist()[0],toks['HMMER_score'].tolist()[0],
				toks['Contig'].tolist()[0],toks['Strand'].tolist()[0], toks['Hit_length'].tolist()[0],toks['Hit_start'].tolist()[0],toks['Hit_stop'].tolist()[0],
				toks['Sequence'].tolist()[0],toks['Reason1'].tolist()[0], toks['Reason2'].tolist()[0]] )) + "\n")
		outfile.close()

	def _write_antitoxin_output(self,network_type,label,curr_nodes,num_copies,domains,scores,toxins_lengths,directions,antitoxins_lengths,deltas):
		utils.assure_path_exists(os.path.join(self.group_dir,network_type+"_clusters"))
		outfile = open(os.path.join(self.group_dir,network_type+"_clusters",label +".txt"),"w")
		outfile.write("##  ID: " + label + "\n")
		outfile.write("##  Num Copies: " + str(num_copies) + "\n")
		outfile.write("##  Domains: " + ",".join(domains)+ "\n")	
		outfile.write("##  Average Hit Length: " +str(np.mean(toxins_lengths))+ "\n")
		outfile.write("##  Average HMMER score: " + str(np.mean(scores)) + "\n")
		if network_type != "hits":
			outfile.write("##  Average Delta: " +str(np.mean(deltas))+ "\n")
			outfile.write("##  Average Partner Length: " +str(np.mean(antitoxins_lengths))+ "\n")
		if self.order == "either":
			outfile.write("##  Order: Upstream: " + str(directions["standard"]) + " Downstream: " + str(directions["reverse"]) + "\n")
		## add information on each of the hits seperately
		outfile.write(self.sep.join(["Strain","Hit_Cluster","Domain","HMMER_score","Order","Contig","Strand","Partner_Length",
			"Delta","Partner_Start","Partner_Stop","Hit_Length","Hit_Start","Hit_Stop","Partner","Sequence","Source"])+"\n")
		for n in curr_nodes:
			hit_toks = n.split("*")
			hit_ID = hit_toks[0]
			curr_hit = self.hits[hit_ID]
			toks = curr_hit.toks
			strain = hit_ID.split("|")[0]
			cluster = curr_hit.clusters["hit"]
			if hit_toks[1] == "upstream" or self.order == "upstream":
				if self.order == "upstream":
					cluster =  curr_hit.clusters["upstream"]
				outfile.write(self.sep.join(map(str,[strain, cluster,toks['Domain'].tolist()[0],toks['HMMER_score'].tolist()[0],"upstream",
					toks['Contig'].tolist()[0],toks['Strand'].tolist()[0],toks['Upstream_length'].tolist()[0],toks['Upstream_delta'].tolist()[0],toks['Upstream_start'].tolist()[0],
					toks['Upstream_stop'].tolist()[0],toks['Hit_length'].tolist()[0],toks['Hit_start'].tolist()[0],toks['Hit_stop'].tolist()[0],toks['Upstream'].tolist()[0],
					toks['Hit'].tolist()[0],toks['Source'].tolist()[0]])) + "\n")
			else:
				if self.order == "downstream":
					cluster =  curr_hit.clusters["downstream"]
				outfile.write(self.sep.join(map(str,[strain, cluster,toks['Domain'].tolist()[0],toks['HMMER_score'].tolist()[0],"downstream",
					toks['Contig'].tolist()[0],toks['Strand'].tolist()[0],toks['Downstream_length'].tolist()[0],toks['Downstream_delta'].tolist()[0],toks['Downstream_start'].tolist()[0],
					toks['Downstream_stop'].tolist()[0],toks['Hit_length'].tolist()[0],toks['Hit_start'].tolist()[0],toks['Hit_stop'].tolist()[0],toks['Downstream'].tolist()[0],
					toks['Hit'].tolist()[0],toks['Source'].tolist()[0]])) + "\n")
		outfile.close()

	def _update_unfit_label(self, unfit_dict, n, label, cluster_id):
		new_clusters = map(str,list(unfit_dict[n]))
		new_clusters = [s + "H" for s in new_clusters]
		
		if ("(") in label: # merge previous clusters with new ones
			prev_clusters = label.split("(")[-1]
			prev_clusters = prev_clusters.split(")")[0]
			prev_clusters = prev_clusters.split("+")
			new_clusters = list(set(new_clusters + prev_clusters))
		
		return str(cluster_id) + "UNFIT(" + "+".join(new_clusters)  + ")"

	def _write_file_per_cluster(self,components,network_type, unfit_dict = {}):
		cluster_id = 1
		for c in components:
			if network_type=="hits":
				label =  str(cluster_id) + "H"
			elif network_type == "partners":
				label =  str(cluster_id) + "P"
			elif network_type == "upstream":
				label =  str(cluster_id) + "U"
			elif network_type == "unfit":
				label = str(cluster_id) + "UNFIT"
			else:
				label =  str(cluster_id) + "D"
			
			curr_nodes = c.nodes()
			## get attributes of current CC
			directions = {"standard":0,"reverse":0}
			domains=set()
			scores=[]
			deltas=[]
			toxins_lengths=[]
			antitoxins_lengths=[]
			reasons = set()
			num_copies = 0
			for n in curr_nodes:
				
				if network_type == "unfit" and n in unfit_dict: ## this node already has an ID from the hits		
					label = self._update_unfit_label(unfit_dict, n, label, cluster_id)

				hit_toks = n.split("*")	
				hit_ID = hit_toks[0]
				hit_type = hit_toks[1]

				if network_type == "unfit":
					hits = self.unfits
				else:
					hits = self.hits

				curr_hit = hits[hit_ID]

				num_copies+=1
				domains.add(curr_hit.toks['Domain'].tolist()[0]) 
				scores.append(float(curr_hit.toks['HMMER_score'].tolist()[0]))
				toxins_lengths.append(int(curr_hit.toks['Hit_length'].tolist()[0]))

				if hit_type == "upstream":
					directions["standard"] += 1
					deltas.append(float(curr_hit.toks['Upstream_delta'].tolist()[0]))
					antitoxins_lengths.append(int(curr_hit.toks['Upstream_length'].tolist()[0]))
				elif hit_type == "downstream":
					directions["reverse"] += 1
					deltas.append(float(curr_hit.toks['Downstream_delta'].tolist()[0]))
					antitoxins_lengths.append(int(curr_hit.toks['Downstream_length'].tolist()[0]))

				if network_type == "unfit":
					reasons.add(curr_hit.toks['Reason1'].tolist()[0])
					
					reason2 = str(curr_hit.toks['Reason2'].tolist()[0]) ## if using only upstream or downstream, the second reason will be "NaN"
					if not reason2=="nan":
						reasons.add(reason2)


			## write all results of this cluster
			if network_type == "hits" and (self.order == "either" or self.order == "both"):
				self._write_toxin_output(label,curr_nodes, num_copies,domains,scores,toxins_lengths)
			elif network_type == "unfit":
				self._write_unfit_output(label,curr_nodes, num_copies, domains, scores, toxins_lengths, reasons)
			else:
				self._write_antitoxin_output(network_type,label,curr_nodes, num_copies,domains,scores,toxins_lengths,directions,antitoxins_lengths,deltas)
			cluster_id += 1

	def _write_matrix_to_file(self,network_type, binary):
		# when completed going over all components, save the matrix
		utils.assure_path_exists(os.path.join(self.group_dir ,network_type + "_clusters"))
		
		with open(os.path.join(self.group_dir, network_type +".csv"), "wb") as f:
			writer = csv.writer(f)
			writer.writerows(binary)
		
		if self.save_to_ITOL:
			
			## get the maximum value of the matrix for ITOL
			max_val = [item for sublist in binary for item in sublist]
			max_val = [x for x in max_val if isinstance(x, int)]
			max_val = max(max_val)

			utils.assure_path_exists(os.path.join(self.group_dir, "ITOL"))
			
			with open(os.path.join(self.group_dir,"ITOL",network_type + ".txt"),"wb") as f:
				f.write("DATASET_HEATMAP\nSEPARATOR COMMA\nDATASET_LABEL,"+network_type+"\nCOLOR,#ff0000\nFIELD_LABELS," + ",".join(binary[0][1:]) + "\nDATA\n")
				for i in range(1,len(binary)):
					f.write(",".join(map(str,binary[i]))+ "\n")
			strains = [row[0] for row in binary][1:]
			for i in range(1,len(binary[0])):
				feature_vec = [row[i] for row in binary]
				ta = feature_vec[0]
				feature_vec = feature_vec[1:]
				self._annotate_system(feature_vec,ta,strains,network_type,max_val)

	def _annotate_system(self,feature_vec,ta,strains,network_type,max_val):

		## preset the colors to be loaded in ITOL
		max_colors = {"complete": "#005900","hits":"#3b1365","upstream":"#000099","partners":"#000099","downstream":"#cf4c0b","unfit":"#990000"}

		outdir = os.path.join(self.group_dir, "ITOL/", network_type + "_clusters")
		utils.assure_path_exists(outdir)
		out = open(os.path.join(outdir, ta  + ".txt"),"w")

		out.write("DATASET_HEATMAP\nSEPARATOR COMMA\nDATASET_LABEL," + ta + "\nCOLOR,#ff0000\nFIELD_LABELS," + ta +
			"\nCOLOR_MIN,#eeeded\nCOLOR_MAX," + max_colors[network_type] + "\nUSER_MIN_VALUE,0\nUSER_MAX_VALUE," +str(max_val) +"\nDATA\n")
		for i in range(0,len(strains)):
			out.write(strains[i]+ "," + str(feature_vec[i]) + "\n") 
		out.close()


	def _write_graph_to_matrix(self,G,network_type, unfit_dict = {}):

		components = list(nx.connected_component_subgraphs(G)) # get connected components, each cc is a cluster	 
		
		binary= [["Strain"] + map(str,range(1,1+ len(components)))] # initate binary annotation matrix

		for strain in self.strains:
			binary.append([strain] + [0] * len(components))

		## go over the connected components
		cluster_id=1
		for c in components: 
			curr_nodes =  c.nodes() # all nodes of current CC
			if network_type=="hits":
				label =  str(cluster_id) + "H"
			elif network_type == "partners":
				label =  str(cluster_id) + "P"
			elif network_type == "upstream":
				label = str(cluster_id) + "U"
			elif network_type == "unfit":
				label = str(cluster_id) + "UNFIT"
			else:
				label = str(cluster_id) + "D"
			
			binary[0][cluster_id] = label

			for n in curr_nodes:
				
				### for the unfit network, check which cluster these hits would have belonged to if they were fit
				if network_type == "unfit" and n in unfit_dict: ## this node already has an ID from the hits			
					label = self._update_unfit_label(unfit_dict, n, label, cluster_id)
				else:
					G[n]["cluster_id"] = cluster_id
				
				hit_toks = n.split("*")
				
				hit_ID = hit_toks[0]
				hit_type = hit_toks[1]
				strain = hit_ID.split("|")[0]
				row = self.strains.index(strain) + 1

				if network_type == "unfit":
					hits = self.unfits
				else:
					hits = self.hits
					
				curr_hit = hits[hit_ID]
				
				curr_hit.clusters[hit_type] = label # label the cluster of this hit

				if network_type == "unfit" and n in unfit_dict:
					clusters = [s + "H" for s in map(str,list(unfit_dict[n]))]
					curr_hit.clusters["hit"] = "+".join(clusters)
					curr_hit.clusters[hit_type] = cluster_id

				
				domain = hits[hit_ID].toks['Domain'].tolist()[0]
				
				if network_type == "unfit":
					binary[0][cluster_id] = label
				
				curr_domains = binary[0][cluster_id]
				if curr_domains.find(domain) < 0:  # add a new domain to this label
					curr_domains = str(curr_domains) + "-" + domain  
					binary[0][cluster_id] = curr_domains
				
				binary[row][cluster_id] += 1 # update matrix in this TA system for this strain
			cluster_id += 1
		

		self._write_matrix_to_file(network_type, binary)
		
		return (components, G)


	def _read_blast_to_network(self,network_type):
		
		results =  os.path.join(self.blast_dir, network_type + "_blast_results" )
		
		G = nx.Graph()

		with open(results) as f: # add the edges to the graph
			for line in f:
				line = line.strip().split()
				id1=line[0]
				id2=line[1]
				node_type_1 = id1.split("*")[1]
				node_type_2 = id2.split("*")[1]
				score = float(line[2])

				if score >= self.min_identity and node_type_1 != "unfit" and node_type_2 != "unfit": # only add hits to the network
					G.add_edge(id1,id2) 


		return self._write_graph_to_matrix(G,network_type)

	def _write_line_to_complete(self,toxin,key,curr_files,line,header,completes,domain):
		key = toxin + "_" + key
		if key not in curr_files:
			completes[key] = []
			curr_files[key] = open(os.path.join(self.group_dir,"complete_clusters",key + ".txt"),"w")
			curr_files[key].write(header)
		completes[key].append(line.split(self.sep)[0])
		curr_files[key].write(line)

	def _write_line_to_both(self,toxin,upstream,downstream,curr_files,line,header,completes,domain):
		if self.order == "either" and upstream == "-":
			key = toxin + "_" + downstream
		elif self.order == "either" and downstream == "-":
			key = upstream + "_" + toxin
		else:
			key = upstream + "_" + toxin + "_" + downstream
		if key not in curr_files:
			completes[key] = []
			curr_files[key] = open(os.path.join(self.group_dir, "complete_clusters/", key + ".txt"),"w")
			curr_files[key].write(header)
		completes[key].append(line.split(self.sep)[0])
		curr_files[key].write(line)

	def _create_complete_files(self):
		completes = {}
		utils.assure_path_exists(os.path.join(self.group_dir ,"complete_clusters"))
		for file in os.listdir(os.path.join(self.group_dir,"hits_clusters")):
			curr_files = {}
			if not file.endswith(".txt"):
				continue
			with open(os.path.join(self.group_dir,"hits_clusters",file)) as f:
				toxin = file.split(".")[0]
				for line in f:
					if line.startswith("#"):
						continue
					if line.startswith("Strain"):
						header = line
						continue
					toks = line.strip().split(self.sep)
					domain = toks[3]
					if self.order == "upstream" or self.order == "downstream":
						self._write_line_to_complete(toxin,toks[1],curr_files,line,header,completes,domain)
					else:
						self._write_line_to_both(toxin,toks[1],toks[2],curr_files,line,header,completes,domain)

			for c_file in curr_files:
				curr_files[c_file].close()

		## create a matrix of the results
		binary= [["Strain"] + map(str,range(1,1+ len(completes)))]
		for strain in self.strains:
			binary.append([strain] + [0] * len(completes))
		curr_column = 1
		for complete in completes:
			binary[0][curr_column] = complete
			for strain in completes[complete]:
				row = self.strains.index(strain) + 1
				binary[row][curr_column] += 1
			curr_column += 1

		self._write_matrix_to_file("complete", binary)


	def _create_unfits_network(self,hits):
		results = os.path.join(self.blast_dir,"hits_blast_results" )

		G = nx.Graph()

		unfit_dict = {} # dict of the clusters of an unfit hit
		with open(results) as f: # add the edges to the graph
			for line in f:
				line = line.strip().split()
				id1=line[0]
				id2=line[1]
				node_type_1 = id1.split("*")[1]
				node_type_2 = id2.split("*")[1]
				score = float(line[2])

				if score < self.min_identity:
					continue
				
				if node_type_1 == "unfit" and node_type_2 == "unfit":
					G.add_edge(id1,id2)
				elif node_type_1 == "unfit" and node_type_2 == "hit":
					if id1 not in unfit_dict:
						unfit_dict[id1] = set()
					unfit_dict[id1].add(hits[id2]["cluster_id"])
				elif node_type_2 == "unfit" and node_type_1 == "hit":
					if id2 not in unfit_dict:
						unfit_dict[id2] = set()
					unfit_dict[id2].add(hits[id1]["cluster_id"])
		return (G, unfit_dict)


	def _report_unfit(self,hits):
		if not self.report_unfit:
			return
		self._parse_hits_files("UNFIT") # read all the unfit hits into a file
		unfit_network, unfit_dict = self._create_unfits_network(hits)
		unfit_components, unfit_network = self._write_graph_to_matrix(unfit_network,"unfit", unfit_dict)
		self._write_file_per_cluster(unfit_components,"unfit", unfit_dict = unfit_dict)


	
	def run(self):
		# create a list of all the strains
		for file in os.listdir(self.filter_dir):
			if file.endswith(".csv"):
				strain = os.path.basename(file)
				strain = strain.replace(".csv","")
				self.strains.append(strain)

		self._parse_hits_files("") # create an object for each hit

		hits, hits_graph = self._read_blast_to_network("hits")
		
		if self.order != "both":
			partners, partners_graph = self._read_blast_to_network("partners")
			self._write_file_per_cluster(partners,"partners")
		else:
			upstream, upstream_graph = self._read_blast_to_network("upstream")
			downstream, downstream_graph = self._read_blast_to_network("downstream")
			self._write_file_per_cluster(upstream,"upstream")
			self._write_file_per_cluster(downstream,"downstream")	
		self._write_file_per_cluster(hits,"hits")	
		self._create_complete_files()

		self._report_unfit(hits_graph)

class Hit():

	def __init__(self,ID,toks):
		self.ID = ID
		self.toks = toks # a pandas table with all the tokens of this hit
		self.clusters = {"hit":"-", "upstream":"-","downstream":"-"}


