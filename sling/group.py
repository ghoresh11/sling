import utils
import os
import subprocess
import warnings
import networkx as nx
import csv
import sys
import numpy as np

class Error (Exception):
    pass


## quick function to get the versions of BLAST used for log files
def get_blast_version(command, name):
    ''' get blast version being used'''
    p = subprocess.Popen([command, "-h"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    rc = p.returncode
    output = output.split()
    for i in range(0,len(output)):
        if "version" == output[i] and name=="makeblastdb":
            return ("BLAST-" + output[i+1])
        elif "version" in output[i] and name == "blastp":
            return ("BLAST-" + output[i+4])
    return ("BLAST version not found\n")

def hits_to_fasta(args, out_dir, filter_dir):
    '''combine all the hits and write them into
    a single fasta file to run blast'''
    hits = open(os.path.join(out_dir,"hits.fasta"),"w")
    if args.order == "both": ## create two seperate fasta files for up and downstream
        downstream = open(os.path.join(out_dir,"downstream.fasta"), "w")
        upstream = open(os.path.join(out_dir ,"upstream.fasta"),"w")
    else: # otherwise treat them the same
        partners = open(os.path.join(out_dir ,"partners.fasta"),"w")

    for file in os.listdir(filter_dir):
        if not file.endswith(".csv"): ## not a results file
            continue

        strain = os.path.basename(file)
        strain = strain.replace(".csv","")

        with open(os.path.join(filter_dir,file)) as f:
            line_num = 0
            for line in f:
                toks = line.strip().split(args.sep)
                if line.startswith("Domain"):
                    hit_index = toks.index("Hit")
                    upstream_index = toks.index("Upstream")
                    downstream_index = toks.index("Downstream")
                    continue
                line_num += 1
                identifier = ">" + strain + "|" + str(line_num)
                ## hits always exist
                hits.write(identifier + "*hit" + "\n" + toks[hit_index] + "\n")

                if args.order == "both": ## with "both" there are 2 output files
                    upstream.write(identifier +"*upstream"+ "\n" + toks[upstream_index] + "\n")
                    downstream.write(identifier + "*downstream" + "\n" + toks[downstream_index] + "\n")
                else:  ## otherwise, write them if found
                    if upstream_index > 0 and toks[upstream_index] != "-":
                        partners.write(identifier +"*upstream"+ "\n" + toks[upstream_index] + "\n")
                    if downstream_index > 0 and toks[downstream_index] != "-":
                        partners.write(identifier +"*downstream"+ "\n" + toks[downstream_index] + "\n")

    #### UNFIT ####
    ## if in the previous step, the unfit were also reported, add them to the "hits" in the network analysis
    if args.report_unfit and not os.path.exists(os.path.join(filter_dir,"UNFIT")):
        warnings.warn("Could not find UNFIT files from SUMMARISE step. To report unfit, turn on --report_unfit / -u flag in FILTER and run again.")
    elif args.report_unfit:
        for file in os.listdir(os.path.join(filter_dir,"UNFIT")):
            if not file.endswith(".csv"):
                continue
            with open(os.path.join(filter_dir,"UNFIT",file)) as f:
                strain = os.path.basename(file)
                strain = strain.replace(".csv","")
                line_num = 0
                for line in f:
                    toks = line.strip().split(args.sep)
                    if line.startswith("Strain") or line.startswith("Domain"):
                        hit_index = toks.index("Sequence")
                        continue
                    line_num += 1
                    strain = strain
                    identifier = ">" + strain + "|" + str(line_num) + "*unfit"
                    hits.write(identifier + "\n" + toks[hit_index] + "\n")
    ## close all the files
    hits.close()
    if args.order == "both":
        upstream.close()
        downstream.close()
    else:
        partners.close()
    return

def call_blast_command(args, out_dir, file_type,):
    ''' run blast on hits, partnerts, downstream or upstream'''
    command = map(str,[args.makeblastdb, "-in",os.path.join(out_dir, file_type + ".fasta"),"-dbtype","prot"] )
    command2 = map(str,[args.blastp, "-db",os.path.join(out_dir, file_type + ".fasta"),
            "-query", os.path.join(out_dir,file_type + ".fasta"),"-out",
            os.path.join(out_dir,file_type + "_blast_results"),
            "-outfmt","6 qseqid sseqid pident length qlen slen evalue bitscore","-evalue",args.min_blast_evalue, "-num_threads", args.cpu])
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
    return

def run_blast(args, group_dir, filter_dir, out_dir):
    ''' run blast on all hits and partners'''
    utils.assure_path_exists(out_dir)
    hits_to_fasta(args, out_dir, filter_dir)

    call_blast_command(args, out_dir, "hits")
    if args.order == "both":
        call_blast_command(args, out_dir, "upstream")
        call_blast_command(args, out_dir, "downstream")
    else:
        call_blast_command(args, out_dir, "partners")
    return


def parse_hits_files(args, results_dir):
    ''' read all the outputs of filter into hits'''
    hit_dict = {}
    for file in os.listdir(results_dir):
        if not file.endswith(".csv"): ## not a hits file
            continue
        strain = os.path.basename(file)
        strain = strain.replace(".csv","")

        with open(os.path.join(results_dir ,file)) as f:
            line_num = 1
            for line in f:
                toks = line.strip().split(args.sep)
                if line.startswith("Domain"):
                    keys = toks # the header of the file are the values
                    continue
                ID = strain + "|" + str(line_num)
                hit_dict[ID] = {}
                for k in range(len(keys)):
                    hit_dict[ID][keys[k]] = toks[k]
                hit_dict[ID]["clusters"] = {"hit":"-", "upstream":"-", "downstream":"-"}
                line_num += 1
    return hit_dict, keys


def label_cluster(cluster_id, network_type):
    ''' label the cluster depending on the network type'''
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
    return label


def generate_single_itol_output(group_dir, feature_vec, ta,strains, network_type, max_val):
    ''' generate a single output file for a cluster to be loaded into ITOL'''
    ## colors to be loaded in ITOL
    max_colors = {"complete": "#005900","hits":"#3b1365","upstream":"#000099","partners":"#000099","downstream":"#cf4c0b","unfit":"#990000"}
    outdir = os.path.join(group_dir, "ITOL/", network_type + "_clusters")
    utils.assure_path_exists(outdir)
    out = open(os.path.join(outdir, ta  + ".txt"),"w")
    out.write("DATASET_HEATMAP\nSEPARATOR COMMA\nDATASET_LABEL," + ta + "\nCOLOR,#ff0000\nFIELD_LABELS," + ta +	"\nCOLOR_MIN,#eeeded\nCOLOR_MAX," + max_colors[network_type] + "\nUSER_MIN_VALUE,0\nUSER_MAX_VALUE," +str(max_val) +"\nDATA\n")
    for i in range(0,len(strains)):
        out.write(strains[i]+ "," + str(feature_vec[i]) + "\n")
    out.close()
    return

def write_matrix_to_file(args, group_dir, network_type, binary):
    ''' write the matrix to a file, both as a CSV
    and for reading into ITOL if asked'''
    with open(os.path.join(group_dir, network_type +".csv"), "wb") as f:
        writer = csv.writer(f)
        writer.writerows(binary)

    if args.itol:
        ## get the maximum value of the matrix for ITOL
        max_val = [item for sublist in binary for item in sublist]
        max_val = [x for x in max_val if isinstance(x, int)]
        max_val = max(max_val)
        utils.assure_path_exists(os.path.join(group_dir, "ITOL"))
        with open(os.path.join(group_dir,"ITOL",network_type + ".txt"),"wb") as f:
            f.write("DATASET_HEATMAP\nSEPARATOR COMMA\nDATASET_LABEL,"+network_type+"\nCOLOR,#ff0000\nFIELD_LABELS," + ",".join(binary[0][1:]) + "\nDATA\n")
            for i in range(1,len(binary)):
                f.write(",".join(map(str,binary[i]))+ "\n")
        strains = [row[0] for row in binary][1:]
        for i in range(1,len(binary[0])):
            feature_vec = [row[i] for row in binary]
            ta = feature_vec[0]
            feature_vec = feature_vec[1:]
            generate_single_itol_output(group_dir, feature_vec, ta, strains, network_type, max_val) ## create a single file for one system
    return

def update_unfit_label(args, unfits, n, label, cluster_id):
    ''' update the label of an UNFIT cluster to include the
    hits it was match with
    return: the label'''
    n = n.split("*")[0] ## get the node ID
    mapped_clusters = unfits[n]["clusters"]["hit"]
    if mapped_clusters == "-":
        return label
    new_clusters = map(str,list(mapped_clusters))
    new_clusters = [s + "H" for s in new_clusters]

    if ("(") in label: # merge previous clusters with new ones
        prev_clusters = label.split("(")[-1]
        prev_clusters = prev_clusters.split(")")[0]
        prev_clusters = prev_clusters.split("+")
        new_clusters = list(set(new_clusters + prev_clusters))
    return str(cluster_id) + "UNFIT(" + "+".join(new_clusters)  + ")"


def convert_graph_to_matrix(args, G, network_type, strains, hits):
    ''' write the networkx graph into a CSV file where tows are
    clusters and columns are the strains
    Return: 1. the connected components of the Graph
            2. A binary presence absence matrix'''
    components = list(sorted(nx.connected_components(G), key = len, reverse=True))
    binary= [["Strain"] + map(str,range(1,1+ len(components)))] # initate binary annotation matrix

    for strain in strains:
        binary.append([strain] + [0] * len(components))

    ## go over the connected components
    cluster_id=0
    for curr_nodes in components:
        cluster_id += 1
        label = label_cluster(cluster_id, network_type)
        binary[0][cluster_id] = label
        for n in curr_nodes:
            ### for the unfit network, check which cluster these hits would have belonged to if they were fit
            if network_type == "unfit":
                label = update_unfit_label(args, hits, n, label, cluster_id)

            G.nodes[n]["cluster_id"] = cluster_id

            hit_toks = n.split("*") ## hits look like this "NC_022083.1|5*hit"
            hit_ID = hit_toks[0] ## NC_022083.1|5
            hit_type = hit_toks[1] ## hit
            strain = hit_ID.split("|")[0] ## NC_022083.1
            row = strains.index(strain) + 1
            curr_hit = hits[hit_ID]
            curr_hit["clusters"][hit_type] = label # label the cluster of this hit

            if network_type == "unfit":
                curr_hit["clusters"][hit_type] = cluster_id
                binary[0][cluster_id] = label

            domain = hits[hit_ID]["Domain"] ## add this domain to this cluster label
            curr_domains = binary[0][cluster_id]
            if curr_domains.find(domain) < 0:  # add a new domain to this label
                curr_domains = str(curr_domains) + "-" + domain
                binary[0][cluster_id] = curr_domains

            binary[row][cluster_id] += 1 # update matrix in this TA system for this strain

    return components, binary



def is_blast_match(toks, args):
    ''' read a line from the blast output
    return True if sequences should have an edge,
    False otherwise'''
    if toks[0] == toks[1]:
        return False
    curr_ident = float(toks[2])
    curr_coverage = float(toks[3]) / float(toks[4])
    if curr_ident >= args.min_identity and curr_coverage >= args.length_coverage:
        return True
    return False


def read_blast_to_network(args, network_type, blast_dir):
    ''' create a network from the blast output files
    using the identity requirements provided by the user'''
    results =  os.path.join(blast_dir, network_type + "_blast_results" )
    G = nx.Graph()
    with open(results) as f: # add the edges to the graph
        for line in f:
            toks = line.strip().split()
            id1=toks[0]
            id2=toks[1]
            node_type_1 = id1.split("*")[1]
            node_type_2 = id2.split("*")[1]
            if is_blast_match(toks, args) and node_type_1 != "unfit" and node_type_2 != "unfit": # only add hits to the network
                G.add_edge(id1,id2)
    return G


def write_single_output(args, network_type, keys, hits, group_dir, label, curr_nodes, num_copies, domains, scores, toxins_lengths, antitoxins_lengths, deltas, directions, reasons):
    ''' create the output file for a toxin cluster
    Slightly messy because the outputs are different,
    depending if its a toxin, antitoxin or both...'''
    utils.assure_path_exists(os.path.join(group_dir,network_type + "_clusters"))
    outfile = open(os.path.join(group_dir,network_type + "_clusters",label +".txt"),"w")

    ### Attributes of this cluster ####
    outfile.write("#  ID: " + label + "\n")
    outfile.write("#  Num_Strains: " + str(num_copies) + "\n")
    outfile.write("#  Domains: " + ",".join(domains)+ "\n")
    outfile.write("#  Average Hit Length: " +str(np.mean(toxins_lengths))+ "\n")
    outfile.write("#  Average HMMER score: " + str(np.mean(scores)) + "\n")
    if network_type in ["partners", "upstream", "downstream"]:
        outfile.write("##  Average Delta: " +str(np.mean(deltas))+ "\n")
        outfile.write("##  Average Partner Length: " +str(np.mean(antitoxins_lengths))+ "\n")
    if args.order == "either" and network_type != "unfit":
        outfile.write("#  Order: Upstream: " + str(directions["standard"]) + " Downstream: " + str(directions["reverse"]) + "\n")
    if network_type == "unfit":
        outfile.write("# Reasons: " + ",".join(map(str,list(reasons))) + "\n")

    ### HEADER ###
    if network_type == "hits":
        if args.order in ["either", "both"]:
            outfile.write(args.sep.join(['Strain','Upstream_Cluster','Downstream_Cluster'] + keys) + "\n")
        else: ## either upstream or downstream
            ## only upstream antitoxin
            outfile.write(args.sep.join(['Strain','Partner_Cluster'] + keys) + "\n")
    elif network_type == "unfit":
        outfile.write(args.sep.join(['Strain', 'Hit_Cluster'] + keys) + "\n")
    else: # network_type in ["partners", "upstream", "downstream"]: ## antitoxins output
        ## only antitoxins
        outfile.write(args.sep.join(["Strain","Hit_Cluster","Domain","HMMER_score","Order","Contig","Strand","Partner_Length",
        "Delta","Partner_Start","Partner_Stop","Hit_Length","Hit_Start","Hit_Stop","Partner","Sequence","Source"]) + "\n")

    ### one line per member of this cluster ###
    for n in curr_nodes:
        hit_toks = n.split("*")
        hit_type = hit_toks[1]
        hit_ID = hit_toks[0]
        curr_hit = hits[hit_ID]
        strain = hit_ID.split("|")[0]
        if hit_type == "hit": ## toxins
            if args.order in ["either", "both"]:
                outfile.write(args.sep.join([strain, curr_hit["clusters"]["upstream"],curr_hit["clusters"]["downstream"]]))
            elif args.order == "upstream":
                outfile.write(args.sep.join([strain, curr_hit["clusters"]["upstream"]]))
            else:
                outfile.write(args.sep.join([strain, curr_hit["clusters"]["downstream"]]))
            for k in keys:
                outfile.write("," + str(curr_hit[k]))
            outfile.write("\n")
        elif hit_type == "unfit":
            hit_clusters = "+".join(list(curr_hit["clusters"]["hit"]))
            outfile.write(args.sep.join([strain, hit_clusters]))
            for k in keys:
                outfile.write("," + str(curr_hit[k]))
            outfile.write("\n")
        else: ## antitoxin
            cluster = curr_hit["clusters"]["hit"]
            if hit_type == "upstream":
                outfile.write(args.sep.join(map(str,[strain, cluster, curr_hit['Domain'], curr_hit['HMMER_score'], "upstream",
                    curr_hit['Contig'], curr_hit['Strand'], curr_hit['Upstream_length'], curr_hit['Upstream_delta'], curr_hit['Upstream_start'],
                    curr_hit['Upstream_stop'],curr_hit['Hit_length'],curr_hit['Hit_start'],curr_hit['Hit_stop'],curr_hit['Upstream'],
                    curr_hit['Hit'],curr_hit['Source']])) + "\n")
            else: ## Downstream
                outfile.write(args.sep.join(map(str,[strain, cluster, curr_hit['Domain'], curr_hit['HMMER_score'], "downstream",
                    curr_hit['Contig'], curr_hit['Strand'], curr_hit['Downstream_length'], curr_hit['Downstream_delta'], curr_hit['Downstream_start'],
                    curr_hit['Downstream_stop'],curr_hit['Hit_length'],curr_hit['Hit_start'],curr_hit['Hit_stop'],curr_hit['Downstream'],
                    curr_hit['Hit'],curr_hit['Source']])) + "\n")
    outfile.close()
    return

def write_file_per_cluster(args, group_dir, keys, hits, components, network_type, unfit_dict = {}):
    ''' create a single output file for one SLING cluster'''
    cluster_id = 0
    for curr_nodes in components:
        cluster_id += 1
        label = label_cluster(cluster_id, network_type)
        ## initiate attributes current CC
        directions = {"standard":0,"reverse":0}
        domains=set()
        scores=[]
        deltas=[]
        toxins_lengths=[]
        antitoxins_lengths=[]
        reasons = set()
        num_copies = 0
        for n in curr_nodes:
            hit_toks = n.split("*") ## hits look like this "NC_022083.1|5*hit"
            hit_ID = hit_toks[0]
            hit_type = hit_toks[1]
            curr_hit = hits[hit_ID]
            num_copies+=1

            domains.add(curr_hit["Domain"])
            scores.append(float(curr_hit["HMMER_score"]))
            toxins_lengths.append(int(curr_hit["Hit_length"]))

            if hit_type == "upstream":
                directions["standard"] += 1
                deltas.append(float(curr_hit['Upstream_delta']))
                antitoxins_lengths.append(int(curr_hit['Upstream_length']))
            elif hit_type == "downstream":
                directions["reverse"] += 1
                deltas.append(float(curr_hit['Downstream_delta']))
                antitoxins_lengths.append(int(curr_hit['Downstream_length']))

            if network_type == "unfit":
                reasons.add(curr_hit['Reasons'])

        ## write all results of this cluster
        write_single_output(args, network_type, keys, hits, group_dir, label,
            curr_nodes, num_copies, domains, scores, toxins_lengths,
            antitoxins_lengths, deltas, directions, reasons)
    return


def summarise_blast_resulsts(args, network_type, group_dir, blast_dir, hits, strains):
    ''' create a network from the blast results,
    find all the connected components
    output a CSV file and all the ITOL files'''
    graph = read_blast_to_network(args, network_type, blast_dir)
    components, matrix = convert_graph_to_matrix(args, graph, network_type, strains, hits)
    write_matrix_to_file(args, group_dir, network_type, matrix)
    return components


def write_line_to_complete(args, group_dir, toxin, upstream, downstream, curr_files, line, header, completes, domain):
    ''' the output of the complete clusters looks the same as hits,
    so it requires to copy all lines that have the same hits and partners
    UPSTREAM AND DOWNSTREAM: operon only with single antitoxin
    EITHER AND BOTH: operon can have more than one partner'''
    if args.order in ["upstream", "downstream"]:
        key = toxin + "_" + upstream ## upstream key is the partner
    if args.order == "either" and upstream == "-":
        key = toxin + "_" + downstream
    elif args.order == "either" and downstream == "-":
        key = upstream + "_" + toxin
    else:
        key = upstream + "_" + toxin + "_" + downstream
    if key not in curr_files: ## new file
        completes[key] = []
        curr_files[key] = open(os.path.join(group_dir,"complete_clusters",key + ".txt"),"w")
        curr_files[key].write(header)
    curr_files[key].write(line)
    ## save tthe strain as a member of this complete cluster
    completes[key].append(line.split(args.sep)[0])
    return

def create_complete_files(args, group_dir, strains):
    ''' aggregate the hit clusters and partners clusters
    to report on the full operons
    Generate all the output files for the complete clusters'''
    completes = {}
    utils.assure_path_exists(os.path.join(group_dir ,"complete_clusters"))
    for filename in os.listdir(os.path.join(group_dir,"hits_clusters")):
        curr_files = {} ## create a new file for each unique operon
        if not filename.endswith(".txt"):
            continue
        with open(os.path.join(group_dir,"hits_clusters",filename)) as f:
            toxin = filename.split(".")[0] ## get the toxin ID
            for line in f:
                if line.startswith("#"):
                    continue
                toks = line.strip().split(args.sep)
                if line.startswith("Strain"):
                    header = line
                    domain_index = toks.index("Domain")
                    continue
                domain = toks[domain_index]
                write_line_to_complete(args, group_dir, toxin, toks[1], toks[2], curr_files, line, header, completes, domain)

    ## close all the cluster files that were opened and copied
    for c_file in curr_files:
        curr_files[c_file].close()
    ## create a matrix of the results
    binary= [["Strain"] + map(str,range(1,1+ len(completes)))]
    for strain in strains:
        binary.append([strain] + [0] * len(completes))
    curr_column = 1
    for complete in completes:
        binary[0][curr_column] = complete
        for strain in completes[complete]:
            row = strains.index(strain) + 1
            binary[row][curr_column] += 1
        curr_column += 1

    ## write the matrix to a file
    write_matrix_to_file(args,group_dir,"complete",binary)
    return

def map_hit_to_unfit(hit_id, unfit_id, unfits, hits):
    ''' given an edge between a hit and unfit,
    add that hit id to the unfit cluster'''
    unfit_id = unfit_id.split("*")[0]
    hit_id = hit_id.split("*")[0]
    if unfits[unfit_id]["clusters"]["hit"] == "-":
        unfits[unfit_id]["clusters"]["hit"] = set()
    unfits[unfit_id]["clusters"]["hit"].add(hits[hit_id]["clusters"]["hit"])
    return

def create_unfits_network(args, unfits, hits, blast_dir):
    ''' build the special network that connects the unfit hits
    to the rest of the hit clusters'''
    results = os.path.join(blast_dir,"hits_blast_results" )
    G = nx.Graph()
    with open(results) as f: # add the edges to the graph
        for line in f:
            toks = line.strip().split()
            id1 = toks[0]
            id2 = toks[1]
            node_type_1 = id1.split("*")[1]
            node_type_2 = id2.split("*")[1]
            if not is_blast_match(toks, args):
                continue
            if node_type_1 == "unfit" and node_type_2 == "unfit":
                G.add_edge(id1, id2) ## two unfits should be connected
            elif node_type_1 == "unfit" and node_type_2 == "hit":
                map_hit_to_unfit(id2, id1, unfits, hits)
            elif node_type_2 == "unfit" and node_type_1 == "hit":
                map_hit_to_unfit(id1, id2, unfits, hits)
	return G


def report_unfit(args, group_dir, filter_dir, blast_dir, hits, strains):
    ''' Generate all the output files for the unfit clusters'''
    if not args.report_unfit:
        return
    ## read all the unfit hits from the UNFIT directory
    unfit_dir = os.path.join(filter_dir,"UNFIT")
    if not os.path.exists(unfit_dir):
        warnings.warn("Could not find the output for UNFITs. Run filter is -u/--report_unfit\n")
        return
    unfits, keys = parse_hits_files(args, unfit_dir)
    unfits_graph = create_unfits_network(args, unfits, hits, blast_dir)

    unfit_components, unfit_matrix = convert_graph_to_matrix(args, unfits_graph, "unfit", strains, unfits)
    write_file_per_cluster(args, group_dir, keys, unfits, unfit_components, "unfit")
    return

def group_operons(args, group_dir, filter_dir, blast_dir):
    ''' use the blast output to group all the operons'''
    strains = []
    ## read all the strains from the output directory
    for file in os.listdir(filter_dir):
        if file.endswith(".csv"):
            strain = os.path.basename(file)
            strain = strain.replace(".csv","")
            strains.append(strain)
    ### parse filter output
    hits, keys = parse_hits_files(args, filter_dir) # create an object for each hit

    ## hits
    hits_components = summarise_blast_resulsts(args, "hits", group_dir, blast_dir, hits, strains)

    ## partners
    if args.order != "both":
        partners = summarise_blast_resulsts(args, "partners", group_dir, blast_dir, hits, strains)
        write_file_per_cluster(args, group_dir, keys, hits, partners, "partners")
    else: # in both, treat "up" and "down" seperately
        upstream = summarise_blast_resulsts(args, "upstream", group_dir, blast_dir, hits, strains)
        downstream = summarise_blast_resulsts(args, "downstream", group_dir, blast_dir, hits, strains)
        write_file_per_cluster(args, group_dir, keys, hits, upstream, "upstream")
        write_file_per_cluster(args, group_dir, keys, hits, downstream, "downstream")

    ### hits can only be written now because of partner clusters in output
    write_file_per_cluster(args, group_dir, keys, hits, hits_components, "hits")

    ## generate the ouputs of the complete operons
    create_complete_files(args, group_dir, strains)
    ## generate all the outputs to report the unfit hits
    report_unfit(args, group_dir, filter_dir, blast_dir, hits, strains)
    return


def run(args):
    if args.filter_id is None:
        args.filter_id = args.id

    args.out_dir = os.path.abspath(args.out_dir)
    filter_dir = os.path.join(args.out_dir, args.filter_id+ "_FILTER")
    group_dir = os.path.join(args.out_dir, args.id+ "_GROUP")
    utils.assure_path_exists(group_dir)
    blast_dir = os.path.join(group_dir,"blast_files")

    if args.order is None: ## user didn't override the order argument
        d = os.path.abspath(os.path.dirname(__file__))
        data_env = os.path.join(d, 'data/')
        args.order = utils.get_order(args.hmm_db,data_env)

    ## write LOG file
    get_blast_version(args.blastp, "blastp")
    get_blast_version(args.makeblastdb, "makeblastdb")
    utils.write_log(os.path.join(group_dir, "LOG"), "STEP 4 : GROUP ORFs", vars(args), "")

    run_blast(args, group_dir, filter_dir, blast_dir)
    group_operons(args, group_dir, filter_dir, blast_dir)
    return
