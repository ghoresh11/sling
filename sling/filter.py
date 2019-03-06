import utils
import os
import pandas
import sys
import copy

class Error (Exception): pass

def get_requirements(args, data_env):
	# path to the requirement document of each built in DB
    reqs = utils.databases
    if args["hmm_db"] in reqs:  # built in collections
        req_file = os.path.join(data_env, args["hmm_db"] + ".txt")
    else:  # not a built in collection, take the default values
        req_file = os.path.join(data_env, "default.txt")

    with open(req_file) as f:  # get all the values from the req file
        for line in f:
            key, val = line.strip().split()
            if key != "order":
                val = int(val)
            if args[key] is None: # otherwise, keep what the user gave
                args[key] = val
    utils.check_reqs(args)
    return


def parse_domains(hmm_db, domains_file, data_env):
    ''' get the expcted profile lengths from the domains file
    return as dictionary'''
    profile_lengths = {}
    if domains_file is None:
        if hmm_db in utils.databases:
            domains_file = os.path.join(data_env, "domains.txt")
        else:
            return profile_lengths
    else:
        domains_file = os.path.abspath(domains_file)
    with open(domains_file) as f:
        for line in f:
            line = line.strip().split()
            profile_lengths[line[0]] = float(line[1])
    return profile_lengths


def parse_domains_to_ignore(hmm_db, ignore_file, data_env):
    ''' get a list of profiles to ignore in the analysis,
    return as list'''
    if ignore_file is None:
        if hmm_db in utils.databases:
            ignore_file = os.path.join(data_env, "domains_to_ignore.txt")
        else:
            return []
    else:
        ignore_file = os.path.abspath(self.args["ignore_file"])
    profiles_to_ignore = [line.rstrip() for line in open(ignore_file)]
    return profiles_to_ignore


def create_job_list(results_dir, prep_dir, scan_dir, args):
    # keeping a text file of all the genomes used
    log_other = "###   INPUT  ### \ngenome\thmmer_result\n"
    jobs = []
    for f in os.listdir(scan_dir):
        if f.endswith(".result"):
            basename = os.path.basename(f)
            basename = basename.replace(".result", "")

            scan_summary = copy.copy(args)

            scan_summary["basename"] = basename
            scan_summary["out_file"] = os.path.join(
                results_dir, basename + ".csv")
            scan_summary["results_file"] = os.path.join(
                scan_dir, basename + ".result")
            scan_summary["orf_locs"] = os.path.join(
                prep_dir, basename + ".bed")
            scan_summary["orfs"] = {}
            scan_summary["unfit_orfs"] = {}

            if scan_summary["report_unfit"]:
                scan_summary["report_unfit"] = os.path.join(
                    scan_summary["report_unfit"], basename + ".csv")

            jobs.append(scan_summary)
            log_other = log_other + basename + "\t" + f + "\n"
    utils.write_log(os.path.join(results_dir, "LOG"), "STEP 3 : FILTER HMMER HITS",
     args, log_other)


    return jobs



def parse_hmmer_results(args):
    ''' parse the output of HMMER to find  hits'''
    with open(args["results_file"]) as f:
        prev_name = ""  # save previous name, to not repeat actions
        for line in f:
            if (line.startswith("#")):
                continue
            curr_line = line.split()
            if len(curr_line) < 8:  # run hasn't completed
                continue

            name = curr_line[0]  # orf name
            domain_name = curr_line[3]
            score = float(curr_line[7])  # score of current hit

            # ignore hits with low scores
            if score < args["min_hmmscan_score"] or domain_name in args["profiles_to_ignore"]:
                continue

            # if using hmmsearch, the name is toks[0] and the domain is in curr_line[3]
            # check this by checking if name starts with source
            if not name.lower().startswith("annotation") and not name.lower().startswith("sixframe"):
                name = curr_line[3]
                domain_name = curr_line[0]

            ## get length requirements for current domain
            min_hit_length = args["min_hit_length"]
            max_hit_length = args["max_hit_length"]

            # if domain length given, recalc the required length
            if domain_name in args["profile_lengths"]:
                min_hit_length = max(
                    args["profile_lengths"][domain_name] - args["max_diff_avg_length"], args["min_hit_length"])
                max_hit_length = args["profile_lengths"][domain_name] + \
                    args["max_diff_avg_length"]

            if name != prev_name:  # prevents looking up the name when it appears in multiple rows
                prev_name = name
                curr_orf = args["orf_locs"][args["orf_locs"][3] == name]
                curr_orf = curr_orf.values.tolist()[0]
                sequence = curr_orf[6] # sequence is in the 7th column of BED file, in nucleotide sequence
                length = len(sequence)/3

                # get all the values from the BED file
                source, contig, strand, start, stop = name.split(
                    "|")[0], curr_orf[0], curr_orf[5], curr_orf[1], curr_orf[2]


            if length < min_hit_length or length > max_hit_length:  # doesn't meet length requirements
                ## Add to dictionary of unfit hits
                args["unfit_orfs"][name] = {
                    "name": name,
                    "domain": domain_name,
                    "score": score,
                    "sequence": sequence,
                    "contig": contig,
                    "strand": strand,
                    "start": start,
                    "stop": stop,
                    "length": length,
                    "source": source,
                    "reason": "hit length",
                    "unfit": True
                }
                continue

            if name in args["orfs"]:  # this ORF already has a hit with another domain
                # update only if current score better
                if score > args["orfs"][name]["score"]:
                    args["orfs"][name]["score"] = score
                    args["orfs"][name]["domain"] = domain_name
            else:
                # Add a new ORF
                args["orfs"][name] = {
                    "name": name,
                    "domain": domain_name,
                    "score": score,
                    "sequence": sequence,
                    "contig": contig,
                    "strand": strand,
                    "start": start,
                    "stop": stop,
                    "length": length,
                    "keep": True,
                    "unfit": False,
                    "reason": "",
                    "source": source}
	return


def make_unfit(orf, args, reason):
    ''' change the status of an ORF to "unfit"'''
    if orf["name"] not in args["unfit_orfs"]: ## add to unfit dict
        args["unfit_orfs"][orf["name"]] = orf
    ## update the reason
    if orf["reason"] != "": # not the only reason
        orf["reason"] += "+"
    orf["reason"] = orf["reason"] + reason
    return


def choose_delta(deltas):
    ''' choose the closest candidate antitoxin
    return the chosen row of the dataframe'''
    deltas = list(deltas)
    chosen = deltas[0]
    for d in deltas:
        if abs(d) < abs(chosen):
            chosen = d
    return chosen

def set_upstream_at(orf, start, stop, sequence, delta):
    orf["upstream_at_start"] = start
    orf["upstream_at_stop"] = stop
    orf["upstream_at_sequence"] = sequence
    orf["delta_up"] = delta
    return

def set_downstream_at(orf, start, stop, sequence, delta):
    orf["downstream_at_start"] = start
    orf["downstream_at_stop"] = stop
    orf["downstream_at_sequence"] = sequence
    orf["delta_down"] = delta
    return

# cases in which we have T - AT (downstream for +, upstream for -)
def find_antitoxin(args, orf, orf_locs, orientation):
    if orf["strand"] == "+" and args["order"] not in [orientation, "either", "both"]: ## not required
        return False
    if orf["strand"] == "-" and args["order"] == orientation:
        return False

    ## filter by distance
    if orientation == "downstream":
        ## AT_start <= T_stop + MAX_distance
        filtered = orf_locs[orf_locs[1] <= orf["stop"] + args["max_distance"]] # antitoxin start location
        ## AT_start >= T_stop - MAX_overlap
        filtered = filtered[filtered[1] >= orf["stop"] - args["max_overlap"]] # antitoxin start location
    else: ## upstream
        # AT stop >= T start - MAX_distance
        filtered = orf_locs[orf_locs[2] >= orf["start"] - args["max_distance"]]
        # AT stop <= T start + max overlap
        filtered = filtered[filtered[2] <= orf["start"] + args["max_overlap"]]

    if filtered.shape[0] == 0:
        if orf["strand"] == "+":
            make_unfit(orf, args, "No adjacent " + orientation + " ORF")
        else:
            if orientation == "upstream":
                make_unfit(orf, args, "No adjacent downstream ORF")
            else:
                make_unfit(orf, args, "No adjacent upstream ORF")
        return False

    ## filter by length
    ### Downstream
    if (orf["strand"] == "+" and orientation == "downstream") or (orf["strand"] == "-" and orientation == "upstream"):
        filtered = filtered[filtered[2] - filtered[1] >=  args["min_downstream_length"]*3]
        filtered = filtered[filtered[2] - filtered[1] <=  args["max_downstream_length"]*3]
    else: # - and downstream or + and upstream
        filtered = filtered[filtered[2] - filtered[1] >=  args["min_upstream_length"]*3]
        filtered = filtered[filtered[2] - filtered[1] <=  args["max_upstream_length"]*3]

    if filtered.shape[0] == 0:
        if orf["strand"] == "+":
            make_unfit(orf, args, orientation + " partner length")
        else:
            if orientation == "downstream":
                make_unfit(orf, args, "Upstream partner length")
            else:
                make_unfit(orf, args, "Downstream partner length")
        return False

    if orientation == "downstream":
        filtered = filtered[filtered[1] - orf["stop"] == choose_delta(filtered[1] - orf["stop"])]
    else:
        filtered = filtered[orf["start"] - filtered[2]  == choose_delta(orf["start"] - filtered[2])]
    filtered = filtered.values.tolist()[0]
    start, stop, sequence = filtered[1], filtered[2], filtered[6]

    if orientation == "upstream":
        delta =  orf["start"] - stop
    else:
        delta = start - orf["stop"]
    if (orf["strand"] == "+" and orientation == "downstream") or (orf["strand"] == "-" and orientation == "upstream"):
        set_downstream_at(orf, start, stop, sequence, delta)
    else:
        set_upstream_at(orf, start, stop, sequence, delta)

    return True

def find_candidate_antitoxins(args):
    ''' look for antitoxins upstream/downstream to valid hits'''
    # remove ORFs that have a large number of "XXX" as potentail antitoxins
    args["orf_locs"] = args["orf_locs"][~args["orf_locs"][6].str.contains("XXXXXXXX")]
    args["orf_locs"]= args["orf_locs"][~args["orf_locs"][6].str.contains("NNNNNNNN")]


    for o in args["orfs"]:
        curr_orf = args["orfs"][o]
        orf_locs = args["orf_locs"][args["orf_locs"][0] == curr_orf["contig"]] # filter by contig
        orf_locs = orf_locs[orf_locs[5] == curr_orf["strand"]] # filter by strand

        up = find_antitoxin(args, curr_orf, orf_locs, "upstream") # returns True if found
        down = find_antitoxin(args, curr_orf, orf_locs, "downstream") # returns True if found

        # For 3 component operons, must find up and down
        if args["order"] == "both" and not (up and down):
            curr_orf["keep"] = False
            curr_orf["unfit"] = True
        # couldn't find a candidate antitoxin
        if not up and not down:
            curr_orf["keep"] = False
            curr_orf["unfit"] = True
    return


def check_orfs_overlap(o1, o2):
    ''' return TRUE if two ORFs overlap
    Consider the ORFs to overlap if their start or stop are ~100bp apart.'''
    if o1["strand"] == o2["strand"] and o1["contig"] == o2["contig"] and o1["start"]-100 <= o2["start"] and o2["start"] <= o1["start"] + 100:
        return True
    if o1["strand"] == o2["strand"] and o1["contig"] ==o2["contig"] and o1["stop"]-100 <= o2["stop"] and o2["stop"] <= o1["stop"] + 100:
        return True
    return False

def change_status(keep, lose):
    ''' change the status of ORFs to keep and lose'''
    lose["keep"] = False
    keep["keep"] = True
    if keep["source"] != lose["source"]:
        keep["source"] = "both"
    return

def merge_hits(args):
    ''' because there are multiple ORFs of the same location
    because of
    1. format of 6-frame sixframe_translation
    2. both annotation and sixframe_translation
    Some of the hits need to be merged and referred to as a single ORF.
    Merging is only applied in this step because merging earlier would change
    the structure of the ORFs and would make us lose some hits
    The "keep" value in each ORF determines int the end if it is kept or removed
    from the outputs'''
    for o1_k in args["orfs"]:
        for o2_k in args["orfs"]:
            o1 = args["orfs"][o1_k]
            o2 = args["orfs"][o2_k]
            if o1["name"] == o2["name"]:
                continue # ignore an ORF against itself

            if o1["keep"] == False or o2["keep"] == False: # have already been checked
                continue

            if check_orfs_overlap(o1, o2):
                # if losing an ORF, prefer to keep both
                ## it could be that a partner was found only for one of the two orfs
                ## In that case, keep the ORF that found a partner
                o1_downstream_found = "downstream_at_sequence" in o1.keys()
                o1_upstream_found = "upstream_at_sequence" in o1.keys()
                o2_downstream_found = "downstream_at_sequence" in o2.keys()
                o2_upstream_found = "upstream_at_sequence" in o2.keys()

                if not o1_downstream_found and o2_downstream_found and o2_upstream_found:
                    change_status(o2, o1)
                elif not  o1_upstream_found and o2_downstream_found and o2_upstream_found:
                    change_status(o2, o1)
                elif not o2_upstream_found and o1_downstream_found and o1_upstream_found:
                    change_status(o1, o2)
                elif not o2_upstream_found and o1_downstream_found and o1_upstream_found:
                    change_status(o1, o2)
                elif o1["source"] != o2["source"] and o2["source"] == "annotation": ## keep annotation ORFs compared to sixframe
                    change_status(o2, o1)
                elif o1["source"] != o2["source"] and o1["source"] == "annotation": ## keep annotation ORFs compared to sixframe
                    change_status(o1, o2)
                elif o1["score"] >= o2["score"]: ## if all else fails go by HMMER score
                    change_status(o1, o2)
                else:
                    change_status(o2, o1)
    return

def orf_to_line(orf, sep, order):
    ''' write all the ORF results into a single file '''
    ### either is used on one of the partners was not found
    if order  == "either" and "downstream_at_sequence" not in orf:
        orf["downstream_at_sequence"] = "-"
        orf["downstream_at_start"] = "-"
        orf["downstream_at_stop"] = "-"
        orf["delta_down"] = "-"
    elif order  == "either" and "upstream_at_sequence" not in orf:
        orf["upstream_at_sequence"] = "-"
        orf["upstream_at_start"] = "-"
        orf["upstream_at_stop"] = "-"
        orf["delta_up"] = "-"

    prefix = [orf["domain"], orf["score"], orf["contig"], orf["strand"], orf["length"], orf["start"], orf["stop"]]
    upstream = []
    downstream = []
    if order in ["upstream", "both", "either"]:
        upstream = [len(orf["upstream_at_sequence"]), orf["delta_up"], orf["upstream_at_start"], orf["upstream_at_stop"]]
    if order in ["downstream", "both", "either"]:
        downstream = [len(orf["downstream_at_sequence"]), orf["delta_down"], orf["downstream_at_start"], orf["downstream_at_stop"]]

    line = sep.join(map(str, prefix + upstream + downstream))
    line += sep + orf["sequence"]
    if order in ["upstream", "both", "either"]:
        line += sep + orf["upstream_at_sequence"]
    if order in ["downstream", "both", "either"]:
        line += sep + orf["downstream_at_sequence"]
    line += sep + orf["source"] + "\n"
    return line


def write_results(args):
    ''' write the results files for this strain'''
    out = open(args["out_file"], "w")
    prefix = ["Domain","HMMER_score","Contig","Strand","Hit_length","Hit_start","Hit_stop"]
    upstream = []
    downstream = []
    if args["order"] in ["upstream", "both", "either"]:
        upstream =  ["Upstream_length","Upstream_delta","Upstream_start","Upstream_stop"]
    if args["order"] in ["downstream", "both", "either"]:
        downstream = ["Downstream_length","Downstream_delta","Downstream_start","Downstream_stop"]

    out.write(args["sep"].join(prefix + upstream + downstream + ["Hit,"]))
    if args["order"] in ["upstream", "both", "either"]:
        out.write("Upstream,")
    if args["order"] in ["downstream", "both", "either"]:
        out.write("Downstream,")
    out.write("Source\n")

    for o in args["orfs"]:
        if args["orfs"][o]["keep"]:
            out.write(orf_to_line(args["orfs"][o], args["sep"], args["order"]))
    out.close()
    return


def merge_unfit(args):
    ''' again there are duplicates of ORFs in the unfit hits,
    need to merge them to report only one ORF'''
    for o1_k in args["unfit_orfs"]:
        for o2_k in args["unfit_orfs"]:
            o1 = args["unfit_orfs"][o1_k]
            o2 = args["unfit_orfs"][o2_k]
        if o1["name"] == o2["name"]:
            continue
        if not o1["unfit"] or not o2["unfit"]: ## one is already removed
            continue
        if check_orfs_overlap(o1, o2): ## keep one at random
            o1["unfit"] = False
            o2["unfit"] = True
    return

def unfit_to_line(orf,sep):
    return sep.join(map(str,[orf["domain"],orf["score"],orf["contig"],orf["strand"],orf["length"],orf["start"],orf["stop"],orf["sequence"],orf["reason"]])) + "\n"

def report_unfit(args):
    ''' generate a report with all the significant hits that were
    discarded because they did not meet the structural requirements'''
    if not args["report_unfit"]: ## user didn't request to report unfit
        return
    merge_unfit(args)
    report = open(args["report_unfit"],"w")
    report.write(args["sep"].join(["Domain","HMMER_score","Contig","Strand","Hit_length","Hit_start","Hit_stop","Sequence","Reasons"])+ "\n")

    for o in args["unfit_orfs"]:
        if args["unfit_orfs"][o]["unfit"]:
            report.write(unfit_to_line(args["unfit_orfs"][o], args["sep"]))
    report.close()
    return

def run_summarise(args):
    # read in the ORF locs file
    args["orf_locs"] = pandas.read_csv(args["orf_locs"], header = None, sep = "\t")
    parse_hmmer_results(args)
    find_candidate_antitoxins(args)
    merge_hits(args)
    write_results(args)
    report_unfit(args)
    return


def run(args):
    args = copy.deepcopy(args)
    # define input and output directories
    args.out_dir = os.path.abspath(args.out_dir)

    if "prep_id" not in vars(args):
        prep_id = args.id
    elif args.prep_id is None:
        prep_id = args.id
    else:
        prep_id = args.prep_id

    if "scan_id" not in vars(args):
        scan_id = args.id
    elif args.scan_id is None:
        scan_id = args.id
    else:
        scan_id = args.scan_id

    prep_dir = os.path.join(args.out_dir, prep_id + "_PREPARE")
    scan_dir = os.path.join(args.out_dir, scan_id + "_SCAN")
    results_dir = os.path.join(args.out_dir, args.id + "_FILTER")

    if args.report_unfit:
        args.report_unfit = os.path.join(results_dir, "UNFIT")
        utils.assure_path_exists(args.report_unfit)

    # create the results directory
    utils.assure_path_exists(results_dir)

    # get the SLING data environment
    d = os.path.abspath(os.path.dirname(__file__))
    data_env = os.path.join(d, 'data/')

    # get all FILTERING requirements
    args = vars(args)
    get_requirements(args, data_env)

    # get the profile lengths and profiles to ignore
    profile_lengths = parse_domains(
        args["hmm_db"], args["domains_file"], data_env)
    profiles_to_ignore = parse_domains_to_ignore(
        args["hmm_db"], args["domains_to_ignore"], data_env)

    # covert args to a dictionary
    args["profiles_to_ignore"] = profiles_to_ignore
    args["profile_lengths"] = profile_lengths

    jobs = create_job_list(results_dir, prep_dir, scan_dir, args)
    ### run all the jobs as a pool
    utils.run_pool(jobs, args, run_summarise)
    return
