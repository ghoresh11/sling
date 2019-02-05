import os
import argparse
import sys

class Error (Exception): pass


def _get_strain_dict(f):
    ''' read a gff file into a strain dict
    return: the dictionary for the strain with contig -> strand -> start -> stop -> CDS'''
    strain_dict = {}
    with open(f) as f_open:
        for line in f_open:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")

            if len(toks) < 7:
                continue

            if toks[2] != "CDS":
                continue

            contig = toks[0]
            if contig not in strain_dict:
                strain_dict[contig] = {}

            strand = toks[6]
            if strand not in strain_dict[contig]:
                strain_dict[contig][strand] = {}

            start = toks[3]
            if start not in strain_dict[contig][strand]:
                strain_dict[contig][strand][start] = {}

            stop = toks[4]
            strain_dict[contig][strand][start][stop] = toks[-1]
    return strain_dict


def read_gffs(gff_dir, gff_suffix):
    ''' read gff files into a dictionary of
    strain -> contig -> strand -> start -> stop -> CDS
    return: the dictionary for all the strains '''

    gff_files = os.listdir(os.path.abspath(gff_dir))

    strains = {}
    for f in gff_files:
        if not f.endswith(gff_suffix):
            continue
        strain = os.path.basename(f)
        strain = strain.replace(gff_suffix,"")
        print("Reading in strain %s" %strain)
        strains[strain] = _get_strain_dict(os.path.join(os.path.abspath(gff_dir), f))
    return strains

def _match_hit_to_cds(clusters_dir, f, strains, out_dir):
    ''' go over a hit file line by line
    If the last column says "both" or "annotation", find the relevant CDS identifier
    in the strains dict.
    Output: new files with first column as the CDS identifier
    Return: nothing
    '''
    out = open(os.path.join(out_dir,f), "w")
    is_filter = False
    with open(os.path.join(clusters_dir,f)) as f_open:
        for line in f_open:
            if line.startswith("#"):
                continue
            if line.startswith("Domain"):
                is_filter = True

            if line.startswith("Strain") or line.startswith("Domain"):
                toks = line.lower().strip().split(",")
                out.write("CDS_identifier," + line)
                contig_index = toks.index("contig")
                start_index = toks.index("hit_start")
                stop_index = toks.index("hit_stop")
                strand_index = toks.index("strand")
                source_index = toks.index("source")
                continue

            toks = line.strip().split(",")
            identifier = "None"

            strain = toks[0]
            if is_filter:
                strain = f.replace(".csv","")

            if toks[source_index] in ["annotation", "both"]:
                identifier = strains[strain][toks[contig_index]][toks[strand_index]][toks[start_index]][toks[stop_index]]

            out.write(identifier + "," + line)
    out.close()
    return

def match_to_hits(strains, clusters_dir):
    ''' Go over all the hit outputs and find their CDS identifier from the GFF files'''
    clusters_dir = os.path.abspath(clusters_dir)
    out_dir =  clusters_dir + "_w_cds"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    hit_files = os.listdir(clusters_dir)
    for f in hit_files:
        if not f.endswith(".txt") and not f.endswith(".csv"):
            continue
        print("Finding identifiers for %s..." %f)
        _match_hit_to_cds(clusters_dir, f, strains, out_dir)
    return


def run(args):
    ''' take the outputs from GROUP, and match
    annotation CDSs to the GFF file'''
    strains = read_gffs(args["gff_dir"], args["gff_suffix"])
    match_to_hits(strains, os.path.join(args["group_dir"], "hits_clusters"))
    match_to_hits(strains, os.path.join(args["group_dir"], "complete_clusters"))

    if args["filter_dir"] != "None":
        match_to_hits(strains, args["filter_dir"])
    return

def init():
    ''' init getting all the arguments from the user'''
    parser = argparse.ArgumentParser(
        description='Find the gff CDS identifier for all the hits in the GROUP step',
        usage='python cds_identifier [options] <gff_dir> <group_dir>')
    parser.add_argument('--gff_suffix', type=str,
                        help='Suffix used for the GFF files [%(default)s]', default=".gff", metavar='STR')
    parser.add_argument('--filter_dir', type=str,
                        help='Directory containing all outputs of sling GROUP [%(default)s]', default="None", metavar='DIR')
    parser.add_argument(
        'gff_dir', type=str, help="Directory containing all the GFF files", metavar='DIR')
    parser.add_argument(
        'group_dir', type=str, help="Directory containing all outputs of sling GROUP", metavar='DIR')
    args = vars(parser.parse_args())
    run(args)
    return

if __name__ == "__main__":
    init()
