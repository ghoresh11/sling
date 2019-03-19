import argparse
import os


def gff_to_fasta(gff_file, fasta_file):
    out = open(fasta_file, "w")
    contigs = {}
    with open(gff_file) as f:
        fasta = False
        for line in f:
            if fasta:
                out.write(line)
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
    out.close()
    return


def run(args):
    gff_files = os.listdir(args["gff_dir"])
    if args["fasta_dir"] is None:
        args["fasta_dir"] = args["gff_dir"]
    for f in gff_files: 
        if not f.endswith(args["gff_suffix"]):
            continue
        outfile = os.path.abspath(os.path.join(args["fasta_dir"], f.replace(args["gff_suffix"],".fasta")))
        gff_to_fasta(os.path.abspath(os.path.join( args["gff_dir"],f)), outfile)
    return

def init():
    ''' init getting all the arguments from the user'''
    parser = argparse.ArgumentParser(
        description='Convert GFF files to FASTA files with same identifier (for SLING input)',
        usage='python gff_to_fasta.py [options] <gff_dir> <gff_suffix>')
    parser.add_argument('--gff_dir', type=str, required = True,
                        help='Directory containing all the GFF files [%(default)s]', default="None", metavar='DIR')
    parser.add_argument('--gff_suffix', type=str,
                        help='Suffix used for the GFF files [%(default)s]', default=".gff", metavar='STR')
    parser.add_argument('--fasta_dir', type = str, 
        help="Directory for output FASTA files", metavar='DIR')
    args = vars(parser.parse_args())
    run(args)
    return

if __name__ == "__main__":
    init()
