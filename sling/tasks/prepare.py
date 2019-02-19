import argparse
import sling

def run():
    parser = argparse.ArgumentParser(
        description = 'Prepare genomes for pipeline',
        usage = 'sling prepare [options] <id> --fasta_dir <fasta_dir> AND/OR --gff_dir <gff_dir>')
    parser.add_argument('-fs','--fasta_suffix', type=str, help='Suffix of FASTA files in <input_dir> [%(default)s]', metavar='STR', default=".fasta")
    parser.add_argument('-gd','--gff_dir', help='Name of directory containing GFF files. [<input_dir>]', metavar='PATH', default=None)
    parser.add_argument('-fd', '--fasta_dir', help='Name of directory containing FASTA files  [%(default)s]', metavar='PATH', default=None)
    parser.add_argument('-gs','--gff_suffix', type=str, help='Suffix of GFF files in <gff_dir> [%(default)s] ', metavar='STR', default=".gff")
    parser.add_argument('-mol','--min_orf_length', type=int, help='Minimun length of an open reading frame [%(default)s]', metavar='INT', default=20)
    parser.add_argument('-c','--cpu', type=int, help='Number of CPUs to be used [%(default)s]', default = 1, metavar="INT")
    parser.add_argument('-ct','--codon_table', type=str, help='Codon table to use in translation [%(default)s]', default = "Standard", metavar="STR")
    parser.add_argument('-sc','--start_codons', type=str, help='Accepted start codons written in hierarchical order of usage [%(default)s]', default="atg,gtg,ttg", metavar='STR')
    parser.add_argument('-o','--out_dir', help='Directory for all the output files', metavar="PATH",default=".")
    parser.add_argument('id',help = 'ID of prepare run, to be used by downstream tasks',metavar="STR")

    options = parser.parse_args()
    sling.prepare.run(options)
