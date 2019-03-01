### run the full search on a multiple genomes

import argparse
import sling

def run():
    parser = argparse.ArgumentParser(
        description = 'Run the full search strategy',
        usage = 'sling run [options] <id> <input_dir> <hmm_db> ')

    ## General
    parser.add_argument('-c','--cpu', type=int, help='Number of CPUs to be used [%(default)s]', default = 8, metavar="INT")
    parser.add_argument('-o','--out_dir', help='Directory for all the output files', metavar="PATH",default=".")
    parser.add_argument('-u','--report_unfit', action='store_true', help='Generate reports for HMMER hits that did not meet requirements  [%(default)s]', default=False)
    parser.add_argument('-s','--sep', type=str, help='Delimiter to use in the output file [%(default)s]', metavar='PATH', default=",")
    parser.add_argument('-it','--itol', action='store_true', help='Generate files that can be loaded into ITOL [%(default)s]', default=False)

    ## PREPARE
    parser.add_argument('-gd','--gff_dir', help='Name of directory containing GFF files, if different from <input_dir>', metavar='PATH', default=None)
    parser.add_argument('-fs','--fasta_suffix', type=str, help='Suffix of FASTA files in <input_dir> [%(default)s]', metavar='STR', default=".fasta")
    parser.add_argument('-gs','--gff_suffix', type=str, help='Suffix of GFF files in <gff_dir> [%(default)s] ', metavar='STR', default=".gff")
    parser.add_argument('-mol','--min_orf_length', type=int, help='Minimun length of an open reading frame [%(default)s]', metavar='INT', default=20)
    parser.add_argument('-ct','--codon_table', type=str, help='Codon table to use in translation [%(default)s]', default = "Standard", metavar="STR")
    parser.add_argument('-sc','--start_codons', type=str, help='Accepted start codons written in hierarchical order of usage [%(default)s]', default="atg,gtg,ttg", metavar='STR')

    ## SCAN
    parser.add_argument('--hmmsearch', help='HMM search executable (set to hmmscan if wish to run scan not search) [Default: hmmsearch]', metavar="STR",default="hmmsearch")
    parser.add_argument('--hmmpress', help='HMM press executable [relevant for hmmscan only] [Default: hmmpress]', metavar="STR",default="hmmpress")

    ## FILTER
    parser.add_argument('-t','--order', type=str, help='Location of partner gene relative to hit. Options: upstream, downstream, either, both [either]', metavar='STR', default=None)
    parser.add_argument('-mhl','--min_hit_length', type=int, help='Minimum length of a hit, if not in DOMAINS file [1]', metavar='INT', default=None)
    parser.add_argument('-Mhl','--max_hit_length', type=int, help='Maximum length of a hit, if not in DOMAINS file [10000000]', metavar='INT', default=None)
    parser.add_argument('-mul','--min_upstream_length', type=int, help='Minimum length of the upstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mul','--max_upstream_length', type=int, help='Maximum length of the upstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-mdl','--min_downstream_length', type=int, help='Minimum length of the downstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mdl','--max_downstream_length', type=int, help='Maximum length of the downstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-Mo','--max_overlap', type=int, help='Maximum overlap between two operon proteins [300]', metavar='INT', default=None)
    parser.add_argument('-Md','--max_distance', type=int, help='Maximum distance between two opern proteins [10000000]', metavar='INT', default=None)
    parser.add_argument('-mhs','--min_hmmscan_score', type=float, help='Minimum HMMER score to use for significant hits [%(default)s]', default=20, metavar='FLOAT')
    parser.add_argument('-di','--domains_to_ignore', type=str, help='File with line delemited hmmer domains to ignore in summary', default = None, metavar="FILE")
    parser.add_argument('-df','--domains_file', type=str, help='Tab delimited file of HMMER domains and the expected length of their hits', metavar='FILE', default = None)
    parser.add_argument('-Mda','--max_diff_avg_length', type=int, help='Maximum difference between hit length and its average length defined in <domains_file> [10000000]', metavar='INT', default=None)

    ## GROUP
    parser.add_argument('-mbe','--min_blast_evalue', type=float, help='Minimum BLAST evalue to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=0.01)
    parser.add_argument('-mi','--min_identity', type=int, help='Minimum BLAST identity to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=30)
    parser.add_argument('-lc','--length_coverage', type=float, help='Minimum alignment coverage when comparing two sequences [%(default)s]', metavar='FLOAT', default=0.75)
    parser.add_argument('--makeblastdb', type=str, help='makeblastdb executable [%(default)s]', default="makeblastdb", metavar='STR')
    parser.add_argument('--blastp', type=str, help='blastp executable [%(default)s]', default="blastp", metavar='STR')

    ## Required
    parser.add_argument('id', type=str, help='ID for this SLING run', metavar='PATH')
    parser.add_argument('input_dir',  help='Name of directory containing GFF and/or FASTA files', metavar='PATH')
    parser.add_argument('hmm_db', help='Name of the predefined HMM database ' + str(sling.utils.databases) + ' OR path to custom HMM file', metavar='STR/FILE')


    options = parser.parse_args()
    sling.prepare.run(options)
    sling.scan.run(options)
    sling.filter.run(options)
    sling.group.run(options)
    return
