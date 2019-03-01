import argparse
import sling

def run():
    parser = argparse.ArgumentParser(
        description = 'Group the operons based on sequence similarity',
        usage = 'sling group [options] <id> <hmm_db>')

    parser.add_argument('--filter_id',type=str,  help="ID of filter run [default: same as <id>]", metavar='STR')
        
    parser.add_argument('-it','--itol', action='store_true', help='Generate files that can be loaded into ITOL [%(default)s]', default=False)
    parser.add_argument('-t','--order', type=str, help = 'Location of partner gene relative to hit if overriding <req_file>. Options: upstream, downstream, either, both [either]', metavar='STR', default=None)
    parser.add_argument('-u','--report_unfit', action='store_true', help='Generate outputs for HMMER hits that did not meet requirements  [%(default)s]', default=False)

    parser.add_argument('-mbe','--min_blast_evalue', type=float, help='Minimum BLAST evalue to use for an edge in the sequence similarity network [%(default)s]', metavar='FLOAT', default=0.01)
    parser.add_argument('-mi','--min_identity', type=int, help='Minimum BLAST identity to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=75)
    parser.add_argument('-lc','--length_coverage', type=float, help='Minimum alignment coverage when comparing two sequences [%(default)s]', metavar='FLOAT', default=0.75)

    parser.add_argument('-c','--cpu', type=int, help='Number of CPUs to use [%(default)s]', default=8, metavar='INT')

    parser.add_argument('-o','--out_dir', type=str, help='Directory for all the output files', metavar="PATH",default=".")
    parser.add_argument('-s','--sep', type=str, help='Seperator for the input and output files, [%(default)s]', default=",", metavar='STR')

    parser.add_argument('--makeblastdb', type=str, help='makeblastdb executable [%(default)s]', default="makeblastdb", metavar='STR')
    parser.add_argument('--blastp', type=str, help='blastp executable [%(default)s]', default="blastp", metavar='STR')

    parser.add_argument('id',type=str,  help="ID of group run [default: same as <id>]", metavar='STR')
    parser.add_argument('hmm_db', help='Name of the predefined HMM database ' + str(sling.utils.databases) + ' OR path to custom HMM file', metavar='STR/FILE')

    options = parser.parse_args()
    sling.group.run(options)
    return
