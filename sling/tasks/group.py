import argparse
import sling
import sys
import os

def run():
    parser = argparse.ArgumentParser(
        description = 'Group the operons based on sequence similarity',
        usage = 'sling group [options] <filter_id> <group_id> <hmm_db>')
    parser.add_argument('-it','--save_to_ITOL', action='store_true', help='Generate files that can be loaded into ITOL [%(default)s]', default=False)
    parser.add_argument('-t','--order', type=str, help = 'Location of partner gene relative to hit if overriding <req_file>. Options: upstream, downstream, either, both [either]', metavar='STR', default=None)
    parser.add_argument('-u','--report_unfit', action='store_true', help='Generate outputs for HMMER hits that did not meet requirements  [%(default)s]', default=False)
    parser.add_argument('-mbe','--min_blast_evalue', type=float, help='Minimum BLAST evalue to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=0.01)
    parser.add_argument('-mi','--min_identity', type=int, help='Minimum BLAST identity to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=30)
    parser.add_argument('-o','--out_dir', type=str, help='Directory for all the output files', metavar="PATH",default=".")
    parser.add_argument('-s','--sep', type=str, help='Seperator for the input and output files, [%(default)s]', default=",", metavar='STR')
    parser.add_argument('filter_id',type=str,  help="ID of filter run", metavar='STR')
    parser.add_argument('group_id',type=str,  help="ID of group run", metavar='STR')
    parser.add_argument('hmm_db', help='Name of the predefined HMM database ' + str(sling.utils.databases) + ' OR path to custom HMM file', metavar='STR/FILE')

    options = parser.parse_args()

    group = sling.group.Group(options.filter_id,
        options.group_id,
        options.hmm_db,
        out_dir = options.out_dir,
        order = options.order,   
        min_identity = options.min_identity,
        min_blast_evalue = options.min_blast_evalue,
        save_to_ITOL = options.save_to_ITOL,
        sep = options.sep,
        report_unfit = options.report_unfit)
    group.run()