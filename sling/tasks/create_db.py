import argparse
import sling
import os

def run():
    
    parser = argparse.ArgumentParser(
        description = 'Create a built-in database for SLING',
        usage = 'sling create_db [options] <name> <hmm_file>')
    parser.add_argument('-mhl','--min_hit_length', type=int, help='Minimum length of a hit, if not in DOMAINS file [1]', metavar='INT', default=None)
    parser.add_argument('-Mhl','--max_hit_length', type=int, help='Maximum length of a hit, if not in DOMAINS file [10000000]', metavar='INT', default=None)
    parser.add_argument('-mul','--min_upstream_length', type=int, help='Minimum length of the upstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mul','--max_upstream_length', type=int, help='Maximum length of the upstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-mdl','--min_downstream_length', type=int, help='Minimum length of the downstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mdl','--max_downstream_length', type=int, help='Maximum length of the downstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-Mo','--max_overlap', type=int, help='Maximum overlap between two operon proteins [300]', metavar='INT', default=None)
    parser.add_argument('-Md','--max_distance', type=int, help='Maximum distance between two opern proteins [10000000]', metavar='INT', default=None)
    parser.add_argument('-Mda','--max_diff_avg_length', type=int, help='Maximum difference between hit length and its average length defined in <domains_file> [10000000]', metavar='INT', default=None)
    parser.add_argument('-t','--order', type=str, help='Location of partner gene relative to hit. Options: upstream, downstream, either, both [either]', metavar='STR', default=None)  
    parser.add_argument('name',type=str, help='Name of the predefined HMM database', metavar='STR')
    parser.add_argument('hmm_file',type=str, help='Path to custom HMM file', metavar='FILE')
    parser.add_argument('sling_dir',type=str, help='Path to SLING download', metavar='FILE')

    

    options = parser.parse_args()
    
    create_db = sling.create_db.CreateDB(options.name,
        options.hmm_file,
        options.sling_dir,
        order = options.order,
        max_diff_avg_length = options.max_diff_avg_length,
        min_hit_length = options.min_hit_length,
        max_hit_length = options.max_hit_length,
        min_downstream_length = options.min_downstream_length,
        max_downstream_length = options.max_downstream_length,
        min_upstream_length = options.min_upstream_length,
        max_upstream_length = options.max_upstream_length,
        max_distance = options.max_distance,
        max_overlap = options.max_overlap)
    create_db.run()