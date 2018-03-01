import argparse
import sling
import os

def run():
    
    parser = argparse.ArgumentParser(
        description = 'Summarise HMMER results and search for operons in each strain',
        usage = 'sling filter [options] <prep_id> <scan_id> <filter_id> <req>')
    parser.add_argument('-di','--domains_to_ignore', type=str, help='File with line delemited hmmer domains to ignore in summary', default = "", metavar="FILE")
    parser.add_argument('-mhs','--min_hmmscan_score', type=float, help='Minimum HMMER score to use for significant hits [%(default)s]', default=20, metavar='FLOAT')
    parser.add_argument('-Mda','--max_diff_avg_length', type=int, help='Maximum difference between hit length and its average length defined in <domains_file> [%(default)s]', metavar='INT', default=None)
    parser.add_argument('-mhl','--min_hit_length', type=int, help='Minimum length of a hit, if not in DOMAINS file', metavar='INT', default=None)
    parser.add_argument('-Mhl','--max_hit_length', type=int, help='Maximum length of a hit, if not in DOMAINS file', metavar='INT', default=None)
    parser.add_argument('-mul','--min_upstream_length', type=int, help='Minimum length of the upstream gene [%(default)s]', metavar='INT', default=None)
    parser.add_argument('-Mul','--max_upstream_length', type=int, help='Maximum length of the upstream gene [%(default)s]', metavar='INT', default=None)
    parser.add_argument('-mdl','--min_downstream_length', type=int, help='Minimum length of the downstream gene [%(default)s]', metavar='INT', default=None)
    parser.add_argument('-Mdl','--max_downstream_length', type=int, help='Maximum length of the downstream gene [%(default)s]', metavar='INT', default=None)
    parser.add_argument('-Mo','--max_overlap', type=int, help='Maximum overlap between two operon proteins [%(default)s]', metavar='INT', default=None)
    parser.add_argument('-Md','--max_distance', type=int, help='Maximum distance between two opern proteins [%(default)s]', metavar='INT', default=None)
    parser.add_argument('-df','--domains_file', type=str, help='Tab delimited file of HMMER domains and the expected length of their hits', metavar='FILE', default="")
    parser.add_argument('-s','--sep', type=str, help='Delimiter to use in the output file [%(default)s]', metavar='PATH', default=",")
    parser.add_argument('-t','--order', type=str, help='Location of partner gene relative to hit. Options: upstream, downstream, either, both [%(default)s]', metavar='STR', default=None)  
    parser.add_argument('-u','--report_unfit', action='store_true', help='Generate reports for HMMER hits that did not meet requirements  [%(default)s]', default=False)
    parser.add_argument('-o','--out_dir', type=str, help='Path to all output files', metavar='PATH',default=".")
    parser.add_argument('prep_id',type=str, help="ID of prepare run", metavar='STR')
    parser.add_argument('scan_id',type=str,  help="ID of scan run", metavar='STR')
    parser.add_argument('filter_id',type=str,  help="ID of filter run", metavar='STR')
    parser.add_argument('req',type=str, help="The name of one of the predefined databases OR a path to a requirement file for custom DBs", metavar='STR/PATH')

    options = parser.parse_args()
    
    summarise = sling.filter.Summarise(options.prep_id,
        options.scan_id,
        options.filter_id,
        options.req,
        out_dir = options.out_dir,
        sep = options.sep,
        min_hmmscan_score = options.min_hmmscan_score,
        order = options.order,
        max_diff_avg_length = options.max_diff_avg_length,
        min_hit_length = options.min_hit_length,
        max_hit_length = options.max_hit_length,
        min_downstream_length = options.min_downstream_length,
        max_downstream_length = options.max_downstream_length,
        min_upstream_length = options.min_upstream_length,
        max_upstream_length = options.max_upstream_length,
        max_distance = options.max_distance,
        max_overlap = options.max_overlap,
        domains_file = options.domains_file,
        ignore_file = options.domains_to_ignore,
        report_unfit = options.report_unfit)
    summarise.run()