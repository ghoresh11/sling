import argparse
import sling

def run():

    parser = argparse.ArgumentParser(
        description = 'Summarise HMMER results and search for operons in each strain',
        usage = 'sling filter [options] <id> <hmm_db>')

    parser.add_argument('--prep_id',type=str, help="ID of prepare run [Default: same as --id]", metavar='STR', default = None)
    parser.add_argument('--scan_id',type=str,  help="ID of scan run [Default: same as --id]", metavar='STR', default = None)

    parser.add_argument('-c','--cpu', type=int, help='Number of CPUs to be used [%(default)s]', default = 8, metavar="INT")

    parser.add_argument('-o','--out_dir', type=str, help='Path to all output files [.]', metavar='PATH',default=".")
    parser.add_argument('-u','--report_unfit', action='store_true', help='Generate reports for HMMER hits that did not meet requirements  [%(default)s]', default=False)
    parser.add_argument('-s','--sep', type=str, help='Delimiter to use in the output file [%(default)s]', metavar='PATH', default=",")

    parser.add_argument('-mhl','--min_hit_length', type=int, help='Minimum length of a hit (if not in domains_file) [1]', metavar='INT', default=None)
    parser.add_argument('-Mhl','--max_hit_length', type=int, help='Maximum length of a hit (if not in domains_file) [10000000]', metavar='INT', default=None)
    parser.add_argument('-mul','--min_upstream_length', type=int, help='Minimum length of the upstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mul','--max_upstream_length', type=int, help='Maximum length of the upstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-mdl','--min_downstream_length', type=int, help='Minimum length of the downstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mdl','--max_downstream_length', type=int, help='Maximum length of the downstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-Mo','--max_overlap', type=int, help='Maximum overlap between two operon proteins [50]', metavar='INT', default=None)
    parser.add_argument('-Md','--max_distance', type=int, help='Maximum distance between two operon proteins [10000000]', metavar='INT', default=None)
    parser.add_argument('-mhs','--min_hmmscan_score', type=float, help='Minimum HMMER score to use for significant hits [%(default)s]', default=20, metavar='FLOAT')
    parser.add_argument('-t','--order', type=str, help='Location of partner gene relative to hit. Options: upstream, downstream, either, both [either]', metavar='STR', default=None)


    parser.add_argument('-di','--domains_to_ignore', type=str, help='File with line delemited hmmer domains to ignore in summary [%(default)s]', default = None , metavar="FILE")
    parser.add_argument('-df','--domains_file', type=str, help='Tab delimited file of HMMER domains and the expected length of their hits [%(default)s]', metavar='FILE', default=None)
    parser.add_argument('-Mda','--max_diff_avg_length', type=int, help='Maximum difference between hit length and its average length defined in <domains_file> [10000000]', metavar='INT', default=None)

    parser.add_argument('id', type=str,  help="ID of this run", metavar='STR')
    parser.add_argument('hmm_db',type=str, help='Name of the predefined HMM database ' + str(sling.utils.databases) + ' OR path to custom HMM file', metavar='STR/FILE')

    options = parser.parse_args()
    sling.filter.run(options)
    return


    return
