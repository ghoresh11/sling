import argparse
import sling
import os



def run():
    
    parser = argparse.ArgumentParser(
        description = 'Run hmmscan to search for hits in the genome',
        usage = 'sling scan [options] <prep_id> <scan_id> <hmm_db>')
    parser.add_argument('-c','--cpu', type=int, help='Number of CPUs to use [%(default)s]', default=2, metavar='INT')
    parser.add_argument('-o','--out_dir', help='Working for all the output files', metavar="PATH",default=".")
    parser.add_argument('prep_id', help='ID of prepare run', metavar='STR')
    parser.add_argument('scan_id', help='ID of scan run', metavar='STR')
    parser.add_argument('hmm_db', help='Name of the predefined HMM database ' + str(sling.utils.databases) +  ' OR path to custom HMM file', metavar='STR/FILE')
    options = parser.parse_args()
    
    scan = sling.scan.Scan(
            options.prep_id,
            options.scan_id,
            options.hmm_db,
            out_dir = options.out_dir,
            cpu = options.cpu,     
        )
    scan.run()