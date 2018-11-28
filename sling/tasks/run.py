### run the full search on a multiple genomes

import argparse
import sling

def run():
    parser = argparse.ArgumentParser(
        description = 'Run the full search strategy',
        usage = 'sling run [options] <run_id> <input_dir> <hmm_db> ')
    ## PREPARE
    parser.add_argument('-fs','--fasta_suffix', type=str, help='Suffix of FASTA files in <input_dir> [%(default)s]', metavar='STR', default=".fasta")
    parser.add_argument('-gd','--gff_dir', help='Name of directory containing GFF files. [<input_dir>]', metavar='PATH', default="")
    parser.add_argument('-gs','--gff_suffix', type=str, help='Suffix of GFF files in <gff_dir> [%(default)s] ', metavar='STR', default=".gff")
    parser.add_argument('-mol','--min_orf_length', type=int, help='Minimun length of an open reading frame [%(default)s]', metavar='INT', default=20)
    parser.add_argument('-c','--cpu', type=int, help='Number of CPUs to be used [%(default)s]', default = 1, metavar="INT")
    parser.add_argument('-ct','--codon_table', type=str, help='Codon table to use in translation [%(default)s]', default = "Standard", metavar="STR")
    parser.add_argument('-sc','--start_codons', type=str, help='Accepted start codons written in hierarchical order of usage [%(default)s]', default="atg,gtg,ttg", metavar='STR')
    parser.add_argument('-o','--out_dir', help='Directory for all the output files', metavar="PATH",default=".")
    
    ## SCAN
    parser.add_argument('--hmmsearch', help='HMM search executable (set to hmmscan if wish to run scan not search) [Default: hmmsearch]', metavar="STR",default="hmmsearch")
    parser.add_argument('--hmmpress', help='HMM press executable [relevant for hmmscan only] [Default: hmmpress]', metavar="STR",default="hmmpress")

    ## FILTER
    parser.add_argument('-u','--report_unfit', action='store_true', help='Generate reports for HMMER hits that did not meet requirements  [%(default)s]', default=False)
    parser.add_argument('-mhl','--min_hit_length', type=int, help='Minimum length of a hit, if not in DOMAINS file [1]', metavar='INT', default=None)
    parser.add_argument('-Mhl','--max_hit_length', type=int, help='Maximum length of a hit, if not in DOMAINS file [10000000]', metavar='INT', default=None)
    parser.add_argument('-mul','--min_upstream_length', type=int, help='Minimum length of the upstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mul','--max_upstream_length', type=int, help='Maximum length of the upstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-mdl','--min_downstream_length', type=int, help='Minimum length of the downstream gene [1]', metavar='INT', default=None)
    parser.add_argument('-Mdl','--max_downstream_length', type=int, help='Maximum length of the downstream gene [10000000]', metavar='INT', default=None)
    parser.add_argument('-Mo','--max_overlap', type=int, help='Maximum overlap between two operon proteins [300]', metavar='INT', default=None)
    parser.add_argument('-Md','--max_distance', type=int, help='Maximum distance between two opern proteins [10000000]', metavar='INT', default=None)
    parser.add_argument('-mhs','--min_hmmscan_score', type=float, help='Minimum HMMER score to use for significant hits [%(default)s]', default=20, metavar='FLOAT')
    parser.add_argument('-t','--order', type=str, help='Location of partner gene relative to hit. Options: upstream, downstream, either, both [either]', metavar='STR', default=None)  
    parser.add_argument('-s','--sep', type=str, help='Delimiter to use in the output file [%(default)s]', metavar='PATH', default=",")
    parser.add_argument('-di','--domains_to_ignore', type=str, help='File with line delemited hmmer domains to ignore in summary', default = "", metavar="FILE")
    parser.add_argument('-df','--domains_file', type=str, help='Tab delimited file of HMMER domains and the expected length of their hits', metavar='FILE', default="")
    parser.add_argument('-Mda','--max_diff_avg_length', type=int, help='Maximum difference between hit length and its average length defined in <domains_file> [10000000]', metavar='INT', default=None)

    ## GROUP
    parser.add_argument('-it','--save_to_ITOL', action='store_true', help='Generate files that can be loaded into ITOL [%(default)s]', default=False)
    parser.add_argument('-mbe','--min_blast_evalue', type=float, help='Minimum BLAST evalue to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=0.01)
    parser.add_argument('-mi','--min_identity', type=int, help='Minimum BLAST identity to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=30)
    parser.add_argument('--makeblastdb', type=str, help='makeblastdb executable [%(default)s]', default="makeblastdb", metavar='STR')
    parser.add_argument('--blastp', type=str, help='blastp executable [%(default)s]', default="blastp", metavar='STR')

    ## Required
    parser.add_argument('run_id', type=str, help='Directory to save all result files', metavar='PATH')
    parser.add_argument('input_dir',  help='Name of directory containing FASTA files', metavar='PATH')
    parser.add_argument('hmm_db', help='Name of the predefined HMM database ' + str(sling.utils.databases) + ' OR path to custom HMM file', metavar='STR/FILE')

    
    options = parser.parse_args()

    prepare = sling.prepare.Prepare(options.input_dir, 
        options.run_id,
        out_dir = options.out_dir,
        fasta_suffix = options.fasta_suffix,
        gff_suffix = options.gff_suffix,
        gff_dir = options.gff_dir,
        min_orf_length = options.min_orf_length,
        start_codons = options.start_codons,
        codon_table = options.codon_table,
        cpu = options.cpu)

    prepare.run()

    scan = sling.scan.Scan(
        options.run_id,
        options.run_id,
        options.hmm_db,
        hmmsearch = options.hmmsearch,
        hmmpress = options.hmmpress,
        out_dir = options.out_dir,
        cpu = options.cpu,     
    )
    scan.run()

    summarise = sling.filter.Summarise(options.run_id,
        options.run_id,
        options.run_id,
        options.hmm_db,
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

    group = sling.group.Group(options.run_id,
        options.run_id,
        options.hmm_db,
        out_dir = options.out_dir,
        order = options.order,   
        min_identity = options.min_identity,
        min_blast_evalue = options.min_blast_evalue,
        save_to_ITOL = options.save_to_ITOL,
        sep = options.sep,
        report_unfit = options.report_unfit,
        cpu = options.cpu,
        makeblastdb = options.makeblastdb,
        blastp = options.blastp)
    group.run()




