from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
from Bio.Seq import translate
import numpy as np
import sys
import os
from utils import *
import re
import multiprocessing
import signal
import copy


class Error (Exception):
    pass


def check_start_codons(start_codons):
    ''' check that the start codons from the user are valid'''
    start_codons = start_codons.split(",")
    if len(start_codons) < 1:
        sys.exit(
            "Must provide at least one start codon! Start codons need to be comma delemited.")

    for i in range(0, len(start_codons)):
        s = start_codons[i]
        start_codons[i] = s.lower()
        if len(s) != 3 or not re.match("^[acgt]*$", s):
            sys.exit(
                "Start codons must be of length 3 and contain only the letters a,c,g or t: " + s)
    return start_codons


def check_basename(f, suffix):
    basename = os.path.basename(f)
    basename = basename.replace(suffix, "")
    if len(set('[~!@$%^&*()+{}":;\']+$').intersection(basename)) > 0:
        sys.exit("Avoid using special characters in file name! \nFile: " + basename)
    return basename


def convert_gff_to_fasta(gff_file, fasta_file):
    ''' generate a FASTA file from the GFF file'''
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
    if not fasta:
        sys.exit("Could not find FASTA sequence for GFF file %s!" % gff_file)
    return


def generate_missing_fastas(files, fasta_dir, fasta_suffix):
    ''' every FASTA file that isn't found, generate it from the GFF file'''
    for basename in files:
        if files[basename]["fasta"] == "Not found":
            out_fasta = os.path.abspath(os.path.join(
                fasta_dir, basename + fasta_suffix))
            convert_gff_to_fasta(files[basename]["gff"], out_fasta)
            files[basename]["fasta"] = out_fasta
    return


def get_all_files(gff_dir, fasta_dir, gff_suffix, fasta_suffix):
    ''' read all the files into a dictionary contaning:
    basename -> "gff": location of gff , "fasta" : location of fasta'''
    files = {}

    ''' generate a list of all the GFF files'''
    if gff_dir is None:
        gff_dir = fasta_dir
    else:
        gff_dir = gff_dir

    # get all the GFF files
    get_file_from_dir(files, gff_dir, gff_suffix, "gff")

    if fasta_dir is not None:  # if a FASTA file is given, get all the FASTA files
        get_file_from_dir(files, fasta_dir, fasta_suffix, "fasta")
    else:
        fasta_dir = os.path.join(gff_dir)

    # go over all the files -> if a FASTA file is missing generate it
    generate_missing_fastas(files, fasta_dir, fasta_suffix)

    # tell the user exactly what files it's running on
    print("Running preparation step on files: ")
    for basename in files:
        print("############## %s  ##############\nFASTA: %s\nGFF: %s\n" %
              (basename, files[basename]["fasta"], files[basename]["gff"]))
    return files


def get_file_from_dir(files, directory, suffix, file_type):
    ''' read all the files from a directory into the dictionary'''
    directory = os.path.abspath(directory)
    filenames = os.listdir(directory)
    for f in filenames:
        if f.endswith(suffix):
            basename = check_basename(f, suffix)
            if basename not in files:
                files[basename] = {"fasta": "Not found", "gff": "Not found"}
            files[basename][file_type] = os.path.join(directory, f)
    return


def create_jobs_dict(files, args):
    ''' create a list of all the jobs that need to be call
    by the pool object'''

    # LOG output
    log_other = "###   INPUT ### \ncnt\tgenome\tfasta_file\tgff_file\n"
    jobs = []

    for basename in files:
        fasta_file = files[basename]["fasta"]
        gff_file = files[basename]["gff"]

        # convert the args class to a dictionary
        genome_prep = copy.copy(vars(args))

        # add missing values
        genome_prep["strain"] = basename
        genome_prep["fasta_file"] = fasta_file
        genome_prep["gff_file"] = gff_file
        # arguments required for the preparation step
        genome_prep["orf_id"] = 0
        genome_prep["contigs"] = {}

        jobs.append(genome_prep)
        log_other += "\t".join(map(str, [len(jobs),
                                         basename, fasta_file, gff_file])) + "\n"

    # write the log file
    write_log(os.path.join(args.out_dir, "LOG"),
              "STEP 1 : GENOME PREP", vars(args), log_other)

    return jobs


def run_pool(jobs, args):
    ''' run the pool of prep workers'''
    # create input handler for using CTRL+C
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = multiprocessing.Pool(processes=args.cpu)
    signal.signal(signal.SIGINT, original_sigint_handler)

    try:
        results = pool.map_async(prepare_one_genome, tuple(jobs))
        results.get(120000)
    except KeyboardInterrupt as e:
        pool.terminate()
        sys.exit("Terminated by user")
    except ValueError as e:
        pool.terminate()
        sys.exit(
            "Names of contigs in the GFF file (column 1) must match name of \
            the descriptors in the FASTA files. If FASTA at the end of GFF, \
            simply call SLING only with --gff_dir and GFFs will be converted \
            (see Wiki).")
    else:
        pool.close()
    pool.join()

    return


def prepare_one_genome(args):
    fasta_out = open(os.path.join(
        args["out_dir"], args["strain"] + ".fasta"), "w")
    orf_locs_out = open(os.path.join(
        args["out_dir"], args["strain"] + ".bed"), "w")
    sixframe_translation(args, fasta_out, orf_locs_out)  # prep FASTA files
    annotated_orf_locs(args, fasta_out, orf_locs_out)  # prep GFF files
    fasta_out.close()
    orf_locs_out.close()
    return


def sixframe_translation(args, sixframe_out, orf_locs_out):
    ''' Perform a 6-frame translation on all the FASTA files '''
    with open(args["fasta_file"]) as handle:
        for values in SimpleFastaParser(handle):
            contig = values[0].split()[0]

            genome = values[1].lower()
            rev_genome = reverse_complement(genome)

            # save the contigs for the annotated ORF locs
            args["contigs"][contig] = genome

            # doesn't assume the contigs are circular
            for i in range(0, 3):
                # standard
                get_frame_orfs(args, translate(
                    genome[i:], table=args["codon_table"]), genome[i:], i, contig, sixframe_out, orf_locs_out, genome)
                # complement
                get_frame_orfs(args, translate(
                    rev_genome[i:], table=args["codon_table"]), rev_genome[i:], i + 3, contig, sixframe_out, orf_locs_out, genome)
    return


def get_frame_orfs(args, protein_frame, nuc_frame, i, contig, sixframe_out, orf_locs_out, genome):
    ''' get all open reading frames in current frame '''
    ORFs = protein_frame.split("*")  # split at stop codons

    stops = map(len, ORFs)  # get the stop codon positions
    stops = [x + 1 for x in stops]
    stops = list(np.cumsum(stops))

    starts = [0] + stops[:-1]  # get the start codons positions

    for index in range(0, len(ORFs)):  # iterate over all possible ORFs

        # get nuc sequence of that ORF
        nuc_seq = nuc_frame[starts[index] * 3: stops[index] * 3]
        # premature stop due to end of contig, don't take this sequence
        if not translate(nuc_seq).endswith("*"):
            continue

        check_orf_conditions(
            args, ORFs[index], nuc_seq, i, contig, sixframe_out, orf_locs_out, genome)
    return


def check_orf_conditions(args, protein_seq, nuc_seq, strand, contig, sixframe_out, orf_locs_out, genome):
    ''' given an ORF, find the relevant start codon and write to file if meets requirements '''
    for i in range(0, len(args["start_codons"])):
        all_res = find_start_codon(
            args, args["start_codons"][i], protein_seq, nuc_seq, strand, genome)
        if all_res != None:  # found an ORF with the higher priority codon
            break

    if all_res == None:  # no start codon with the permitted start codons
        return

    for res in all_res:
        orf = res[0]  # the actual ORF
        strand_symbol = res[1]  # strand
        start = res[2]  # start
        stop = res[3]  # stop

        # give the ORF a name that can later be useful
        ORF_name = "Sixframe|Strain:" + args["strain"] + "|ORF:" + str(args["orf_id"]) + "|Contig:" + str(
            contig) + "|Strand:" + strand_symbol + "|Start:" + str(start) + "|Stop:" + str(stop)

        # write the results to the fasta file
        sixframe_out.write(">" + ORF_name + "\n" + orf + "\n")

        # write the results to the ORF locations file
        orf_locs_out.write("\t".join([contig, str(start), str(
            stop), ORF_name, "0", strand_symbol, orf]) + "\n")

        # increase ORF id by one for next ORF
        args["orf_id"] += 1
    return


def find_start_codon(args, codon, protein_seq, nuc_seq, strand, genome):
    ''' given an ORF, look for a start codon in the candidate ORF '''
    # would return V, M or L depending on the start codon
    aa = translate(codon, args["codon_table"])
    aa_index = protein_seq.find(aa)  # find the relevant aa of the start codon

    if aa_index > 0 and len(protein_seq[aa_index:]) >= args["min_orf_length"] and nuc_seq[aa_index * 3:aa_index * 3 + 3] == codon:
        # use M
        protein_seq = protein_seq[aa_index:]  # take the ORF from the start
        strand_symbol = "+"
        nuc_seq = nuc_seq[aa_index * 3:]

        if strand > 2:
            nuc_seq = reverse_complement(nuc_seq)
            strand_symbol = "-"

        starts = [m.start() for m in re.finditer(nuc_seq, genome)]
        res = []

        for start in starts:
            start = genome.find(nuc_seq) + 1
            stop = start + len(nuc_seq) - 1

            # always put an M at start of protein sequence
            protein_seq = "M" + protein_seq[1:]
            res.append([protein_seq, strand_symbol, start, stop])

        return res  # return all values
    # couldn't find ORF
    return None


def annotated_orf_locs(args, fasta_out, orf_locs_out):
    ''' get all the open reading frames from the GFF file, save as BED and as FASTA '''
    if args["gff_file"] == "Not found":
        return

    with open(args["gff_file"]) as f:
        for line in f:
            if line.startswith("##FASTA"):  # end of file
                break

            if line.startswith("#"):  # comments
                continue

            line = line.strip().split("\t")

            if len(line) < 7:  # not full line, ignore
                continue

            if (line[2] == "CDS"):  # only coding sequence
                strand = line[6]
                start = int(line[3])
                stop = int(line[4])
                contig = line[0]

                if contig not in args["contigs"].keys():
                    raise ValueError(
                        "Error: Name of contigs in FASTA file must match contigs in GFF file")

                orf_name = "Annotation|Strain:" + args["strain"] + "|ORF:" + str(args["orf_id"]) + "|Contig:" + str(
                    contig) + "|Strand:" + strand + "|Start:" + str(start) + "|Stop:" + str(stop)

                if strand == "-":
                    sequence = args["contigs"][contig][start:stop]
                    sequence = reverse_complement(sequence)
                else:
                    sequence = args["contigs"][contig][start - 1:stop - 1]

                sequence = translate(sequence, table=args["codon_table"])

                # write the results to the fasta file
                fasta_out.write(">" + orf_name + "\n" + sequence + "\n")

                # write the results to the ORF locations file
                orf_locs_out.write("\t".join([contig, str(start), str(
                    stop), orf_name, "0", strand, sequence]) + "\n")

                args["orf_id"] += 1
    return


def run(args):

    if args.gff_dir is None and args.fasta_dir is None:
        sys.exit("Error: please provide either <gff_dir> or <fasta_dir>")

    # check that the start codons are valid
    starts_codons = check_start_codons(args.start_codons)
    # input and output dir
    args.out_dir = os.path.join(os.path.abspath(
        args.out_dir), args.prep_id + "_PREPARE")
    assure_path_exists(args.out_dir)
    files = get_all_files(args.gff_dir, args.fasta_dir,
                          args.gff_suffix, args.fasta_suffix)
    jobs = create_jobs_dict(files, args)
    run_pool(jobs, args)

    return
