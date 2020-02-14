# sling
A tool to search for linked gene arrays in bacterial datasets

For details on how to use SLING, please see the [SLING wiki page](https://github.com/ghoresh11/sling/wiki).

To cite SLING: 

Horesh G, Harms A, Fino C, Parts L, Gerdes K, Heinz E, et al. SLING: a tool to search for linked genes in bacterial datasets. Nucleic Acids Res. 2018; [doi:10.1093/nar/gky738](https://doi.org/10.1093/nar/gky738)

### February 14, 2020:  Version 2.0.1

Added option to provide gzipped fasta files.

### March 6, 2019:  Version 2.0

Major changes:

1. Input: GFF format can be provided alone if the FASTA sequence is at the end of the file (for instance, PROKKA output GFF files).

2. No need to provide ID for all previous steps (prepare, scan, filter) when running any other task. SLING will assume the ID is the same as the current ID. This can be modified if running steps separately with multiple IDs (see wiki). 

3. GROUP step now uses a length coverage cutoff for the alignment (how much of the query length does the alignment cover), for two proteins to be considered the same in the sequence similarity network. The default length coverage is 0.75.

4. Default values changed for maximum overlap (from 300 to 50) and minimum blast identity (from 30 to 75).

5. Sixframe ORFs and annotation ORFs are now in a single file at the end of the PREPARE step and are treated in a single file along the entire program.

6. Outputs now return the sequences in nucleotide rather than protein sequence. It is easier to convert nucleotide to protein. The headers have slightly changed. Please refer to the wiki.

7. Major code clean-up and restructuring.

Minor changes:

1. Fixed a bug with `create_db` which was relying on a deprecated version of SLING using a configuration file.

2. `filter` now works with a pool object and is more efficient.

3. Added more checks for input of structural requirements by the user.

4. Nucleotide sequences which have more than 5% unknown bases (Ns or Xs), are removed already in the preparation step.

5. ORFs from the annotation file (GFF) which are shorter than the stated `min_orf_length` in the preparation step are removed from further analysis.

6. Installation requirements set for packages (fixed networkx bug)
