# Helper Scripts

Some helper scripts for SLING inputs/outputs


### Convert GFF file with appended FASTA to FASTA format

`gff_to_fasta`		Convert GFF files to FASTA files with same identifier (for SLING input)

#### Usage: 

`python gff_to_fasta.py [options] <gff_dir> <gff_suffix>`




#### Required arguments:
  
  `--gff_dir DIR `    Directory containing all the GFF files

#### Optional arguments:
  
  `-h, --help `       show this help message and exit

`  --gff_suffix STR ` Suffix used for the GFF files [.gff]

 ` --fasta_dir DIR `  Directory for output FASTA files


### Get GFF CDS identifiers

`cds_identifier.py`    Add all the CDS unique identifiers from the GFF files to SLING outputs.

#### Usage:

`python cds_identifier.py [options] <gff_dir> <group_dir>`

#### Positinal arguments

`<gff_dir> DIR`    Directory containing all the GFF files

` <group_dir> DIR`   Directory containing all outputs of SLING GROUP

#### Optional arguments:
 ` -h, --help  `      show this help message and exit
 
 ` --gff_suffix STR`  Suffix used for the GFF files `[.gff]`
 
`  --filter_dir DIR ` Directory containing all outputs of sling FILTER, if wanted `[None]`
