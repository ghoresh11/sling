# Helper Scripts

Some helper scripts for parsing through SLING outputs

`cds_identifier.py`    Add all the CDS unique identifiers from the GFF files to SLING outputs.

### Usage:
python cds_identifier.py [options] <gff_dir> <group_dir>
<gff_dir> DIR Directory containing all the GFF files
 <group_dir> DIR Directory containing all outputs of SLING GROUP

optional arguments:
  -h, --help        show this help message and exit
  --gff_suffix STR  Suffix used for the GFF files [.gff]
  --filter_dir DIR  Directory containing all outputs of sling FILTER, if wanted [None]
