# Helper Scripts

Some helper scripts for parsing through SLING outputs

### Get GFF CDS identifiers
`cds_identifier.py`    Add all the CDS unique identifiers from the GFF files to SLING outputs.

#### Usage:
`python cds_identifier.py [options] <gff_dir> <group_dir>`

#### Positinal arguments
`<gff_dir> DIR`    Directory containing all the GFF files

` <group_dir> DIR`   Directory containing all outputs of SLING GROUP

#### Pptional arguments:
 ` -h, --help  `      show this help message and exit
 
 ` --gff_suffix STR`  Suffix used for the GFF files `[.gff]`
 
`  --filter_dir DIR ` Directory containing all outputs of sling FILTER, if wanted `[None]`
