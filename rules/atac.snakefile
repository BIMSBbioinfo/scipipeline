import os
# here we might need to set the paths
from src.data.utils import download_folder
from src.data.utils import unpack_labels_to_csv

FIRST_MATE = 
SECOND_MATE = 
TRIM_PATTERN = 
FIRST_MATE_TRIMMED = TRIM_PATTERN + '_1.fastq'
SECOND_MATE_TRIMMED = TRIM_PATTERN + '_2.fastq'
ADAPTERS = /data/ohler/Scott/Asli_Scripts/sciAtacAdapters.fa

FLY_MAPPING =
FLY_INDEX
PSEUDOGENOME_MAPPING = 
PSEUDOGENOME

# ------------------------- #
# Adapter trimming

rule adapter_trimming:
    """Trim adapters using flexbar"""
    input: 
        first=FIRST_MATE, 
        second=SECOND_MATE
    output: 
        first=FIRST_MATE_TRIMMED, 
        second=SECOND_MATE_TRIMMED
    script:
        "/gnu/var/guix/profiles/custom/bimsb/bin/flexbar "
	+ "-r {input.first} -p {input.second} -t " + TRIM_PATTERN
        + " -f i1.8 -u 10 -ae RIGHT -at 1.0"

INPUT_ALL.append(rules.adapter_trimming.output)

# ------------------------- #
# Mapping to fly genome

rule read_mapping_fly:
    input: 
        first=FIRST_MATE_TRIMMED, 
        second=SECOND_MATE_TRIMMED,
        bowtie_index=FLY_INDEX
    output: FLY_MAPPING
    script:
	"bowtie2 -p 4 -X 1500 --no-mixed "
	"--no-discordant -x {input.bowtie_index} "
        "-1 {input.first} -2 {input.second} | samtools view -bS - > "
        "{output}"

INPUT_ALL.append(rules.read_mapping_fly.output)

# ------------------------- #
# Construct a pseudo genome in fasta format

# ------------------------- #
# Create bowtie2 index for pseudo genome

# ------------------------- #
# Mapping to pseudo genome


