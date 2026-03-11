# APOGEE Nanopore Metabarcoding Pipeline

## DESCRIPTION

This pipeline is used for Nanopore metabarcoding. It allows fastq files (post-basecalling)
to be processed into tables ready for import into R. This pipeline was adapted from the work
of [Latorre-Pérez et al. 2021](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2021.768240/full), with the following modifications:

The current version is still using [Porechop](https://github.com/rrwick/Porechop) for adapter 
removal, but we have added an optional primer trimming step using [bbmap's](https://bbmap.org/tools/removesmartbell) removesmartbell tool. [Nanofilt](https://github.com/wdecoster/nanofilt) is used to filter reads based on quality (≥12) and length (≥350 bp). We use [vsearch](https://github.com/torognes/vsearch) --uchime_denovo for chimera filtering and for an optional clustering step. [seqtk](https://github.com/lh3/seqtk) is used for FASTQ to FASTA conversion (when clustering is disabled). Finally, [minimap2](https://github.com/lh3/minimap2) is used for read mapping, while OTU and taxonomy tables are generated using customized Python scripts.

![plot](https://github.com/hvanderheyden/APOGEE_nanopore_metabarcoding/blob/main/APOGEE.png)


## Installation

It is recommended to install the pipeline within a conda environment using the provided yml file. See [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) for instructions to install conda. 


1. Clone repository:
```commandline
# Make sure you have "git" installed in your base environment:
conda install git

# Clone repo to your preferred location:
git clone https://github.com/hvanderheyden/APOGEE_nanopore_metabarcoding
```
2. Create virtual environment:
```commandline
# create conda environment
conda env create -f APOGEE.yml 

# Activate environment
conda activate Apogee-pipeline

# Make the file executable
chmod +x */APOGEE_V03.sh
```


## Usage
```commandline
conda activate Apogee-pipeline

/path/to/APOGEE_V03.sh -i <input_dir> -o <output_dir> -r <ref_db> -t <threads> -T <taxonomy_db> -F <filter_with_confidence.py> -S <taxonomyTable.py> [-b <enable_read_splitting>] [-R <reverse_primer>] [-P <forward_primer>] [-c <enable_clustering>] [-x <clustering_identity>] [-w <wordlength>] [-C <mapping_confidence>] [-v <query_coverage>] [-k <identity_threshold>]
```

## REQUIRED ARGUMENTS
```
-i <input_dir>
     Path to directory containing input FASTQ files

-o <output_dir>
     Path to directory for output files (will be created if it doesn't exist)

-r <ref_db>
     Path to the minimap2 index file (.mmi) for reference database
     Note: Create with: minimap2 -d output.mmi sequences.fasta

-T <taxonomy_db>
     Path to the taxonomy file (.tsv) for reference database

-F <filter_with_confidence.py>
     Path to filter_with_confidence.py script for PAF filtering

-S <taxonomyTable.py>
     Path to directory containing taxonomyTable.py script

-t <threads>
     Number of threads for parallel processing (positive integer)
```

## OPTIONAL ARGUMENTS
```
-b <enable_read_splitting>
     Enable primer trimming using bbmap's removesmartbell.sh (default: false)

-R <reverse_primer>
     Sequence of the reverse primer used in the primer trimming step

-P <forward_primer>
     Sequence of the forward primer used in the primer trimming step

-c <enable_clustering>
     Enable sequence clustering with VSEARCH (true/false, default: false)

-x <identity>
     Identity threshold for clustering (0.0-1.0, default: 0.97)

-w <wordlength>
     Word length for clustering k-mer matching (default: 10)

-C <confidence>
     Mapping confidence threshold for filtering (0.0-1.0, default: 0.3)

-v <query_coverage>
     Query coverage threshold for filtering (0.0-1.0, default: 0.0)

-k <identity_threshold> 
     Identity threshold for filtering (0.0-1.0, default: 0.0)
```

## OUTPUT FILES

The pipeline generates the following output files in the specified output directory:

- **otu_table.csv** - OTU abundance table with mapping statistics (confidence scores, query coverage, identity)
- **phyloseq_taxonomy.csv** - Taxonomy table with progressive filtering based on confidence, query coverage, and identity thresholds

Intermediate files are organized in subdirectories:
- `porechop/` - Adapter-removed FASTQ files
- `trimmed/` - Primer-trimmed FASTQ files (if enabled)
- `nanofilt/` - Quality and length filtered FASTQ files
- `clusters/` - Clustered FASTA centroids (if clustering enabled)
- `chimera/` - Chimera-filtered FASTA files
- `mapped/` - PAF mapping files
- `filteredPAFs/` - Filtered OTU table (TSV format)

