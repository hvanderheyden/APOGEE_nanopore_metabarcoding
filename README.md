
### DESCRIPTION

This pipeline is used for Nanopore metabarcoding. It allows fastq files (post-basecalling)
to be processed into tables ready for import into R. This pipeline was adapted from the work
of [Latorre-Pérez et al. 2021](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2021.768240/full), with the following modifications:

The curent version is still using [Porechop](https://github.com/rrwick/Porechop) for adapter 
removal, but we have added an optional read splitting step using [bbmap's](https://bbmap.org/tools/removesmartbell) removesmartbell tool. [Nanofilt](https://github.com/wdecoster/nanofilt) is still used to filter reads based on length and quality. We use [vsearch](https://github.com/torognes/vsearch) --uchime_denovo for chimera filtering and for an optional clustering step. Finally, [minimap2](https://github.com/lh3/minimap2) is used for read mapping, while OTU and taxonomy tables are generated using customized Python scripts. 

![plot](https://github.com/hvanderheyden/APOGEE_nanopore_metabarcoding/blob/main/APOGEE.png)


## Installation

It is recommanded to install the pipeline within a conda environment using the provided yml file. See [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) for instructions to install conda. 


1. Clone repository:
```commandline
# Make sure you have "git" install in your base enviroment:
conda install git

# Clone repo to your prefered location:
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

### Usage
conda activate Apogee-pipeline

/APOGEE_V03.sh -i <input_dir> -o <output_dir> -r <ref_db> -T <taxonomy_db> -F <filter_with_confidence.py> -S <taxonomyTable.py> -b <enable_read_splitting> -R <reverse_primer> -P <forward_primer> -c <enable_clustering> -x <clustering_identity> -w <wordlength> -C <mapping_confidence> -v <query_coverage> -k <identity_threshold> 


### REQUIRED ARGUMENTS
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

### OPTIONAL ARGUMENTS
```
  -b <enable_read_splitting>
     Enable read splitting using bbmap's removesmartbell.sh

  -R <reverse_primer>
     Sequence of the reverse primer used in the read splitting step

  -P <forward_primer>
     Sequence of the forward primer used in the read splitting step

  -c <enable_clustering>
     Enable sequence clustering with VSEARCH (true/false, default: false)

  -x <identity>
     Identity threshold for clustering (0.0-1.0, default: 0.97)

  -w <wordlength>
     Word length for clustering k-mer matching (default: 10)

  -C <confidence>
     Mapping confidence threshold for filtering

  -v <query_coverage>
     Query coverage threshold for filtering

  -k <identity_threshold> 
     Identity threshold for filtering

```
### OUTPUT FILES
```
  otu_table.csv
    Format: Accession,sample1,sample2,...,TotalCount,MappingConfidence

  taxonomy_table.csv
    Format: #OTU ID,Domain,Phylum,Class,Order,Family,Genus,Species
```






