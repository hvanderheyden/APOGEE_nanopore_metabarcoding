
### DESCRIPTION

This pipeline is used for Nanopore metabarcoding. It allows fastq files (post-basecalling)
to be processed into tables ready for import into R. This pipeline was adapted from the work
of [Latorre-Pérez et al. 2021](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2021.768240/full), with the following modifications:

The curent version is still using [Porechop](https://github.com/rrwick/Porechop) for adapter 
removal, but we have added an optional read splitting step using [bbmap's](https://bbmap.org/tools/removesmartbell) removesmartbell tool. [Nanofilt](https://github.com/wdecoster/nanofilt) is still used to filter reads based on length and quality. We use [vsearch](https://github.com/torognes/vsearch) --uchime_denovo for chimera filtering and for an optional clustering step. Finally, [minimap2](https://github.com/lh3/minimap2) is used for read mapping, while OTU and taxonomy tables are generated using customized Python scripts. 

![plot](https://github.com/hvanderheyden/APOGEE_nanopore_metabarcoding/blob/main/APOGEE.png)

### Dependencies
```
Required software (install via conda):
  - porechop: Adapter removal
  - NanoFilt: Length/quality filtering
  - vsearch: Sequence clustering
  - minimap2: Read mapping and alignment
  - yacrd: Chimera detection
  - Python 3 with pandas: Data processing
```
## Installation
```
Install all dependencies:
  conda env create -f APOGEE.yml
  conda activate Apogee-pipeline

chmod +x */APOGEE.sh
```
### Basic usage
```
conda activate Apogee-pipeline

/minimap2_v02.sh \
  -i <input_dir> \
  -o <output_dir> \
  -r <ref_db> \
  -t <threads> \
  -T <taxonomy_db> \
  -F <filter_script> \
  -S <taxonomy_script> \
  [-c <enable_clustering>] \
  [-x <identity>] \
  [-w <wordlength>] \
  [-C <confidence>]
```
### REQUIRED ARGUMENTS
```
-i <input_dir>
     Path to directory containing input FASTQ files (gzipped or uncompressed)

  -o <output_dir>
     Path to directory for output files (will be created if it doesn't exist)

  -r <ref_db>
     Path to minimap2 index file (.mmi) for reference database
     Note: Create with: minimap2 -d output.mmi sequences.fasta

  -s <script_dir>
     Path to directory containing taxonomyTable.py script

  -t <threads>
     Number of threads for parallel processing (positive integer)

  -T <taxonomy_db>
     Path to taxonomy database file (TSV format: accession<tab>taxonomy)

  -F <filter_script>
     Path to filter_with_confidence.py script for PAF filtering
```
### OPTIONAL ARGUMENTS
```
  -c <enable_clustering>
     Enable sequence clustering with VSEARCH (true/false, default: false)

  -x <identity>
     Identity threshold for clustering (0.0-1.0, default: 0.97)

  -w <wordlength>
     Word length for clustering k-mer matching (default: 10)

  -C <confidence>
     Minimum confidence threshold for filtering (default: 1)

```
### OUTPUT FILES
```
  otu_table.csv
    Format: Accession,sample1,sample2,...,TotalCount,MappingConfidence

  taxonomy_table.csv
    Format: #OTU ID,Domain,Phylum,Class,Order,Family,Genus,Species
```






