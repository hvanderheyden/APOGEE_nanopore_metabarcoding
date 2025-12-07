
### Installation

```
conda env create -f */environment.yaml
conda activate Nanopore_apogee

chmod +x */APOGEE.sh
```
### Usage

```
/minimap2_v02.sh \
 -i /media/herve/10TB/Apogee/6_mock/6_minimap2/reads \
 -o /media/herve/10TB/Apogee/6_mock/11_minimap2_clustered \
 -r /media/herve/10TB/Apogee/6_mock/6_minimap2/ITS-RefDB_V02.mmi \
 -s /media/herve/10TB/Apogee/6_mock/6_minimap2 \
 -t 32 \
 -c true \
 -x 0.98 \
 -T /media/herve/10TB/Apogee/6_mock/6_minimap2/taxonomy.tsv \
 -F /media/herve/10TB/Apogee/5_Scripts/filter_with_confidence.py \
 -C 1
```
```
  -i input_dir 
   -o  output_dir
   -r  ref_db
   -s  script_dir
   -t  threads
   -c  enable_clustering
   -x  identity
   -w  wordlength
   -T  taxonomy_db
   -F  filter_script
   -C  confidence
```

### Arguments:
```
  -i         Directory containing input FASTQ files to analyze
  -o    CSV file containing target sequences with their metadata
  -r            Output directory for results and visualization files
  -s
  -t
  -c true \
  -x 0.98 \
  -T /media/herve/10TB/Apogee/6_mock/6_minimap2/taxonomy.tsv \
  -F /media/herve/10TB/Apogee/5_Scripts/filter_with_confidence.py \
  -C 1

```
### Example 

```
/media/herve/10TB/Apogee/6_mock/6_minimap2/minimap2_v02.sh \
 -i /media/herve/10TB/Apogee/6_mock/6_minimap2/reads \
 -o /media/herve/10TB/Apogee/6_mock/11_minimap2_clustered \
 -r /media/herve/10TB/Apogee/6_mock/6_minimap2/ITS-RefDB_V02.mmi \
 -s /media/herve/10TB/Apogee/6_mock/6_minimap2 \
 -t 32 \
 -c true \
 -x 0.98 \
 -T /media/herve/10TB/Apogee/6_mock/6_minimap2/taxonomy.tsv \
 -F /media/herve/10TB/Apogee/5_Scripts/filter_with_confidence.py \
 -C 1
```
   
