# APOGEE Nanopore Metabarcoding Pipeline — Version 05

## DESCRIPTION

Version 05 builds on previous versions by replacing the simple confidence-based filtering step
with an **advanced LCA (Lowest Common Ancestor) approach** and adding a **species clustering** feature.

Key changes compared to V03/V04:

- **Weighted LCA** replaces single-best-hit filtering: all hits above thresholds are combined,
  weighted by mapping confidence × identity, and the taxonomy is resolved to the deepest node
  supported by the majority of evidence.
- **Species clustering** (`-K`): closely related species that cannot be reliably distinguished
  by Nanopore ITS reads are merged into named clusters (e.g. *Botrytis cinerea* cluster),
  preventing systematic over-splitting of species counts.
- **LCA type** is selectable (`-L`): `weighted` (recommended), `strict`, `majority-rule`, or `bootstrap`.
- Intermediate steps (adapter removal, primer trimming, length filtering, chimera removal, mapping)
  are unchanged and cached — re-running with different LCA/filter parameters only re-executes
  the fast filtering steps.

The pre-processing steps remain: [Porechop](https://github.com/rrwick/Porechop) for adapter
removal, optional primer trimming with [bbmap's](https://bbmap.org/tools/removesmartbell)
removesmartbell, [NanoFilt](https://github.com/wdecoster/nanofilt) for quality/length filtering
(≥Q12, ≥350 bp), [vsearch](https://github.com/torognes/vsearch) `--uchime_denovo` for chimera
removal, [seqtk](https://github.com/lh3/seqtk) for FASTQ→FASTA conversion (when clustering is
disabled), and [minimap2](https://github.com/lh3/minimap2) for read mapping.

## Installation

1. Clone the repository (if not already done):
```bash
conda install git   # if git is not available in base
git clone https://github.com/hvanderheyden/APOGEE_nanopore_metabarcoding
```

2. Create the conda environment:
```bash
conda env create -f Version_05/APOGEE_V05.yml
conda activate Apogee-pipeline-v05
```

3. Make the script executable:
```bash
chmod +x Version_05/APOGEE_V05.sh
```

## Usage

```bash
/path/to/APOGEE_V05.sh \
  -i <input_dir> -o <output_dir> \
  -r <ref_db.mmi> -T <taxonomy.tsv> \
  -F <filter_with_advanced_lca.py> \
  -S <generate_species_proportions.py> \
  -G <collapse_otu_by_taxonomy.py> \
  -t <threads> \
  [-K <species_clusters.tsv>] \
  [-L <lca_type>] [-H <top_hits>] \
  [-v <query_coverage>] [-m <min_length>] [-k <identity>] [-C <confidence>] \
  [-b <enable_trimming>] [-R <reverse_primer>] [-P <forward_primer>] \
  [-c <enable_clustering>] [-x <clustering_identity>] [-w <wordlength>]
```

### Recommended parameters (optimized on mock communities)

```bash
-v 0.90 -m 450 -k 0.90 -C 0.3 -L weighted -H 13
```

## REQUIRED ARGUMENTS

```
-i <input_dir>
     Path to directory containing input FASTQ files

-o <output_dir>
     Path to output directory (created if absent)

-r <ref_db>
     Path to minimap2 index file (.mmi)
     Create with: minimap2 -d output.mmi sequences.fasta

-T <taxonomy_db>
     Path to taxonomy TSV file for the reference database

-F <filter_with_advanced_lca.py>
     Path to the advanced LCA filtering script

-S <generate_species_proportions.py>
     Path to the species proportions summary script

-G <collapse_otu_by_taxonomy.py>
     Path to the OTU collapse script

-t <threads>
     Number of threads (positive integer)
```

## OPTIONAL ARGUMENTS

```
-K <species_clusters.tsv>
     Path to species clusters file (TSV). When provided, species that belong to
     the same cluster are merged into a single cluster OTU after LCA assignment.
     Omit to skip species clustering (V04 behavior).

-L <lca_type>
     LCA algorithm: weighted (default, recommended), strict, majority-rule, bootstrap

-H <top_hits>
     Number of top hits per read considered by the LCA (default: 5)
     Higher values capture more ambiguous hits; 13 recommended for ITS.

-v <query_coverage>
     Minimum fraction of the query read covered by the alignment (default: 0.85)
     Recommended: 0.90

-m <min_qlen>
     Minimum query read length in bp (default: 0)
     Recommended: 450

-k <identity_threshold>
     Minimum alignment identity (1 - divergence) for a hit to be retained (default: 0.9)

-C <confidence>
     Minimum mapping confidence threshold, derived from MapQ: 1 - 10^(MapQ/-10)
     (default: 0.5). Has little effect when minimap2 is run with --secondary=no.

-b <enable_trimming>
     Enable primer trimming with removesmartbell (true/false, default: false)

-R <reverse_primer>
     Reverse primer sequence (required when -b true)

-P <forward_primer>
     Forward primer sequence (required when -b true)

-c <enable_clustering>
     Enable VSEARCH read clustering before mapping (true/false, default: false)

-x <identity>
     Identity threshold for VSEARCH clustering (default: 0.97)

-w <wordlength>
     Word length for VSEARCH k-mer matching (default: 10)
```

## OUTPUT FILES

```
<output_dir>/
├── summary.csv                              Top species per sample (proportions)
├── collapsed_otu_table.csv                  OTU abundance table (collapsed by taxonomy)
├── collapsed_taxonomy.csv                   Taxonomy for each collapsed OTU
├── filteredPAFs/
│   ├── filtered_otu.tsv                     Raw filtered OTU table (pre-clustering)
│   ├── phyloseq_taxonomy_filtered_otu.csv   LCA taxonomy assignments
│   ├── clustered_otu.tsv                    OTU table after species clustering (if -K used)
│   ├── clustered_otu_taxonomy.tsv           Taxonomy after clustering
│   └── clustered_otu_clustering_changes.txt Log of clustering modifications
├── mapped/                                  PAF files from minimap2
├── chimera/                                 Chimera-filtered FASTA files
├── fasta/                                   FASTQ→FASTA converted files
├── nanofilt/                                Quality/length filtered FASTQ files
├── trimmed/                                 Primer-trimmed FASTQ files (if -b true)
└── porechop/                                Adapter-removed FASTQ files
```

## Species clusters file format

The species clusters TSV (`-K`) has no header and two columns:

```
<species_name>\t<cluster_id>
```

Species sharing the same `cluster_id` integer are merged. The cluster name used in the output
is the species name of the representative (first encountered) member.

Example:
```
Botrytis_cinerea	1
Botrytis_pseudocinerea	1
Botrytis_cinerea_var_cinerea	1
Alternaria_alternata	101
Alternaria_arborescens	101
```

## Parameter sweep utility

`param_sweep.sh` tests a grid of `-k`, `-v`, `-m`, `-C` values and reports the unclassified rate
per sample for each combination. It reuses existing PAF files (only re-runs the fast filtering
steps). Results are written to `/tmp/apogee_param_sweep_<date>.tsv`.

```bash
# Edit the parameter grids at the top of the script, then:
bash param_sweep.sh
```

## To do

- Replace Porechop (unmaintained) with `dorado trim` for adapter removal.
- Add support for multiple reference databases.
