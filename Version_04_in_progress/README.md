# APOGEE V04 - Nanopore Metabarcoding Pipeline

## DESCRIPTION

APOGEE V04 is an enhanced version of V03 with advanced LCA (Lowest Common Ancestor) filtering and OTU collapsing functionality.

### Major Changes from V03:
- **Advanced LCA Filtering**: Replaces simple confidence-based filtering with sophisticated LCA strategies (strict, majority-rule, weighted)
- **OTU Collapsing & Species Summary**: New step to collapse OTU by complete taxonomic classification and generate species-level abundance tables
- **Bug Fix**: Resolves OTU ID collision where collapsed OTU are now prefixed as COTU_ (distinct from filtered OTU_)
- **Improved Taxonomy Resolution**: Multiple strategies for handling taxonomic conflicts in read assignments

The pipeline takes Nanopore fastq files (post-basecalling) and processes them into species-level abundance tables ready for import into R or other statistical analysis software.

![plot](https://github.com/hvanderheyden/APOGEE_nanopore_metabarcoding/blob/main/APOGEE.png)


## Installation

It is recommended to install the pipeline within a conda environment using the provided yml file. See [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) for instructions to install conda.

```commandline
# Make sure you have "git" installed in your base environment:
conda install git

# Clone repo to your preferred location:
git clone https://github.com/hvanderheyden/APOGEE_nanopore_metabarcoding

# Create conda environment
conda env create -f APOGEE.yml 

# Activate environment
conda activate Apogee-pipeline

# Make the file executable
chmod +x */APOGEE_V04.sh
```


## Usage

```commandline
conda activate Apogee-pipeline

/path/to/APOGEE_V04.sh \
  -i <input_dir> \
  -o <output_dir> \
  -r <ref_db> \
  -t <threads> \
  -T <taxonomy_db> \
  -F <filter_with_advanced_lca.py> \
  -S <generate_species_proportions.py> \
  -G <collapse_otu_by_taxonomy.py> \
  [-b <enable_read_splitting>] \
  [-R <reverse_primer>] \
  [-P <forward_primer>] \
  [-c <enable_clustering>] \
  [-x <clustering_identity>] \
  [-w <wordlength>] \
  [-C <mapping_confidence>] \
  [-v <query_coverage>] \
  [-m <length_filter>] \
  [-k <identity_threshold>] \
  [-L <lca_strategy>]
```


## REQUIRED ARGUMENTS

```
-i <input_dir>
     Path to directory containing input FASTQ files

-o <output_dir>
     Path to directory for output files (will be created if it doesn't exist)
     NOTE: Use absolute paths, not relative paths!

-r <ref_db>
     Path to the minimap2 index file (.mmi) for reference database
     Create with: minimap2 -d output.mmi sequences.fasta

-T <taxonomy_db>
     Path to the taxonomy file (.tsv) for reference database

-F <filter_with_advanced_lca.py>
     Path to filter_with_advanced_lca.py script for advanced PAF filtering

-S <generate_species_proportions.py>
     Path to generate_species_proportions.py script for species-level summary generation

-G <collapse_otu_by_taxonomy.py>
     Path to collapse_otu_by_taxonomy.py script for OTU collapsing

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
     Query coverage threshold for filtering (0.0-1.0, default: 0.85)

-m <length_filter>
     Minimum read length filter in bp (default: 450)
     Set to 0 to disable length filtering

-k <identity_threshold> 
     Identity threshold for filtering (0.0-1.0, default: 0.9)

-L <lca_strategy>
     LCA resolution strategy: 'strict', 'majority-rule', 'weighted', or 'bootstrap' (default: 'strict')
     - strict: Stops at first disagreement in taxonomy (most conservative)
     - majority-rule: Continues if >50% of alignments agree (balanced)
     - weighted: Uses weighted voting based on alignment confidence scores (more permissive)
     - bootstrap: Resamples alignments N=100 times per read to calculate taxonomy support (most robust)
       Use --bootstrap-n to modify the number of replicates
```


## DIFFERENCES WITH V03

### 1. Advanced LCA (Lowest Common Ancestor) Filtering (Major Change)
**V03**: Used simple confidence-based filtering with single LCA strategy
```
V03 filtering criteria:
- Mapping confidence >= 0.3
- Query coverage >= 0.0 (no filtering by default)
- Identity >= 0.0 (no filtering by default)
```

**V04**: Implements sophisticated LCA strategies with multiple resolution methods
```
V04 LCA strategies (-L parameter):
- strict: Stops at first disagreement in taxonomy (most conservative)
- majority-rule: Continues if >50% of alignments agree (balanced)
- weighted: Uses weighted voting based on confidence scores (more permissive)
- bootstrap: Resamples alignments to calculate taxonomy support + confidence (most robust)
```

**Impact**: Better handling of multi-mapping reads that align to multiple species. Previous version could arbitrarily assign reads; V04 chooses based on explicit LCA resolution strategy with confidence metrics.


## LCA (Lowest Common Ancestor) STRATEGIES EXPLAINED

When a read aligns to multiple species in the reference database, LCA determines the most specific taxonomic level that can be confidently assigned. V04 offers 4 strategies:

### STRICT
- Stops at the first taxonomic level where alignments disagree
- Requires 100% agreement to proceed to next level
- **Use when**: High-quality data, conservative analysis preferred
- **Example**: If read maps to [A. alternata, A. chlamydospora, A. infectoria] → assigns Alternaria genus (first level where 100% agree)

### MAJORITY-RULE
- Requires >50% of alignments to agree at each taxonomic level
- Continues until majority agreement fails
- **Use when**: Balanced approach, data with some noise
- **Example**: With 5 alignments, if 3 agree on species → allows species-level assignment

### WEIGHTED
- Uses confidence/identity scores to weight voting
- Alignments with higher scores count more strongly
- **Use when**: Data with variable quality, want maximum taxonomic resolution
- **Example**: Strong hit (conf: 0.99) vs weak hit (conf: 0.50) are weighted differently

### BOOTSTRAP
- Resamples alignments N times and recalculates LCA for each sample
- Reports final taxonomy + support percentage (% of resamples supporting this assignment)
- **Default N**: 100 replicates per read (can be modified with `--bootstrap-n` parameter)
- **Use when**: Data is noisy, need confidence metric, want robustness assessment
- **Example**: If bootstrapped taxonomy shows 87% support out of 100 replicates, assignment is robust


## 2. OTU Collapsing & Species-Level Summary (New Feature)
**V03**: Stopped at OTU level - no collapsing or species-level aggregation
**V04**: New 2-step aggregation:
1. **Step 8b (collapse_otu_by_taxonomy.py)**: Collapses 2.2M+ OTU → ~600 unique taxonomic profiles (COTU_)
2. **Step 8c (generate_species_proportions.py)**: Aggregates COTU → species level, generates summary.csv

**Output difference**:
- V03: `otu_table.csv`, `phyloseq_taxonomy.csv` (no species-level summary)
- V04: `collapsed_otu_table.csv`, `collapsed_taxonomy.csv`, **`summary.csv`** (species-level abundances)

### 3. OTU ID Naming (Bug Fix)
**V03**: Collapsed OTU named `OTU_1` through `OTU_570`
- Problem: Collision with filtered OTU IDs in taxonomy file → wrong species assignment

**V04**: Collapsed OTU named `COTU_1` through `COTU_570`  
- Solution: Unique prefix prevents ID collision
- Impact: Species abundances now mathematically correct

### 4. New Required Scripts
**V03**: `filter_with_confidence.py`
**V04**: 
- `filter_with_advanced_lca.py` (replaces filter_with_confidence.py - uses LCA strategies)
- `generate_species_proportions.py` (NEW - species aggregation)
- `collapse_otu_by_taxonomy.py` (NEW - OTU collapsing)


## OUTPUT FILES

The pipeline generates the following output files in the specified output directory:

### Main Results
- **summary.csv** - Species-level abundance table with read counts and proportions for all samples
- **collapsed_otu_table.csv** - Collapsed OTU table with COTU_ prefixed IDs and aggregated abundances
- **collapsed_taxonomy.csv** - Taxonomy assignments for collapsed OTU

### Intermediate Files (organized in subdirectories)
- `porechop/` - Adapter-removed FASTQ files
- `trimmed/` - Primer-trimmed FASTQ files (if enabled)
- `nanofilt/` - Quality and length filtered FASTQ files
- `clusters/` - Clustered FASTA centroids (if clustering enabled)
- `chimera/` - Chimera-filtered FASTA files
- `mapped/` - PAF mapping files from minimap2
- `filteredPAFs/` - Filtered OTU table and taxonomy from advanced LCA filtering


## EXAMPLE USAGE

### Default Parameters (Recommended)
```commandline
/path/to/APOGEE_V04.sh \
  -i /path/to/fastq/files \
  -o /path/to/output \
  -r /path/to/reference.mmi \
  -t 24 \
  -c false \
  -b false \
  -T /path/to/taxonomy.tsv \
  -F /path/to/filter_with_advanced_lca.py \
  -S /path/to/generate_species_proportions.py \
  -G /path/to/collapse_otu_by_taxonomy.py \
  -C 0.3 -v 0.85 -m 450 -k 0.9 -L strict
```

**Parameter explanation**:
- `-C 0.3` - Mapping confidence threshold
- `-v 0.85` - Query coverage threshold  
- `-m 450` - Minimum read length (bp)
- `-k 0.9` - Identity threshold
- `-L strict` - Use strict LCA strategy (recommended for most cases)


## VERSION INFORMATION

### V04 (Current) 
- Advanced LCA filtering strategies (strict, majority-rule, weighted)
- OTU collapsing by complete taxonomy
- Species-level abundance summary (summary.csv)
- Bug fix: Collapsed OTU use COTU_ prefix (prevents ID collision)

### V03 (Previous)
- Simple confidence-based filtering
- No OTU collapsing or species-level aggregation
- Known issue: OTU naming collision → incorrect species abundances

## SCRIPT COMPONENTS

### New in V04 (Replacements/Additions from V03)
- **filter_with_advanced_lca.py** - REPLACES V03's filter_with_confidence.py
  - Implements LCA strategies instead of simple confidence filtering
  - Handles multi-mapping reads with taxonomic voting
  
- **collapse_otu_by_taxonomy.py** - NEW
  - Collapses OTU with identical taxonomies
  - Renames to COTU_ prefix to avoid ID collisions
  
- **generate_species_proportions.py** - NEW
  - Aggregates collapsed OTU to species level
  - Generates summary.csv with species abundances

### Unchanged from V03
- **APOGEE_V04.sh** - Main pipeline orchestrator (updated to call new scripts)
- Porechop (adapter removal)
- Nanofilt (quality filtering)
- VSEARCH (chimera detection, optional clustering)
- minimap2 (sequence mapping)

## ACKNOWLEDGMENTS

### Development Tools
- V04 development benefited from AI-assisted code generation and debugging
- LCA algorithm implementations and documentation were developed with AI assistance
- All code has been reviewed, tested, and validated on mock community data


## TROUBLESHOOTING

### Common Issues

**"Using V03 results - species abundances seem wrong"**
- V03 had an OTU ID collision bug. Always use V04 for accurate species-level analysis.
- If you have V03 results, re-run with V04 to get correct proportions.

**"summary.csv missing or incomplete"**
- Verify `collapsed_taxonomy.csv` exists in output directory (created by collapse_otu_by_taxonomy.py)
- Check that `-G` parameter points to correct collapse_otu_by_taxonomy.py script
- Ensure output directory uses absolute path, not relative path

**"Wrong LCA strategy applied"**
- Specify `-L` parameter explicitly: `-L strict`, `-L majority-rule`, or `-L weighted`
- Default is `strict` if not specified
- Check pipeline logs to confirm which LCA method was used

**Output directory contains only porechop/ nanofilt/ etc. but no final results**
- Pipeline may still be running or may have failed after mapping
- Check for errors in mapping/filtering steps
- verify all Python scripts are executable: `chmod +x *.py` 
- Ensure all paths (especially `-o`) use absolute paths


## CITATION

If you use this pipeline, please cite:
- Original APOGEE: [Latorre-Pérez et al. 2021](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2021.768240/full)
- This version (V04) with modifications and bug fixes

