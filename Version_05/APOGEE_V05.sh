#!/bin/bash
set -euo pipefail

trap 'echo "ERROR: Script failed at line $LINENO" >&2' EXIT

###################
## Example usage ##
###################

#chmod +x path/to/APOGEE_V05.sh

# Complete example with CLUSTERING (V05 NEW FEATURE)
# /path/to/APOGEE_V05.sh \
# -i /input/fastq \
# -o /output \
# -r /ref/ITS-RefDB_V03_fixed.mmi \
# -T /ref/ITSRefDB_V03_taxonomy.tsv \
# -F /scripts/filter_with_advanced_lca.py \
# -S /scripts/generate_species_proportions.py \
# -G /scripts/collapse_otu_by_taxonomy.py \
# -K /ref/species_clusters.tsv \
# -t 24 \
# -x 0.99 -w 8 \
# -v 0.85 -m 450 -k 0.9 \
# -C 0.5 \
# -L weighted \
# -H 5

# Example without clustering (original V04 mode)
# /path/to/APOGEE_V05.sh \
# -i /input/fastq \
# -o /output \
# -r /ref/ITS-RefDB_V03_fixed.mmi \
# -T /ref/ITSRefDB_V03_taxonomy.tsv \
# -F /scripts/filter_with_advanced_lca.py \
# -S /scripts/generate_species_proportions.py \
# -G /scripts/collapse_otu_by_taxonomy.py \
# -t 24 \
# -x 0.99 -w 8 \
# -v 0.85 -m 450 -k 0.9 \
# -L weighted \
# -H 5
# (omit -K to skip clustering)

#############
# ARGUMENTS #
#############

# Parse command-line arguments
while getopts ":i:o:r:t:c:x:w:T:F:S:C:v:m:k:b:R:P:L:G:H:K:" opt; do
  case ${opt} in
    i ) input_dir="$OPTARG" ;;
    o ) output_dir="$OPTARG" ;;
    r ) ref_db="$OPTARG" ;;
    t ) threads="$OPTARG" ;;
    c ) enable_clustering="$OPTARG" ;;
    x ) identity="$OPTARG" ;;
    w ) wordlength="$OPTARG" ;;
    T ) taxonomy_db="$OPTARG" ;;
    F ) filter_script="$OPTARG" ;;
    S ) taxonomy_script="$OPTARG" ;;
    C ) confidence="$OPTARG" ;;
    v ) qcov="$OPTARG" ;;
    m ) min_qlen="$OPTARG" ;;
    k ) identity_threshold="$OPTARG" ;;
    b ) enable_trimming="$OPTARG" ;;
    R ) reverse_primer="$OPTARG" ;;
    P ) forward_primer="$OPTARG" ;;
    L ) lca_type="$OPTARG" ;;
    G ) collapse_script="$OPTARG" ;;
    H ) top_species_count="$OPTARG" ;;
    K ) clusters_file="$OPTARG" ;;
    \? ) echo "ERROR: Invalid option: -$OPTARG" 1>&2; exit 1 ;;
    : ) echo "ERROR: Invalid option: -$OPTARG requires an argument" 1>&2; exit 1 ;;
  esac
done
shift $((OPTIND -1))

########################
# Validate Arguments   #
########################

if [[ -z "${input_dir:-}" ]] || [[ -z "${output_dir:-}" ]] || [[ -z "${ref_db:-}" ]] || [[ -z "${threads:-}" ]] || [[ -z "${taxonomy_db:-}" ]] || [[ -z "${filter_script:-}" ]] || [[ -z "${taxonomy_script:-}" ]] || [[ -z "${collapse_script:-}" ]]; then
  echo "ERROR: Missing required arguments" >&2
  echo "Usage: $0 -i <input_dir> -o <output_dir> -r <ref_db> -t <threads> -T <taxonomy_db> -F <filter_script> -S <taxonomy_script> -G <collapse_script> [-c <enable_clustering>] [-x <identity>] [-w <wordlength>] [-C <confidence>] [-v <qcov>] [-m <min_qlen>] [-k <identity_threshold>] [-b <enable_trimming>] [-R <reverse_primer>] [-P <forward_primer>] [-L <lca_type>]" >&2
  exit 1
fi

# Validate that input directory and reference DB exist
if [[ ! -d "${input_dir}" ]]; then
  echo "ERROR: Input directory does not exist: ${input_dir}" >&2
  exit 1
fi

if [[ ! -f "${ref_db}" ]]; then
  echo "ERROR: Reference database does not exist: ${ref_db}" >&2
  exit 1
fi

if [[ ! -f "${taxonomy_db}" ]]; then
  echo "ERROR: Taxonomy database does not exist: ${taxonomy_db}" >&2
  exit 1
fi

if [[ ! -f "${filter_script}" ]]; then
  echo "ERROR: Filter script does not exist: ${filter_script}" >&2
  exit 1
fi

if [[ ! -f "${taxonomy_script}" ]]; then
  echo "ERROR: Taxonomy script does not exist: ${taxonomy_script}" >&2
  exit 1
fi

if [[ ! -f "${collapse_script}" ]]; then
  echo "ERROR: Collapse script does not exist: ${collapse_script}" >&2
  exit 1
fi

# Check clusters file if provided (V05 species clustering feature)
if [[ -n "${clusters_file:-}" ]] && [[ ! -f "${clusters_file}" ]]; then
  echo "WARNING: Clusters file not found: ${clusters_file}" >&2
  echo "         Species clustering will be skipped" >&2
  clusters_file=""
fi

# Set defaults if not provided
identity="${identity:-0.97}"
wordlength="${wordlength:-10}"
enable_clustering="${enable_clustering:-false}"
enable_trimming="${enable_trimming:-false}"
confidence="${confidence:-1.0}"
qcov="${qcov:-0.85}"
min_qlen="${min_qlen:-0}"
identity_threshold="${identity_threshold:-0.9}"
reverse_primer="${reverse_primer:-}"
forward_primer="${forward_primer:-}"
lca_type="${lca_type:-weighted}"
top_species_count="${top_species_count:-5}"
clusters_file="${clusters_file:-}"

# Validate numeric parameters
if ! [[ "${threads}" =~ ^[0-9]+$ ]] || [[ "${threads}" -lt 1 ]]; then
  echo "ERROR: threads must be a positive integer" >&2
  exit 1
fi

if ! [[ "${identity}" =~ ^0\.[0-9]+$ ]] || (( $(echo "${identity} < 0 || ${identity} > 1" | bc -l) )); then
  echo "ERROR: identity must be between 0 and 1" >&2
  exit 1
fi

if ! [[ "${confidence}" =~ ^[0-9]+\.?[0-9]*$ ]] || (( $(echo "${confidence} < 0 || ${confidence} > 1" | bc -l) )); then
  echo "ERROR: confidence must be between 0 and 1" >&2
  exit 1
fi

if ! [[ "${qcov}" =~ ^[0-9]+\.?[0-9]*$ ]] || (( $(echo "${qcov} < 0 || ${qcov} > 1" | bc -l) )); then
  echo "ERROR: qcov must be between 0 and 1" >&2
  exit 1
fi

if ! [[ "${min_qlen}" =~ ^[0-9]+$ ]] || [[ "${min_qlen}" -lt 0 ]]; then
  echo "ERROR: min_qlen must be a non-negative integer" >&2
  exit 1
fi

if ! [[ "${identity_threshold}" =~ ^[0-9]+\.?[0-9]*$ ]] || (( $(echo "${identity_threshold} < 0 || ${identity_threshold} > 1" | bc -l) )); then
  echo "ERROR: identity_threshold must be between 0 and 1" >&2
  exit 1
fi

# Validate lca_type parameter
valid_lca_types="strict weighted bootstrap majority-rule"
if ! echo "${valid_lca_types}" | grep -q "${lca_type}"; then
  echo "ERROR: lca_type must be one of: ${valid_lca_types}" >&2
  exit 1
fi

# Validate top_species_count parameter
if ! [[ "${top_species_count}" =~ ^[0-9]+$ ]] || [[ "${top_species_count}" -lt 1 ]]; then
  echo "ERROR: top_species_count must be a positive integer" >&2
  exit 1
fi

# Validate trimming parameters
if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
  if [[ -z "${reverse_primer}" ]] || [[ -z "${forward_primer}" ]]; then
    echo "ERROR: When trimming is enabled, both -R and -P primers must be provided" >&2
    exit 1
  fi
fi

echo "==============================================="
echo "Metabarcoding Pipeline v05 - With Species Clustering - Multi-filter with Optional Trimming"
echo "==============================================="
echo "Input dir:           ${input_dir}"
echo "Output dir:          ${output_dir}"
echo "Reference DB:        ${ref_db}"
echo "Threads:             ${threads}"
echo "Trimming:            ${enable_trimming}"
if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
  echo "  Forward primer:    ${forward_primer}"
  echo "  Reverse primer:    ${reverse_primer}"
fi
echo "Clustering:          ${enable_clustering}"
echo "Identity (clustering): ${identity}"
echo "Word length:         ${wordlength}"
echo "Taxonomy DB:         ${taxonomy_db}"
echo "Filter script:       ${filter_script}"
echo "Taxonomy script:     ${taxonomy_script}"
echo "Mapping confidence:  ${confidence}"
echo "Query coverage:      ${qcov}"
echo "Min query length:    ${min_qlen} bp"
echo "Identity threshold:  ${identity_threshold}"
echo "LCA type:            ${lca_type}"
echo "Top species count:   ${top_species_count} per sample"
echo "==============================================="
echo ""

#############################
# Create Output Directories #
#############################

mkdir -p "${output_dir}/porechop"
if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
  mkdir -p "${output_dir}/trimmed"
fi
mkdir -p "${output_dir}/nanofilt"
mkdir -p "${output_dir}/clusters"
mkdir -p "${output_dir}/chimera"
mkdir -p "${output_dir}/mapped"
mkdir -p "${output_dir}/filteredPAFs"

echo "Output directories created."
echo ""

########################
# Read Pre-processing  #
########################

echo "Starting analysis"
echo "Pre-processing the reads"
echo "------------------------"
echo ""

###################
# Adapter Removal #
###################

if [[ -n "$(find "${output_dir}/porechop" -maxdepth 1 -type f -name '*.fastq' 2>/dev/null)" ]]; then
    echo "✓ Step 1: Skipping adapter removal (output already exists)"
else
    echo "► Step 1: Adapter removal"
    echo "  Processing files from: ${input_dir}"
    input_count=$(find "${input_dir}" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) 2>/dev/null | wc -l)
    if [[ ${input_count} -eq 0 ]]; then
      echo "ERROR: No FASTQ files found in ${input_dir}" >&2
      exit 1
    fi
    echo "  Found ${input_count} input file(s)"
    
    time find "${input_dir}" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) -print0 | \
      xargs -0 -I {} -P "${threads}" bash -c 'porechop -t 1 -i "$1" -o "${2}/porechop/$(basename "$1" | sed "s/\(\.fastq\(\.gz\)\?\)$/-porechop.fastq/")" || { echo "ERROR: Adapter removal failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
    echo "  ✓ Adapter removal completed"
fi
echo ""

###################################################
# Primer Trimming (Optional) with removesmartbell #
###################################################

if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
    if [[ -n "$(find "${output_dir}/trimmed" -maxdepth 1 -type f -name '*.fastq' 2>/dev/null)" ]]; then
        echo "✓ Step 2: Skipping primer trimming (output already exists)"
        trimming_input="${output_dir}/trimmed"
    else
        echo "► Step 2: Primer trimming with removesmartbell"
        echo "  Forward primer: ${forward_primer}"
        echo "  Reverse primer: ${reverse_primer}"
        
        time find "${output_dir}/porechop" -maxdepth 1 -type f -name '*.fastq' -print0 | \
          xargs -0 -I {} -P "${threads}" bash -c 'removesmartbell.sh in="$1" out="${2}/trimmed/$(basename "$1" .fastq).fastq" adapter="'"${forward_primer}"'",RC_"'"${reverse_primer}"'" qin=33 2>/dev/null || { echo "ERROR: Primer trimming failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
        echo "  ✓ Primer trimming completed"
        trimming_input="${output_dir}/trimmed"
    fi
    nanofilt_input="${trimming_input}"
else
    echo "✓ Step 2: Primer trimming skipped (enable_trimming=false)"
    nanofilt_input="${output_dir}/porechop"
fi
echo ""

#################################
# Length Filtering              #
#################################

if [[ -n "$(find "${output_dir}/nanofilt" -maxdepth 1 -type f -name '*.fastq' 2>/dev/null)" ]]; then
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 3: Skipping length filtering (output already exists)"
    else
        echo "✓ Step 2: Skipping length filtering (output already exists)"
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "► Step 3: Length filtering"
    else
        echo "► Step 2: Length filtering"
    fi
    echo "  Filtering: quality >= 12, length: 350 bp"
    
    time find "${nanofilt_input}" -maxdepth 1 -type f -name '*.fastq' -print0 | \
      xargs -0 -I {} -P "${threads}" bash -c 'NanoFilt -q 12 -l 350 < "$1" > "${2}/nanofilt/$(basename "$1" | sed "s/porechop/nanofilt/; s/trimmed/nanofilt/")" || { echo "ERROR: Length filtering failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
    echo "  ✓ Length filtering completed"
fi
echo ""

######################################
# Clustering (Optional) with VSEARCH #
######################################

if [[ "${enable_clustering}" == "true" ]] || [[ "${enable_clustering}" == "True" ]] || [[ "${enable_clustering}" == "TRUE" ]]; then
    if [[ -n "$(find "${output_dir}/clusters" -maxdepth 1 -type f \( -name '*.fasta' -o -name '*.fastq' \) 2>/dev/null)" ]]; then
        if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
            echo "✓ Step 4: Skipping clustering (output already exists)"
        else
            echo "✓ Step 3: Skipping clustering (output already exists)"
        fi
        clustering_input="${output_dir}/clusters"
    else
        if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
            echo "► Step 4: Clustering with VSEARCH"
        else
            echo "► Step 3: Clustering with VSEARCH"
        fi
        echo "  Identity threshold: ${identity}"
        echo "  Word length: ${wordlength}"
        
        time find "${output_dir}/nanofilt" -maxdepth 1 -type f -name '*.fastq' -print0 | \
          xargs -0 -I {} -P "${threads}" bash -c 'vsearch --cluster_fast "$1" --id '"${identity}"' --wordlength '"${wordlength}"' --centroids "${2}/clusters/$(basename "$1" .fastq)-centroids.fasta" 2>/dev/null || { echo "ERROR: Clustering failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
        echo "  ✓ Clustering completed"
        clustering_input="${output_dir}/clusters"
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 4: Clustering skipped (enable_clustering=false)"
    else
        echo "✓ Step 3: Clustering skipped (enable_clustering=false)"
    fi
    clustering_input="${output_dir}/nanofilt"
fi
echo ""

#############################################
# Convert FASTQ to FASTA (if no clustering) #
#############################################

if [[ "${enable_clustering}" != "true" ]] && [[ "${enable_clustering}" != "True" ]] && [[ "${enable_clustering}" != "TRUE" ]]; then
    if [[ -n "$(find "${output_dir}/fasta" -maxdepth 1 -type f -name '*.fasta' 2>/dev/null)" ]]; then
        if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
            echo "✓ Step 5: Skipping FASTQ to FASTA conversion (output already exists)"
        else
            echo "✓ Step 4: Skipping FASTQ to FASTA conversion (output already exists)"
        fi
        fasta_input="${output_dir}/fasta"
    else
        if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
            echo "► Step 5: Converting FASTQ to FASTA"
        else
            echo "► Step 4: Converting FASTQ to FASTA"
        fi
        
        mkdir -p "${output_dir}/fasta"
        time find "${output_dir}/nanofilt" -maxdepth 1 -type f -name '*.fastq' -print0 | \
          xargs -0 -I {} bash -c 'seqtk seq -A "$1" > "${2}/fasta/$(basename "$1" .fastq).fasta" || { echo "ERROR: FASTQ to FASTA conversion failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
        echo "  ✓ FASTQ to FASTA conversion completed"
        fasta_input="${output_dir}/fasta"
    fi
    echo ""
else
    fasta_input="${output_dir}/clusters"
fi
echo ""

###################
# Chimera Removal #
###################

if [[ -n "$(find "${output_dir}/chimera" -maxdepth 1 -type f -name '*.fasta' 2>/dev/null)" ]]; then
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 6: Skipping chimera removal (output already exists)"
    else
        echo "✓ Step 5: Skipping chimera removal (output already exists)"
    fi
    mapping_input="${output_dir}/chimera"
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "► Step 6: Chimera removal"
    else
        echo "► Step 5: Chimera removal"
    fi
    echo "  Using VSEARCH uchime_denovo for chimera detection"
    
    mkdir -p "${output_dir}/chimera"
    
    # Calculate number of FASTA files for parallel processing
    # uchime_denovo is single-threaded, so parallelize by file
    num_files=$(find "${fasta_input}" -maxdepth 1 -type f -name '*.fasta' | wc -l)
    
    # Limit parallel jobs to avoid I/O bottleneck (max 6)
    parallel_jobs=$((threads / 4))
    if [[ ${parallel_jobs} -lt 1 ]]; then
      parallel_jobs=1
    fi
    if [[ ${parallel_jobs} -gt 6 ]]; then
      parallel_jobs=6
    fi
    if [[ ${parallel_jobs} -gt ${num_files} ]]; then
      parallel_jobs=${num_files}
    fi
    
    echo "  Processing ${num_files} FASTA file(s) with ${parallel_jobs} parallel jobs"
    
    time find "${fasta_input}" -maxdepth 1 -type f -name '*.fasta' -print0 | \
      xargs -0 -I {} -P "${parallel_jobs}" bash -c 'vsearch --uchime_denovo "$1" --nonchimeras "${2}/chimera/$(basename "$1" .fasta).fasta" 2>/dev/null || { echo "WARNING: VSEARCH chimera detection failed for $1, copying original" >&2; cp "$1" "${2}/chimera/$(basename "$1")"; }' _ {} "${output_dir}" || true
    
    echo "  ✓ Chimera removal completed"
    mapping_input="${output_dir}/chimera"
fi
echo ""

##########################
# Read Mapping (minimap2)#
##########################

if [[ -n "$(find "${output_dir}/mapped" -maxdepth 1 -type f -name '*.paf' 2>/dev/null)" ]]; then
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 7: Skipping mapping (output already exists)"
    else
        echo "✓ Step 6: Skipping mapping (output already exists)"
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "► Step 7: Read mapping with minimap2"
    else
        echo "► Step 6: Read mapping with minimap2"
    fi
    echo "  Preset: map-ont (optimized for long reads ~800bp metabarcodes)"
    echo "  Reference: ${ref_db}"
    
    time find "${mapping_input}" -maxdepth 1 -type f -name '*.fasta' -print0 | \
      xargs -0 -I {} bash -c 'minimap2 -x map-ont -Q -t "$3" --secondary=no -K 10M "$2" "$1" > "${4}/mapped/$(basename "$1" .fasta).paf" 2>/dev/null || { echo "ERROR: Mapping failed for $1" >&2; exit 1; }' _ {} "${ref_db}" "${threads}" "${output_dir}" || exit 1
    echo "  ✓ Mapping completed"
fi
echo ""

#############
# Filtering #
#############

if [[ -f "${output_dir}/filteredPAFs/filtered_otu.tsv" ]]; then
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 8: Skipping PAF filtering with multi-filter (output already exists)"
    else
        echo "✓ Step 7: Skipping PAF filtering with multi-filter (output already exists)"
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "► Step 8: Filtering PAF files with multi-filter"
    else
        echo "► Step 7: Filtering PAF files with multi-filter"
    fi
    echo "  Script: ${filter_script}"
    echo "  Parameters: -C ${confidence} -v ${qcov} -m ${min_qlen} -k ${identity_threshold} -L ${lca_type}"
    
    if [[ ! -f "${filter_script}" ]]; then
      echo "ERROR: Filter script not found at ${filter_script}" >&2
      exit 1
    fi
    
    mkdir -p "${output_dir}/filteredPAFs"
    time "${filter_script}" -i "${output_dir}/mapped" -T "${taxonomy_db}" -o "${output_dir}/filteredPAFs/filtered_otu.tsv" -C "${confidence}" -v "${qcov}" -m "${min_qlen}" -k "${identity_threshold}" --lca-type "${lca_type}" || { echo "ERROR: PAF filtering with multi-filter failed" >&2; exit 1; }
    echo "  ✓ PAF filtering completed"
fi
echo ""

#####################################################
# Apply Species Clustering (NEW V05 FEATURE)       #
#####################################################

if [[ -n "${clusters_file}" ]] && [[ -f "${clusters_file}" ]]; then
    if [[ -f "${output_dir}/filteredPAFs/clustered_otu.tsv" ]]; then
        if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
            echo "✓ Step 8a: Skipping species clustering (output already exists)"
        else
            echo "✓ Step 8: Skipping species clustering (output already exists)"
        fi
        clustering_otu_file="${output_dir}/filteredPAFs/clustered_otu.tsv"
    else
        if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
            echo "► Step 8a: Applying species clustering"
        else
            echo "► Step 8: Applying species clustering"
        fi
        echo "  Clusters database: ${clusters_file}"
        echo "  Confidence thresholds: high=${confidence} (≥0.85), medium=0.5"
        
        if [[ ! -f "${output_dir}/filteredPAFs/filtered_otu.tsv" ]]; then
            echo "ERROR: Filtered OTU file not found: ${output_dir}/filteredPAFs/filtered_otu.tsv" >&2
            exit 1
        fi
        
        # Find the apply_species_clustering.py script
        # Assume it's in the same directory as APOGEE_V05.sh or check environment
        clustering_script="apply_species_clustering.py"
        if [[ ! -f "${clustering_script}" ]]; then
            clustering_script="$(dirname "${BASH_SOURCE[0]}")/apply_species_clustering.py"
        fi
        
        if [[ ! -f "${clustering_script}" ]]; then
            echo "ERROR: apply_species_clustering.py not found. Cannot apply species clustering." >&2
            echo "  Searched: ./apply_species_clustering.py and $(dirname "${BASH_SOURCE[0]}")/apply_species_clustering.py" >&2
            echo "  Proceeding without clustering (using filtered OTU table)" >&2
            clustering_otu_file="${output_dir}/filteredPAFs/filtered_otu.tsv"
        else
            # Merge taxonomy with OTU table for clustering
            # (clustering script needs Species column + confidence)
            echo "  Preparing taxonomy+confidence file for clustering..."
            
            merged_tax_file="${output_dir}/filteredPAFs/taxonomy_with_confidence.csv"
            otu_input="${output_dir}/filteredPAFs/filtered_otu.tsv"
            tax_input="${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu.csv"
            
            OTU_INPUT="${otu_input}" TAX_INPUT="${tax_input}" MERGED_TAX_FILE="${merged_tax_file}" python3 << 'EOFPY'
import pandas as pd
import sys
import os

try:
    otu_input = os.environ['OTU_INPUT']
    tax_input = os.environ['TAX_INPUT']
    merged_tax_file = os.environ['MERGED_TAX_FILE']
    
    # Read OTU table with confidences
    otu_table = pd.read_csv(otu_input, sep='\t')
    
    # Read taxonomy file (skip comment line)
    with open(tax_input) as f:
        header = f.readline().strip().replace('#', '').split(',')
    
    taxonomy = pd.read_csv(
        tax_input,
        header=None,
        skiprows=1
    )
    taxonomy.columns = header
    
    # Get confidence column (handle different possible names)
    conf_col = 'MappingConfidence' if 'MappingConfidence' in otu_table.columns else 'AvgConfidence'
    
    # Rename confidence column for consistency
    otu_table_copy = otu_table.copy()
    otu_table_copy['confidence'] = otu_table_copy[conf_col]
    
    # Normalize confidence to 0-1 if needed (in case it's 0-100)
    if otu_table_copy['confidence'].max() > 1:
        otu_table_copy['confidence'] = otu_table_copy['confidence'] / 100.0
    
    # Merge OTU table with taxonomy
    merged = otu_table_copy.merge(
        taxonomy,
        left_on='otu',
        right_on='OTU ID',
        how='inner'
    )
    
    # Keep only essential columns:
    # - OTU ID (otu or OTU ID from merge)
    # - Sample columns (mock3-nanofilt, mock4-nanofilt, mock5-nanofilt, etc)
    # - confidence (for clustering logic)
    # - Taxonomy columns (Kingdom, Phylum, Class, Order, Family, Genus, Species)
    
    # Identify sample columns (all numeric columns except otu, OTU ID, and confidence)
    sample_cols = []
    for col in otu_table.columns:
        if col not in ['otu', 'OTU ID', 'confidence', conf_col] and col != 'OTU ID':
            # Check if this looks like a sample name (not metadata like TotalCount, NumHits, etc)
            if col.lower() not in ['totalcount', 'numhits', 'avgconfidence', 'medianqcov', 'avgidentity', 'medianconfidence']:
                sample_cols.append(col)
    
    # Build list of columns to keep
    keep_cols = ['otu'] + sample_cols + ['confidence'] + [c for c in taxonomy.columns if c not in ['OTU ID', 'Method', 'LCADepth']]
    
    # Select only those columns
    merged = merged[keep_cols]
    
    # Save (OTU as first column, index=False to keep it as column)
    merged.to_csv(merged_tax_file, index=False)
    print(f"Merged {len(merged)} OTUs with taxonomies and confidence scores")
    
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
EOFPY
            
            if [[ ! -f "${merged_tax_file}" ]]; then
                echo "ERROR: Failed to create merged taxonomy file" >&2
                exit 1
            fi
            
            # Apply species clustering
            time python3 "${clustering_script}" \
                -i "${merged_tax_file}" \
                -c "${clusters_file}" \
                -o "${output_dir}/filteredPAFs/clustered_otu.tsv" \
                --high-conf 0.85 \
                --medium-conf 0.5 \
                --confidence-col confidence \
                --species-col Species \
                --log-changes || { echo "ERROR: Species clustering failed" >&2; exit 1; }
            echo "  ✓ Species clustering completed"
            clustering_otu_file="${output_dir}/filteredPAFs/clustered_otu.tsv"
            
            # Show clustering statistics
            if [[ -f "${output_dir}/filteredPAFs/clustered_otu_clustering_changes.txt" ]]; then
                changes=$(wc -l < "${output_dir}/filteredPAFs/clustered_otu_clustering_changes.txt")
                echo "  Assignments modified: $((changes - 1)) OTUs"
            fi
            
            # Use the taxonomy file created by apply_species_clustering.py
            # which already contains the clustered species assignments
            if [[ -f "${output_dir}/filteredPAFs/clustered_otu_taxonomy.tsv" ]]; then
                echo "  Converting clustered taxonomy to CSV format..."
                
                TSV_FILE="${output_dir}/filteredPAFs/clustered_otu_taxonomy.tsv" \
                CSV_FILE="${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu_clustered.csv" \
                python3 << 'EOFCONVERT'
import pandas as pd
import os
import sys

try:
    # Read TSV file
    tsv_file = os.environ['TSV_FILE']
    csv_file = os.environ['CSV_FILE']
    
    df = pd.read_csv(tsv_file, sep='\t', index_col=0)
    df.to_csv(csv_file)
    print(f"Converted {len(df)} OTUs to CSV format", file=sys.stderr)
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
EOFCONVERT
                
                if [[ -f "${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu_clustered.csv" ]]; then
                    echo "  ✓ Clustered taxonomy file created"
                else
                    echo "WARNING: Failed to create clustered taxonomy file"
                fi
            else
                echo "WARNING: clustered_otu_taxonomy.tsv not found from apply_species_clustering.py"
            fi
        fi
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 8a: Species clustering skipped (-K not provided or file not found)"
    else
        echo "✓ Step 8: Species clustering skipped (-K not provided or file not found)"
    fi
    clustering_otu_file="${output_dir}/filteredPAFs/filtered_otu.tsv"
fi
echo ""

#################################################
# Collapse OTU Table by Taxonomic Assignment   #
#################################################

if [[ -f "${output_dir}/collapsed_otu_table.csv" ]]; then
    echo "✓ Step 8b: Skipping OTU collapse (output already exists)"
else
    echo "► Step 8b: Collapsing OTU by taxonomic assignment"
    echo "  Script: ${collapse_script}"
    echo "  Using input: $(basename "${clustering_otu_file}")"
    echo "  Grouping identical taxonomies into single OTU"
    
    if [[ ! -f "${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu.csv" ]]; then
      echo "ERROR: Taxonomy file not found at ${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu.csv" >&2
      exit 1
    fi
    
    if [[ ! -f "${clustering_otu_file}" ]]; then
      echo "ERROR: OTU file not found at ${clustering_otu_file}" >&2
      exit 1
    fi
    
    # Determine which taxonomy file to use
    # If species clustering was applied, use the clustered taxonomy file instead of original
    taxonomy_file_for_collapse="${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu.csv"
    
    if [[ "${clustering_otu_file}" == "${output_dir}/filteredPAFs/clustered_otu.tsv" ]] && [[ -f "${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu_clustered.csv" ]]; then
        echo "  Using clustered taxonomy file (reflecting species clustering changes)"
        taxonomy_file_for_collapse="${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu_clustered.csv"
    fi
    
    time "${collapse_script}" \
        -i "${clustering_otu_file}" \
        -t "${taxonomy_file_for_collapse}" \
        -o "${output_dir}/collapsed_otu_table.csv" \
        -x "${output_dir}/collapsed_taxonomy.csv" || { echo "ERROR: OTU collapse failed" >&2; exit 1; }
    echo "  ✓ OTU collapse completed"
    
    # Get collapse statistics
    original_count=$(wc -l < "${clustering_otu_file}")
    collapsed_count=$(wc -l < "${output_dir}/collapsed_otu_table.csv")
    echo "  Original OTU: $((original_count - 1))"
    echo "  Collapsed OTU: $((collapsed_count - 1))"
fi
echo ""



#####################################
# Create Summary Taxonomy Table     #
#####################################

if [[ -f "${output_dir}/summary.csv" ]]; then
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 9: Skipping summary table creation (output already exists)"
    else
        echo "✓ Step 8c: Skipping summary table creation (output already exists)"
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "► Step 9: Creating summary table from collapsed OTU"
    else
        echo "► Step 8c: Creating summary table from collapsed OTU"
    fi
    echo "  Script: ${taxonomy_script}"
    echo "  Parameters: -C ${confidence} -v ${qcov} -m ${min_qlen} -k ${identity_threshold} -H ${top_species_count}"
    
    # Use collapsed OTU table for species proportions
    if [[ ! -f "${output_dir}/collapsed_otu_table.csv" ]]; then
      echo "ERROR: Collapsed OTU table not found at ${output_dir}/collapsed_otu_table.csv" >&2
      exit 1
    fi
    
    # Determine which taxonomy database to use
    # If species clustering was applied and clustered taxonomy exists, use that for consistency
    summary_taxonomy_db="${taxonomy_db}"
    
    if [[ -f "${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu_clustered.csv" ]]; then
        echo "  Using clustered taxonomy for summary (reflecting species clustering changes)"
        summary_taxonomy_db="${output_dir}/filteredPAFs/phyloseq_taxonomy_filtered_otu_clustered.csv"
    fi
    
    "${taxonomy_script}" -i "${output_dir}/collapsed_otu_table.csv" -T "${summary_taxonomy_db}" -C "${confidence}" -v "${qcov}" -m "${min_qlen}" -k "${identity_threshold}" -H "${top_species_count}" > "${output_dir}/summary.csv" || { echo "ERROR: Summary table creation failed" >&2; exit 1; }
    echo "  ✓ Summary table created: ${output_dir}/summary.csv"
fi
echo ""
echo "==============================================="
echo "Pipeline completed successfully!"
echo "==============================================="
echo "Output files:"
echo "  Collapsed OTU table:   ${output_dir}/collapsed_otu_table.csv"
echo "  Collapsed taxonomy:    ${output_dir}/collapsed_taxonomy.csv"
echo "  Summary (species):     ${output_dir}/summary.csv"
echo "==============================================="

# Disable trap and exit successfully
trap - EXIT
exit 0
