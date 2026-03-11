#!/bin/bash
set -euo pipefail

trap 'echo "ERROR: Script failed at line $LINENO" >&2' EXIT

###################
## Example usage ##
###################

#chmod +x path/to/APOGEE_V03.sh

# /path/to/APOGEE_V03.sh \
# -i /path/to/fastq \
# -o /path/to/outputfolder \
# -r /path/to/ITS-RefDB_V02_2_fixed.mmi \
# -t 24 \
# -b true \
# -R TCCGTAGGTGAACCTGCGG \
# -P TCCTCCGCTTATTGATATGC \
# -c false \
# -x 0.99 \
# -w 8 \
# -T /path/to/taxonomy.tsv \
# -F /path/to/filter_with_confidence.py \
# -S /path/to/taxonomyTable.py \
# -C 0.5 \
# -v 0.7 \
# -k 0.9

#############
# ARGUMENTS #
#############

# Parse command-line arguments
while getopts ":i:o:r:t:c:x:w:T:F:S:C:v:k:b:R:P:" opt; do
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
    k ) identity_threshold="$OPTARG" ;;
    b ) enable_trimming="$OPTARG" ;;
    R ) reverse_primer="$OPTARG" ;;
    P ) forward_primer="$OPTARG" ;;
    \? ) echo "ERROR: Invalid option: -$OPTARG" 1>&2; exit 1 ;;
    : ) echo "ERROR: Invalid option: -$OPTARG requires an argument" 1>&2; exit 1 ;;
  esac
done
shift $((OPTIND -1))

########################
# Validate Arguments   #
########################

if [[ -z "${input_dir:-}" ]] || [[ -z "${output_dir:-}" ]] || [[ -z "${ref_db:-}" ]] || [[ -z "${threads:-}" ]] || [[ -z "${taxonomy_db:-}" ]] || [[ -z "${filter_script:-}" ]] || [[ -z "${taxonomy_script:-}" ]]; then
  echo "ERROR: Missing required arguments" >&2
  echo "Usage: $0 -i <input_dir> -o <output_dir> -r <ref_db> -t <threads> -T <taxonomy_db> -F <filter_script> -S <taxonomy_script> [-c <enable_clustering>] [-x <identity>] [-w <wordlength>] [-C <confidence>] [-v <qcov>] [-k <identity_threshold>] [-b <enable_trimming>] [-R <reverse_primer>] [-P <forward_primer>]" >&2
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

# Set defaults if not provided
identity="${identity:-0.97}"
wordlength="${wordlength:-10}"
enable_clustering="${enable_clustering:-false}"
enable_trimming="${enable_trimming:-false}"
confidence="${confidence:-0.3}"
qcov="${qcov:-0.0}"
identity_threshold="${identity_threshold:-0.0}"
reverse_primer="${reverse_primer:-}"
forward_primer="${forward_primer:-}"

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

if ! [[ "${identity_threshold}" =~ ^[0-9]+\.?[0-9]*$ ]] || (( $(echo "${identity_threshold} < 0 || ${identity_threshold} > 1" | bc -l) )); then
  echo "ERROR: identity_threshold must be between 0 and 1" >&2
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
echo "Metabarcoding Pipeline v03 - Multi-filter with Optional Trimming"
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
echo "Identity threshold:  ${identity_threshold}"
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
    echo "  Parameters: -C ${confidence} -v ${qcov} -k ${identity_threshold}"
    
    if [[ ! -f "${filter_script}" ]]; then
      echo "ERROR: Filter script not found at ${filter_script}" >&2
      exit 1
    fi
    
    mkdir -p "${output_dir}/filteredPAFs"
    time "${filter_script}" -i "${output_dir}/mapped" -T "${taxonomy_db}" -o "${output_dir}/filteredPAFs/filtered_otu.tsv" -C "${confidence}" -v "${qcov}" -k "${identity_threshold}" || { echo "ERROR: PAF filtering with multi-filter failed" >&2; exit 1; }
    echo "  ✓ PAF filtering completed"
fi
echo ""

#######################################
# Create OTU Table from Filtered PAFs #
#######################################

if [[ -f "${output_dir}/otu_table.csv" ]]; then
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 9: Skipping OTU table creation (output already exists)"
    else
        echo "✓ Step 8: Skipping OTU table creation (output already exists)"
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "► Step 9: Creating OTU table from filtered PAFs"
    else
        echo "► Step 8: Creating OTU table from filtered PAFs"
    fi
    echo "  Input: ${output_dir}/filteredPAFs/filtered_otu.tsv"
    
    if [[ ! -f "${output_dir}/filteredPAFs/filtered_otu.tsv" ]]; then
      echo "ERROR: Filtered OTU file not found at ${output_dir}/filteredPAFs/filtered_otu.tsv" >&2
      exit 1
    fi
    
    # Convert TSV to CSV format - keep all columns including confidence scores
    awk 'BEGIN {FS="\t"; OFS=","} {$1=$1; print}' "${output_dir}/filteredPAFs/filtered_otu.tsv" > "${output_dir}/otu_table.csv" 2>/dev/null || { echo "ERROR: OTU table conversion failed" >&2; exit 1; }
    echo "  ✓ OTU table created: ${output_dir}/otu_table.csv"
fi
echo ""

#####################################
# Create Phyloseq Taxonomy Table    #
#####################################

if [[ -f "${output_dir}/phyloseq_taxonomy.csv" ]]; then
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "✓ Step 10: Skipping taxonomy table creation (output already exists)"
    else
        echo "✓ Step 9: Skipping taxonomy table creation (output already exists)"
    fi
else
    if [[ "${enable_trimming}" == "true" ]] || [[ "${enable_trimming}" == "True" ]] || [[ "${enable_trimming}" == "TRUE" ]]; then
        echo "► Step 10: Creating taxonomy table with progressive filtering"
    else
        echo "► Step 9: Creating taxonomy table with progressive filtering"
    fi
    echo "  Script: ${taxonomy_script}"
    echo "  Parameters: -C ${confidence} -v ${qcov} -k ${identity_threshold}"
    
    "${taxonomy_script}" -i "${output_dir}/otu_table.csv" -T "${taxonomy_db}" -C "${confidence}" -v "${qcov}" -k "${identity_threshold}" > "${output_dir}/phyloseq_taxonomy.csv" 2>/dev/null || { echo "ERROR: Taxonomy table creation failed" >&2; exit 1; }
    echo "  ✓ Taxonomy table created: ${output_dir}/phyloseq_taxonomy.csv"
fi
echo ""
echo "==============================================="
echo "Pipeline completed successfully!"
echo "==============================================="
echo "Output files:"
echo "  OTU table:        ${output_dir}/otu_table.csv"
echo "  Taxonomy table:   ${output_dir}/phyloseq_taxonomy.csv"
echo "==============================================="

# Disable trap and exit successfully
trap - EXIT
exit 0
