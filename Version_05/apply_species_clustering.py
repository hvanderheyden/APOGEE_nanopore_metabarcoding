#!/usr/bin/env python3
"""
Apply species clustering based on LCA confidence levels
Implements progressive taxonomy assignment:
  - confidence >= 0.85: Use exact LCA species
  - 0.5 <= confidence < 0.85: Use species cluster if match found
  - confidence < 0.5: Use LCA fallback to Genus/Family

Part of APOGEE V05 pipeline
"""

import sys
import os
import argparse
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(description='Apply species clustering to taxonomy assignments')
parser.add_argument('-i', '--input_file', required=True, help='Input taxonomy table (from filter_with_advanced_lca)')
parser.add_argument('-c', '--clusters_db', required=True, help='Species clusters database (TSV)')
parser.add_argument('-o', '--output_file', required=True, help='Output taxonomy table with clusters applied')
parser.add_argument('--confidence-col', default='confidence', help='Column name for confidence scores')
parser.add_argument('--species-col', default='Species', help='Column name for species assignments')
parser.add_argument('--high-conf', type=float, default=0.85, help='High confidence threshold (species)')
parser.add_argument('--medium-conf', type=float, default=0.5, help='Medium confidence threshold (clusters)')
parser.add_argument('--log-changes', action='store_true', help='Log which assignments were changed')

args = parser.parse_args()

# Load clusters database
print("Loading species clusters database...", file=sys.stderr)
clusters_map = defaultdict(lambda: None)  # species -> cluster_name

try:
    with open(args.clusters_db, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                cluster_num = int(parts[0])
                species_a = parts[1].strip().replace(' ', '_')  # Normalize spaces to underscores
                species_b = parts[2].strip().replace(' ', '_')  # Normalize spaces to underscores
                cluster_name = parts[3].strip().replace(' ', '_')  # Also normalize cluster name
                
                # Map both species to this cluster
                clusters_map[species_a] = cluster_name
                clusters_map[species_b] = cluster_name
except Exception as e:
    print(f"ERROR: Cannot read clusters database: {e}", file=sys.stderr)
    sys.exit(1)

print(f"Loaded {len(clusters_map)} unique species in clustering database", file=sys.stderr)

# Load taxonomy table
print(f"Loading taxonomy table: {args.input_file}", file=sys.stderr)
try:
    tax = pd.read_csv(args.input_file, sep=',')
    # Move first column to index if it's OTU ID or otu
    if tax.columns[0] == 'otu' or tax.columns[0] == 'OTU ID':
        tax = tax.set_index(tax.columns[0])
except Exception as e:
    print(f"ERROR: Cannot read taxonomy file: {e}", file=sys.stderr)
    sys.exit(1)

# Apply clustering strategy
print("Applying clustering strategy...", file=sys.stderr)

# Check if required columns exist
if args.species_col not in tax.columns:
    print(f"ERROR: Column '{args.species_col}' not found in input file", file=sys.stderr)
    sys.exit(1)

if args.confidence_col not in tax.columns:
    print(f"WARNING: Column '{args.confidence_col}' not found, using default confidence 0.9", file=sys.stderr)
    tax[args.confidence_col] = 0.9

# Create new column for clustered species
tax['Species_clustered'] = tax[args.species_col]
changes_log = []

for idx, row in tax.iterrows():
    confidence = float(row[args.confidence_col])
    current_species = str(row[args.species_col]).replace(' ', '_')  # Normalize spaces to underscores
    
    # Skip if already unclassified
    if 'unclassified' in current_species:
        continue
    
    # Decision logic
    if confidence >= args.high_conf:
        # Keep exact species
        action = "HIGH_CONF_EXACT"
    elif args.medium_conf <= confidence < args.high_conf:
        # Try to find cluster
        if current_species in clusters_map and clusters_map[current_species]:
            clustered_name = clusters_map[current_species]
            tax.at[idx, 'Species_clustered'] = clustered_name
            action = "MEDIUM_CONF_CLUSTER"
            changes_log.append({
                'OTU': idx,
                'original': str(row[args.species_col]),
                'clustered': clustered_name,
                'confidence': confidence,
                'action': action
            })
        else:
            # No cluster found, keep original
            action = "MEDIUM_CONF_NO_CLUSTER"
    else:
        # Low confidence - LCA fallback already in place from filter_with_advanced_lca
        action = "LOW_CONF_LCA_FALLBACK"

# SECOND PASS: Re-apply clustering to ensure consistency
# This guarantees all cluster members converge to the same cluster name
print(f"Applying second clustering pass for consistency...", file=sys.stderr)
second_pass_count = 0
for idx, row in tax.iterrows():
    # Check if the current Species_clustered value can be mapped to a cluster
    current_clustered = str(tax.at[idx, 'Species_clustered']).replace(' ', '_')
    
    if 'unclassified' in current_clustered:
        continue
    
    # If this value is also in the clusters map, update it
    if current_clustered in clusters_map and clusters_map[current_clustered]:
        new_cluster = clusters_map[current_clustered]
        if new_cluster != current_clustered:  # Only log if it actually changed
            tax.at[idx, 'Species_clustered'] = new_cluster
            second_pass_count += 1

if second_pass_count > 0:
    print(f"Second pass modified {second_pass_count} additional OTUs", file=sys.stderr)

# Output results
print(f"Applied clustering to {len(changes_log)} OTUs", file=sys.stderr)
print(f"Writing output to: {args.output_file}", file=sys.stderr)

# For output, keep only:
# - Sample columns (numeric columns like mock3-nanofilt, mock4-nanofilt, etc.)
# - Species_clustered
# Remove metadata and confidence columns

# Identify sample columns (numeric columns, not taxonomy/metadata)
sample_cols = []
for col in tax.columns:
    # Skip these metadata columns
    if col in ['Species', 'Species_clustered', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'LCADepth', 'Method', 'confidence', 'AvgConfidence', 'AvgIdentity', 'MedianQCov', 'NumHits', 'TotalCount', 'OTU ID']:
        continue
    # Keep sample columns (they're numeric and have sample names)
    sample_cols.append(col)

# Build output with only sample cols + Species_clustered
output_cols = sample_cols + ['Species_clustered']
output_df = tax[output_cols].copy()

# Determine separator based on output filename
sep = '\t' if args.output_file.endswith('.tsv') else ','
output_df.to_csv(args.output_file, sep=sep, index=True)

# ALSO create a clustered taxonomy file for downstream processing
# This replaces Species with Species_clustered and preserves all other taxonomy columns
clustered_taxonomy = tax[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species_clustered']].copy()
clustered_taxonomy.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# Get the basename and create corresponding taxonomy filename
if args.output_file.endswith('.tsv'):
    taxonomy_file = args.output_file.replace('.tsv', '_taxonomy.tsv')
else:
    taxonomy_file = args.output_file.replace('.csv', '_taxonomy.csv')

print(f"Writing clustered taxonomy to: {taxonomy_file}", file=sys.stderr)
clustered_taxonomy.to_csv(taxonomy_file, sep=sep, index=True, index_label='OTU ID')

# Log changes if requested
if args.log_changes and changes_log:
    # Replace .csv or .tsv with _clustering_changes.txt
    if args.output_file.endswith('.tsv'):
        log_file = args.output_file.replace('.tsv', '_clustering_changes.txt')
    else:
        log_file = args.output_file.replace('.csv', '_clustering_changes.txt')
    print(f"Writing change log to: {log_file}", file=sys.stderr)
    with open(log_file, 'w') as f:
        f.write("OTU_ID\tOriginal_Species\tClustered_Species\tConfidence\n")
        for change in changes_log:
            f.write(f"{change['OTU']}\t{change['original']}\t{change['clustered']}\t{change['confidence']}\n")

print("Done!", file=sys.stderr)
