#!/usr/bin/env python3
"""
Collapse OTU table by taxonomic assignment to make it comparable to Kraken2.

Instead of 1 OTU per sequence, groups all sequences with identical taxonomic
assignment into a single OTU. Creates two output files:
1. Collapsed OTU table (read counts per taxonomic assignment)
2. Collapsed taxonomy file (same structure as phyloseq_taxonomy format)

Usage:
    python collapse_otu_by_taxonomy.py \
        -i otu_table.csv \
        -t phyloseq_taxonomy.csv \
        -o collapsed_otu_table.csv \
        -x collapsed_taxonomy.csv
"""

import argparse
import os
import sys
from collections import defaultdict
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input_otu', required=True, 
                        help='Input OTU table (CSV from pipeline)')
    parser.add_argument('-t', '--taxonomy', required=True,
                        help='Taxonomy file (phyloseq_taxonomy.csv with full taxonomy columns)')
    parser.add_argument('-o', '--output_otu', required=True,
                        help='Output collapsed OTU table (CSV)')
    parser.add_argument('-x', '--output_tax', required=True,
                        help='Output collapsed taxonomy file (CSV, same structure as input)')
    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.input_otu):
        print(f"ERROR: Input OTU file not found: {args.input_otu}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.taxonomy):
        print(f"ERROR: Taxonomy file not found: {args.taxonomy}", file=sys.stderr)
        sys.exit(1)

    try:
        # Read taxonomy file (with full structure: #OTU ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, LCADepth, Method)
        print(f"Loading taxonomy from {args.taxonomy}...", file=sys.stderr)
        
        # Read first line to get column names (handle #OTU ID format)
        with open(args.taxonomy, 'r') as f:
            header_line = f.readline().strip()
        
        # Remove # if present at start of column name
        if header_line.startswith('#'):
            header_line = header_line[1:]
        
        col_names = [c.strip() for c in header_line.split(',')]
        print(f"Column names: {col_names}", file=sys.stderr)
        
        # Read the CSV with proper column names
        tax_df = pd.read_csv(args.taxonomy, comment='#', names=col_names)
        
        # Get taxonomy columns (all except first which is OTU ID)
        otu_id_col = col_names[0]
        tax_cols = col_names[1:]
        
        print(f"Taxonomy columns found: {tax_cols}", file=sys.stderr)
        
        # Create mapping: OTU_ID -> full taxonomy dict
        otu_to_tax_dict = {}
        for _, row in tax_df.iterrows():
            otu_id = str(row[otu_id_col])
            tax_dict = {col: row[col] for col in tax_cols}
            otu_to_tax_dict[otu_id] = tax_dict
        
        print(f"Loaded {len(otu_to_tax_dict)} taxonomic assignments", file=sys.stderr)

        # Read OTU table (auto-detect separator: TSV or CSV)
        print(f"Loading OTU table from {args.input_otu}...", file=sys.stderr)
        # Auto-detect separator based on file extension, default to comma
        sep = '\t' if args.input_otu.endswith('.tsv') else ','
        otu_df = pd.read_csv(args.input_otu, sep=sep)
        
        # Get sample columns (all except first/ID column)
        id_col = otu_df.columns[0]
        sample_cols = [col for col in otu_df.columns[1:] if col != id_col]
        
        print(f"Found {len(otu_df)} OTU with {len(sample_cols)} samples", file=sys.stderr)

        # Group by taxonomy
        # Key: tuple of taxonomy values, Value: {samples: counts, tax_dict: taxonomy_info}
        tax_grouped = defaultdict(lambda: {'counts': {col: 0 for col in sample_cols}, 'tax_dict': None})
        
        for _, row in otu_df.iterrows():
            otu_id = str(row[id_col])
            
            # Get taxonomy for this OTU
            if otu_id in otu_to_tax_dict:
                tax_dict = otu_to_tax_dict[otu_id]
                # Create taxonomy tuple as key (use all taxonomy values except Method/LCADepth for matching)
                tax_tuple = tuple((col, tax_dict.get(col, 'unassigned')) for col in [c for c in tax_cols if c not in ['Method', 'LCADepth']])
            else:
                # If no taxonomy, create unassigned entry
                tax_dict = {col: 'unassigned' for col in tax_cols}
                tax_tuple = tuple((col, 'unassigned') for col in tax_cols if col not in ['Method', 'LCADepth'])
            
            # Store taxonomy dict ONLY on first occurrence (avoid overwriting)
            if tax_grouped[tax_tuple]['tax_dict'] is None:
                tax_grouped[tax_tuple]['tax_dict'] = tax_dict
            
            # Add read counts to this taxonomy group
            for col in sample_cols:
                try:
                    tax_grouped[tax_tuple]['counts'][col] += int(row[col]) if pd.notna(row[col]) else 0
                except (ValueError, TypeError):
                    pass  # Skip non-numeric values
        
        print(f"Collapsed to {len(tax_grouped)} unique taxonomic assignments", file=sys.stderr)

        # Create output dataframes - ENSURE ONE ROW PER UNIQUE TAXONOMY
        otu_data_list = []
        tax_data_list = []
        
        # Process each unique taxonomy group ONCE
        # Use COTU_ prefix to avoid collision with filtered OTU IDs
        for otu_num, (tax_tuple, data) in enumerate(sorted(tax_grouped.items()), 1):
            new_otu_id = f"COTU_{otu_num}"
            
            # Build OTU table row with counts
            otu_row = {id_col: new_otu_id}
            for sample_col in sample_cols:
                otu_row[sample_col] = data['counts'].get(sample_col, 0)
            otu_data_list.append(otu_row)
            
            # Build taxonomy table row
            tax_row = {otu_id_col: new_otu_id}
            if data['tax_dict']:
                for tax_col in tax_cols:
                    tax_row[tax_col] = data['tax_dict'].get(tax_col, 'unassigned')
            else:
                for tax_col in tax_cols:
                    tax_row[tax_col] = 'unassigned'
            tax_data_list.append(tax_row)
        
        # Create dataframes with proper column order
        col_order = [id_col] + sample_cols
        otu_result_df = pd.DataFrame(otu_data_list)[col_order]
        
        # Preserve taxonomy column order from input
        tax_col_order = [otu_id_col] + tax_cols
        tax_result_df = pd.DataFrame(tax_data_list)[tax_col_order]
        
        # Save outputs
        otu_result_df.to_csv(args.output_otu, index=False)
        print(f"Saved collapsed OTU table to {args.output_otu}", file=sys.stderr)
        
        tax_result_df.to_csv(args.output_tax, index=False)
        print(f"Saved collapsed taxonomy to {args.output_tax}", file=sys.stderr)
        
        # Print summary statistics
        total_reads = otu_result_df[sample_cols].sum().sum()
        print(f"\nCollapse Statistics:", file=sys.stderr)
        print(f"  Original OTU count: {len(otu_df)}", file=sys.stderr)
        print(f"  Collapsed to: {len(otu_result_df)}", file=sys.stderr)
        print(f"  Reduction: {len(otu_df) - len(otu_result_df)} OTU ({100 * (1 - len(otu_result_df)/len(otu_df)):.1f}%)", file=sys.stderr)
        print(f"  Total reads: {int(total_reads)}", file=sys.stderr)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
