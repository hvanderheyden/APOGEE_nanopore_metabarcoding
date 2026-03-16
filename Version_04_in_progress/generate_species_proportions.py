#!/usr/bin/env python3
"""
Generate species proportions from OTU table + taxonomy
"""
import pandas as pd
import sys
import os
from pathlib import Path

def generate_proportions(otu_file, tax_file, output_file=None, sep='\t'):
    """
    Generate species proportions file from OTU and taxonomy tables
    
    Parameters:
    - otu_file: OTU table (samples as columns, OTU as rows) - CSV format
    - tax_file: Taxonomy database TSV (QIIME format - legacy, may not be used)
    - output_file: Output proportions file (None = stdout)
    - sep: Separator char (not used - format auto-detected)
    
    NOTE: This function looks first for a phyloseq_taxonomy file from filter_with_advanced_lca
    which contains pre-computed species assignments organized by rank.
    """
    
    print(f"Loading OTU table: {otu_file}", file=sys.stderr)
    otu = pd.read_csv(otu_file, sep=',', index_col=0)  # OTU table is CSV
    
    # Get sample columns (exclude metadata columns)
    sample_cols = [c for c in otu.columns if 'nanofilt' in c]
    print(f"Found sample columns: {sample_cols}", file=sys.stderr)
    
    # Try to load pre-collapsed taxonomy from collapse_otu_by_taxonomy if available
    # Otherwise fall back to pre-filtered phyloseq taxonomy from filter_with_advanced_lca
    otu_dir = os.path.dirname(otu_file)
    collapsed_tax_file = os.path.join(otu_dir, 'collapsed_taxonomy.csv')
    phyloseq_file = os.path.join(otu_dir, 'filteredPAFs', 'phyloseq_taxonomy_filtered_otu.csv')
    
    if not os.path.exists(phyloseq_file):
        # Try alternate location (directly in otu_dir)
        phyloseq_file = os.path.join(otu_dir, 'phyloseq_taxonomy_filtered_otu.csv')
    
    # PRIORITY 1: Use collapsed taxonomy if available (matches collapsed OTU table)
    if os.path.exists(collapsed_tax_file):
        print(f"Loading collapsed taxonomy: {collapsed_tax_file}", file=sys.stderr)
        collapsed_tax = pd.read_csv(collapsed_tax_file, sep=',', index_col=0)
        # Remove the '#' from the OTU ID column name if present
        if collapsed_tax.index.name and collapsed_tax.index.name.startswith('#'):
            collapsed_tax.index.name = collapsed_tax.index.name.lstrip('#')
        otu_taxa = collapsed_tax[['Species']].copy()
    # PRIORITY 2: Use pre-filtered phyloseq taxonomy (for backward compat with filtered OTU)
    elif os.path.exists(phyloseq_file):
        print(f"Loading pre-filtered phyloseq taxonomy: {phyloseq_file}", file=sys.stderr)
        # This file contains: #OTU ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, LCADepth, Method
        # The # at the start indicates comment, but it's the header, so we skip it
        phyloseq = pd.read_csv(phyloseq_file, sep=',', index_col=0)
        # Remove the '#' from the OTU ID column name if present
        if phyloseq.index.name and phyloseq.index.name.startswith('#'):
            phyloseq.index.name = phyloseq.index.name.lstrip('#')
        otu_taxa = phyloseq[['Species']].copy()
    else:
        print(f"Pre-filtered phyloseq not found. Loading raw taxonomy database: {tax_file}", file=sys.stderr)
        
        # Load taxonomy database
        tax = pd.read_csv(tax_file, sep='\t', index_col=0)
        
        # Extract species-level classification from QIIME taxonomy format
        # Format: k__Eukaryota;p__Ascomycota;...;s__species_name
        def extract_species(taxon_str):
            """Extract species from QIIME taxonomy string"""
            if pd.isna(taxon_str):
                return 'unclassified'
            parts = str(taxon_str).split(';')
            for part in reversed(parts):  # Start from the most specific level
                part = part.strip()
                if part.startswith('s__'):
                    species = part.replace('s__', '').strip('[]')
                    return species if species else 'unclassified'
            # If no species found, use genus
            for part in reversed(parts):
                part = part.strip()
                if part.startswith('g__'):
                    genus = part.replace('g__', '').strip('[]')
                    return genus if genus else 'unclassified'
            return 'unclassified'
        
        # Add Species column
        tax['Species'] = tax['Taxon'].apply(extract_species)
        otu_taxa = tax[['Species']].copy()
    
    # Merge OTU with taxonomy
    otu_with_species = otu[sample_cols].reset_index()
    otu_with_species.columns = ['OTU'] + sample_cols
    
    # Join with taxonomy
    otu_with_species = otu_with_species.set_index('OTU').join(otu_taxa, how='left')
    otu_with_species = otu_with_species.reset_index()
    otu_with_species['Species'] = otu_with_species['Species'].fillna('unclassified')
    
    # Group by species and sum reads
    grouped = otu_with_species.groupby('Species')[sample_cols].sum()
    
    # Calculate total reads per sample
    totals = grouped.sum()
    
    # Calculate percentages
    proportions = grouped.copy()
    for col in sample_cols:
        proportions[f'{col}_pct'] = 100 * grouped[col] / totals[col]
    
    # Build output: Species | Mock_Reads | Mock_Percentage | etc
    output = pd.DataFrame()
    output['Species'] = grouped.index
    
    for col in sample_cols:
        # Rename column to match output format
        col_name = col.replace('-nanofilt', '').upper()
        output[f'{col_name}_Reads'] = grouped[col].values
        output[f'{col_name}_Percentage'] = proportions[f'{col}_pct'].values
    
    # Sort by total reads descending
    total_reads = output[[c for c in output.columns if 'Reads' in c]].sum(axis=1)
    output = output.iloc[total_reads.argsort()[::-1]]
    
    # Save or print to stdout
    if output_file:
        print(f"Writing to: {output_file}", file=sys.stderr)
        output.to_csv(output_file, index=False)
    else:
        # Write to stdout
        output.to_csv(sys.stdout, index=False)
    
    # Print summary
    print(f"\nSummary:", file=sys.stderr)
    for col in sample_cols:
        col_name = col.replace('-nanofilt', '').upper()
        total = output[f'{col_name}_Reads'].sum()
        unclass = output[output['Species'] == 'unclassified'][f'{col_name}_Reads'].values
        unclass_pct = output[output['Species'] == 'unclassified'][f'{col_name}_Percentage'].values
        
        print(f"  {col}:", file=sys.stderr)
        print(f"    Total reads: {int(total):,d}", file=sys.stderr)
        if len(unclass) > 0:
            print(f"    Unclassified: {int(unclass[0]):,d} ({unclass_pct[0]:.1f}%)", file=sys.stderr)
        print(f"    Species count: {len(output[output['Species'] != 'unclassified'])}", file=sys.stderr)
    
    return output

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate species proportions from OTU table')
    parser.add_argument('-i', '--otu-file', required=True, help='OTU table file (CSV)')
    parser.add_argument('-T', '--tax-file', required=True, help='Taxonomy database file (TSV)')
    parser.add_argument('-o', '--output', required=False, help='Output file (optional, defaults to stdout)')
    parser.add_argument('-C', '--confidence', type=float, default=1.0, help='Confidence threshold (not used, for compatibility)')
    parser.add_argument('-v', '--qcov', type=float, default=0.85, help='Query coverage (not used, for compatibility)')
    parser.add_argument('-m', '--min-qlen', type=int, default=0, help='Min query length (not used, for compatibility)')
    parser.add_argument('-k', '--identity', type=float, default=0.9, help='Identity threshold (not used, for compatibility)')
    
    args = parser.parse_args()
    
    otu_file = args.otu_file
    tax_file = args.tax_file
    output_file = args.output
    
    generate_proportions(otu_file, tax_file, output_file)
