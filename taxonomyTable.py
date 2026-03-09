#!/usr/bin/env python3
import pandas as pd
from argparse import ArgumentParser
import sys

def parse_arguments():
    """Parse command-line arguments"""
    parser = ArgumentParser(description="Generate phyloseq taxonomy table from OTU table and taxonomy database")
    parser.add_argument('-i', '--input', dest='otu_file', required=True,
                        help='Path to OTU table (CSV format)')
    parser.add_argument('-t', '--taxonomy', dest='tax_file', required=True,
                        help='Path to taxonomy database (TSV format: ID\\ttaxonomy_string)')
    parser.add_argument('-C', '--confidence', dest='min_conf', type=float, default=1.0,
                        help='Minimum mapping confidence threshold (0-1, default: 1.0)')
    parser.add_argument('-v', '--qcov', dest='min_qcov', type=float, default=0.0,
                        help='Minimum query coverage threshold (0-1, default: 0.0)')
    parser.add_argument('-p', '--identity', dest='min_identity', type=float, default=0.0,
                        help='Minimum identity threshold (0-1, default: 0.0)')
    
    args = parser.parse_args()
    return args

args = parse_arguments()
otu_file = args.otu_file
tax_file = args.tax_file
min_conf = args.min_conf
min_qcov = args.min_qcov
min_identity = args.min_identity

# Read OTU table (CSV format)
try:
    otu_df = pd.read_csv(otu_file, sep=',')
except Exception as e:
    print(f"ERROR: Cannot read OTU table: {e}", file=sys.stderr)
    sys.exit(1)

# Read taxonomy database (TSV format: ID\ttaxonomy_string)
tax_db = {}
try:
    with open(tax_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                seq_id = parts[0]
                taxonomy_str = parts[1]
                tax_db[seq_id] = taxonomy_str
except Exception as e:
    print(f"ERROR: Cannot read taxonomy database: {e}", file=sys.stderr)
    sys.exit(1)

def parse_taxonomy_string(tax_string):
    """
    Parse taxonomy string from semicolon-delimited format
    e.g. k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
    Returns dict with keys: kingdom, phylum, class, order, family, genus, species
    """
    tax_dict = {
        'kingdom': 'unknown',
        'phylum': 'unknown',
        'class': 'unknown',
        'order': 'unknown',
        'family': 'unknown',
        'genus': 'unknown',
        'species': 'unknown'
    }
    
    if not tax_string:
        return tax_dict
    
    # Split by semicolon
    parts = tax_string.split(';')
    for part in parts:
        if '__' not in part:
            continue
        prefix, name = part.split('__', 1)
        name = name.strip()
        
        if prefix == 'k':
            tax_dict['kingdom'] = name if name else 'unknown'
        elif prefix == 'p':
            tax_dict['phylum'] = name if name else 'unknown'
        elif prefix == 'c':
            tax_dict['class'] = name if name else 'unknown'
        elif prefix == 'o':
            tax_dict['order'] = name if name else 'unknown'
        elif prefix == 'f':
            tax_dict['family'] = name if name else 'unknown'
        elif prefix == 'g':
            tax_dict['genus'] = name if name else 'unknown'
        elif prefix == 's':
            tax_dict['species'] = name if name else 'unknown'
    
    return tax_dict

def get_filtered_taxonomy(accession, tax_dict, conf, qcov, identity):
    """
    Apply progressive filtering:
    - If all filters pass (conf >= min AND qcov >= min AND identity >= min) 
      → return full taxonomy (genus, species)
    - If any filter fails 
      → return genus as-is, species as 'unknown_<genus>'
    """
    all_pass = (conf >= min_conf and qcov >= min_qcov and identity >= min_identity)
    
    if all_pass:
        # All filters passed: return full taxonomy
        return [
            tax_dict['kingdom'],
            tax_dict['phylum'],
            tax_dict['class'],
            tax_dict['order'],
            tax_dict['family'],
            tax_dict['genus'],
            tax_dict['species']
        ]
    else:
        # At least one filter failed: keep genus, put unknown_<genus> in species
        genus = tax_dict['genus']
        return [
            tax_dict['kingdom'],
            tax_dict['phylum'],
            tax_dict['class'],
            tax_dict['order'],
            tax_dict['family'],
            genus,
            'unknown_' + genus if genus != 'unknown' else 'unknown'
        ]

# Prepare output
output_data = []
output_header = ['#OTU ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

try:
    for idx, row in otu_df.iterrows():
        accession = row['Accession']
        conf = row['MappingConfidence']
        qcov = row['MedianQCov']
        identity = row['AvgIdentity']
        
        # Get taxonomy from database
        if accession not in tax_db:
            # OTU not in database: unknown
            tax_dict = {
                'kingdom': 'unknown',
                'phylum': 'unknown',
                'class': 'unknown',
                'order': 'unknown',
                'family': 'unknown',
                'genus': 'unknown',
                'species': 'unknown'
            }
        else:
            tax_string = tax_db[accession]
            tax_dict = parse_taxonomy_string(tax_string)
        
        # Get filtered taxonomy based on thresholds
        filtered_tax = get_filtered_taxonomy(accession, tax_dict, conf, qcov, identity)
        
        # Build output row
        out_row = [accession] + filtered_tax
        output_data.append(out_row)
except Exception as e:
    print(f"ERROR: Processing OTU data: {e}", file=sys.stderr)
    sys.exit(1)

# Write to stdout (CSV format, comma-separated)
output_df = pd.DataFrame(output_data, columns=output_header)
try:
    output_df.to_csv(sys.stdout, sep=',', index=False)
except Exception as e:
    print(f"ERROR: Writing output: {e}", file=sys.stderr)
    sys.exit(1)
