#!/usr/bin/env python3
"""
V03 OTU-level filtering with threshold-based taxonomy assignment
FIXED: Apply thresholds to taxonomy (genus-level for failed OTUs)
"""
import sys
import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Filter OTU table with confidence level')
parser.add_argument('-i', '--paf_folder', required=True, help='Input PAF folder')
parser.add_argument('-T', '--taxonomy', required=True, help='Taxonomy database')
parser.add_argument('-o', '--output_file', required=True, help='Output OTU table')
parser.add_argument('-C', '--conf', type=float, default=0.5, help='Min MappingConfidence')
parser.add_argument('-v', '--qcov', type=float, default=0.7, help='Min query coverage')
parser.add_argument('-p', '--identity', type=float, default=0.9, help='Min identity')
args = parser.parse_args()

paf_folder = args.paf_folder
out_file = args.output_file
min_conf = args.conf
min_qcov = args.qcov
min_identity = args.identity

# Load taxonomy
tax_db = {}
try:
    with open(args.taxonomy, 'r') as f:
        header = next(f)
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
    print(f"ERROR: Cannot read taxonomy: {e}", file=sys.stderr)
    sys.exit(1)

def parse_taxonomy_string(tax_string):
    """Parse QIIME format"""
    tax_dict = {
        'kingdom': 'unknown', 'phylum': 'unknown', 'class': 'unknown',
        'order': 'unknown', 'family': 'unknown', 'genus': 'unknown', 'species': 'unknown'
    }
    rank_map = {'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order',
                'f': 'family', 'g': 'genus', 's': 'species'}
    
    parts = tax_string.split(';')
    for part in parts:
        part = part.strip()
        if '__' in part:
            rank_code, value = part.split('__', 1)
            rank_code = rank_code.strip().lower()
            value = value.strip()
            if rank_code in rank_map:
                rank_name = rank_map[rank_code]
                tax_dict[rank_name] = value
    return tax_dict

def extract_optional_tag(tags_list, tag_name):
    for tag in tags_list:
        if tag.startswith(tag_name + ':'):
            parts = tag.split(':')
            return float(parts[-1])
    return None

def ave_prob(map_qual_list):
    if map_qual_list:
        return sum([10**(int(q) / -10) for q in map_qual_list]) / len(map_qual_list)
    else:
        return None

def prob_2_conf(prob):
    if prob is None:
        return 0.0
    return round(1.0 - prob, 4)

def median(values):
    if not values:
        return None
    sorted_values = sorted(values)
    n = len(sorted_values)
    if n % 2 == 0:
        return (sorted_values[n//2 - 1] + sorted_values[n//2]) / 2
    else:
        return sorted_values[n//2]

# Get PAF files and sample names
paf_list = sorted([os.path.join(paf_folder, f) for f in os.listdir(paf_folder) if f.endswith('.paf')])
sample_list = [os.path.basename(f).replace('.paf', '') for f in paf_list]

# Process PAF files
db_dict = defaultdict(dict)

for paf in paf_list:
    sample_name = os.path.basename(paf).replace('.paf', '')
    
    with open(paf, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            info_list = line.split('\t')
            
            query_name = info_list[0]
            query_len = int(info_list[1])
            query_start = int(info_list[2])
            query_end = int(info_list[3])
            target_name = info_list[5]
            map_qual = int(info_list[11])
            
            qcov = (query_end - query_start) / query_len
            dv_f = extract_optional_tag(info_list[12:], 'dv:f')
            identity = (1 - dv_f) if dv_f is not None else 1.0
            
            # Store for aggregation
            if target_name not in db_dict:
                db_dict[target_name] = dict()
                db_dict[target_name]['qual_list'] = [(map_qual, qcov, identity)]
                db_dict[target_name]['n'] = 1
            else:
                db_dict[target_name]['qual_list'].append((map_qual, qcov, identity))
                db_dict[target_name]['n'] += 1
            
            if sample_name not in db_dict[target_name]:
                db_dict[target_name][sample_name] = 1
            else:
                db_dict[target_name][sample_name] += 1

# Create output with taxonomy assignment based on thresholds
output_dir = os.path.dirname(out_file)
os.makedirs(output_dir, exist_ok=True)

# OTU table
with open(out_file, 'w') as f:
    f.write('Accession\t{}\tTotalCount\tMappingConfidence\tMedianQCov\tAvgIdentity\n'.format(
        '\t'.join(sample_list)))
    
    for target_name, my_dict in db_dict.items():
        qual_tuples = my_dict['qual_list']
        n_reads = my_dict['n']
        
        map_qual_list = [q[0] for q in qual_tuples]
        qcov_list = [q[1] for q in qual_tuples]
        identity_list = [q[2] for q in qual_tuples]
        
        prob = ave_prob(map_qual_list)
        conf = prob_2_conf(prob)
        median_qcov = median(qcov_list) if qcov_list else 0.0
        avg_identity = sum(identity_list) / len(identity_list) if identity_list else 0.0
        
        count_list = [my_dict.get(sample_name, 0) for sample_name in sample_list]
        
        f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            target_name,
            '\t'.join([str(x) for x in count_list]),
            n_reads,
            round(conf, 4),
            round(median_qcov, 4),
            round(avg_identity, 4)))

print(f"✓ OTU table written: {out_file}")

# Taxonomy file with threshold-based assignment
tax_file = os.path.join(output_dir, 'phyloseq_taxonomy_v03.csv')
with open(tax_file, 'w') as f:
    f.write('#OTU ID,Kingdom,Phylum,Class,Order,Family,Genus,Species\n')
    
    for target_name, my_dict in db_dict.items():
        qual_tuples = my_dict['qual_list']
        
        map_qual_list = [q[0] for q in qual_tuples]
        qcov_list = [q[1] for q in qual_tuples]
        identity_list = [q[2] for q in qual_tuples]
        
        prob = ave_prob(map_qual_list)
        conf = prob_2_conf(prob)
        median_qcov = median(qcov_list) if qcov_list else 0.0
        avg_identity = sum(identity_list) / len(identity_list) if identity_list else 0.0
        
        # Check if passes ALL thresholds
        passes_filters = (conf >= min_conf and median_qcov >= min_qcov and avg_identity >= min_identity)
        
        # Get taxonomy
        if target_name not in tax_db:
            tax_dict = {'kingdom': 'unknown', 'phylum': 'unknown', 'class': 'unknown',
                       'order': 'unknown', 'family': 'unknown', 'genus': 'unknown', 'species': 'unknown'}
        else:
            tax_dict = parse_taxonomy_string(tax_db[target_name])
        
        # Apply threshold-based taxonomy
        if not passes_filters:
            # OTU fails: set species to unknown_Genus
            genus_name = tax_dict.get('genus', 'unknown')
            tax_dict['species'] = f'unknown_{genus_name}'
        
        f.write(f"{target_name},{tax_dict['kingdom']},{tax_dict['phylum']},{tax_dict['class']},{tax_dict['order']},{tax_dict['family']},{tax_dict['genus']},{tax_dict['species']}\n")

print(f"✓ Taxonomy written: {tax_file}")
