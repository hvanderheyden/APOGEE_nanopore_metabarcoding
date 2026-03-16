#!/usr/bin/env python3
"""
Advanced LCA filtering with multiple strategies
- Weighted LCA (by confidence/identity)
- Bootstrap-style approach (random subsampling)
- Phylogenetic distance-based filtering
"""
import sys
import os
import argparse
from collections import defaultdict, Counter
import random

parser = argparse.ArgumentParser(description='Advanced LCA taxonomy assignment with multiple strategies')
parser.add_argument('-i', '--paf_folder', required=True, help='Input PAF folder')
parser.add_argument('-T', '--taxonomy', required=True, help='Taxonomy database')
parser.add_argument('-o', '--output_file', required=True, help='Output OTU table')
parser.add_argument('-C', '--conf', type=float, default=0.5, help='Min confidence per hit')
parser.add_argument('-v', '--qcov', type=float, default=0.85, help='Min query coverage (fraction, 0-1)')
parser.add_argument('-m', '--min-qlen', type=int, default=0, help='Min query length in bp (absolute, overrides --qcov if larger)')
parser.add_argument('-k', '--identity', type=float, default=0.9, help='Min identity')
parser.add_argument('--lca-type', choices=['strict', 'weighted', 'bootstrap', 'majority-rule'], 
                    default='weighted', help='LCA computation method')
parser.add_argument('--log-filtering', action='store_true', help='Log hits that were filtered out')
parser.add_argument('--bootstrap-n', type=int, default=100, help='Number of bootstrap replicates')
parser.add_argument('--bootstrap-threshold', type=float, default=0.5, help='Bootstrap consensus threshold')
parser.add_argument('--divergence-penalty', type=float, default=0.1, help='Penalty for disagreement at each rank')
args = parser.parse_args()

paf_folder = args.paf_folder
out_file = args.output_file
min_conf = args.conf
min_qcov = args.qcov
min_qlen = args.min_qlen
min_identity = args.identity
lca_type = args.lca_type
log_filtering = args.log_filtering

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

RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
RANK_MAP = {'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order',
            'f': 'family', 'g': 'genus', 's': 'species'}

def parse_taxonomy_string(tax_string):
    """Parse QIIME format taxonomy string"""
    tax_dict = {rank: 'unclassified' for rank in RANKS}
    
    if not tax_string:
        return tax_dict
    
    parts = tax_string.split(';')
    for part in parts:
        part = part.strip()
        if '__' not in part:
            continue
        rank_code, value = part.split('__', 1)
        rank_code = rank_code.strip().lower()
        value = value.strip()
        if rank_code in RANK_MAP and value:
            rank_name = RANK_MAP[rank_code]
            tax_dict[rank_name] = value
    
    return tax_dict

def extract_optional_tag(tags_list, tag_name):
    """Extract optional tag from PAF"""
    for tag in tags_list:
        if tag.startswith(tag_name + ':'):
            parts = tag.split(':')
            return float(parts[-1])
    return None

def median(values):
    """Calculate median"""
    if not values:
        return None
    sorted_values = sorted(values)
    n = len(sorted_values)
    if n % 2 == 0:
        return (sorted_values[n//2 - 1] + sorted_values[n//2]) / 2
    else:
        return sorted_values[n//2]

def compute_lca_strict(taxon_list):
    """Simple strict LCA: all must agree"""
    if not taxon_list:
        return {rank: 'unclassified' for rank in RANKS}, -1
    
    if len(taxon_list) == 1:
        return taxon_list[0], len(RANKS) - 1
    
    lca_tax = {rank: 'unclassified' for rank in RANKS}
    deepest_rank = -1
    
    for rank_idx, rank in enumerate(RANKS):
        values = [tax[rank] for tax in taxon_list if tax[rank] != 'unclassified']
        
        if not values:
            break
        
        if len(set(values)) == 1:
            lca_tax[rank] = values[0]
            deepest_rank = rank_idx
        else:
            break
    
    return lca_tax, deepest_rank

def compute_lca_weighted(taxon_list, weight_list):
    """
    Weighted LCA: considers confidence/identity scores
    Higher weights = more trust in that taxonomy
    """
    if not taxon_list:
        return {rank: 'unclassified' for rank in RANKS}, -1
    
    if len(taxon_list) == 1:
        return taxon_list[0], len(RANKS) - 1
    
    # Normalize weights
    total_weight = sum(weight_list)
    norm_weights = [w / total_weight for w in weight_list]
    
    lca_tax = {rank: 'unclassified' for rank in RANKS}
    deepest_rank = -1
    
    for rank_idx, rank in enumerate(RANKS):
        # Get weighted vote
        weighted_votes = defaultdict(float)
        
        for tax, weight in zip(taxon_list, norm_weights):
            value = tax[rank]
            if value != 'unclassified':
                weighted_votes[value] += weight
        
        if not weighted_votes:
            break
        
        # Get most likely value
        best_value = max(weighted_votes.items(), key=lambda x: x[1])
        value, vote_weight = best_value
        
        # Accept if >60% weight agrees
        if vote_weight >= 0.6:
            lca_tax[rank] = value
            deepest_rank = rank_idx
        else:
            break
    
    return lca_tax, deepest_rank

def compute_lca_bootstrap(taxon_list, n_replicates=100, threshold=0.5):
    """
    Bootstrap consensus approach:
    1. Randomly subsample hits
    2. Compute LCA on subsample
    3. Track consensus across replicates
    """
    if not taxon_list:
        return {rank: 'unclassified' for rank in RANKS}, -1
    
    if len(taxon_list) == 1:
        return taxon_list[0], len(RANKS) - 1
    
    # Run bootstrap replicates
    rank_votes = [{} for _ in range(len(RANKS))]  # votes per rank
    
    for rep in range(n_replicates):
        sample_size = max(1, len(taxon_list) - 1)
        sample = random.sample(taxon_list, sample_size)
        
        # Compute strict LCA on this sample
        for rank_idx, rank in enumerate(RANKS):
            values = [tax[rank] for tax in sample if tax[rank] != 'unclassified']
            
            if not values:
                continue
            
            if len(set(values)) == 1:
                value = values[0]
                if value not in rank_votes[rank_idx]:
                    rank_votes[rank_idx][value] = 0
                rank_votes[rank_idx][value] += 1
    
    # Build consensus
    lca_tax = {rank: 'unclassified' for rank in RANKS}
    deepest_rank = -1
    
    for rank_idx, rank in enumerate(RANKS):
        if not rank_votes[rank_idx]:
            break
        
        # Get most common value
        best_value = max(rank_votes[rank_idx].items(), key=lambda x: x[1])
        value, count = best_value
        confidence = count / n_replicates
        
        # Accept if meets threshold
        if confidence >= threshold:
            lca_tax[rank] = value
            deepest_rank = rank_idx
        else:
            break
    
    return lca_tax, deepest_rank

def compute_lca_majority_rule(taxon_list):
    """Majority-rule LCA: >50% must agree at each rank"""
    if not taxon_list:
        return {rank: 'unclassified' for rank in RANKS}, -1
    
    if len(taxon_list) == 1:
        return taxon_list[0], len(RANKS) - 1
    
    lca_tax = {rank: 'unclassified' for rank in RANKS}
    deepest_rank = -1
    
    for rank_idx, rank in enumerate(RANKS):
        values = [tax[rank] for tax in taxon_list if tax[rank] != 'unclassified']
        
        if not values:
            break
        
        # Find majority
        counts = Counter(values)
        most_common, count = counts.most_common(1)[0]
        
        if count > len(values) * 0.5:
            lca_tax[rank] = most_common
            deepest_rank = rank_idx
        else:
            break
    
    return lca_tax, deepest_rank

# Get PAF files
paf_list = sorted([os.path.join(paf_folder, f) for f in os.listdir(paf_folder) if f.endswith('.paf')])
sample_list = [os.path.basename(f).replace('.paf', '') for f in paf_list]

# Process PAF files
query_hits = defaultdict(list)  # Good hits (pass filtering)
all_hits = defaultdict(list)     # All hits (for fallback if needed)
filtered_count = 0  # Track filtered hits

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
            qlen_bp = query_end - query_start
            
            dv_f = extract_optional_tag(info_list[12:], 'dv:f')
            identity = (1 - dv_f) if dv_f is not None else 1.0
            
            prob = 10**(int(map_qual) / -10)
            conf = 1.0 - prob
            
            weight = conf * identity  # Combined confidence weight
            
            # Track ALL hits for fallback
            hit_data = {
                'target': target_name,
                'conf': conf,
                'qcov': qcov,
                'identity': identity,
                'weight': weight,
                'qlen_bp': qlen_bp,
                'sample': sample_name
            }
            all_hits[query_name].append(hit_data)
            
            # Compute effective minimum query length (whichever is MORE restrictive)
            effective_min_qlen = max(min_qlen, int(min_qcov * query_len))
            
            # Filter at hit level
            if conf >= min_conf and qcov >= min_qcov and qlen_bp >= effective_min_qlen and identity >= min_identity:
                query_hits[query_name].append(hit_data)
            else:
                filtered_count += 1
                if log_filtering:
                    reason = []
                    if conf < min_conf:
                        reason.append(f"conf={conf:.3f}<{min_conf}")
                    if qcov < min_qcov:
                        reason.append(f"qcov={qcov:.3f}<{min_qcov}")
                    if qlen_bp < effective_min_qlen:
                        reason.append(f"qlen={qlen_bp}bp<{effective_min_qlen}bp")
                    if identity < min_identity:
                        reason.append(f"ident={identity:.3f}<{min_identity}")
                    print(f"FILTERED: {query_name} -> {target_name}: {', '.join(reason)}", file=sys.stderr)

# Fallback: for reads with NO good hits after filtering, keep the best original hit
# This ensures we never lose reads, only filter out secondary fragments
fallback_count = 0
for query_name in all_hits.keys():
    if query_name not in query_hits or len(query_hits[query_name]) == 0:
        # Get best hit from all_hits by weight
        best_hit = max(all_hits[query_name], key=lambda h: h['weight'])
        query_hits[query_name] = [best_hit]
        fallback_count += 1
        if log_filtering:
            print(f"FALLBACK: {query_name} -> {best_hit['target']}: All hits filtered, using best hit", file=sys.stderr)

# Output
output_dir = os.path.dirname(out_file)
os.makedirs(output_dir, exist_ok=True)

print(f"✓ Filtering complete: Kept {sum(len(h) for h in query_hits.values())} hits, filtered {filtered_count} hits, {fallback_count} reads with fallback", file=sys.stderr)

# OTU table
with open(out_file, 'w') as f:
    f.write('otu\t{}\tTotalCount\tNumHits\tAvgConfidence\tMedianQCov\tAvgIdentity\n'.format(
        '\t'.join(sample_list)))
    
    query_id = 1
    for query_name in sorted(query_hits.keys()):
        hits = query_hits[query_name]
        
        conf_list = [h['conf'] for h in hits]
        qcov_list = [h['qcov'] for h in hits]
        identity_list = [h['identity'] for h in hits]
        
        avg_conf = sum(conf_list) / len(conf_list) if conf_list else 0.0
        median_qcov = median(qcov_list) if qcov_list else 0.0
        avg_identity = sum(identity_list) / len(identity_list) if identity_list else 0.0
        
        sample_counts = defaultdict(int)
        for hit in hits:
            sample_counts[hit['sample']] += 1
        
        count_list = [sample_counts.get(s, 0) for s in sample_list]
        total_count = sum(count_list)
        
        otu_name = f"OTU_{query_id}"
        query_id += 1
        
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            otu_name, '\t'.join([str(x) for x in count_list]), total_count,
            len(hits), round(avg_conf, 4), round(median_qcov, 4), round(avg_identity, 4)))

print(f"✓ OTU table written: {out_file}")

# Print filtering summary
print(f"\n=== FILTERING PARAMETERS ===", file=sys.stderr)
print(f"  Min confidence: {min_conf}", file=sys.stderr)
print(f"  Min query coverage (fraction): {min_qcov}", file=sys.stderr)
print(f"  Min query length (absolute): {min_qlen} bp", file=sys.stderr)
print(f"  Min identity: {min_identity}", file=sys.stderr)

# Taxonomy with LCA
# Use the output file name to derive the tax file name
base_name = os.path.basename(out_file)
# Remove extension (.tsv, .csv, etc.) from the base name
tax_base = os.path.splitext(base_name)[0]
tax_file = os.path.join(output_dir, f'phyloseq_taxonomy_{tax_base}.csv')
with open(tax_file, 'w') as f:
    f.write(f'#OTU ID,Kingdom,Phylum,Class,Order,Family,Genus,Species,LCADepth,Method\n')
    
    query_id = 1
    for query_name in sorted(query_hits.keys()):
        hits = query_hits[query_name]
        
        # Get LCA method
        if lca_type == 'strict':
            hit_taxes = [parse_taxonomy_string(tax_db.get(h['target'], '')) for h in hits]
            lca_tax, deepest_rank = compute_lca_strict(hit_taxes)
        
        elif lca_type == 'weighted':
            hit_taxes = [parse_taxonomy_string(tax_db.get(h['target'], '')) for h in hits]
            weights = [h['weight'] for h in hits]
            lca_tax, deepest_rank = compute_lca_weighted(hit_taxes, weights)
        
        elif lca_type == 'bootstrap':
            hit_taxes = [parse_taxonomy_string(tax_db.get(h['target'], '')) for h in hits]
            lca_tax, deepest_rank = compute_lca_bootstrap(hit_taxes, 
                                                          n_replicates=args.bootstrap_n,
                                                          threshold=args.bootstrap_threshold)
        
        else:  # majority-rule
            hit_taxes = [parse_taxonomy_string(tax_db.get(h['target'], '')) for h in hits]
            lca_tax, deepest_rank = compute_lca_majority_rule(hit_taxes)
        
        otu_name = f"OTU_{query_id}"
        query_id += 1
        
        f.write(f"{otu_name},{lca_tax['kingdom']},{lca_tax['phylum']},{lca_tax['class']},{lca_tax['order']},{lca_tax['family']},{lca_tax['genus']},{lca_tax['species']},{deepest_rank},{lca_type}\n")

print(f"✓ Advanced taxonomy with LCA written: {tax_file}")
print(f"✓ Used method: {lca_type}")
