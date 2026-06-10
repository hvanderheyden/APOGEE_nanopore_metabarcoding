#!/bin/bash
# =============================================================================
# APOGEE_V05 Parameter Sweep Script
# Tests combinations of -k, -v, -m, -C to optimize unclassified rate
# -H is fixed (not swept)
#
# NOTE: This script assumes the mapped/*.paf files already exist.
#       Only filtering/LCA steps are re-run (fast: ~2-5 min per run).
# =============================================================================

set -uo pipefail

# ─── Paths ────────────────────────────────────────────────────────────────────
WORK_DIR="/media/herve/10TB/Apogee/6_mock/9_new_V05"
SCRIPTS="/media/herve/10TB/Apogee/git/APOGEE_nanopore_metabarcoding/Version_05"
DB="/media/herve/10TB/Apogee/git/APOGEE_nanopore_metabarcoding/Data_base"
INPUT="/media/herve/10TB/Apogee/3_Fastq/FASTQ_mock/Fastq_pooled/fastq_mock_unzipped"
RESULTS_FILE="/tmp/apogee_param_sweep_$(date +%Y%m%d_%H%M).tsv"

# ─── Fixed parameters (not swept) ────────────────────────────────────────────
FIXED_OPTS="-c false -x 0.99 -w 8 -t 24 -b true -R TCCGTAGGTGAACCTGCGG -P TCCTCCGCTTATTGATATGC -L weighted"
H_VALUE=13    # top hits for LCA (fixed)

# ─── Parameter grid ───────────────────────────────────────────────────────────
# -k : identity threshold for LCA
K_VALUES=(0.70 0.75 0.80 0.85 0.90 0.95)
# -v : minimum query coverage
V_VALUES=(0.70 0.75 0.80 0.85 0.90)
# -m : minimum read length (bp)
M_VALUES=(250 350 450)
# -C : minimum mapping confidence
C_VALUES=(0.3 0.4 0.5)

# Total combinations
TOTAL=$(( ${#K_VALUES[@]} * ${#V_VALUES[@]} * ${#M_VALUES[@]} * ${#C_VALUES[@]} ))

# ─── Check prerequisites ──────────────────────────────────────────────────────
if [[ ! -d "${WORK_DIR}/mapped" ]] || [[ -z "$(ls "${WORK_DIR}/mapped"/*.paf 2>/dev/null)" ]]; then
    echo "ERROR: No PAF files found in ${WORK_DIR}/mapped/" >&2
    echo "       Run the full pipeline at least once before the sweep." >&2
    exit 1
fi

PAF_COUNT=$(ls "${WORK_DIR}/mapped"/*.paf 2>/dev/null | wc -l)
echo "==============================================="
echo "APOGEE_V05 Parameter Sweep"
echo "==============================================="
echo "Work dir:   ${WORK_DIR}"
echo "PAF files:  ${PAF_COUNT} files found (mapping will be skipped)"
echo ""
echo "Parameter grid:"
echo "  -k (identity threshold):  ${K_VALUES[*]}"
echo "  -v (query coverage):      ${V_VALUES[*]}"
echo "  -m (min read length):     ${M_VALUES[*]}"
echo "  -C (confidence):          ${C_VALUES[*]}"
echo "  -H (top hits for LCA):    ${H_VALUE} [fixed]"
echo "  -L (LCA type):            weighted [fixed]"
echo ""
echo "Total combinations: ${TOTAL}"
echo "Results file: ${RESULTS_FILE}"
echo "==============================================="
echo ""

# ─── Write results header ─────────────────────────────────────────────────────
printf "%-6s %-6s %-5s %-5s | %-12s %-12s %-12s %-12s %-12s\n" \
    "k" "v" "m" "C" \
    "MOCK5_total" "MOCK5_unclass" "MOCK5_pct" "MOCK3_pct" "MOCK4_pct" \
    | tee "${RESULTS_FILE%.tsv}_display.txt"
echo "-------|-------|------|------|-----------|-----------|-----------|-----------|-----------" \
    | tee -a "${RESULTS_FILE%.tsv}_display.txt"

# TSV header
echo -e "k\tv\tm\tC\tMOCK5_total\tMOCK5_unclassified\tMOCK5_pct\tMOCK3_pct\tMOCK4_pct" > "${RESULTS_FILE}"

# ─── Sweep loop ───────────────────────────────────────────────────────────────
run_count=0

for k in "${K_VALUES[@]}"; do
for v in "${V_VALUES[@]}"; do
for m in "${M_VALUES[@]}"; do
for C in "${C_VALUES[@]}"; do

    run_count=$(( run_count + 1 ))
    echo ""
    echo "─── Run ${run_count}/${TOTAL}: k=${k}  v=${v}  m=${m}  C=${C}  H=${H_VALUE} ───"

    # Clean only fast-step outputs (keep mapped/*.paf and everything upstream)
    rm -rf "${WORK_DIR}/filteredPAFs"
    rm -f  "${WORK_DIR}/collapsed_otu_table.csv" \
           "${WORK_DIR}/collapsed_taxonomy.csv" \
           "${WORK_DIR}/summary.csv"

    # Run pipeline (will skip all steps up to and including mapping)
    "${SCRIPTS}/APOGEE_V05.sh" \
        -i "${INPUT}" \
        -o "${WORK_DIR}" \
        -r "${DB}/ITS-RefDB_V03_fixed.mmi" \
        -T "${DB}/ITSRefDB_V03_taxonomy.tsv" \
        -F "${SCRIPTS}/filter_with_advanced_lca.py" \
        -S "${SCRIPTS}/generate_species_proportions.py" \
        -G "${SCRIPTS}/collapse_otu_by_taxonomy.py" \
        -K "${SCRIPTS}/species_clusters.tsv" \
        ${FIXED_OPTS} \
        -v "${v}" -m "${m}" -k "${k}" -C "${C}" -H "${H_VALUE}" \
        2>&1 | grep -E "^(✓|►|ERROR|WARNING|MOCK|mock|Step|Pipeline|unclass)" | tail -15

    # ─── Extract and record stats ──────────────────────────────────────────
    K_VAL="${k}" V_VAL="${v}" M_VAL="${m}" C_VAL="${C}" \
    RESULTS_TSV="${RESULTS_FILE}" \
    RESULTS_DISPLAY="${RESULTS_FILE%.tsv}_display.txt" \
    SUMMARY_FILE="${WORK_DIR}/summary.csv" \
    python3 << 'PYEOF'
import pandas as pd
import os, sys

k = os.environ['K_VAL']
v = os.environ['V_VAL']
m = os.environ['M_VAL']
C = os.environ['C_VAL']
results_tsv  = os.environ['RESULTS_TSV']
results_disp = os.environ['RESULTS_DISPLAY']
summary_file = os.environ['SUMMARY_FILE']

def get_unclass(df, sample_keyword):
    """Return (total_reads, unclassified_reads, pct_unclassified) for a sample."""
    rows = df[df['Sample'].str.contains(sample_keyword, case=False, na=False)]
    if rows.empty:
        return 0, 0, float('nan')
    total   = rows['Reads'].sum()
    unclass = rows[rows['Species'].str.lower().str.contains('unclassif', na=False)]['Reads'].sum()
    pct     = 100.0 * unclass / total if total > 0 else float('nan')
    return int(total), int(unclass), pct

try:
    df = pd.read_csv(summary_file)
    mock5_tot, mock5_unc, mock5_pct = get_unclass(df, 'mock5')
    mock3_tot, mock3_unc, mock3_pct = get_unclass(df, 'mock3')
    mock4_tot, mock4_unc, mock4_pct = get_unclass(df, 'mock4')

    # Console output
    print(f"  MOCK5: {mock5_tot:>9,} reads | {mock5_unc:>8,} unclassified ({mock5_pct:5.2f}%)")
    print(f"  MOCK3: {mock3_pct:5.2f}%  MOCK4: {mock4_pct:5.2f}%")

    # TSV row
    with open(results_tsv, 'a') as f:
        f.write(f"{k}\t{v}\t{m}\t{C}\t{mock5_tot}\t{mock5_unc}\t{mock5_pct:.3f}\t{mock3_pct:.3f}\t{mock4_pct:.3f}\n")

    # Display row
    with open(results_disp, 'a') as f:
        f.write(f"{k:<6} {v:<6} {m:<5} {C:<5} | {mock5_tot:<12} {mock5_unc:<12} {mock5_pct:<12.2f} {mock3_pct:<12.2f} {mock4_pct:<12.2f}\n")

except FileNotFoundError:
    msg = "MISSING_SUMMARY"
    print(f"  WARNING: summary.csv not found", file=sys.stderr)
    with open(results_tsv, 'a') as f:
        f.write(f"{k}\t{v}\t{m}\t{C}\t{msg}\t{msg}\t{msg}\t{msg}\t{msg}\n")
except Exception as e:
    print(f"  ERROR parsing summary: {e}", file=sys.stderr)
    with open(results_tsv, 'a') as f:
        f.write(f"{k}\t{v}\t{m}\t{C}\tERROR\tERROR\tERROR\tERROR\tERROR\n")
PYEOF

done  # C
done  # m
done  # v
done  # k

# ─── Final summary ────────────────────────────────────────────────────────────
echo ""
echo "==============================================="
echo "SWEEP COMPLETE — ${TOTAL} combinations tested"
echo "==============================================="
echo ""
echo "Full results table:"
cat "${RESULTS_FILE%.tsv}_display.txt"
echo ""

# Find best combination
RESULTS_FILE="${RESULTS_FILE}" python3 << 'PYEOF'
import pandas as pd

results_file = os.environ['RESULTS_FILE'] if 'RESULTS_FILE' in __import__('os').environ else None

import os
results_file = os.environ.get('RESULTS_FILE', '/tmp/apogee_sweep_results.tsv')

try:
    df = pd.read_csv(results_file, sep='\t')
    df_ok = df[pd.to_numeric(df['MOCK5_pct'], errors='coerce').notna()].copy()
    df_ok['MOCK5_pct'] = df_ok['MOCK5_pct'].astype(float)
    df_ok['MOCK3_pct'] = df_ok['MOCK3_pct'].astype(float)
    df_ok['MOCK4_pct'] = df_ok['MOCK4_pct'].astype(float)
    df_ok['avg_pct']   = (df_ok['MOCK5_pct'] + df_ok['MOCK3_pct'] + df_ok['MOCK4_pct']) / 3

    best_mock5 = df_ok.loc[df_ok['MOCK5_pct'].idxmin()]
    best_avg   = df_ok.loc[df_ok['avg_pct'].idxmin()]

    print("Best for MOCK5 only:")
    print(f"  -k {best_mock5['k']}  -v {best_mock5['v']}  -m {int(best_mock5['m'])}  -C {best_mock5['C']}")
    print(f"  → MOCK5: {best_mock5['MOCK5_pct']:.2f}%  MOCK3: {best_mock5['MOCK3_pct']:.2f}%  MOCK4: {best_mock5['MOCK4_pct']:.2f}%")
    print()
    print("Best average across all 3 samples:")
    print(f"  -k {best_avg['k']}  -v {best_avg['v']}  -m {int(best_avg['m'])}  -C {best_avg['C']}")
    print(f"  → MOCK5: {best_avg['MOCK5_pct']:.2f}%  MOCK3: {best_avg['MOCK3_pct']:.2f}%  MOCK4: {best_avg['MOCK4_pct']:.2f}%  avg: {best_avg['avg_pct']:.2f}%")

    print()
    print("Top 10 combinations by MOCK5 unclassified %:")
    top10 = df_ok.nsmallest(10, 'MOCK5_pct')[['k','v','m','C','MOCK5_pct','MOCK3_pct','MOCK4_pct','avg_pct']]
    print(top10.to_string(index=False))
except Exception as e:
    print(f"Could not analyze results: {e}")
PYEOF

echo ""
echo "TSV results: ${RESULTS_FILE}"
echo "Display:     ${RESULTS_FILE%.tsv}_display.txt"
