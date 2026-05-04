#!/bin/bash
# Quick summary of dbCAN CAZyme annotations across genomes.
# Run interactively from the project root (no SLURM job needed).

set -euo pipefail

# ---- PATHS ----
DBCAN_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/10_annotation_cazymes_ref/"

# ---- sanity checks ----
echo "Checking inputs..."
[[ -d "$DBCAN_DIR" ]] || { echo "❌ Missing directory: $DBCAN_DIR"; exit 1; }
shopt -s nullglob
overviews=("$DBCAN_DIR"/*/overview.tsv)
[[ "${#overviews[@]}" -gt 0 ]] || { echo "❌ No overview.tsv files found in $DBCAN_DIR"; exit 1; }
echo "✅ Inputs look good (${#overviews[@]} genomes)"

# ---- print summary ----
echo "=========================================="
echo "Summarising dbCAN results..."
echo "Input: $DBCAN_DIR"
echo "=========================================="

# 1. GH family hits per genome
echo ""
echo "=== GH hits ==="
for f in "${overviews[@]}"; do
    genome=$(basename "$(dirname "$f")")
    gh=$(awk -F'\t' 'NR>1 {print $7}' "$f" | grep -oE 'GH[0-9]+' | wc -l)
    printf "%s\t%s\n" "$genome" "$gh"
done

# 2. PL family hits per genome
echo ""
echo "=== PL hits ==="
for f in "${overviews[@]}"; do
    genome=$(basename "$(dirname "$f")")
    pl=$(awk -F'\t' 'NR>1 {print $7}' "$f" | grep -oE 'PL[0-9]+' | wc -l)
    printf "%s\t%s\n" "$genome" "$pl"
done

# 3. Total recommendation counts per genome
echo ""
echo "=== Total recommendations ==="
for f in "${overviews[@]}"; do
    genome=$(basename "$(dirname "$f")")
    total=$(awk -F'\t' 'NR>1 && $7!="-" {print $7}' "$f" | wc -l)
    printf "%s\t%s\n" "$genome" "$total"
done

# 4. Unique predicted substrates per genome
echo ""
echo "=== Unique substrates ==="
for f in "${overviews[@]}"; do
    genome=$(basename "$(dirname "$f")")
    subs=$(awk -F'\t' 'NR>1 && $8!="-" {print $8}' "$f" | tr ',' '\n' | sort -u | wc -l)
    printf "%s\t%s\n" "$genome" "$subs"
done

echo ""
echo "Done."
