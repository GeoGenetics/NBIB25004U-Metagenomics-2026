#!/bin/bash
# Quick summary of dbCAN CAZyme annotations across the 3 in-class reference genomes.
# Run interactively from a login node with: bash 13_2_cazyme_summary.sh
#
# Output:
#   - Status messages stay on stdout.
#   - Summary tables are written to $SUMMARY_FILE (a TSV alongside the dbCAN
#     output folders) and echoed back to stdout at the end for visibility.

set -euo pipefail

# ---- PATHS ----
DBCAN_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/10_annotation_cazymes_ref/"
SUMMARY_FILE="$DBCAN_DIR/dbcan_summary_ref.tsv"

# ---- counters (awk gsub: returns 0 cleanly when no matches) ----
count_family() {
    local file="$1" pat="$2"
    awk -F'\t' -v pat="$pat" '
        NR>1 { total += gsub(pat, "&", $7) }
        END  { print total+0 }
    ' "$file"
}

count_recommendations() {
    awk -F'\t' 'NR>1 && $7!="-" { n++ } END { print n+0 }' "$1"
}

count_unique_substrates() {
    awk -F'\t' 'NR>1 && $8!="-" { print $8 }' "$1" \
        | tr ',' '\n' \
        | awk 'NF { sub(/^[ \t]+/, ""); sub(/[ \t]+$/, ""); print }' \
        | sort -u \
        | awk 'END { print NR+0 }'
}

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
echo "Input:  $DBCAN_DIR"
echo "Output: $SUMMARY_FILE"
echo "=========================================="

# ---- write summary file ----
{
    echo "# dbCAN summary (in-class reference genomes)"
    echo "# Input:    $DBCAN_DIR"
    echo "# Genomes:  ${#overviews[@]}"
    echo "# Generated: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo

    # 1. Combined per-genome table
    echo "## Per-genome table (genome | total | GH | PL | unique_substrates)"
    echo -e "genome\ttotal\tGH\tPL\tunique_substrates"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        total=$(count_recommendations "$f")
        gh=$(count_family "$f" "GH[0-9]+")
        pl=$(count_family "$f" "PL[0-9]+")
        subs=$(count_unique_substrates "$f")
        printf "%s\t%s\t%s\t%s\t%s\n" "$genome" "$total" "$gh" "$pl" "$subs"
    done | sort

    echo

    # 2. GH family hits per genome
    echo "## GH hits (genome | count)"
    echo -e "genome\tGH_count"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        gh=$(count_family "$f" "GH[0-9]+")
        printf "%s\t%s\n" "$genome" "$gh"
    done | sort

    echo

    # 3. PL family hits per genome
    echo "## PL hits (genome | count)"
    echo -e "genome\tPL_count"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        pl=$(count_family "$f" "PL[0-9]+")
        printf "%s\t%s\n" "$genome" "$pl"
    done | sort

    echo

    # 4. Total recommendation counts per genome
    echo "## Total recommendations (genome | count)"
    echo -e "genome\ttotal_recs"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        total=$(count_recommendations "$f")
        printf "%s\t%s\n" "$genome" "$total"
    done | sort

    echo

    # 5. Unique predicted substrates per genome
    echo "## Unique substrates (genome | count)"
    echo -e "genome\tunique_substrates"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        subs=$(count_unique_substrates "$f")
        printf "%s\t%s\n" "$genome" "$subs"
    done | sort
} > "$SUMMARY_FILE"

# ---- echo file back to stdout for visibility ----
echo
cat "$SUMMARY_FILE"
echo
echo "✅ Summary written to: $SUMMARY_FILE  ($(wc -l < "$SUMMARY_FILE") lines)"
echo "Done."
