#!/bin/bash
#SBATCH --job-name=dbcan
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_summary_%j.out
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_summary_%j.err
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=6G
#SBATCH --time=10:00:00
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching#!/bin/bash
#SBATCH --job-name=dbcan
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_summary_%j.out
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_summary_%j.err
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=6G
#SBATCH --time=10:00:00
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

# Per-MAG dbCAN summary, joined with GTDB-Tk taxonomy.
# Submit with: sbatch 13_4_cazyme_summary_mags.sh
#
# Output:
#   - SLURM stdout (.out file) carries status messages.
#   - All summary tables go to $SUMMARY_FILE as one TSV with section headers.
#
# Bodysite grouping: each MAG name is expected to contain its sample accession
# (e.g. "ERR2641635_bin.3"). Reference genomes that don't match any accession
# get bodysite = "n/a".
#
# Taxonomy: looked up against the shared GTDB-Tk bac120 summary by exact match
# on the MAG name (column 1, user_genome). Missing matches → taxonomy = "n/a".

# ---- PATHS ----
DBCAN_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/10_annotation_cazymes/"
GTDBTK_FILE="//maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/08_taxa_gtdbtk/gtdbtk.bac120.summary.tsv"
SUMMARY_FILE="$DBCAN_DIR/dbcan_summary.tsv"

# ---- bodysite lookup ----
get_bodysite() {
    local mag="$1"
    case "$mag" in
        *ERR2641635*|*ERR2641677*|*ERR2641733*) echo "gut_adult"  ;;
        *SRR8692206*|*SRR8692207*|*SRR8692213*) echo "gut_infant" ;;
        *SRR059458*|*SRR059459*|*SRR513791*)    echo "vaginal"    ;;
        *)                                       echo "n/a"        ;;
    esac
}

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
echo "✅ dbCAN: found ${#overviews[@]} genomes"

# ---- load GTDB-Tk taxonomy lookup ----
declare -A TAXA
if [[ -s "$GTDBTK_FILE" ]]; then
    while IFS=$'\t' read -r genome classification _rest; do
        [[ -n "$genome" && "$genome" != "user_genome" ]] || continue
        TAXA["$genome"]="$classification"
    done < "$GTDBTK_FILE"
    echo "✅ GTDB-Tk: loaded ${#TAXA[@]} taxonomy entries from $GTDBTK_FILE"
else
    echo "⚠️  GTDB-Tk file not found at $GTDBTK_FILE — taxonomy column will be 'n/a'"
fi

# ---- write summary file ----
echo "=========================================="
echo "Summarising dbCAN results..."
echo "Input:  $DBCAN_DIR"
echo "Output: $SUMMARY_FILE"
echo "=========================================="

{
    echo "# dbCAN summary"
    echo "# Input:    $DBCAN_DIR"
    echo "# Taxonomy: $GTDBTK_FILE"
    echo "# Genomes:  ${#overviews[@]}"
    echo "# Generated: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo

    # ---- 1. MAGs per bodysite ----
    echo "## MAGs per bodysite"
    echo -e "bodysite\tn_MAGs"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        get_bodysite "$genome"
    done | sort | uniq -c | sort -rn | awk '{printf "%s\t%s\n", $2, $1}'
    echo

    # ---- 2. Combined per-MAG table ----
    # Columns: bodysite | genome | total | GH | PL | unique_substrates | taxonomy
    echo "## Per-MAG table (bodysite | genome | total | GH | PL | unique_substrates | taxonomy)"
    echo -e "bodysite\tgenome\ttotal\tGH\tPL\tunique_substrates\ttaxonomy"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        bodysite=$(get_bodysite "$genome")
        total=$(count_recommendations "$f")
        gh=$(count_family "$f" "GH[0-9]+")
        pl=$(count_family "$f" "PL[0-9]+")
        subs=$(count_unique_substrates "$f")
        taxa="${TAXA[$genome]:-n/a}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$bodysite" "$genome" "$total" "$gh" "$pl" "$subs" "$taxa"
    done | sort
    echo

    # ---- 3. Aggregate totals per bodysite ----
    echo "## Totals per bodysite (bodysite | n_MAGs | total_GH | total_PL)"
    echo -e "bodysite\tn_MAGs\ttotal_GH\ttotal_PL"
    declare -A N_MAGS GH_TOT PL_TOT
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        bodysite=$(get_bodysite "$genome")
        gh=$(count_family "$f" "GH[0-9]+")
        pl=$(count_family "$f" "PL[0-9]+")
        N_MAGS[$bodysite]=$(( ${N_MAGS[$bodysite]:-0} + 1 ))
        GH_TOT[$bodysite]=$(( ${GH_TOT[$bodysite]:-0} + gh ))
        PL_TOT[$bodysite]=$(( ${PL_TOT[$bodysite]:-0} + pl ))
    done
    for site in "${!N_MAGS[@]}"; do
        printf "%s\t%s\t%s\t%s\n" "$site" "${N_MAGS[$site]}" "${GH_TOT[$site]}" "${PL_TOT[$site]}"
    done | sort
} > "$SUMMARY_FILE"

# ---- final status ----
echo "✅ Summary written to: $SUMMARY_FILE"
echo "   Sections: 3 (MAGs per bodysite, per-MAG table with taxonomy, totals per bodysite)"
echo "   Lines:    $(wc -l < "$SUMMARY_FILE")"
echo "Done."

# Per-MAG dbCAN summary, joined with GTDB-Tk taxonomy.
# Submit with: sbatch dbcan_summary.sh
#
# Output:
#   - SLURM stdout (.out file) carries status messages.
#   - All summary tables go to $SUMMARY_FILE as one TSV with section headers.
#
# Bodysite grouping: each MAG name is expected to contain its sample accession
# (e.g. "ERR2641635_bin.3"). Reference genomes that don't match any accession
# get bodysite = "n/a".
#
# Taxonomy: looked up against the shared GTDB-Tk bac120 summary by exact match
# on the MAG name (column 1, user_genome). Missing matches → taxonomy = "n/a".

# ---- PATHS ----
DBCAN_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/10_annotation_cazymes/"
GTDBTK_FILE="//maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/08_taxa_gtdbtk/gtdbtk.bac120.summary.tsv"
SUMMARY_FILE="$DBCAN_DIR/dbcan_summary.tsv"

# ---- bodysite lookup ----
get_bodysite() {
    local mag="$1"
    case "$mag" in
        *ERR2641635*|*ERR2641677*|*ERR2641733*) echo "gut_adult"  ;;
        *SRR8692206*|*SRR8692207*|*SRR8692213*) echo "gut_infant" ;;
        *SRR059458*|*SRR059459*|*SRR513791*)    echo "vaginal"    ;;
        *)                                       echo "n/a"        ;;
    esac
}

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
echo "✅ dbCAN: found ${#overviews[@]} genomes"

# ---- load GTDB-Tk taxonomy lookup ----
declare -A TAXA
if [[ -s "$GTDBTK_FILE" ]]; then
    while IFS=$'\t' read -r genome classification _rest; do
        [[ -n "$genome" && "$genome" != "user_genome" ]] || continue
        TAXA["$genome"]="$classification"
    done < "$GTDBTK_FILE"
    echo "✅ GTDB-Tk: loaded ${#TAXA[@]} taxonomy entries from $GTDBTK_FILE"
else
    echo "⚠️  GTDB-Tk file not found at $GTDBTK_FILE — taxonomy column will be 'n/a'"
fi

# ---- write summary file ----
echo "=========================================="
echo "Summarising dbCAN results..."
echo "Input:  $DBCAN_DIR"
echo "Output: $SUMMARY_FILE"
echo "=========================================="

{
    echo "# dbCAN summary"
    echo "# Input:    $DBCAN_DIR"
    echo "# Taxonomy: $GTDBTK_FILE"
    echo "# Genomes:  ${#overviews[@]}"
    echo "# Generated: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo

    # ---- 1. MAGs per bodysite ----
    echo "## MAGs per bodysite"
    echo -e "bodysite\tn_MAGs"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        get_bodysite "$genome"
    done | sort | uniq -c | sort -rn | awk '{printf "%s\t%s\n", $2, $1}'
    echo

    # ---- 2. Combined per-MAG table ----
    # Columns: bodysite | genome | total | GH | PL | unique_substrates | taxonomy
    echo "## Per-MAG table (bodysite | genome | total | GH | PL | unique_substrates | taxonomy)"
    echo -e "bodysite\tgenome\ttotal\tGH\tPL\tunique_substrates\ttaxonomy"
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        bodysite=$(get_bodysite "$genome")
        total=$(count_recommendations "$f")
        gh=$(count_family "$f" "GH[0-9]+")
        pl=$(count_family "$f" "PL[0-9]+")
        subs=$(count_unique_substrates "$f")
        taxa="${TAXA[$genome]:-n/a}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$bodysite" "$genome" "$total" "$gh" "$pl" "$subs" "$taxa"
    done | sort
    echo

    # ---- 3. Aggregate totals per bodysite ----
    echo "## Totals per bodysite (bodysite | n_MAGs | total_GH | total_PL)"
    echo -e "bodysite\tn_MAGs\ttotal_GH\ttotal_PL"
    declare -A N_MAGS GH_TOT PL_TOT
    for f in "${overviews[@]}"; do
        genome=$(basename "$(dirname "$f")")
        bodysite=$(get_bodysite "$genome")
        gh=$(count_family "$f" "GH[0-9]+")
        pl=$(count_family "$f" "PL[0-9]+")
        N_MAGS[$bodysite]=$(( ${N_MAGS[$bodysite]:-0} + 1 ))
        GH_TOT[$bodysite]=$(( ${GH_TOT[$bodysite]:-0} + gh ))
        PL_TOT[$bodysite]=$(( ${PL_TOT[$bodysite]:-0} + pl ))
    done
    for site in "${!N_MAGS[@]}"; do
        printf "%s\t%s\t%s\t%s\n" "$site" "${N_MAGS[$site]}" "${GH_TOT[$site]}" "${PL_TOT[$site]}"
    done | sort
} > "$SUMMARY_FILE"

# ---- final status ----
echo "✅ Summary written to: $SUMMARY_FILE"
echo "   Sections: 3 (MAGs per bodysite, per-MAG table with taxonomy, totals per bodysite)"
echo "   Lines:    $(wc -l < "$SUMMARY_FILE")"
echo "Done."
