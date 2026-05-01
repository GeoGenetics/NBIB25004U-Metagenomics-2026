#!/bin/bash
#SBATCH --job-name=extract_hq_mags
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/logs/extract_hq_mags_%x_%j.out   # stdout log
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/extract_hq_mags_%x_%j.err    # stderr log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail

# ---- INPUTS ----
   #### The MAGs and CheckM2 results come from Jonas’ output directory.
   #### This ensures everyone works with the same input data, giving consistent and comparable results across the class.
MAGS_DIR="/maps/projects/course_1/people/fvb335/07_all_hq_mags/"
CHECKM2_file="/maps/projects/course_1/people/fvb335/06_quality_checkm2/quality_report.tsv"
HQ_DIR="/maps/projects/course_1/scratch/<group#>/<group-project-group-#>/07_1_hq_mags/"
EXTENSION="fa"

# ---- sanity checks ----
echo "Checking inputs..."
[[ -d "$MAGS_DIR" ]] || { echo "❌ Missing directory: $MAGS_DIR"; exit 1; }
[[ -f "$CHECKM2_file" ]] || { echo "❌ Missing file: $CHECKM2_file"; exit 1; }
echo "✅ Inputs look good"

# ---- make output dir ----
mkdir -p "$HQ_DIR"

# ---- print summary ----
echo "=========================================="
echo "Starting HQ MAG extraction"
echo "Input dir:  $MAGS_DIR"
echo "Report:     $CHECKM2_file"
echo "Output dir: $HQ_DIR"
echo "Criteria:   Completeness > 90 % AND Contamination < 5 %"
echo "=========================================="

count=0

# ---- main loop ----
while read -r mag; do
    src="${MAGS_DIR}/${mag}.${EXTENSION}"
    HQ_MAGS="${HQ_DIR}/${mag}.${EXTENSION}"

    echo "Processing: $mag"

    # ---- resume logic ----
    if [[ -f "$HQ_MAGS" ]]; then
        echo "  -> Already exists, skipping"
        continue
    fi

    echo "  -> Copying: $src"
    cp "$src" "$HQ_MAGS"
    ((++count))

    echo "------------------------------------------"

done < <(awk -F '\t' 'NR>1 && $2 > 90 && $3 < 5 {print $1}' "$CHECKM2_file")

# ---- summary ----
echo "=========================================="
echo "Finished."
echo "Total MAGs copied: $count"
echo "Output folder:     $HQ_DIR"
echo "=========================================="
