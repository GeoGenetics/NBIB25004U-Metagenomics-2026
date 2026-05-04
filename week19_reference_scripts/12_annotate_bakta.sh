#!/bin/bash
#SBATCH --job-name=bakta
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/logs/bakta_%j.out
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/bakta_%j.err
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6G
#SBATCH --time=10:00:00
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

set -euo pipefail

# ---- PATHS ----
INPUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/07_1_hq_mags/"
OUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/09_annotation_bakta/"
DB="/maps/projects/course_1/data/bakta_db_light"

# ---- setup for Bakta ----
export PATH=/opt/shared_software/shared_envmodules/conda/bakta-1.11.3/bin:$PATH

# ---- sanity checks ----
bakta --version
command -v bakta >/dev/null 2>&1 || {
    echo "ERROR: 'bakta' not found in PATH." >&2
    exit 1
}

echo "Checking inputs..."
[[ -d "$INPUT_DIR" ]] || { echo "❌ Missing directory: $INPUT_DIR"; exit 1; }
[[ -d "$DB" ]]        || { echo "❌ Missing directory: $DB"; exit 1; }
echo "✅ Inputs look good"

# ---- make output dir ----
mkdir -p "$OUT_DIR"

# ---- collect input bins (.fa) ----
mapfile -t BINS < <(ls "$INPUT_DIR"/*.fa 2>/dev/null | sort)
[[ "${#BINS[@]}" -gt 0 ]] || { echo "❌ No .fa bin files found in $INPUT_DIR"; exit 1; }

# ---- print summary ----
echo "=========================================="
echo "Running Bakta on HQ MAGs..."
echo "Input:    $INPUT_DIR"
echo "Output:   $OUT_DIR"
echo "Database: $DB"
echo "Bins:     ${#BINS[@]}"
echo "=========================================="

# ---- run Bakta ----
for bin in "${BINS[@]}"; do
    sample=$(basename "$bin" .fa)
    sample_out="$OUT_DIR/$sample"

    # Skip if already annotated
    if [[ -s "$sample_out/$sample.gff3" ]]; then
        echo "[$sample] Already annotated — skipping."
        continue
    fi

    mkdir -p "$sample_out"
    echo "[$sample] Annotating..."

    bakta \
        --threads 6 \
        --db "$DB" \
        --compliant \
        --verbose \
        --force \
        --prefix "$sample" \
        --output "$sample_out" \
        "$bin"

    echo "[$sample] Done"
done

echo "Done."
