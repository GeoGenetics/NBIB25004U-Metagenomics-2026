#!/bin/bash
#SBATCH --job-name=dbcan
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_%j.out
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_%j.err
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=6G
#SBATCH --time=10:00:00
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

set -euo pipefail

# ---- PATHS ----
INPUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/09_annotation_bakta/"
OUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/10_annotation_cazymes/"
DB="/maps/projects/course_1/data/dbcan_db"

# ---- setup for dbCAN ----
export PATH=/opt/shared_software/shared_envmodules/conda/dbcan-5.2.8/bin:$PATH

# ---- sanity checks ----
run_dbcan version
command -v run_dbcan >/dev/null 2>&1 || {
    echo "ERROR: 'run_dbcan' not found in PATH." >&2
    exit 1
}

echo "Checking inputs..."
[[ -d "$INPUT_DIR" ]] || { echo "❌ Missing directory: $INPUT_DIR"; exit 1; }
[[ -d "$DB" ]]        || { echo "❌ Missing directory: $DB"; exit 1; }
echo "✅ Inputs look good"

# ---- make output dir ----
mkdir -p "$OUT_DIR"

# ---- collect input proteomes (.faa from Bakta) ----
mapfile -t FAAS < <(ls "$INPUT_DIR"/*/*.faa 2>/dev/null | sort)
[[ "${#FAAS[@]}" -gt 0 ]] || { echo "❌ No .faa files found in $INPUT_DIR"; exit 1; }

# ---- print summary ----
echo "=========================================="
echo "Running dbCAN on Bakta proteomes..."
echo "Input:    $INPUT_DIR"
echo "Output:   $OUT_DIR"
echo "Database: $DB"
echo "Genomes:  ${#FAAS[@]}"
echo "=========================================="

# ---- run dbCAN ----
for faa in "${FAAS[@]}"; do
    sample=$(basename "$faa" .faa)
    sample_out="$OUT_DIR/$sample"

    # Skip if already annotated
    if [[ -s "$sample_out/overview.tsv" ]]; then
        echo "[$sample] Already annotated — skipping."
        continue
    fi

    mkdir -p "$sample_out"
    echo "[$sample] Annotating CAZymes from $faa"

    run_dbcan CAZyme_annotation \
        --input_raw_data "$faa" \
        --mode protein \
        --output_dir "$sample_out" \
        --db_dir "$DB" \
        --threads 15

    echo "[$sample] Done"
done

echo "Done."
