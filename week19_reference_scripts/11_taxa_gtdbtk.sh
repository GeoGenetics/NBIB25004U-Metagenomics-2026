#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/logs/gtdbtk_%j.out
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/gtdbtk_%j.err
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=6G         # GTDB-Tk R232 requires ≥140 GB RAM (30*6=180 GB)
#SBATCH --time=10:00:00

set -euo pipefail

# ---- PATHS ----
INPUT_DIR="/maps/projects/course_1/scratch/<group#>/<group-project-group-#>/07_1_hq_mags/"
OUT_DIR="/maps/projects/course_1/scratch/<group#>/<group-project-group-#>/08_gtdbtk/"
DB="/maps/projects/course_1/data/gtdb232/release232"

# ---- setup for GTDB-Tk ----
export PATH=/opt/shared_software/shared_envmodules/conda/gtdbtk-2.7.1/bin:$PATH
export GTDBTK_DATA_PATH="$DB"

# ---- sanity checks ----             
command -v gtdbtk >/dev/null 2>&1 || {
    echo "ERROR: 'gtdbtk' not found in PATH." >&2
    exit 1
}
gtdbtk --version    

echo "Checking inputs..."
[[ -d "$INPUT_DIR" ]] || { echo "❌ Missing directory: $INPUT_DIR"; exit 1; }
[[ -d "$DB" ]] || { echo "❌ Missing directory: $DB"; exit 1; }
echo "✅ Inputs look good"

# ---- make output dir ----
mkdir -p "$OUT_DIR"

# ---- print summary ----
echo "=========================================="
echo "Running GTDB-Tk on HQ MAGs..."
echo "Input:  $INPUT_DIR"
echo "Output: $OUT_DIR"
echo "=========================================="

# ---- run GTDB-Tk ----
gtdbtk classify_wf \
    --genome_dir "$INPUT_DIR" \
    --out_dir "$OUT_DIR" \
    --cpus 30 \
    --extension fa \
    --place_species

echo "Done."
