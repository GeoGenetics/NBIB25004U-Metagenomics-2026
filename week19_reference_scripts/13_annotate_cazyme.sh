#!/bin/bash
#SBATCH --job-name=dbcan
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_%j.out
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/dbcan_%j.err
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=6G
#SBATCH --time=10:00:00
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

# ---- PATHS ----
# In class: run on Bakta-annotated genomes from the three shared reference genomes.
# After class: swap INPUT_DIR / OUT_DIR to drop the _ref suffix:
#   INPUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/09_annotation_bakta/"
#   OUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/10_annotation_cazymes/"
INPUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/09_annotation_bakta_ref/"
OUT_DIR="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/10_annotation_cazymes_ref/"
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

# ---- enumerate Bakta-annotated genomes via .fna ----
# Why .fna for enumeration? Bakta writes both <sample>.faa and
# <sample>_hypotheticals.faa per genome, so globbing *.faa would list each
# genome twice. The .fna file is unique per genome, so we use it as the
# enumeration anchor and then look up the matching <sample>.faa for the
# actual CAZyme annotation below (--mode protein).
mapfile -t bins < <(ls "$INPUT_DIR"/*/*.fna 2>/dev/null | sort)
[[ "${#bins[@]}" -gt 0 ]] || { echo "❌ No .fna files found in $INPUT_DIR"; exit 1; }

# ---- print summary ----
echo "=========================================="
echo "Running dbCAN on Bakta-annotated genomes..."
echo "Input:    $INPUT_DIR"
echo "Output:   $OUT_DIR"
echo "Database: $DB"
echo "Genomes:  ${#bins[@]}"
echo "=========================================="

# ---- run dbCAN ----
for bin in "${bins[@]}"; do
    sample=$(basename "$bin" .fna)
    sample_out="$OUT_DIR/$sample"
    faa_file="$INPUT_DIR/$sample/$sample.faa"

    # Make sure the protein file actually exists
    if [[ ! -s "$faa_file" ]]; then
        echo "[$sample] No .faa file found at $faa_file — skipping."
        continue
    fi

    # Skip if already annotated (overview.tsv is the canonical dbCAN output)
    if [[ -s "$sample_out/overview.tsv" ]]; then
        echo "[$sample] Already annotated — skipping."
        continue
    fi

    mkdir -p "$sample_out"
    echo "[$sample] Running dbCAN on $faa_file"

    run_dbcan CAZyme_annotation \
        --input_raw_data "$faa_file" \
        --mode protein \
        --output_dir "$sample_out" \
        --db_dir "$DB" \
        --threads 15

    echo "[$sample] Done"
done

echo "Done."
echo "=== Job end: $(date) ==="
