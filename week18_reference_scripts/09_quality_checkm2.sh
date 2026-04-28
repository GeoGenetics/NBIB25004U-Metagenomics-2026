#!/bin/bash
#SBATCH --job-name=checkm2_quality
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=4G         # memory per cpu-core
#SBATCH --time=05:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching
#
# run_checkm2.sh — Run CheckM2 'predict' on a folder of MAGs / bins.
#
# Usage:
#   run_checkm2.sh -i <folder_with_bins> -o <output_folder> [options]
#
# Required:
#   -i, --input DIR          Folder containing the bin/MAG FASTA files
#   -o, --output DIR         Output directory for CheckM2 results
#
# Optional:
#   -d, --database PATH      Path to CheckM2 diamond database
#                            (default: /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd)
#   -t, --threads INT        Number of threads                 (default: 25)
#   -x, --extension STR      File extension of bins            (default: fa)
#   -e, --extra "STR"        Extra raw arguments forwarded to checkm2 (quoted)
#   -f, --force              Overwrite existing output directory
#   -h, --help               Show this help and exit
# ---------------------------------------------------------------------------

# ---- 0. Load Modules -------------------------------------------------------

export PATH=/opt/shared_software/shared_envmodules/conda/checkm2-1.0.2/bin:$PATH
checkm2 --version    # sanity check

set -euo pipefail

# ---------- defaults --------------------------------------------------------
THREADS=25
EXTENSION="fa"
DATABASE="/maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd"
EXTRA_ARGS=""
FORCE=0

INPUT_DIR=""
OUT_DIR=""

# ---------- usage -----------------------------------------------------------
usage() {
    cat <<EOF
Usage:
  $(basename "$0") -i <folder_with_bins> -o <output_folder> [options]

Required:
  -i, --input DIR          Folder containing the bin/MAG FASTA files
  -o, --output DIR         Output directory for CheckM2 results

Optional:
  -d, --database PATH      Path to CheckM2 diamond database
                           (default: ${DATABASE})
  -t, --threads INT        Number of threads                 (default: ${THREADS})
  -x, --extension STR      File extension of bins            (default: ${EXTENSION})
  -e, --extra "STR"        Extra raw arguments forwarded to checkm2 (quoted)
  -f, --force              Overwrite existing output directory
  -h, --help               Show this help and exit
EOF
}

# ---------- arg parsing -----------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)       INPUT_DIR="$2";  shift 2 ;;
        -o|--output)      OUT_DIR="$2";    shift 2 ;;
        -d|--database)    DATABASE="$2";   shift 2 ;;
        -t|--threads)     THREADS="$2";    shift 2 ;;
        -x|--extension)   EXTENSION="$2";  shift 2 ;;
        -e|--extra)       EXTRA_ARGS="$2"; shift 2 ;;
        -f|--force)       FORCE=1;         shift   ;;
        -h|--help)        usage; exit 0 ;;
        *)  echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ---------- sanity checks ---------------------------------------------------
command -v checkm2 >/dev/null 2>&1 || {
    echo "ERROR: 'checkm2' not found in PATH." >&2; exit 1;
}

if [[ -z "$INPUT_DIR" || -z "$OUT_DIR" ]]; then
    echo "ERROR: -i/--input and -o/--output are required." >&2
    usage; exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory not found: $INPUT_DIR" >&2; exit 1
fi

if [[ ! -f "$DATABASE" ]]; then
    echo "ERROR: database file not found: $DATABASE" >&2; exit 1
fi

# Count input bins as a sanity check
shopt -s nullglob
BINS=( "${INPUT_DIR%/}"/*."${EXTENSION}" )
shopt -u nullglob
if [[ ${#BINS[@]} -eq 0 ]]; then
    echo "ERROR: no '*.${EXTENSION}' files found in: $INPUT_DIR" >&2
    echo "       (use -x/--extension if your bins use a different suffix)" >&2
    exit 1
fi

# ---------- prepare output directory ---------------------------------------
mkdir -p "$(dirname "$OUT_DIR")"

# ---------- build & run command --------------------------------------------
cmd=( checkm2 predict
      --threads "$THREADS"
      --input "$INPUT_DIR"
      --output-directory "$OUT_DIR"
      --database_path "$DATABASE"
      -x "$EXTENSION" )

[[ "$FORCE" -eq 1 ]] && cmd+=( --force )
# shellcheck disable=SC2206
[[ -n "$EXTRA_ARGS" ]] && cmd+=( $EXTRA_ARGS )

echo "======================================================================"
echo "[$(date '+%F %T')] Running CheckM2 predict"
echo "  Input:    $INPUT_DIR  (${#BINS[@]} bin(s) with .${EXTENSION})"
echo "  Output:   $OUT_DIR"
echo "  Database: $DATABASE"
echo "  Threads:  $THREADS"
echo "  Cmd:      ${cmd[*]}"
echo "======================================================================"

"${cmd[@]}"

echo "[$(date '+%F %T')] Done."
