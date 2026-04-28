#!/bin/bash
#SBATCH --job-name=MetaphlanRun
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem-per-cpu=8G    # memory per cpu-core
#SBATCH --time=03:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching
#
# run_metaphlan.sh — Run MetaPhlAn taxonomic profiling on cleaned paired-end reads.
#
# Two modes of operation:
#
#   (A) Batch mode — profile every accession of a given bodysite:
#       run_metaphlan.sh -p <preprocessing_dir> -b <bodysite> -o <metaphlan_dir> [options]
#
#   (B) Single-sample mode — profile one pair of reads directly:
#       run_metaphlan.sh -1 <R1.fq.gz> -2 <R2.fq.gz> -a <accession> -o <out_dir> [options]
#
# In batch mode the expected input layout is:
#   <preprocessing_dir>/<bodysite>/<accession>/cleaned/<accession>_1_clean.fq.gz
#   <preprocessing_dir>/<bodysite>/<accession>/cleaned/<accession>_2_clean.fq.gz
#
# And the output layout produced is:
#   <metaphlan_dir>/<bodysite>/<accession>/<accession>_metaphlan.txt
#   <metaphlan_dir>/<bodysite>/<accession>/<accession>.mapout.bz2
# ---------------------------------------------------------------------------

# ---- 0. Load Modules -------------------------------------------------------

module load metaphlan/4.1.1 bowtie2

set -euo pipefail

# ---------- defaults --------------------------------------------------------
NPROC=11
INPUT_TYPE="fastq"
EXTRA_ARGS=""

PREPROC_DIR=""
BODYSITE=""
R1=""
R2=""
ACCESSION=""
OUT_DIR=""
METAPHLAN_DB_DIR="/maps/datasets/globe_databases/metaphlan/20241118"
INDEX_DB="mpa_vJun23_CHOCOPhlAnSGB_202403"

# ---------- usage -----------------------------------------------------------
usage() {
    cat <<EOF
Usage:
  Batch mode:
    $(basename "$0") -p <preprocessing_dir> -b <bodysite> -o <metaphlan_dir> [options]

  Single-sample mode:
    $(basename "$0") -1 <R1.fq.gz> -2 <R2.fq.gz> -a <accession> -o <out_dir> [options]

Required (batch mode):
  -p, --preproc DIR        Preprocessing directory (contains <bodysite>/<acc>/cleaned/...)
  -b, --bodysite NAME      Bodysite subdirectory (e.g. gut_adult, gut_infant, vaginal)
  -o, --output DIR         Root metaphlan output directory

Required (single-sample mode):
  -1, --reads1 FILE        Forward reads (fastq[.gz])
  -2, --reads2 FILE        Reverse reads (fastq[.gz])
  -a, --accession NAME     Accession ID (used to name output files)
  -o, --output DIR         Output directory for this accession

MetaPhlAn options:
  -t, --nproc INT          Number of threads                (default: ${NPROC})
  -i, --input-type STR     MetaPhlAn --input_type           (default: ${INPUT_TYPE})
  -e, --extra "STR"        Extra raw arguments forwarded to metaphlan (quoted)

  -h, --help               Show this help and exit
EOF
}

# ---------- arg parsing -----------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--preproc)     PREPROC_DIR="$2"; shift 2 ;;
        -b|--bodysite)    BODYSITE="$2";    shift 2 ;;
        -1|--reads1)      R1="$2";          shift 2 ;;
        -2|--reads2)      R2="$2";          shift 2 ;;
        -a|--accession)   ACCESSION="$2";   shift 2 ;;
        -o|--output)      OUT_DIR="$2";     shift 2 ;;
        -t|--nproc)       NPROC="$2";       shift 2 ;;
        -i|--input-type)  INPUT_TYPE="$2";  shift 2 ;;
        -e|--extra)       EXTRA_ARGS="$2";  shift 2 ;;
        -h|--help)        usage; exit 0 ;;
        *)  echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ---------- sanity checks ---------------------------------------------------
command -v metaphlan >/dev/null 2>&1 || {
    echo "ERROR: 'metaphlan' not found in PATH." >&2; exit 1;
}

if [[ -z "$OUT_DIR" ]]; then
    echo "ERROR: -o/--output is required." >&2; usage; exit 1
fi

# ---------- helper: run a single profiling ---------------------------------
run_one_metaphlan() {
    local r1="$1"
    local r2="$2"
    local outdir="$3"
    local acc="$4"

    if [[ ! -f "$r1" ]]; then echo "ERROR: R1 not found: $r1" >&2; return 1; fi
    if [[ ! -f "$r2" ]]; then echo "ERROR: R2 not found: $r2" >&2; return 1; fi

    mkdir -p "$outdir"

    local profile_out="${outdir}/${acc}_metaphlan.txt"
    local mapout="${outdir}/${acc}.mapout.bz2"

    if [[ -f "$profile_out" ]]; then
        echo "WARN: profile already exists — skipping: $profile_out" >&2
        return 0
    fi

    # MetaPhlAn refuses to overwrite an existing mapout file; remove stale ones.
    [[ -f "$mapout" ]] && rm -f "$mapout"

    # MetaPhlAn does not accept gzipped input — decompress to a temp dir if needed.
    local tmpdir=""
    local r1_in="$r1"
    local r2_in="$r2"
    if [[ "$r1" == *.gz || "$r2" == *.gz ]]; then
        tmpdir="$(mktemp -d "${outdir}/tmp_decompressed_XXXXXX")"
        # Ensure cleanup of decompressed reads even on error/interrupt.
        trap 'rm -rf "$tmpdir"' RETURN
        if [[ "$r1" == *.gz ]]; then
            r1_in="${tmpdir}/${acc}_1.fq"
            echo "  decompressing R1 -> $r1_in"
            gunzip -c "$r1" > "$r1_in"
        fi
        if [[ "$r2" == *.gz ]]; then
            r2_in="${tmpdir}/${acc}_2.fq"
            echo "  decompressing R2 -> $r2_in"
            gunzip -c "$r2" > "$r2_in"
        fi
    fi

    local cmd=( metaphlan
                "${r1_in},${r2_in}"
                --bowtie2out "$mapout"
               	--bowtie2db "$METAPHLAN_DB_DIR"
                --index "$INDEX_DB"
                --nproc "$NPROC"
                --input_type "$INPUT_TYPE"
                -o "$profile_out" )

    # shellcheck disable=SC2206
    [[ -n "$EXTRA_ARGS" ]] && cmd+=( $EXTRA_ARGS )

    echo "======================================================================"
    echo "[$(date '+%F %T')] Profiling: $acc"
    echo "  R1:     $r1"
    echo "  R2:     $r2"
    echo "  Out:    $profile_out"
    echo "  Cmd:    ${cmd[*]}"
    echo "======================================================================"

    "${cmd[@]}"
}

# ---------- dispatch: single vs batch mode ----------------------------------
if [[ -n "$R1" || -n "$R2" || -n "$ACCESSION" ]]; then
    # Single-sample mode
    if [[ -z "$R1" || -z "$R2" || -z "$ACCESSION" ]]; then
        echo "ERROR: single-sample mode needs -1, -2 and -a." >&2; exit 1
    fi
    if [[ -n "$PREPROC_DIR" || -n "$BODYSITE" ]]; then
        echo "ERROR: do not mix -1/-2/-a with -p/-b." >&2; exit 1
    fi
    run_one_metaphlan "$R1" "$R2" "$OUT_DIR" "$ACCESSION"
    echo "[$(date '+%F %T')] Done."
    exit 0
fi

# Batch mode from here on
if [[ -z "$PREPROC_DIR" || -z "$BODYSITE" ]]; then
    echo "ERROR: batch mode needs -p and -b (or use -1/-2/-a for single-sample)." >&2
    usage; exit 1
fi

BODYSITE_DIR="${PREPROC_DIR%/}/${BODYSITE}"
if [[ ! -d "$BODYSITE_DIR" ]]; then
    echo "ERROR: bodysite directory not found: $BODYSITE_DIR" >&2; exit 1
fi

METAPHLAN_BODYSITE_DIR="${OUT_DIR%/}/${BODYSITE}"
mkdir -p "$METAPHLAN_BODYSITE_DIR"

# Collect accessions = immediate subdirs that contain a 'cleaned' folder.
mapfile -t ACCESSIONS < <(
    find "$BODYSITE_DIR" -mindepth 1 -maxdepth 1 -type d \
        -exec test -d '{}/cleaned' \; -print | sort
)

if [[ ${#ACCESSIONS[@]} -eq 0 ]]; then
    echo "ERROR: no accession directories with a 'cleaned/' subfolder found in $BODYSITE_DIR" >&2
    exit 1
fi

echo "Found ${#ACCESSIONS[@]} accession(s) under $BODYSITE_DIR"

FAILED=()
for acc_path in "${ACCESSIONS[@]}"; do
    acc="$(basename "$acc_path")"
    r1="${acc_path}/cleaned/${acc}_1_clean.fq.gz"
    r2="${acc_path}/cleaned/${acc}_2_clean.fq.gz"
    outdir="${METAPHLAN_BODYSITE_DIR}/${acc}"

    if ! run_one_metaphlan "$r1" "$r2" "$outdir" "$acc"; then
        echo "WARN: metaphlan failed for ${acc}" >&2
        FAILED+=( "$acc" )
    fi
done

echo "======================================================================"
echo "[$(date '+%F %T')] Batch finished."
echo "  Bodysite:   $BODYSITE"
echo "  Succeeded:  $(( ${#ACCESSIONS[@]} - ${#FAILED[@]} )) / ${#ACCESSIONS[@]}"
if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "  Failed:     ${FAILED[*]}"
    exit 2
fi
