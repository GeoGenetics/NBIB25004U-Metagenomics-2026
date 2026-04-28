#!/bin/bash
#SBATCH --job-name=assembly_megahit
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem-per-cpu=10G    # memory per cpu-core
#SBATCH --time=05:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching
#
# run_megahit.sh — Run MEGAHIT assembly on metagenomic paired-end reads.
#
# Two modes of operation:
#
#   (A) Batch mode — assemble every accession of a given bodysite:
#       run_megahit.sh -p <preprocessing_dir> -b <bodysite> -o <assembly_dir> [megahit opts]
#
#   (B) Single-sample mode — assemble one pair of reads directly:
#       run_megahit.sh -1 <R1.fq.gz> -2 <R2.fq.gz> -o <out_dir> [megahit opts]
#
# In batch mode the expected input layout is:
#   <preprocessing_dir>/<bodysite>/<accession>/cleaned/<accession>_1_clean.fq.gz
#   <preprocessing_dir>/<bodysite>/<accession>/cleaned/<accession>_2_clean.fq.gz
#
# And the output layout produced is:
#   <assembly_dir>/<bodysite>/<accession>/   (this is MEGAHIT's -o directory)
#
# MEGAHIT options exposed:
#   -t / --threads            number of threads           (default: 8)
#   -l / --min-contig-len     minimum contig length       (default: 1500)
#   -x / --preset             MEGAHIT --presets value
#                             (meta | meta-sensitive | meta-large | bulk | single-cell)
#   -e / --extra              extra raw args passed to megahit (quoted string)
# ---------------------------------------------------------------------------

# ---- 0. Load Modules -------------------------------------------------------

export PATH=/opt/shared_software/shared_envmodules/conda/megahit-1.2.9/bin:$PATH
megahit --version    # sanity check

set -euo pipefail

# ---------- defaults --------------------------------------------------------
THREADS=11
MIN_CONTIG_LEN=1500
PRESET=""
EXTRA_ARGS=""

PREPROC_DIR=""
BODYSITE=""
R1=""
R2=""
OUT_DIR=""

# ---------- usage -----------------------------------------------------------
usage() {
    cat <<EOF
Usage:
  Batch mode:
    $(basename "$0") -p <preprocessing_dir> -b <bodysite> -o <assembly_dir> [options]

  Single-sample mode:
    $(basename "$0") -1 <R1.fq.gz> -2 <R2.fq.gz> -o <out_dir> [options]

Required (batch mode):
  -p, --preproc DIR        Preprocessing directory (contains <bodysite>/<acc>/cleaned/...)
  -b, --bodysite NAME      Bodysite subdirectory (e.g. gut_adult, gut_infant, vaginal)
  -o, --output DIR         Root assembly output directory

Required (single-sample mode):
  -1, --reads1 FILE        Forward reads (fastq[.gz])
  -2, --reads2 FILE        Reverse reads (fastq[.gz])
  -o, --output DIR         MEGAHIT output directory (must not exist)

MEGAHIT options:
  -t, --threads INT        Threads                         (default: ${THREADS})
  -l, --min-contig-len INT Minimum contig length           (default: ${MIN_CONTIG_LEN})
  -x, --preset STR         MEGAHIT --presets value:
                             meta | meta-sensitive | meta-large | bulk | single-cell
  -e, --extra "STR"        Extra raw arguments forwarded to megahit (quoted)

  -h, --help               Show this help and exit
EOF
}

# ---------- arg parsing -----------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--preproc)         PREPROC_DIR="$2"; shift 2 ;;
        -b|--bodysite)        BODYSITE="$2";    shift 2 ;;
        -1|--reads1)          R1="$2";          shift 2 ;;
        -2|--reads2)          R2="$2";          shift 2 ;;
        -o|--output)          OUT_DIR="$2";     shift 2 ;;
        -t|--threads)         THREADS="$2";     shift 2 ;;
        -l|--min-contig-len)  MIN_CONTIG_LEN="$2"; shift 2 ;;
        -x|--preset)          PRESET="$2";      shift 2 ;;
        -e|--extra)           EXTRA_ARGS="$2";  shift 2 ;;
        -h|--help)            usage; exit 0 ;;
        *)  echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ---------- sanity checks ---------------------------------------------------
command -v megahit >/dev/null 2>&1 || {
    echo "ERROR: 'megahit' not found in PATH." >&2; exit 1;
}

if [[ -n "$PRESET" ]]; then
    case "$PRESET" in
        meta|meta-sensitive|meta-large|bulk|single-cell) ;;
        *) echo "ERROR: invalid --preset '$PRESET'." >&2
           echo "       choose one of: meta, meta-sensitive, meta-large, bulk, single-cell" >&2
           exit 1 ;;
    esac
fi

if [[ -z "$OUT_DIR" ]]; then
    echo "ERROR: -o/--output is required." >&2; usage; exit 1
fi

# ---------- helper: run a single assembly -----------------------------------
run_one_assembly() {
    local r1="$1"
    local r2="$2"
    local outdir="$3"
    local label="$4"   # purely cosmetic, for logging

    if [[ ! -f "$r1" ]]; then echo "ERROR: R1 not found: $r1" >&2; return 1; fi
    if [[ ! -f "$r2" ]]; then echo "ERROR: R2 not found: $r2" >&2; return 1; fi

    # MEGAHIT refuses to run if -o already exists; clean up empty dirs we created.
    if [[ -e "$outdir" ]]; then
        if [[ -d "$outdir" && -z "$(ls -A "$outdir" 2>/dev/null || true)" ]]; then
            rmdir "$outdir"
        else
            echo "WARN: output dir already exists and is not empty — skipping: $outdir" >&2
            return 0
        fi
    fi

    # Parent must exist; MEGAHIT creates the final -o dir itself.
    mkdir -p "$(dirname "$outdir")"

    local cmd=( megahit
                -t "$THREADS"
                --verbose
                --min-contig-len "$MIN_CONTIG_LEN"
                -1 "$r1"
                -2 "$r2"
                -o "$outdir" )

    [[ -n "$PRESET" ]] && cmd+=( --presets "$PRESET" )
    # shellcheck disable=SC2206
    [[ -n "$EXTRA_ARGS" ]] && cmd+=( $EXTRA_ARGS )

    echo "======================================================================"
    echo "[$(date '+%F %T')] Assembling: $label"
    echo "  R1:     $r1"
    echo "  R2:     $r2"
    echo "  Out:    $outdir"
    echo "  Cmd:    ${cmd[*]}"
    echo "======================================================================"

    "${cmd[@]}"
}

# ---------- dispatch: single vs batch mode ----------------------------------
if [[ -n "$R1" || -n "$R2" ]]; then
    # Single-sample mode
    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "ERROR: single-sample mode needs both -1 and -2." >&2; exit 1
    fi
    if [[ -n "$PREPROC_DIR" || -n "$BODYSITE" ]]; then
        echo "ERROR: do not mix -1/-2 with -p/-b." >&2; exit 1
    fi
    run_one_assembly "$R1" "$R2" "$OUT_DIR" "$(basename "$OUT_DIR")"
    echo "[$(date '+%F %T')] Done."
    exit 0
fi

# Batch mode from here on
if [[ -z "$PREPROC_DIR" || -z "$BODYSITE" ]]; then
    echo "ERROR: batch mode needs -p and -b (or use -1/-2 for single-sample)." >&2
    usage; exit 1
fi

BODYSITE_DIR="${PREPROC_DIR%/}/${BODYSITE}"
if [[ ! -d "$BODYSITE_DIR" ]]; then
    echo "ERROR: bodysite directory not found: $BODYSITE_DIR" >&2; exit 1
fi

ASSEMBLY_BODYSITE_DIR="${OUT_DIR%/}/${BODYSITE}"
mkdir -p "$ASSEMBLY_BODYSITE_DIR"

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
    outdir="${ASSEMBLY_BODYSITE_DIR}/${acc}"

    if ! run_one_assembly "$r1" "$r2" "$outdir" "${BODYSITE}/${acc}"; then
        echo "WARN: assembly failed for ${acc}" >&2
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
