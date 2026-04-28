#!/bin/bash
#SBATCH --job-name=binning_metabat2
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=10G    # memory per cpu-core
#SBATCH --time=05:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching
#
# run_metabat2.sh — Bin a MEGAHIT assembly with MetaBAT2.
#
# For each sample, this script runs FOUR steps:
#   1. assembly_index     — build a bowtie2 index of the assembly
#   2. assembly_map       — map the cleaned reads back to the assembly
#   3. assembly_map_depth — compute per-contig depth (metabat + maxbin formats)
#   4. metabat2           — bin the contigs with MetaBAT2
#
# ---------------------------------------------------------------------------
# Expected INPUT layout
# ---------------------------------------------------------------------------
#   <preprocessing_dir>/<bodysite>/<acc>/cleaned/<acc>_1_clean.fq.gz
#   <preprocessing_dir>/<bodysite>/<acc>/cleaned/<acc>_2_clean.fq.gz
#   <assembly_dir>/<bodysite>/<acc>/final.contigs.fa     (from MEGAHIT)
#
# ---------------------------------------------------------------------------
# OUTPUT layout (produced under <output_dir>/<bodysite>/<acc>/)
# ---------------------------------------------------------------------------
#   bowtie2_index/<acc>.*.bt2        bowtie2 index files
#   <acc>.bam                        sorted BAM of reads vs assembly
#   <acc>_metabat.depth              depth file for metabat2
#   <acc>_maxbin.depth               depth file for maxbin2 (bonus)
#   bins/<acc>_bin.<N>.fa            one fasta per bin
#   <acc>.tsv                        list of all bin fastas (summary)
#
# ---------------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------------
#   Batch mode (all accessions of a bodysite):
#     run_metabat2.sh \
#         -p <preprocessing_dir> \
#         -a <assembly_dir> \
#         -b <bodysite> \
#         -o <output_dir> \
#         [-t THREADS] [-m MIN_CONTIG]
#
#   Single-sample mode:
#     run_metabat2.sh \
#         -f <assembly.fa> \
#         -1 <R1.fq.gz> -2 <R2.fq.gz> \
#         -n <accession> \
#         -o <out_dir_for_this_sample> \
#         [-t THREADS] [-m MIN_CONTIG]
# ---------------------------------------------------------------------------

# ---- 0. Load Modules -------------------------------------------------------

module load bowtie2/2.4.2 samtools/1.21
export PATH=/opt/shared_software/shared_envmodules/conda/metabat2-2.17/bin:$PATH

set -euo pipefail

# ---------- defaults --------------------------------------------------------
THREADS=15
MIN_CONTIG=1500
ASSEMBLY_FILENAME="final.contigs.fa"   # MEGAHIT default

PREPROC_DIR=""
ASSEMBLY_DIR=""
BODYSITE=""

ASSEMBLY_FA=""
R1=""
R2=""
ACCESSION=""
OUT_DIR=""

# ---------- usage -----------------------------------------------------------
usage() {
    cat <<EOF
Usage:
  Batch mode:
    $(basename "$0") -p <preproc_dir> -a <assembly_dir> -b <bodysite> -o <output_dir> [options]

  Single-sample mode:
    $(basename "$0") -f <assembly.fa> -1 <R1.fq.gz> -2 <R2.fq.gz> -n <accession> -o <out_dir> [options]

Batch-mode inputs:
  -p, --preproc DIR      Preprocessing dir with <bodysite>/<acc>/cleaned/*_clean.fq.gz
  -a, --assembly DIR     Assembly dir with <bodysite>/<acc>/${ASSEMBLY_FILENAME}
  -b, --bodysite NAME    Bodysite subfolder (gut_adult | gut_infant | vaginal | ...)
  -o, --output DIR       Root output dir (metabat2_binning)

Single-sample inputs:
  -f, --fasta FILE       Assembly FASTA
  -1, --reads1 FILE      Forward reads (fastq[.gz])
  -2, --reads2 FILE      Reverse reads (fastq[.gz])
  -n, --accession NAME   Accession ID (used to name all outputs)
  -o, --output DIR       Output dir for this sample

Options:
  -t, --threads INT      Threads for bowtie2 / samtools   (default: ${THREADS})
  -m, --min-contig INT   MetaBAT2 minimum contig length   (default: ${MIN_CONTIG})
  -h, --help             Show this help and exit
EOF
}

# ---------- arg parsing -----------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--preproc)    PREPROC_DIR="$2";  shift 2 ;;
        -a|--assembly)   ASSEMBLY_DIR="$2"; shift 2 ;;
        -b|--bodysite)   BODYSITE="$2";     shift 2 ;;
        -f|--fasta)      ASSEMBLY_FA="$2";  shift 2 ;;
        -1|--reads1)     R1="$2";           shift 2 ;;
        -2|--reads2)     R2="$2";           shift 2 ;;
        -n|--accession)  ACCESSION="$2";    shift 2 ;;
        -o|--output)     OUT_DIR="$2";      shift 2 ;;
        -t|--threads)    THREADS="$2";      shift 2 ;;
        -m|--min-contig) MIN_CONTIG="$2";   shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ---------- sanity checks ---------------------------------------------------
for tool in bowtie2-build bowtie2 samtools jgi_summarize_bam_contig_depths metabat2; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: '$tool' not found in PATH." >&2
        echo "       Make sure bowtie2, samtools and metabat2 modules are loaded." >&2
        exit 1
    }
done

if [[ -z "$OUT_DIR" ]]; then
    echo "ERROR: -o/--output is required." >&2; usage; exit 1
fi

# ===========================================================================
# CORE FUNCTION: process one sample through all four steps.
# ===========================================================================
run_one_sample() {
    local assembly_fa="$1"
    local r1="$2"
    local r2="$3"
    local acc="$4"
    local sample_out="$5"   # final dir for this accession

    echo "======================================================================"
    echo "[$(date '+%F %T')] Sample: $acc"
    echo "  assembly: $assembly_fa"
    echo "  R1:       $r1"
    echo "  R2:       $r2"
    echo "  outdir:   $sample_out"
    echo "======================================================================"

    # --- Input checks ------------------------------------------------------
    [[ -f "$assembly_fa" ]] || { echo "ERROR: assembly not found: $assembly_fa" >&2; return 1; }
    [[ -f "$r1" ]]          || { echo "ERROR: R1 not found: $r1" >&2; return 1; }
    [[ -f "$r2" ]]          || { echo "ERROR: R2 not found: $r2" >&2; return 1; }

    # --- Paths -------------------------------------------------------------
    local index_dir="${sample_out}/bowtie2_index"
    local index_base="${index_dir}/${acc}"
    local bam="${sample_out}/${acc}.bam"
    local depth_metabat="${sample_out}/${acc}_metabat.depth"
    local depth_maxbin="${sample_out}/${acc}_maxbin.depth"
    local bins_dir="${sample_out}/bins"
    local bin_base="${bins_dir}/${acc}_bin"
    local summary_tsv="${sample_out}/${acc}.tsv"

    mkdir -p "$index_dir" "$bins_dir"

    # ------------------------------------------------------------------
    # STEP 1 — assembly_index: build bowtie2 index of the assembly
    # ------------------------------------------------------------------
    if [[ -f "${index_base}.rev.2.bt2" ]]; then
        echo "[step 1] index already exists — skipping bowtie2-build"
    else
        echo "[step 1] bowtie2-build ${assembly_fa}"
        bowtie2-build --threads "$THREADS" "$assembly_fa" "$index_base"
    fi

    # ------------------------------------------------------------------
    # STEP 2 — assembly_map: map cleaned reads back to the assembly
    # ------------------------------------------------------------------
    if [[ -f "$bam" ]]; then
        echo "[step 2] BAM already exists — skipping mapping"
    else
        echo "[step 2] bowtie2 mapping -> sorted BAM"
        bowtie2 -p "$THREADS" -x "$index_base" -1 "$r1" -2 "$r2" \
            | samtools view -bS - \
            | samtools sort -@ "$THREADS" -o "$bam" -
        samtools index "$bam"
    fi

    # ------------------------------------------------------------------
    # STEP 3 — assembly_map_depth: per-contig depth for metabat + maxbin
    # ------------------------------------------------------------------
    if [[ -f "$depth_metabat" && -f "$depth_maxbin" ]]; then
        echo "[step 3] depth files already exist — skipping"
    else
        echo "[step 3] jgi_summarize_bam_contig_depths"
        jgi_summarize_bam_contig_depths \
            --outputDepth "$depth_metabat" "$bam"
        # maxbin wants only contig name + mean depth, no header
        cut -f1,3 "$depth_metabat" | tail -n +2 > "$depth_maxbin"
    fi

    # ------------------------------------------------------------------
    # STEP 4 — metabat2: bin the contigs
    # ------------------------------------------------------------------
    echo "[step 4] metabat2 binning"
    # Clean any stale bins from a previous partial run
    rm -f "${bin_base}."*.fa
    metabat2 \
        -i "$assembly_fa" \
        -a "$depth_metabat" \
        -o "$bin_base" \
        -m "$MIN_CONTIG" \
        -t "$THREADS" \
        --saveCls

    # Summary file: list every bin fasta we just produced
    find "$bins_dir" -maxdepth 1 -type f -name "${acc}_bin.*.fa" | sort > "$summary_tsv"

    local n_bins
    n_bins=$(wc -l < "$summary_tsv")
    echo "[done] $acc — produced ${n_bins} bin(s)"
}

# ===========================================================================
# Dispatch: single-sample vs batch mode
# ===========================================================================
if [[ -n "$ASSEMBLY_FA" || -n "$R1" || -n "$R2" || -n "$ACCESSION" ]]; then
    # --- single-sample mode ---
    if [[ -z "$ASSEMBLY_FA" || -z "$R1" || -z "$R2" || -z "$ACCESSION" ]]; then
        echo "ERROR: single-sample mode needs -f, -1, -2 and -n." >&2; exit 1
    fi
    if [[ -n "$PREPROC_DIR" || -n "$ASSEMBLY_DIR" || -n "$BODYSITE" ]]; then
        echo "ERROR: do not mix single-sample flags (-f/-1/-2/-n) with batch flags (-p/-a/-b)." >&2
        exit 1
    fi
    mkdir -p "$OUT_DIR"
    run_one_sample "$ASSEMBLY_FA" "$R1" "$R2" "$ACCESSION" "$OUT_DIR"
    exit 0
fi

# --- batch mode ---
if [[ -z "$PREPROC_DIR" || -z "$ASSEMBLY_DIR" || -z "$BODYSITE" ]]; then
    echo "ERROR: batch mode needs -p, -a and -b (or use single-sample flags)." >&2
    usage; exit 1
fi

PREPROC_BODYSITE="${PREPROC_DIR%/}/${BODYSITE}"
ASSEMBLY_BODYSITE="${ASSEMBLY_DIR%/}/${BODYSITE}"
OUT_BODYSITE="${OUT_DIR%/}/${BODYSITE}"

[[ -d "$PREPROC_BODYSITE" ]]  || { echo "ERROR: missing $PREPROC_BODYSITE"  >&2; exit 1; }
[[ -d "$ASSEMBLY_BODYSITE" ]] || { echo "ERROR: missing $ASSEMBLY_BODYSITE" >&2; exit 1; }
mkdir -p "$OUT_BODYSITE"

# Discover accessions from the preprocessing directory (same convention as
# the earlier scripts: an accession is a dir that has a 'cleaned/' subfolder).
mapfile -t ACCESSIONS < <(
    find "$PREPROC_BODYSITE" -mindepth 1 -maxdepth 1 -type d \
        -exec test -d '{}/cleaned' \; -print | sort
)

if [[ ${#ACCESSIONS[@]} -eq 0 ]]; then
    echo "ERROR: no accessions with cleaned/ found under $PREPROC_BODYSITE" >&2
    exit 1
fi

echo "Found ${#ACCESSIONS[@]} accession(s) under $PREPROC_BODYSITE"

FAILED=()
for acc_path in "${ACCESSIONS[@]}"; do
    acc="$(basename "$acc_path")"

    r1="${acc_path}/cleaned/${acc}_1_clean.fq.gz"
    r2="${acc_path}/cleaned/${acc}_2_clean.fq.gz"
    assembly_fa="${ASSEMBLY_BODYSITE}/${acc}/${ASSEMBLY_FILENAME}"
    sample_out="${OUT_BODYSITE}/${acc}"

    mkdir -p "$sample_out"

    if ! run_one_sample "$assembly_fa" "$r1" "$r2" "$acc" "$sample_out"; then
        echo "WARN: sample failed: $acc" >&2
        FAILED+=( "$acc" )
    fi
done

echo "======================================================================"
echo "[$(date '+%F %T')] Batch finished."
echo "  Bodysite:  $BODYSITE"
echo "  Succeeded: $(( ${#ACCESSIONS[@]} - ${#FAILED[@]} )) / ${#ACCESSIONS[@]}"
if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "  Failed:    ${FAILED[*]}"
    exit 2
fi
