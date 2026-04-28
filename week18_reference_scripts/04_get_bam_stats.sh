#!/bin/bash
#SBATCH --job-name=StatsBam
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G    # memory per cpu-core
#SBATCH --time=03:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

# =============================================================================
# BAM Stats Script (samtools)
# -----------------------------------------------------------------------------
# Runs samtools index/flagstat/idxstats/stats on the host BAM produced by
# map_host.sh and writes the reports next to the BAM.
#
# Expected input (from map_host.sh):
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>.host.bam
#
# Output (written next to the BAM):
#   <ACCESSION>.host.bam.bai     BAM index
#   <ACCESSION>.flagstat.txt     quick mapped/unmapped/paired counts
#   <ACCESSION>.idxstats.txt     per-reference read counts
#   <ACCESSION>.stats.txt        detailed alignment statistics
#
# Usage:
#   Single: ./bam_stats.sh -a <ACCESSION> -b <BODYSITE> [-t N] [-p PREP_DIR]
#   Batch : ./bam_stats.sh -A [-b <BODYSITE>] [-t N] [-p PREP_DIR]
# =============================================================================

# ---- 0. Load Modules -------------------------------------------------------

module load samtools/1.21

# ---- 1. Defaults ------------------------------------------------------------
ACCESSION=""
BODYSITE=""
BATCH_MODE=0
THREADS=8
PREP_DIR=""

usage() {
    echo "Usage:"
    echo "  Single: $0 -a <ACCESSION> -b <BODYSITE> [-t N] [-p PREP_DIR]"
    echo "  Batch : $0 -A [-b <BODYSITE>] [-t N] [-p PREP_DIR]"
    exit 1
}

# ---- 2. Parse flags ---------------------------------------------------------
while getopts ":a:b:At:p:h" opt; do
    case "${opt}" in
        a) ACCESSION="${OPTARG}" ;;
        b) BODYSITE="${OPTARG}"  ;;
        A) BATCH_MODE=1          ;;
        t) THREADS="${OPTARG}"   ;;
        p) PREP_DIR="${OPTARG}"  ;;
        h) usage ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; usage ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; usage ;;
    esac
done

set -e
set -o pipefail

if [ -z "${PREP_DIR}" ]; then
    echo "Error: -p <PREP_DIR> is required."
    usage
fi
if [ ! -d "${PREP_DIR}" ]; then
    echo "Error: base preprocessing directory '${PREP_DIR}' does not exist."
    exit 1
fi

PREP_DIR="${PREP_DIR%/}"

# ---- 3. Core function: compute stats for ONE accession --------------------
stats_one() {
    local bodysite="$1"
    local accession="$2"

    local sample_dir="${PREP_DIR}/${bodysite}/${accession}"
    local in_bam="${sample_dir}/${accession}.host.bam"

    # Four output files, all next to the BAM.
    local out_bai="${in_bam}.bai"
    local out_flagstat="${sample_dir}/${accession}.flagstat.txt"
    local out_idxstats="${sample_dir}/${accession}.idxstats.txt"
    local out_stats="${sample_dir}/${accession}.stats.txt"

    if [ ! -f "${in_bam}" ]; then
        echo "WARNING: BAM not found for ${bodysite}/${accession}, skipping."
        return 0
    fi

    echo "=============================================="
    echo "Accession : ${accession}  (${bodysite})"
    echo "BAM       : ${in_bam}"
    echo "=============================================="

    # 1) Build a .bai index so idxstats (and random-access viewers) work.
    echo "[1/4] samtools index"
    samtools index -@ "${THREADS}" "${in_bam}" "${out_bai}"

    # 2) flagstat -> quick counts of mapped/unmapped/paired/etc.
    echo "[2/4] samtools flagstat"
    samtools flagstat -@ "${THREADS}" "${in_bam}" > "${out_flagstat}"

    # 3) idxstats -> per-reference (per-chromosome) read counts. Uses the .bai.
    echo "[3/4] samtools idxstats"
    samtools idxstats -@ "${THREADS}" "${in_bam}" > "${out_idxstats}"

    # 4) stats -> detailed alignment statistics (insert size, GC, error rate).
    echo "[4/4] samtools stats"
    samtools stats -@ "${THREADS}" "${in_bam}" > "${out_stats}"

    echo "Done: stats written to ${sample_dir}"
}

# ---- 4. Dispatch: single vs batch ------------------------------------------
if [ "${BATCH_MODE}" -eq 1 ]; then
    if [ -n "${BODYSITE}" ]; then
        bodysite_dirs=( "${PREP_DIR}/${BODYSITE}" )
    else
        bodysite_dirs=( "${PREP_DIR}"/*/ )
    fi

    for bs_dir in "${bodysite_dirs[@]}"; do
        [ -d "${bs_dir}" ] || continue
        bs_name="$(basename "${bs_dir}")"
        for acc_dir in "${bs_dir%/}"/*/; do
            [ -d "${acc_dir}" ] || continue
            stats_one "${bs_name}" "$(basename "${acc_dir}")"
        done
    done
else
    if [ -z "${ACCESSION}" ] || [ -z "${BODYSITE}" ]; then
        echo "Error: single mode requires -a <ACCESSION> and -b <BODYSITE>."
        usage
    fi
    stats_one "${BODYSITE}" "${ACCESSION}"
fi
