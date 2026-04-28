#!/bin/bash
#SBATCH --job-name=HostReadSplit
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem-per-cpu=6G    # memory per cpu-core
#SBATCH --time=03:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

# =============================================================================
# Host-Microbiota Split Script (samtools)
# -----------------------------------------------------------------------------
# Splits the host BAM from map_host.sh into:
#   - microbial reads (both mates unmapped) -> paired gzipped FASTQs
#   - host reads      (at least one mapped) -> sorted BAM
# and records read/base counts for each class. All outputs are written to a
# new "cleaned/" subdirectory inside the accession folder.
#
# Expected input (from map_host.sh):
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>.host.bam
#
# Output (written to <ACCESSION>/cleaned/):
#   <ACCESSION>_1_clean.fq.gz        microbial R1
#   <ACCESSION>_2_clean.fq.gz        microbial R2
#   <ACCESSION>.metareads            count of microbial reads
#   <ACCESSION>.metabases            total microbial bases
#   <ACCESSION>.host.sorted.bam      sorted BAM of host-mapped reads
#   <ACCESSION>.hostreads            count of host reads
#   <ACCESSION>.hostbases            total host bases
#
# Usage:
#   Single: ./split_reads.sh -a <ACCESSION> -b <BODYSITE> [-t N] [-p PREP_DIR]
#   Batch : ./split_reads.sh -A [-b <BODYSITE>] [-t N] [-p PREP_DIR]
# =============================================================================

# ---- 0. Load Modules -------------------------------------------------------

module load bowtie2/2.4.2 samtools/1.21

# ---- 1. Defaults ------------------------------------------------------------
ACCESSION=""
BODYSITE=""
BATCH_MODE=0
THREADS=11
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

# ---- 3. Core function: split ONE BAM into host + microbial -----------------
split_one() {
    local bodysite="$1"
    local accession="$2"

    local sample_dir="${PREP_DIR}/${bodysite}/${accession}"
    local in_bam="${sample_dir}/${accession}.host.bam"

    # New: dedicated subfolder for the final cleaned outputs.
    local clean_dir="${sample_dir}/cleaned"

    # Microbial outputs (note the "_clean" tag in the FASTQ names).
    local out_r1="${clean_dir}/${accession}_1_clean.fq.gz"
    local out_r2="${clean_dir}/${accession}_2_clean.fq.gz"
    local out_metareads="${clean_dir}/${accession}.metareads"
    local out_metabases="${clean_dir}/${accession}.metabases"

    # Host outputs.
    local out_host_bam="${clean_dir}/${accession}.host.sorted.bam"
    local out_hostreads="${clean_dir}/${accession}.hostreads"
    local out_hostbases="${clean_dir}/${accession}.hostbases"

    if [ ! -f "${in_bam}" ]; then
        echo "WARNING: BAM not found for ${bodysite}/${accession}, skipping."
        return 0
    fi

    # Create the cleaned/ subfolder for this sample.
    mkdir -p "${clean_dir}"

    echo "=============================================="
    echo "Accession : ${accession}  (${bodysite})"
    echo "BAM       : ${in_bam}"
    echo "Clean dir : ${clean_dir}"
    echo "=============================================="

    # -------- Microbial side: both mates unmapped (-f 12) -------------------
    # Extract unmapped pairs, name-sort them (required for proper mate
    # pairing in samtools fastq), and write gzipped paired FASTQs.
    echo "[1/3] Extracting microbial (unmapped) read pairs..."
    samtools view -b -f 12 -@ "${THREADS}" "${in_bam}" \
        | samtools sort -n -@ "${THREADS}" - \
        | samtools fastq -@ "${THREADS}" \
            -1 "${out_r1}" -2 "${out_r2}" \
            -0 /dev/null -s /dev/null -n -

    # Count microbial reads and bases.
    echo "[2/3] Counting microbial reads and bases..."
    samtools view -c -f 12 -@ "${THREADS}" "${in_bam}" > "${out_metareads}"
    samtools view    -f 12 -@ "${THREADS}" "${in_bam}" \
        | awk '{sum += length($10)} END {print sum+0}' > "${out_metabases}"

    # -------- Host side: NOT both unmapped (-F 12) --------------------------
    # Keep reads where at least one mate mapped, then coordinate-sort.
    echo "[3/3] Extracting host-mapped reads..."
    samtools view -b -F 12 -@ "${THREADS}" "${in_bam}" \
        | samtools sort -@ "${THREADS}" -o "${out_host_bam}" -

    # Count host reads and bases.
    samtools view -c -F 12 -@ "${THREADS}" "${in_bam}" > "${out_hostreads}"
    samtools view    -F 12 -@ "${THREADS}" "${in_bam}" \
        | awk '{sum += length($10)} END {print sum+0}' > "${out_hostbases}"

    # Small summary printout.
    echo "  microbial reads : $(cat "${out_metareads}")"
    echo "  microbial bases : $(cat "${out_metabases}")"
    echo "  host reads      : $(cat "${out_hostreads}")"
    echo "  host bases      : $(cat "${out_hostbases}")"
    echo "Done: ${clean_dir}"
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
            split_one "${bs_name}" "$(basename "${acc_dir}")"
        done
    done
else
    if [ -z "${ACCESSION}" ] || [ -z "${BODYSITE}" ]; then
        echo "Error: single mode requires -a <ACCESSION> and -b <BODYSITE>."
        usage
    fi
    split_one "${BODYSITE}" "${ACCESSION}"
fi
