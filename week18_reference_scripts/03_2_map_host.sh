#!/bin/bash
#SBATCH --job-name=mapHost
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10G    # memory per cpu-core
#SBATCH --time=03:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

# =============================================================================
# Host Genome Mapping Script (bowtie2)
# -----------------------------------------------------------------------------
# Maps fastp-cleaned paired-end reads against a bowtie2 index and writes
# a sorted BAM next to the cleaned reads.
#
# Expected input layout (produced by run_fastp.sh):
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_1.trimmed.fastq.gz
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_2.trimmed.fastq.gz
#
# Output (written next to the cleaned reads):
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>.host.bam
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>.bowtie2.log
#
# Usage:
#   Single: ./map_host.sh -a <ACCESSION> -b <BODYSITE> -i <INDEX_DIR> [-t N] [-p PREP_DIR]
#   Batch : ./map_host.sh -A -i <INDEX_DIR> [-b <BODYSITE>] [-t N] [-p PREP_DIR]
#
# Flags:
#   -a  run accession (required in single mode)
#   -b  bodysite folder name (required in single mode, optional in batch)
#   -A  batch mode: process every accession dir under PREP_DIR
#   -i  directory containing the bowtie2 index (required)
#   -t  threads for bowtie2/samtools     [default: 20]
#   -p  root dir with fastp-cleaned reads (required)
#   -h  show this help message
# =============================================================================

# ---- 0. Load Modules -------------------------------------------------------

module load bowtie2/2.4.2 samtools/1.21

# ---- 1. Defaults ------------------------------------------------------------
ACCESSION=""
BODYSITE=""
BATCH_MODE=0
THREADS=20
PREP_DIR=""
INDEX_DIR=""

usage() {
    echo "Usage:"
    echo "  Single: $0 -a <ACCESSION> -b <BODYSITE> -i <INDEX_DIR> [-t N] [-p PREP_DIR]"
    echo "  Batch : $0 -A -i <INDEX_DIR> [-b <BODYSITE>] [-t N] [-p PREP_DIR]"
    exit 1
}

# ---- 2. Parse flags ---------------------------------------------------------
while getopts ":a:b:Ai:t:p:h" opt; do
    case "${opt}" in
        a) ACCESSION="${OPTARG}" ;;
        b) BODYSITE="${OPTARG}"  ;;
        A) BATCH_MODE=1          ;;
        i) INDEX_DIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}"   ;;
        p) PREP_DIR="${OPTARG}"  ;;
        h) usage ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; usage ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; usage ;;
    esac
done

set -e
set -o pipefail

PREP_DIR="${PREP_DIR%/}"
INDEX_DIR="${INDEX_DIR%/}"

# ---- 3. Resolve the bowtie2 index basename ---------------------------------
# The user gives a directory; bowtie2 needs a basename (path without the
# .1.bt2 extension). We find the .1.bt2 (or .1.bt2l) file inside the
# directory and strip the suffix to build the basename.
if [ -z "${INDEX_DIR}" ]; then
    echo "Error: -i <INDEX_DIR> is required."
    usage
fi
if [ ! -d "${INDEX_DIR}" ]; then
    echo "Error: index directory '${INDEX_DIR}' does not exist."
    exit 1
fi

if [ -z "${PREP_DIR}" ]; then
    echo "Error: -p <PREP_DIR> is required."
    usage
fi
if [ ! -d "${PREP_DIR}" ]; then
    echo "Error: base preprocessing directory '${PREP_DIR}' does not exist."
    exit 1
fi

INDEX_FILE=$(ls "${INDEX_DIR}"/*.1.bt2 "${INDEX_DIR}"/*.1.bt2l 2>/dev/null | head -n1 || true)
if [ -z "${INDEX_FILE}" ]; then
    echo "Error: no bowtie2 index found in '${INDEX_DIR}' (expected *.1.bt2)."
    exit 1
fi
# Strip the .1.bt2 or .1.bt2l suffix to get the basename bowtie2 expects.
INDEX="${INDEX_FILE%.1.bt2}"
INDEX="${INDEX%.1.bt2l}"

echo "Using bowtie2 index: ${INDEX}"

# ---- 4. Core function: map ONE accession -----------------------------------
map_host_one() {
    local bodysite="$1"
    local accession="$2"

    local sample_dir="${PREP_DIR}/${bodysite}/${accession}"
    local in_r1="${sample_dir}/${accession}_1.trimmed.fastq.gz"
    local in_r2="${sample_dir}/${accession}_2.trimmed.fastq.gz"
    local out_bam="${sample_dir}/${accession}.host.bam"
    local bt2_log="${sample_dir}/${accession}.bowtie2.log"

    if [ ! -f "${in_r1}" ] || [ ! -f "${in_r2}" ]; then
        echo "WARNING: missing cleaned FASTQ for ${bodysite}/${accession}, skipping."
        return 0
    fi

    echo "=============================================="
    echo "Accession : ${accession}  (${bodysite})"
    echo "Threads   : ${THREADS}"
    echo "Output    : ${out_bam}"
    echo "=============================================="

    # bowtie2 -> SAM to BAM -> sort.
    # Stderr goes to the per-sample bowtie2 log (alignment summary).
    bowtie2 -x "${INDEX}" \
            -1 "${in_r1}" -2 "${in_r2}" \
            -p "${THREADS}" \
            2> "${bt2_log}" \
        | samtools view -bS -@ "${THREADS}" - \
        | samtools sort -@ "${THREADS}" -o "${out_bam}" -

    echo "Done: ${out_bam}"
}

# ---- 5. Dispatch: single vs batch ------------------------------------------
if [ "${BATCH_MODE}" -eq 1 ]; then
    # Batch: walk every <bodysite>/<accession>/ under PREP_DIR.
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
            map_host_one "${bs_name}" "$(basename "${acc_dir}")"
        done
    done
else
    if [ -z "${ACCESSION}" ] || [ -z "${BODYSITE}" ]; then
        echo "Error: single mode requires -a <ACCESSION> and -b <BODYSITE>."
        usage
    fi
    map_host_one "${BODYSITE}" "${ACCESSION}"
fi
