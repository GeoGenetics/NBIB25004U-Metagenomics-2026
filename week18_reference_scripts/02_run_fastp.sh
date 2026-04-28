#!/bin/bash
#SBATCH --job-name=fastpQC
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G    # memory per cpu-core
#SBATCH --time=03:00:00
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

# =============================================================================
# fastp Preprocessing Script for Metagenomics Class
# -----------------------------------------------------------------------------
# Runs fastp on paired-end FASTQ files produced by download_sra.sh and
# organises the outputs by bodysite/accession inside a preprocessing folder.
#
# Two modes:
#   SINGLE mode : process one accession (-a + -b)
#   BATCH  mode : process every accession directory found under RAW_DIR (-A).
#                 If -b is also given, only that bodysite is scanned;
#                 otherwise all bodysites under RAW_DIR are scanned.
#
# Expected input layout (produced by download_sra.sh):
#   <RAW_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_1.fastq
#   <RAW_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_2.fastq
#
# Output layout:
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_1.trimmed.fastq.gz
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_2.trimmed.fastq.gz
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_fastp.html
#   <PREP_DIR>/<BODYSITE>/<ACCESSION>/<ACCESSION>_fastp.json
#   <PREP_DIR>/fastp_runs.log   <-- appended parameters for every run
#
# Usage:
#   Single:  ./run_fastp.sh -a <ACCESSION> -b <BODYSITE> [-t N] [-r RAW_DIR] [-p PREP_DIR]
#   Batch:   ./run_fastp.sh -A [-b <BODYSITE>] [-t N] [-r RAW_DIR] [-p PREP_DIR]
#
# Flags:
#   -a  run accession (required in single mode)
#   -b  bodysite folder name (required in single mode, optional in batch)
#   -A  batch mode: process every accession dir under RAW_DIR (or one bodysite)
#   -t  threads for fastp                                  [default: 8]
#   -r  root directory with raw reads                      [default: ./raw_reads]
#   -p  root directory for preprocessing output            [default: ./preprocessing]
#   -h  show this help message
#
# Examples:
#   ./run_fastp.sh -a SRR513791 -b vaginal
#   ./run_fastp.sh -A                           # every bodysite, every accession
#   ./run_fastp.sh -A -b gut_infant -t 16       # every accession in one bodysite
# =============================================================================

# ---- 0. Load Modules -------------------------------------------------------

module load fastp
fastp -v
# ---- 1. Default values ------------------------------------------------------
ACCESSION=""
BODYSITE=""
BATCH_MODE=0
THREADS=8
RAW_DIR="./raw_reads"
PREP_DIR="./preprocessing"

# ---- 2. Helper: print usage and exit ---------------------------------------
usage() {
    echo "Usage:"
    echo "  Single: $0 -a <ACCESSION> -b <BODYSITE> [-t N] [-r RAW_DIR] [-p PREP_DIR]"
    echo "  Batch : $0 -A [-b <BODYSITE>] [-t N] [-r RAW_DIR] [-p PREP_DIR]"
    echo ""
    echo "  -a  run accession (required in single mode)"
    echo "  -b  bodysite folder name (required in single mode, optional in batch)"
    echo "  -A  batch mode: process every accession dir under RAW_DIR"
    echo "  -t  number of threads for fastp (default: 8)"
    echo "  -r  root dir with raw reads (default: ./raw_reads)"
    echo "  -p  root dir for preprocessing output (default: ./preprocessing)"
    echo "  -h  show this help message"
    exit 1
}

# ---- 3. Parse command-line flags -------------------------------------------
# `-A` has no argument (batch flag), so no ":" after it in the optstring.
while getopts ":a:b:At:r:p:h" opt; do
    case "${opt}" in
        a) ACCESSION="${OPTARG}" ;;
        b) BODYSITE="${OPTARG}"  ;;
        A) BATCH_MODE=1          ;;
        t) THREADS="${OPTARG}"   ;;
        r) RAW_DIR="${OPTARG}"   ;;
        p) PREP_DIR="${OPTARG}"  ;;
        h) usage ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; usage ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; usage ;;
    esac
done

# Exit immediately if any command fails.
set -e

# Strip any trailing slashes so the paths we build stay clean.
RAW_DIR="${RAW_DIR%/}"
PREP_DIR="${PREP_DIR%/}"

# Shared log file in the top-level preprocessing folder (created on demand).
mkdir -p "${PREP_DIR}"
RUN_LOG="${PREP_DIR}/fastp_runs.log"

# ---- 4. Core function: run fastp on ONE accession --------------------------
# Takes two arguments: bodysite and accession. Builds paths, runs fastp,
# and appends an entry to the shared run log.
run_fastp_one() {
    local bodysite="$1"
    local accession="$2"

    local in_r1="${RAW_DIR}/${bodysite}/${accession}/${accession}_1.fastq"
    local in_r2="${RAW_DIR}/${bodysite}/${accession}/${accession}_2.fastq"
    local out_dir="${PREP_DIR}/${bodysite}/${accession}"
    local out_r1="${out_dir}/${accession}_1.trimmed.fastq.gz"
    local out_r2="${out_dir}/${accession}_2.trimmed.fastq.gz"
    local out_html="${out_dir}/${accession}_fastp_report.html"
    local out_json="${out_dir}/${accession}_fastp_report.json"

    # Sanity check: both paired-end files must exist. In batch mode we warn
    # and skip rather than aborting, so one bad sample doesn't kill the run.
    if [ ! -f "${in_r1}" ] || [ ! -f "${in_r2}" ]; then
        echo "WARNING: missing paired-end FASTQ for ${bodysite}/${accession}, skipping."
        echo "  expected: ${in_r1}"

        echo "  expected: ${in_r2}"
        return 0
    fi

    mkdir -p "${out_dir}"

    echo "=============================================="
    echo "Accession : ${accession}"
    echo "Bodysite  : ${bodysite}"
    echo "Threads   : ${THREADS}"
    echo "Output dir: ${out_dir}"
    echo "=============================================="

    # Run fastp with the class's standard parameters.
    fastp \
        --in1 "${in_r1}" --in2 "${in_r2}" \
        --out1 "${out_r1}" --out2 "${out_r2}" \
        --trim_poly_g \
        --trim_poly_x \
        --low_complexity_filter \
        --n_base_limit 5 \
        --qualified_quality_phred 20 \
        --length_required 60 \
        --thread "${THREADS}" \
        --html "${out_html}" \
        --json "${out_json}"


    # Append parameters to the shared run log (audit trail).
    {
        echo "----------------------------------------------"
        echo "Timestamp : $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Accession : ${accession}"
        echo "Bodysite  : ${bodysite}"
        echo "Input R1  : ${in_r1}"
        echo "Input R2  : ${in_r2}"
        echo "Output R1 : ${out_r1}"
        echo "Output R2 : ${out_r2}"
        echo "HTML      : ${out_html}"
        echo "JSON      : ${out_json}"
        echo "Threads   : ${THREADS}"
        echo "Parameters:"
        echo "  --trim_poly_g"
        echo "  --trim_poly_x"
        echo "  --low_complexity_filter"
        echo "  --n_base_limit 5"
        echo "  --qualified_quality_phred 20"
        echo "  --length_required 60"
        echo "fastp version: $(fastp --version 2>&1)"
        echo "----------------------------------------------"
    } >> "${RUN_LOG}"

    echo "Done: ${out_dir}"
}

# ---- 5. Dispatch: single vs batch mode -------------------------------------
if [ "${BATCH_MODE}" -eq 1 ]; then
    # -------- BATCH MODE ----------------------------------------------------
    # Decide which bodysite directories to scan:
    #   - if -b was given, scan only that one
    #   - otherwise, scan every subdirectory of RAW_DIR
    if [ -n "${BODYSITE}" ]; then
        bodysite_dirs=( "${RAW_DIR}/${BODYSITE}" )
    else
        bodysite_dirs=( "${RAW_DIR}"/*/ )
    fi

    # Counters so we can print a summary at the end.
    total=0
    processed=0

    for bs_dir in "${bodysite_dirs[@]}"; do
        # Skip anything that isn't a real directory (handles empty globs).
        [ -d "${bs_dir}" ] || continue

        # Extract the bodysite name from the directory path.
        bs_name="$(basename "${bs_dir}")"

        # Iterate every accession subdirectory inside this bodysite.
        for acc_dir in "${bs_dir%/}"/*/; do
            [ -d "${acc_dir}" ] || continue
            acc_name="$(basename "${acc_dir}")"
            total=$((total + 1))
            echo ""
            echo ">>> [${total}] ${bs_name}/${acc_name}"
            run_fastp_one "${bs_name}" "${acc_name}"
            processed=$((processed + 1))
        done
    done

    echo ""
    echo "=============================================="
    echo "Batch complete. Accession dirs seen: ${total}"
    echo "Run log: ${RUN_LOG}"
    echo "=============================================="

else
    # -------- SINGLE MODE ---------------------------------------------------
    if [ -z "${ACCESSION}" ] || [ -z "${BODYSITE}" ]; then
        echo "Error: single mode requires -a <ACCESSION> and -b <BODYSITE>."
        echo "       (Use -A for batch mode.)"
        usage
    fi
    run_fastp_one "${BODYSITE}" "${ACCESSION}"
fi
