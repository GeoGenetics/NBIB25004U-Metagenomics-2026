#!/bin/bash
#SBATCH --job-name=bowtie2_build
#SBATCH --output=/maps/projects/course_1/people/fvb335/logs/%x_%j.out
#SBATCH --error=/maps/projects/course_1/people/fvb335/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=31
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=jonas.kasmanas@bio.ku.dk
#SBATCH --reservation=NBIB25004U
#SBATCH --account=teaching

set -euo pipefail

usage() {
    echo "Usage: sbatch $0 -i <path/to/genome.fa> [-o <output_dir>]"
    echo "  -i   Input FASTA file (required)"
    echo "  -o   Output directory for index (default: same as input file)"
    exit 1
}

FA_INPUT=""
OUT_DIR=""

while getopts ":i:o:h" opt; do
    case "$opt" in
        i) FA_INPUT="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

[[ -z "$FA_INPUT" ]] && { echo "Error: -i is required"; usage; }
[[ ! -f "$FA_INPUT" ]] && { echo "Error: input file not found: $FA_INPUT"; exit 1; }

# Default output dir = directory of input file
if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR="$(dirname "$(readlink -f "$FA_INPUT")")"
fi
mkdir -p "$OUT_DIR"

# Build index prefix from basename minus extension
BASENAME="$(basename "$FA_INPUT")"
PREFIX="${BASENAME%.fa}"
PREFIX="${PREFIX%.fasta}"
PREFIX="${PREFIX%.fna}"
INDEX_PREFIX="${OUT_DIR}/${PREFIX}"

module load bowtie2/2.4.2

echo "Input FASTA : $FA_INPUT"
echo "Output dir  : $OUT_DIR"
echo "Index prefix: $INDEX_PREFIX"
echo "Threads     : ${SLURM_CPUS_PER_TASK}"

bowtie2-build "$FA_INPUT" "$INDEX_PREFIX" --threads "${SLURM_CPUS_PER_TASK}"
