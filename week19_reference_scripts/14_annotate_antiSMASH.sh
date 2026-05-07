#!/bin/bash
#SBATCH --job-name=antiSMASH
#SBATCH --output=/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/logs/antiSMASH_%j.out
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/logs/antiSMASH_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=10:00:00

# antiSMASH BGC annotation per Bakta-annotated genome.
# Submit with: sbatch 14_annotate_antismash.sh

module purge
module load antismash/8.0.1

##########################
###### first genome ######
##########################

# 1. Pick the sample (change this to your genome's basename)
sample1="Lactobacillus_crispatus"

# 2. Set paths
INPUT1="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/09_annotation_bakta_ref/${sample1}/${sample1}.gbff"
OUTDIR1="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/11_annotation_BGC_ref/${sample1}"

# 3. Create output directory
mkdir -p "${OUTDIR1}"

# 4. print status
echo "Starting antiSMASH annotation for ${sample1}"

# 5. Run antiSMASH
antismash --genefinding-tool none --cpus 8 --cb-knownclusters --cb-subclusters --asf --rre --tfbs --output-dir "${OUTDIR1}" "${INPUT1}"

#6. print status
echo =========================================================================
echo "antiSMASH annotation completed for ${sample1}. Output in ${OUTDIR1}"
echo =========================================================================


##########################
###### second genome #####
##########################

sample2="Bifidobacterium_infantis"
INPUT2="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/09_annotation_bakta_ref/${sample2}/${sample2}.gbff"
OUTDIR2="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/11_annotation_BGC_ref/${sample2}"
mkdir -p "${OUTDIR2}"
echo "Starting antiSMASH annotation for ${sample2}"

antismash --genefinding-tool none --cpus 8 --cb-knownclusters --cb-subclusters --asf --rre --tfbs --output-dir "${OUTDIR2}" "${INPUT2}"

echo =========================================================================
echo "antiSMASH annotation completed for ${sample2}. Output in ${OUTDIR2}"
echo =========================================================================


##########################
###### third genome ######
##########################

sample3="Bacteroides_thetaiotaomicron"
INPUT3="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/09_annotation_bakta_ref/${sample3}/${sample3}.gbff"
OUTDIR3="/maps/projects/course_1/scratch/<group_#>/<group-project-group-#>/11_annotation_BGC_ref/${sample3}"
mkdir -p "${OUTDIR3}"
echo "Starting antiSMASH annotation for ${sample3}"

antismash --genefinding-tool none --cpus 8 --cb-knownclusters --cb-subclusters --asf --rre --tfbs --output-dir "${OUTDIR3}" "${INPUT3}"

echo =========================================================================
echo "antiSMASH annotation completed for ${sample3}. Output in ${OUTDIR3}"
echo =========================================================================
