#!/bin/bash

# --- dsub Job Submission Script (To be run on the VM/Client machine) ---

# NOTE: This script reads gene regions from a file named 'chrom_gene_region.txt'.

# Define Environment Variables (CENSORED PLACEHOLDERS)
# Collaborator must update all variables marked with <...>
DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"
GOOGLE_PROJECT="<YOUR_GOOGLE_PROJECT_ID>"
AOU_NETWORK="global/networks/<NETWORK_NAME>"
AOU_SUBNETWORK="regions/us-central1/subnetworks/<SUBNETWORK_NAME>"
MACHINE_TYPE="n2-standard-16"
WORKSPACE_BUCKET="gs://<YOUR_WORKSPACE_BUCKET>"

# Define GCS Paths for Scripts and Tools (CENSORED)
PIPELINE_SCRIPT="gs://<YOUR_SECURE_BUCKET>/CMP_Gene_Data/Scripts/plink_subsetter.sh"
PLINK_TOOL="gs://<PLINK_TOOL_BUCKET>/plink2"

# Define GCS Paths for Input Datasets (Controlled Tier Plink files - CENSORED BASE PATH)
BASE_INPUT_BED_ACAF="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed/chr"
BASE_INPUT_BED_EXOME="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/plink_bed/exome.chr"
BASE_INPUT_BED_CLINVAR="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/clinvar/plink_bed/chr"

# Define GCS Paths for Output (Your Secure Bucket - CENSORED BASE PATH)
BASE_OUTPUT_PATH="gs://<YOUR_SECURE_BUCKET>/CMP_Gene_Data/Plink_Format"
OUTPUT_PATH1="${BASE_OUTPUT_PATH}/ACAF_threshold/"
OUTPUT_PATH2="${BASE_OUTPUT_PATH}/EXOME_Subset/"
OUTPUT_PATH3="${BASE_OUTPUT_PATH}/Clinvar_Subset/"


while IFS=$'\t' read -r chrom gene region; do
  JOB_NAME="proc_${gene}_chr${chrom}"

  echo "Submitting job for ${gene} on chr${chrom}"
  
  dsub \
    --provider google-batch \
    --project "${GOOGLE_PROJECT}" \
    --user-project "${GOOGLE_PROJECT}" \
    --image "marketplace.gcr.io/google/ubuntu1804:latest" \
    --network "${AOU_NETWORK}" \
    --subnetwork "${AOU_SUBNETWORK}" \
    --service-account "$(gcloud config get-value account)" \
    --use-private-address \
    --user "${DSUB_USER_NAME}" \
    --regions us-central1 \
    --logging "${WORKSPACE_BUCKET}/CMP_Gene_Data/Logs/$(date +'%Y%m%d/%H%M%S')/{job-id}-{task-id}-{task-attempt}.log" \
    --boot-disk-size 3000 \
    --disk-size 1200 \
    --machine-type "${MACHINE_TYPE}" \
    --name "${JOB_NAME}" \
    --env chrom="${chrom}" \
    --env gene="${gene}" \
    --env region="${region}" \
    \
    --input plink="${PLINK_TOOL}" \
    \
    --input input_bed1="${BASE_INPUT_BED_ACAF}${chrom}.bed" \
    --input input_bim1="${BASE_INPUT_BED_ACAF}${chrom}.bim" \
    --input input_fam1="${BASE_INPUT_BED_ACAF}${chrom}.fam" \
    --output-recursive OUTPUT_PATH1="${OUTPUT_PATH1}" \
    \
    --input input_bed2="${BASE_INPUT_BED_EXOME}${chrom}.bed" \
    --input input_bim2="${BASE_INPUT_BED_EXOME}${chrom}.bim" \
    --input input_fam2="${BASE_INPUT_BED_EXOME}${chrom}.fam" \
    --output-recursive OUTPUT_PATH2="${OUTPUT_PATH2}" \
    \
    --input input_bed3="${BASE_INPUT_BED_CLINVAR}${chrom}.bed" \
    --input input_bim3="${BASE_INPUT_BED_CLINVAR}${chrom}.bim" \
    --input input_fam3="${BASE_INPUT_BED_CLINVAR}${chrom}.fam" \
    --output-recursive OUTPUT_PATH3="${OUTPUT_PATH3}" \
    \
    --script "${PIPELINE_SCRIPT}"
done < chrom_gene_region.txt