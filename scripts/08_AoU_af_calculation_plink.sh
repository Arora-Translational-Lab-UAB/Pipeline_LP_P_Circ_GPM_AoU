#!/bin/bash

# --- Ancestry Filtering and AF Calculation using PLINK2 ---
# This script downloads Plink files, subsets samples by ancestry, and calculates AFs.

set -o pipefail
set -o errexit
set -u

echo "--- Starting PopMax AF Calculation Setup ---"

# --- 1. Configuration (CENSORED PLACEHOLDERS) ---
# NOTE: Collaborator MUST replace these placeholders with actual paths/IDs.
GOOGLE_PROJECT="<YOUR_GOOGLE_PROJECT_ID>"
GCS_INPUT_BUCKET="gs://<YOUR_SECURE_BUCKET>/CMP_Gene_Data/Plink_Format"
LOCAL_PLINK_DIR="geno"
PLINK_EXE="plink2" # Assumes plink2 is in PATH, otherwise use full path.

LOCAL_EXOME_SUBSET="${LOCAL_PLINK_DIR}/exome_subset"
LOCAL_CLINVAR_SUBSET="${LOCAL_PLINK_DIR}/clinvar_subset"

LOCAL_FREQ_EXOME_DIR="${LOCAL_PLINK_DIR}/Freq_Exomes"
LOCAL_FREQ_CLINVAR_DIR="${LOCAL_PLINK_DIR}/Freq_Clinvar"

# --- 2. Environment Setup and Data Download ---
mkdir -p "$LOCAL_PLINK_DIR"
mkdir -p "$LOCAL_EXOME_SUBSET" "$LOCAL_CLINVAR_SUBSET"
mkdir -p "$LOCAL_FREQ_EXOME_DIR" "$LOCAL_FREQ_CLINVAR_DIR"

echo "Downloading Plink subsets from GCS..."
# Download Exome and Clinvar subsets generated in Script 2
gsutil -m cp -r "${GCS_INPUT_BUCKET}/EXOME_Subset/*" "$LOCAL_PLINK_DIR"
gsutil -m cp -r "${GCS_INPUT_BUCKET}/Clinvar_Subset/*" "$LOCAL_PLINK_DIR"

# Move files into segregated directories (matching your notebook logic)
echo "Organizing downloaded files..."
mv "$LOCAL_PLINK_DIR"/*_2.* "$LOCAL_EXOME_SUBSET"
mv "$LOCAL_PLINK_DIR"/*_3.* "$LOCAL_CLINVAR_SUBSET"

# --- 3. Ancestry Sample File Generation (R) ---
echo "Running R to prepare ancestry sample lists..."
# NOTE: This R script relies on the 'ancestry_preds.tsv' file being in the execution directory.
R_SCRIPT_PATH="${BASE_DIR}/scripts/06b_ancestry_prep.R"

# We must ensure the R script to generate the .samples files is executed.
# This part is complex and best handled by calling R directly.
# This R script will create european_all.samples, african_all.samples, etc.

# Create a temporary R script for this task:
cat > "${BASE_DIR}/ancestry_prep_temp.R" <<- 'END'
require(data.table)
library(dplyr)
# NOTE: Assumes 'ancestry_preds.tsv' is in the working directory
data = fread("ancestry_preds.tsv")
data$FID = 0
data$IID = data$research_id

data_eur = data %>% filter(ancestry_pred %in% c("eur")) %>% select(FID,IID)
data_afr = data %>% filter(ancestry_pred %in% c("afr")) %>% select(FID,IID)
data_amr = data %>% filter(ancestry_pred %in% c("amr")) %>% select(FID,IID)
data_eas = data %>% filter(ancestry_pred %in% c("eas")) %>% select(FID,IID)
data_sas = data %>% filter(ancestry_pred %in% c("sas")) %>% select(FID,IID)
data_mid = data %>% filter(ancestry_pred %in% c("mid")) %>% select(FID,IID)

write.table(data_eur,file="european_all.samples",row.names=F,col.names=F,sep="\t",dec=".",quote=F)
write.table(data_afr,file="african_all.samples",row.names=F,col.names=F,sep="\t",dec=".",quote=F)
write.table(data_amr,file="amr_all.samples",row.names=F,col.names=F,sep="\t",dec=".",quote=F)
write.table(data_eas,file="eas_all.samples",row.names=F,col.names=F,sep="\t",dec=".",quote=F)
write.table(data_sas,file="sas_all.samples",row.names=F,col.names=F,sep="\t",dec=".",quote=F)
write.table(data_mid,file="mid_all.samples",row.names=F,col.names=F,sep="\t",dec=".",quote=F)
END
Rscript "${BASE_DIR}/ancestry_prep_temp.R"
rm -f "${BASE_DIR}/ancestry_prep_temp.R"


# --- 4. PLINK AF Calculation - Exome Subset ---
echo "Calculating AFs for Exome subsets..."
list_exome=$(ls -1 ${LOCAL_EXOME_SUBSET}/*.bed | sed 's/\.bed$//' | xargs -n 1 basename)

for bfile_prefix in $list_exome; do
    echo "Processing Exome: $bfile_prefix"
    # Execute PLINK AF calculation for all six ancestries (matching original logic)
    "${PLINK_EXE}" --bfile "${LOCAL_EXOME_SUBSET}/${bfile_prefix}" --keep european_all.samples --freq --out "${LOCAL_FREQ_EXOME_DIR}/eur_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_EXOME_SUBSET}/${bfile_prefix}" --keep african_all.samples --freq --out "${LOCAL_FREQ_EXOME_DIR}/afr_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_EXOME_SUBSET}/${bfile_prefix}" --keep amr_all.samples --freq --out "${LOCAL_FREQ_EXOME_DIR}/amr_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_EXOME_SUBSET}/${bfile_prefix}" --keep eas_all.samples --freq --out "${LOCAL_FREQ_EXOME_DIR}/eas_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_EXOME_SUBSET}/${bfile_prefix}" --keep sas_all.samples --freq --out "${LOCAL_FREQ_EXOME_DIR}/sas_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_EXOME_SUBSET}/${bfile_prefix}" --keep mid_all.samples --freq --out "${LOCAL_FREQ_EXOME_DIR}/mid_${bfile_prefix}"
done


# --- 5. PLINK AF Calculation - ClinVar Subset ---
echo "Calculating AFs for ClinVar subsets..."
list_clinvar=$(ls -1 ${LOCAL_CLINVAR_SUBSET}/*.bed | sed 's/\.bed$//' | xargs -n 1 basename)

for bfile_prefix in $list_clinvar; do
    echo "Processing ClinVar: $bfile_prefix"
    # Execute PLINK AF calculation for all six ancestries (matching original logic)
    "${PLINK_EXE}" --bfile "${LOCAL_CLINVAR_SUBSET}/${bfile_prefix}" --keep european_all.samples --freq --out "${LOCAL_FREQ_CLINVAR_DIR}/eur_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_CLINVAR_SUBSET}/${bfile_prefix}" --keep african_all.samples --freq --out "${LOCAL_FREQ_CLINVAR_DIR}/afr_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_CLINVAR_SUBSET}/${bfile_prefix}" --keep amr_all.samples --freq --out "${LOCAL_FREQ_CLINVAR_DIR}/amr_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_CLINVAR_SUBSET}/${bfile_prefix}" --keep eas_all.samples --freq --out "${LOCAL_FREQ_CLINVAR_DIR}/eas_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_CLINVAR_SUBSET}/${bfile_prefix}" --keep sas_all.samples --freq --out "${LOCAL_FREQ_CLINVAR_DIR}/sas_${bfile_prefix}"
    "${PLINK_EXE}" --bfile "${LOCAL_CLINVAR_SUBSET}/${bfile_prefix}" --keep mid_all.samples --freq --out "${LOCAL_FREQ_CLINVAR_DIR}/mid_${bfile_prefix}"
done

echo "PLINK AF calculation complete. AFREQ files are ready for R aggregation."