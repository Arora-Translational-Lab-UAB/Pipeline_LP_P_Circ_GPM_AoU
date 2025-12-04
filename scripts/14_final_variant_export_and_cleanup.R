# --- Final Variant List Export and Cleanup ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")
FINAL_EXPORT_DIR <- file.path(BASE_DIR, "final_variant_lists") # Local output directory for CSVs

# GCS path for final output (CENSORED)
GCS_FINAL_EXPORT_BUCKET="gs://<YOUR_SECURE_BUCKET>/CMP_Gene_Data/Annotations/final_variants_list/"
PROJECT_ID="<YOUR_PROJECT_ID>" # For gsutil command

dir.create(FINAL_EXPORT_DIR, showWarnings = FALSE)

require(data.table)
library(dplyr)

# Load final filtered variant dataframes (from Script 10)
load(file.path(OUTPUT_DIR, "final_filtered_variants_by_disease_v8_12032025.rda"))

# --- 2. Define Gene Lists (Needed for final ClinVar categorization, defined in Script 10) ---
genes_for_dcm1 = c('TTN','MYH7','LMNA','FLNC','DSP','DES','TNNT2','SCN5A','ACTC1','PLN','TNNC1','BAG3','TNNI3','TPM1')
genes_for_hcm1 = c('MYH7','MYBPC3','TNNT2','TPM1','TNNI3','MYL2','MYL3','ACTC1','TNNC1') 
genes_for_arvc1 = c('DSC2','DSG2','DSP','PKP2','TMEM43','JUP','DES','PLN')

# --- 3. Re-Filter ClinVar P/LP by Gene (Matching the logic in the provided code) ---
# NOTE: This step ensures the ClinVar variants (which passed 2-star filter) are subsetted 
# based on the target disease gene lists (genes_for_dcm1, etc.)
# It uses the pre-PopMax filtered variants (clinvar_final_exomes/wgs).

# Helper function for ClinVar subsetting:
subset_clinvar <- function(df, disease_genes) {
    # Assuming GENEINFO_1 exists and contains the gene symbol (from Script 10 cleaning)
    if ("GENEINFO_1" %in% colnames(df)) {
        return(df[which(df$GENEINFO_1 %in% disease_genes),])
    } else {
        # Fallback to GENE column if GENEINFO_1 is not available/clean
        return(df[which(df$GENE %in% disease_genes),])
    }
}

# Apply final gene filtering for ClinVar P/LP lists:
clinvar_dcm_exomes = subset_clinvar(clinvar_final_exomes, genes_for_dcm1)
clinvar_hcm_exomes = subset_clinvar(clinvar_final_exomes, genes_for_hcm1)
clinvar_arvc_exomes = subset_clinvar(clinvar_final_exomes, genes_for_arvc1)

clinvar_dcm_wgs = subset_clinvar(clinvar_final_wgs, genes_for_dcm1)
clinvar_hcm_wgs = subset_clinvar(clinvar_final_wgs, genes_for_hcm1)
clinvar_arvc_wgs = subset_clinvar(clinvar_final_wgs, genes_for_arvc1)

print("Final ClinVar P/LP lists subsetted by disease gene.")

# --- 4. Export Final Variant Lists (CSVs) ---
# Export files use the naming convention: [IMPACT]_[DISEASE]_[AF_FILTER]_[SOURCE]_[VERSION].csv

export_csv <- function(df, filename) {
    write.table(unique(df$ID), 
                file=file.path(FINAL_EXPORT_DIR, filename),
                row.names=F, col.names=T, sep="\t", dec=".", quote=F)
    print(paste("Exported:", filename))
}

# 4.1 ClinVar P/LP Variants (using date suffix 11042025 and 11132025)
# Using the more conservative (final) 11132025 date for the main output files:
export_csv(clinvar_dcm_exomes, "lp_p_dcm_aou_gnomad_popmax_exomes_v8_11132025.csv")
export_csv(clinvar_hcm_exomes, "lp_p_hcm_aou_gnomad_popmax_exomes_v8_11132025.csv")
export_csv(clinvar_arvc_exomes, "lp_p_arvc_aou_gnomad_popmax_exomes_v8_11132025.csv")

# Exporting the full set (Exomes & WGS) with 11042025 suffix:
export_csv(clinvar_dcm_exomes, "lp_p_dcm_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(clinvar_hcm_exomes, "lp_p_hcm_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(clinvar_arvc_exomes, "lp_p_arvc_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(clinvar_dcm_wgs, "lp_p_dcm_aou_gnomad_popmax_wgs_v8_11042025.csv")
export_csv(clinvar_hcm_wgs, "lp_p_hcm_aou_gnomad_popmax_wgs_v8_11042025.csv")
export_csv(clinvar_arvc_wgs, "lp_p_arvc_aou_gnomad_popmax_wgs_v8_11042025.csv")

# 4.2 LoF Variants (Filtered by Known Mechanism/pLI and FaF)
export_csv(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm, "lof_dcm_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm, "lof_hcm_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc, "lof_arvc_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm_wgs, "lof_dcm_aou_gnomad_popmax_wgs_v8_11042025.csv")
export_csv(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm_wgs, "lof_hcm_aou_gnomad_popmax_wgs_v8_11042025.csv")
export_csv(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc_wgs, "lof_arvc_aou_gnomad_popmax_wgs_v8_11042025.csv")

# 4.3 Missense Variants (Filtered by Known Mechanism/Z-score and FaF)
export_csv(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm, "missense_dcm_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm, "missense_hcm_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc, "missense_arvc_aou_gnomad_popmax_exomes_v8_11042025.csv")
export_csv(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm_wgs, "missense_dcm_aou_gnomad_popmax_wgs_v8_11042025.csv")
export_csv(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm_wgs, "missense_hcm_aou_gnomad_popmax_wgs_v8_11042025.csv")
export_csv(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc_wgs, "missense_arvc_aou_gnomad_popmax_wgs_v8_11042025.csv")

# --- 5. GCS Upload and R Cleanup ---
# Upload all generated CSV files to the final GCS directory
system(paste0("gsutil cp -r ", FINAL_EXPORT_DIR, "/*.csv ", GCS_FINAL_EXPORT_BUCKET))

# Clear specific variables (Code used exactly as provided)
rm(list = ls()[grep("faf",ls())])

# Save the final R workspace image (Code used exactly as provided)
save.image("backup_variants_lp_p_v8_11042025.rda")

# Upload the final R workspace image to GCS (Code used exactly as provided)
system(paste0("gsutil cp -r backup_variants_lp_p_v8_11042025.rda ", GCS_FINAL_EXPORT_BUCKET))

print("Final variant lists exported and R workspace saved.")