# --- Final Tiered Filtering: AF, Functional Impact, and Disease Subtype ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")

require(data.table)
library(dplyr)
library(splitstackshape)
library(rlang)

# Load previously filtered data
load(file.path(OUTPUT_DIR, "verified_initial_filtered_v8_12032025.rda"))
load(file.path(OUTPUT_DIR, "aou_gnomad_popmax_filters_v8_12032025.rda")) 

print("Loaded functionally filtered variants and PopMax data.")

# --- 2. Define Gene Lists and AF Cutoffs ---
# Note: These lists are defined exactly as provided.

## P/LP Category Genes
genes_for_dcm1 = c('TTN','MYH7','LMNA','FLNC','DSP','DES','TNNT2','SCN5A','ACTC1','PLN','TNNC1','BAG3','TNNI3','TPM1')
genes_for_hcm1 = c('MYH7','MYBPC3','TNNT2','TPM1','TNNI3','MYL2','MYL3','ACTC1','TNNC1') 
genes_for_arvc1 = c('DSC2','DSG2','DSP','PKP2','TMEM43','JUP','DES','PLN')

## LoF/Missense Genes for final PopMax/FaF analysis
genes_for_dcm = c("BAG3","DES","FLNC","LMNA","MYH7","PLN","RBM20","SCN5A","TNNC1","TNNT2","TTN","DSP","ACTC1","ACTN2")
genes_for_hcm = c("MYBPC3","MYH7","TNNT2","TNNI3","TPM1","ACTC1","MYL3","MYL2") 
genes_for_arvc = c("PKP2","DSP","DSG2","DSC2","JUP","TMEM43")

## Known Mechanism/pLI/Z-score lists (Master filtering lists)
known_mechanism_plof_pli_0.90  = c("PKP2","DSP","DSC2","JUP","DSG2","BAG3","FLNC","TTN","DSP","PLN","LMNA","RBM20","MYBPC3","PLN","FLNC","DSP","ACTN2","ACTN2")
known_mechanism_misse_zscore_1.64  = c("MYBPC3","MYH7","TNNT2","TNNI3","TPM1","ACTC1","MYL3","MYL2","PLN","DES","ACTC1","DSG2","DES","MYH7","SCN5A","TNNC1","TNNT2","PLN","LMNA","RBM20")

# FaF Cutoffs (Absolute values based on your code)
FAF_DCM_CUTOFF = 8.4e-05
FAF_ARVC_CUTOFF = 9.2e-05
FAF_HCM_CUTOFF = 4e-05


# --- 3. WGS Data Pre-processing (Handling Multi-allelics) ---
# NOTE: This complex step corrects the WGS ALT allele field for multi-allelic variants
# by expanding the ALT column and re-calculating the ID based on the ALLELE field.

# Step 3.1: ClinVar WGS (Requires splitting ALT field)
clinvar_final_wgs = cSplit(clinvar_2star_wgs,"ALT",",")

# Loop through ALT columns (ALT_1, ALT_2, etc.) and pick the one that matches ALLELE
alt_cols_cln <- sprintf("ALT_%02d", 1:ncol(clinvar_final_wgs %>% select(starts_with("ALT_")))) 
cases_list_cln <- lapply(alt_cols_cln, function(col) expr(.data$ALLELE == .data[[!!col]] ~ .data[[!!col]]))

clinvar_final_wgs <- clinvar_final_wgs %>%
  mutate(
    ALT_corrected = case_when(!!!cases_list_cln, TRUE ~ NA_character_)
  )
# Recreate ID using the corrected ALT allele for WGS
clinvar_final_wgs$ID = paste0(clinvar_final_wgs$CHROM,":",clinvar_final_wgs$POS,":",clinvar_final_wgs$REF_1,":",clinvar_final_wgs$ALLELE)


# Step 3.2: LoF WGS (Requires splitting ALT field)
lof_final_wgs = cSplit(df_lof_wgs_ttn_hipsi2, "ALT", ",")
alt_cols_lof <- sprintf("ALT_%02d", 1:ncol(lof_final_wgs %>% select(starts_with("ALT_")))) 
cases_list_lof <- lapply(alt_cols_lof, function(col) expr(.data$ALLELE == .data[[!!col]] ~ .data[[!!col]]))

lof_final_wgs <- lof_final_wgs %>%
  mutate(
    ALT_corrected = case_when(!!!cases_list_lof, TRUE ~ NA_character_)
  )
lof_final_wgs$ID = paste0(lof_final_wgs$CHROM,":",lof_final_wgs$POS,":",lof_final_wgs$REF_1,":",lof_final_wgs$ALT_corrected)

# Step 3.3: Missense WGS (Requires splitting ALT field)
missense_final_wgs = cSplit(df_missense_wgs4, "ALT", ",")
# Note: The subsequent missense split was duplicated in the original notebook, using only the first split here.


# --- 4. Apply PopMax (AF <= 0.001) Filter ---
# Note: af_popmax_aou_gnomad_filter is the list of IDs passing the AF <= 0.001 cutoff (AoU OR gnomAD).

# 4.1 ClinVar PopMax Filter
clinvar_final_exomes = cSplit(clinvar_2star_exomes,"GENEINFO",":") # Re-split GENEINFO before PopMax filter
clinvar_final_exomes = cSplit(clinvar_final_exomes,"GENEINFO_2","|") 

clinvar_popmax_aou_gnomad_exomes = clinvar_final_exomes[which(clinvar_final_exomes$ID %in% af_popmax_aou_gnomad_filter$ID),]
clinvar_popmax_aou_gnomad_wgs = clinvar_final_wgs[which(clinvar_final_wgs$ID %in% af_popmax_aou_gnomad_filter$ID),]

# 4.2 LoF/Missense PopMax Filter
lof_popmax_aou_gnomad_exomes = df_lof_exomes_ttn_hipsi2[which(df_lof_exomes_ttn_hipsi2$ID %in% af_popmax_aou_gnomad_filter$ID),]
missense_popmax_aou_gnomad_exomes = df_missense_exomes4[which(df_missense_exomes4$ID %in% af_popmax_aou_gnomad_filter$ID),]

lof_popmax_aou_gnomad_wgs = lof_final_wgs[which(lof_final_wgs$ID %in% af_popmax_aou_gnomad_filter$ID),]
missense_popmax_aou_gnomad_wgs = missense_final_wgs[which(missense_final_wgs$ID %in% af_popmax_aou_gnomad_filter$ID),]


# --- 5. Apply Known Mechanism (pLI/Z-score) Filter ---

## LoF Known Mechanism Filter
lof_popmax_aou_gnomad_wgs_filter = lof_popmax_aou_gnomad_wgs[which(lof_popmax_aou_gnomad_wgs$GENE %in% known_mechanism_plof_pli_0.90),]
lof_popmax_aou_gnomad_exomes_filter = lof_popmax_aou_gnomad_exomes[which(lof_popmax_aou_gnomad_exomes$GENE %in% known_mechanism_plof_pli_0.90),]

## Missense Known Mechanism Filter
missense_popmax_aou_gnomad_wgs_filter = missense_popmax_aou_gnomad_wgs[which(missense_popmax_aou_gnomad_wgs$GENE %in% known_mechanism_misse_zscore_1.64),]
missense_popmax_aou_gnomad_exomes_filter = missense_popmax_aou_gnomad_exomes[which(missense_popmax_aou_gnomad_exomes$GENE %in% known_mechanism_misse_zscore_1.64),]


# --- 6. Apply FaF (Final Filtered Allele Frequency) Cutoffs ---

# Filter AoU PopMax data (df_comb2) and gnomAD data (gnomad1) by disease-specific FaF cutoffs
dcm_aou_faf = df_comb2[which(df_comb2$max_value  <= FAF_DCM_CUTOFF),]
arvc_aou_faf = df_comb2[which(df_comb2$max_value <= FAF_ARVC_CUTOFF),]
hcm_aou_faf = df_comb2[which(df_comb2$max_value  <= FAF_HCM_CUTOFF),]

dcm_gnomad_faf = gnomad1[which(gnomad1$joint_freq_popmax  <= FAF_DCM_CUTOFF),]
arvc_gnomad_faf = gnomad1[which(gnomad1$joint_freq_popmax  <= FAF_ARVC_CUTOFF),]
hcm_gnomad_faf = gnomad1[which(gnomad1$joint_freq_popmax   <= FAF_HCM_CUTOFF),]


# --- 7. Final Tiered Filtering by Disease Subtype (DCM, HCM, ARVC) ---

## 7.1 DCM Filtering
# LoF DCM
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm = lof_popmax_aou_gnomad_exomes_filter[which(lof_popmax_aou_gnomad_exomes_filter$ID %in% dcm_aou_faf$ID | lof_popmax_aou_gnomad_exomes_filter$ID %in% dcm_gnomad_faf$ID),] 
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm = known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm[which(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm$GENE %in% genes_for_dcm),]
# Missense DCM
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm = missense_popmax_aou_gnomad_exomes_filter[which(missense_popmax_aou_gnomad_exomes_filter$ID %in% dcm_aou_faf$ID | missense_popmax_aou_gnomad_exomes_filter$ID %in% dcm_gnomad_faf$ID),] 
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm  = known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm[which(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm$GENE %in% genes_for_dcm),]
# WGS DCM
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm_wgs = lof_popmax_aou_gnomad_wgs_filter[which(lof_popmax_aou_gnomad_wgs_filter$ID %in% dcm_aou_faf$ID | lof_popmax_aou_gnomad_wgs_filter$ID %in% dcm_gnomad_faf$ID),] 
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm_wgs = known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm_wgs[which(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm_wgs$GENE %in% genes_for_dcm),]
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm_wgs = missense_popmax_aou_gnomad_wgs_filter[which(missense_popmax_aou_gnomad_wgs_filter$ID %in% dcm_aou_faf$ID | missense_popmax_aou_gnomad_wgs_filter$ID %in% dcm_gnomad_faf$ID),] 
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm_wgs  = known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm_wgs[which(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm_wgs$GENE %in% genes_for_dcm),]

## 7.2 HCM Filtering
# LoF HCM
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm = lof_popmax_aou_gnomad_exomes_filter[which(lof_popmax_aou_gnomad_exomes_filter$ID %in% hcm_aou_faf$ID | lof_popmax_aou_gnomad_exomes_filter$ID %in% hcm_gnomad_faf$ID),] 
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm = known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm[which(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm$GENE %in% genes_for_hcm),]
# Missense HCM
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm = missense_popmax_aou_gnomad_exomes_filter[which(missense_popmax_aou_gnomad_exomes_filter$ID %in% hcm_aou_faf$ID | missense_popmax_aou_gnomad_exomes_filter$ID %in% hcm_gnomad_faf$ID),] 
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm  = known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm[which(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm$GENE %in% genes_for_hcm),]
# WGS HCM
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm_wgs = lof_popmax_aou_gnomad_wgs_filter[which(lof_popmax_aou_gnomad_wgs_filter$ID %in% hcm_aou_faf$ID | lof_popmax_aou_gnomad_wgs_filter$ID %in% hcm_gnomad_faf$ID),] 
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm_wgs = known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm_wgs[which(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm_wgs$GENE %in% genes_for_hcm),]
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm_wgs = missense_popmax_aou_gnomad_wgs_filter[which(missense_popmax_aou_gnomad_wgs_filter$ID %in% hcm_aou_faf$ID | missense_popmax_aou_gnomad_wgs_filter$ID %in% hcm_gnomad_faf$ID),] 
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm_wgs  = known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm_wgs[which(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm_wgs$GENE %in% genes_for_hcm),]

## 7.3 ARVC Filtering
# LoF ARVC
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc = lof_popmax_aou_gnomad_exomes_filter[which(lof_popmax_aou_gnomad_exomes_filter$ID %in% arvc_aou_faf$ID | lof_popmax_aou_gnomad_exomes_filter$ID %in% arvc_gnomad_faf$ID),] 
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc = known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc[which(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc$GENE %in% genes_for_arvc),]
# Missense ARVC
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc = missense_popmax_aou_gnomad_exomes_filter[which(missense_popmax_aou_gnomad_exomes_filter$ID %in% arvc_aou_faf$ID | missense_popmax_aou_gnomad_exomes_filter$ID %in% arvc_gnomad_faf$ID),] 
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc  = known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc[which(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc$GENE %in% genes_for_arvc),]
# WGS ARVC
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc_wgs = lof_popmax_aou_gnomad_wgs_filter[which(lof_popmax_aou_gnomad_wgs_filter$ID %in% arvc_aou_faf$ID | lof_popmax_aou_gnomad_wgs_filter$ID %in% arvc_gnomad_faf$ID),] 
known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc_wgs = known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc_wgs[which(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc_wgs$GENE %in% genes_for_arvc),]
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc_wgs = missense_popmax_aou_gnomad_wgs_filter[which(missense_popmax_aou_gnomad_wgs_filter$ID %in% arvc_aou_faf$ID | missense_popmax_aou_gnomad_wgs_filter$ID %in% arvc_gnomad_faf$ID),] 
known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc_wgs  = known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc_wgs[which(known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc_wgs$GENE %in% genes_for_arvc),]


# --- 8. Combine Final Filtered Lists (For Diagnostics and Export) ---
# Combined LoF/Missense list for Exomes (v)
v = rbind(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc,
          known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm,
          known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm, known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm)

# Combined ClinVar P/LP list for Exomes (v1)
v1 = rbind(clinvar_dcm_exomes, clinvar_hcm_exomes, clinvar_arvc_exomes)


# --- 9. Save Final Lists for Export ---
save(known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm,
     known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm,
     known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc,
     known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_dcm_wgs, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_dcm_wgs,
     known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_hcm_wgs, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_hcm_wgs,
     known_mechanism_pli_lof_pli_popmax_aou_gnomad_filter_arvc_wgs, known_mechanism_missense_zscore_popmax_aou_gnomad_filter_arvc_wgs,
     v, v1,
     file=file.path(OUTPUT_DIR, "final_filtered_variants_by_disease_v8_12032025.rda"))

print("Final tiered filtering (Script 10) complete. Variants separated by disease and impact.")