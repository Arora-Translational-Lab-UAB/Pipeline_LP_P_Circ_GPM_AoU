# --- Final Variant Counting and Summary Tables ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")
OUTPUT_TABLES_DIR <- file.path(BASE_DIR, "final_summary_tables")

dir.create(OUTPUT_TABLES_DIR, showWarnings = FALSE)

require(data.table)
library(dplyr)

# Load input data from previous steps
# NOTE: df_comb2 (AoU PopMax filtered) is assumed to be loaded from Script 7/10
# NOTE: v (LoF/Missense), v1 (ClinVar P/LP) are combined lists from Script 10
load(file.path(OUTPUT_DIR, "aou_gnomad_popmax_filters_v8_12032025.rda")) 
load(file.path(OUTPUT_DIR, "final_filtered_variants_by_disease_v8_12032025.rda")) 

print("Loaded AoU PopMax data (df_comb2) and final filtered lists (v, v1).")

# --- 2. Define Master Gene List for Counting ---
# Note: These lists are defined exactly as provided in Script 10.
genes_for_dcm = c("BAG3","DES","FLNC","LMNA","MYH7","PLN","RBM20","SCN5A","TNNC1","TNNT2","TTN","DSP","ACTC1","ACTN2")
genes_for_hcm = c("MYBPC3","MYH7","TNNT2","TNNI3","TPM1","ACTC1","MYL3","MYL2") 
genes_for_arvc = c("PKP2","DSP","DSG2","DSC2","JUP","TMEM43")
genes_master = c(genes_for_arvc, genes_for_dcm, genes_for_hcm)


# --- 3. Clean Gene Assignments and Filter Combined Variants ---

# Extract and clean gene names from AF data (Code used exactly as provided)
df_comb2$gene_eur1 <- sub("eur_([^_]+)_.*", "\\1", df_comb2$gene_eur)
df_comb2$gene_afr1 <- sub("afr_([^_]+)_.*", "\\1", df_comb2$gene_afr)
df_comb2$gene_amr1 <- sub("amr_([^_]+)_.*", "\\1", df_comb2$gene_amr)
df_comb2$gene_eas1 <- sub("eas_([^_]+)_.*", "\\1", df_comb2$gene_eas)
df_comb2$gene_sas1 <- sub("sas_([^_]+)_.*", "\\1", df_comb2$gene_sas)
df_comb2$gene_mid1 <- sub("mid_([^_]+)_.*", "\\1", df_comb2$gene_mid)

# Filter AoU PopMax data (df_comb2) to only include variants present in the final filtered lists (v or v1)
df_comb3 = df_comb2[which(df_comb2$ID %in% v$ID | df_comb2$ID %in% v1$ID),]
print(paste("Initial combined variant list (df_comb3) dimension:", paste(dim(df_comb3), collapse = " x ")))

# Deduplicate variants (Code used exactly as provided)
df_comb_aou_gnomad_popmax1 = df_comb3[!duplicated(df_comb3$ID),]

# Assign a single final gene name to each variant (Code used exactly as provided)
df_comb_aou_gnomad_popmax1$genes_final <- apply(df_comb_aou_gnomad_popmax1[, c("gene_eur1", "gene_amr1", "gene_afr1", "gene_eas1","gene_sas1","gene_mid1")], 1, function(x) {
  return(na.omit(x)[1])  # Picks the first non-NA value
})


# --- 4. Flag Variant Categories (LP/P vs. LoF/Missense) ---

# Create binary flags for categories based on inclusion in filtered lists (Code used exactly as provided)
df_comb_aou_gnomad_popmax1$lp_p = ifelse(df_comb_aou_gnomad_popmax1$ID %in% v1$ID,1,0)
df_comb_aou_gnomad_popmax1$lof_missense = ifelse(df_comb_aou_gnomad_popmax1$ID %in% v$ID,1,0) # Note: 'v' is LoF/Missense combined


# --- 5. Final Filtering to Master Gene List and Cross-Tabulation ---

# Filter the final variant set to only include the genes targeted by the disease criteria (Code used exactly as provided)
df_comb_aou_gnomad_popmax2 = df_comb_aou_gnomad_popmax1 %>% filter(genes_final %in% genes_master)

# Generate and save cross-tabulation table (Gene vs. LP/P status)
print("--- Final Variant Counts by Gene and P/LP Status (df_comb_aou_gnomad_popmax2) ---")
variant_counts_lp_p <- ftable(df_comb_aou_gnomad_popmax2$genes_final, df_comb_aou_gnomad_popmax2$lp_p)
print(variant_counts_lp_p)

# Save the final table and dataframe
write.table(variant_counts_lp_p, file.path(OUTPUT_TABLES_DIR, "final_variant_counts_lp_p_ftable.tsv"), sep="\t", quote=F, row.names=T)
save(df_comb_aou_gnomad_popmax2, file=file.path(OUTPUT_TABLES_DIR, "final_counted_variants_df.rda"))

print("Variant counting and summary table generation (Script 12) complete.")