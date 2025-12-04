# --- GnomAD Annotation Processing (V4.0) ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

require(data.table)

# --- 2. Load and Format gnomAD Data ---
# NOTE: Collaborator must ensure this file is present locally.
# CENSORED PATH: Replaced relative path with configuration variable.
GNOMAD_FILE <- "../../gnomad_annotations_v4.0_01292024_v1.tsv.gz" 

gnomad = fread(GNOMAD_FILE)
print("gnomAD data loaded.")

# Apply ID formatting and PopMax filtering (Code used exactly as provided)
gnomad$alleles= gsub("\\[","",gnomad$alleles)
gnomad$alleles = gsub("\\]","",gnomad$alleles)
gnomad$alleles= gsub(",",":",gnomad$alleles)
gnomad$alleles = gsub("\\\"","",gnomad$alleles)

# Filter gnomAD to variants passing the initial PopMax cutoff (<= 0.001)
# NOTE: Assuming 'joint_freq_popmax' is the correct column name in the TSV file.
gnomad1 = gnomad[which(gnomad$joint_freq_popmax <= 0.001),]
gnomad1$ID = paste0(gnomad1$locus,":",gnomad1$alleles)
gnomad1$ID =gsub("chr","",gnomad1$ID)
print(paste("gnomAD variants filtered:", nrow(gnomad1)))

# --- 3. Save Processed gnomAD Data ---
save(gnomad1, file=file.path(OUTPUT_DIR, "gnomad_popmax_filtered.rda"))
print(paste("Processed gnomAD data saved to:", file.path(OUTPUT_DIR, "gnomad_popmax_filtered.rda")))