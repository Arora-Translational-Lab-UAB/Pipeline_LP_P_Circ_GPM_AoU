# --- Annotation Aggregation and Cleaning ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")
ANNOT_BASE_DIR <- file.path(BASE_DIR, "annot_snpeff")

# Directories for input files (from Script 4/5 output)
LOF_DIR <- file.path(ANNOT_BASE_DIR, "lof")
MISSENSE_DIR <- file.path(ANNOT_BASE_DIR, "missense")
CLINVAR_DIR <- file.path(ANNOT_BASE_DIR, "clinvar")

require(data.table)
library(splitstackshape)
library(dplyr) # Required for later filtering/piping

# Full list of 80 genes (Code used exactly as provided)
genes = c('OBSCN','PLEKHM2','PRDM16','PSEN2','TNNI3K','LMNA','TNNT2','ACTN2','RIT1','NEXN','RYR2','TRIM63','ANKRD1','MYPN','NEBL','COX15','LDB3','RBM20','VCL','BAG3','CTNNA3','ILK','KCNQ1','MYBPC3','CRYAB','CSRP3','MYL2','ABCC9','CACNA1C','PTPN11','PKP2','MYH6','MYH7','TGFB3','TPM1','ALPK3','ACTC1','TJP1','CTF1','TCAP','GAA','JUP','DTNA','TTR','DSC2','DSG2','CDH2','MYOM1','TNNI3','CALR3','TTN','DES','JPH2','MYLK2','SCN5A','MYL3','CAV3','RAF1','TMEM43','TNNC1','SLC25A4','MYOZ2','PDLIM3','NKX2-5','SGCD','EYA4','LAMA4','DSP','MYO6','PLN','KCNH2','GATAD1','TBX20','FLNC','PRKAG2','KLF10','FXN','GLA','FHL1','LAMP2')
print("Gene list defined.")


# --- 2. List Files for Aggregation ---
# (Code used exactly as provided, adjusting directory structure to use variables)
lof = list.files(LOF_DIR, full.names = FALSE)
lof_wgs = lof[grep("wgs",lof)]
lof_exome = lof[grep("exome",lof)]

missense = list.files(MISSENSE_DIR, full.names = FALSE)
missense_wgs = missense[grep("wgs",missense)]
missense_exome = missense[grep("exome",missense)]

clinvar = list.files(CLINVAR_DIR, full.names = FALSE)
clinvar_wgs = clinvar[grep("wgs",clinvar)]
clinvar_exome = clinvar[grep("exome",clinvar)]
print("File lists compiled.")


# --- 3. Load and Combine Categorized TSVs (Initialization and Loops) ---

# Initialize empty dataframes (Code used exactly as provided)
df_lof_wgs = data.frame()
df_lof_exome = data.frame()

df_missense_wgs = data.frame()
df_missense_exome = data.frame()

df_clinvar_wgs = data.frame()
df_clinvar_exome = data.frame()

# Read LoF WGS and Exome separately (Code used exactly as provided)
for (f in lof_wgs) {
  tmp = fread(file.path(LOF_DIR, f));
  df_lof_wgs = rbind(df_lof_wgs, tmp)
}
for (f in lof_exome) {
  tmp = fread(file.path(LOF_DIR, f));
  df_lof_exome = rbind(df_lof_exome, tmp)
}

# Read missense WGS and Exome separately (Code used exactly as provided)
for (f in missense_wgs) {
  tmp = fread(file.path(MISSENSE_DIR, f));
  df_missense_wgs = rbind(df_missense_wgs, tmp)
}
for (f in missense_exome) {
  tmp = fread(file.path(MISSENSE_DIR, f));
  df_missense_exome = rbind(df_missense_exome, tmp)
}

# Read clinvar WGS and Exome separately (Code used exactly as provided)
for (f in clinvar_wgs) {
  tmp = fread(file.path(CLINVAR_DIR, f));
  df_clinvar_wgs = rbind(df_clinvar_wgs, tmp)
}
for (f in clinvar_exome) {
  tmp = fread(file.path(CLINVAR_DIR, f));
  df_clinvar_exome = rbind(df_clinvar_exome, tmp)
}
print("Categorized dataframes aggregated.")


# --- 4. Filter by Master Gene List ---

# Filter Exome dataframes (Code used exactly as provided)
df_lof_exomes1 = df_lof_exome[which(df_lof_exome$GENE %in% genes),]
df_missense_exomes1 = df_missense_exome[which(df_missense_exome$GENE %in% genes),]
df_clinvar_exomes1 = df_clinvar_exome[which(df_clinvar_exome$GENE %in% genes),]

# Filter WGS dataframes (Code used exactly as provided)
df_lof_wgs1 =df_lof_wgs[which(df_lof_wgs$GENE %in% genes),]
df_missense_wgs1 =  df_missense_wgs[which(df_missense_wgs$GENE %in% genes),]
df_clinvar_wgs1 = df_clinvar_wgs[which(df_clinvar_wgs$GENE %in% genes),]
print("Dataframes filtered by 80 target genes.")


# --- 5. Create Unique Variant ID (ID) ---
# (Code used exactly as provided)
df_lof_exomes1$ID = paste0(df_lof_exomes1$CHROM,":",df_lof_exomes1$POS,":",df_lof_exomes1$REF,":",df_lof_exomes1$ALT)
df_lof_wgs1$ID = paste0(df_lof_wgs1$CHROM,":",df_lof_wgs1$POS,":",df_lof_wgs1$REF,":",df_lof_wgs1$ALT)
df_missense_exomes1$ID = paste0(df_missense_exomes1$CHROM,":",df_missense_exomes1$POS,":",df_missense_exomes1$REF,":",df_missense_exomes1$ALT)
df_missense_wgs1$ID = paste0(df_missense_wgs1$CHROM,":",df_missense_wgs1$POS,":",df_missense_wgs1$REF,":",df_missense_wgs1$ALT)
df_clinvar_exomes1$ID = paste0(df_clinvar_exomes1$CHROM,":",df_clinvar_exomes1$POS,":",df_clinvar_exomes1$REF,":",df_clinvar_exomes1$ALT)
df_clinvar_wgs1$ID = paste0(df_clinvar_wgs1$CHROM,":",df_clinvar_wgs1$POS,":",df_clinvar_wgs1$REF,":",df_clinvar_wgs1$ALT)
print("Unique variant IDs created.")


# --- 6. Initial Data Cleaning (GENEINFO splitting) ---

# Apply splitstackshape transformation to ClinVar Exomes (Code used exactly as provided)
# Note: This is an important data structure transformation needed for downstream filtering.
df_clinvar_exomes1 = cSplit(df_clinvar_exomes1,"GENEINFO",":")
df_clinvar_exomes1 = cSplit(df_clinvar_exomes1,"GENEINFO_2","|")
print("GENEINFO field split in ClinVar Exomes data.")


# --- 7. Save Intermediate Filtered Data ---
# This saves the cleaned, categorized dataframes for the next filtering step (Step 9).
save(df_lof_exomes1, df_missense_exomes1, df_clinvar_exomes1, 
     df_lof_wgs1, df_missense_wgs1, df_clinvar_wgs1, 
     file=file.path(OUTPUT_DIR, "variants_categorized_and_cleaned_12032025.rda"))

print("Annotation Aggregation and Cleaning (Script 8) complete.")