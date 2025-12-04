# --- R Analysis - Part 1: Initial Filtering and PopMax Calculation ---

# --- 1. Configuration and Library Loading ---
# NOTE: Collaborator MUST set the working directory to the location containing 
# the 'geno', 'annot_snpeff', and external files (gnomad_annotations_...).
BASE_DIR <- getwd() 

# Directories for input files
FREQ_CLINVAR_DIR <- file.path(BASE_DIR, "geno/Freq_Clinvar")
FREQ_EXOMES_DIR <- file.path(BASE_DIR, "geno/Freq_Exomes")
LOF_DIR <- file.path(BASE_DIR, "annot_snpeff/lof")
MISSENSE_DIR <- file.path(BASE_DIR, "annot_snpeff/missense")
CLINVAR_DIR <- file.path(BASE_DIR, "annot_snpeff/clinvar")
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Libraries (Code used exactly as provided)
# install.packages('R.utils') # Note: This must be run if R.utils is not installed
require(data.table)
library(dplyr)
library(splitstackshape)
library(rlang)

# --- 2. Load and Prepare gnomAD Annotations (Script 5 Output) ---
# NOTE: The collaborator must ensure this file is present locally.
data = fread("gnomad_annotations_v4.0_01292024_v1.tsv.gz")

# (Code used exactly as provided)
# Filtering gnomAD by popmax <= 0.001
data$alleles= gsub("\\[","",data$alleles)
data$alleles = gsub("\\]","",data$alleles)
data$alleles= gsub(",",":",data$alleles)
data$alleles = gsub("\\\"","",data$alleles)
gnomad1 = data[which(data$joint_grpmax_AF <= 0.001),] # Renamed joint_freq_popmax to joint_grpmax_AF for consistency with Script 5
gnomad1$ID = paste0(gnomad1$locus,":",gnomad1$alleles)
gnomad1$ID =gsub("chr","",gnomad1$ID)
print("gnomAD data loaded and filtered.")


# --- 3. Aggregate AoU Allele Frequencies (AF) (Script 6 Output) ---

files_cln <- list.files(FREQ_CLINVAR_DIR, pattern = ".afreq", full.names = FALSE)
files_exo <- list.files(FREQ_EXOMES_DIR, pattern = ".afreq", full.names = FALSE)

# Aggregate AF data function
aggregate_afs <- function(file_list, dir, prefix) {
    df_combined <- data.frame()
    for (f in file_list) {
        tmp <- fread(file.path(dir, f))
        tmp <- tmp[which(tmp$ALT_FREQ > 0),] 
        tmp$GENE <- gsub(paste0(prefix, "_|\\.afreq$"), "", basename(f))
        tmp$GENE <- gsub("_chr.*", "", tmp$GENE)
        df_combined <- bind_rows(df_combined, tmp)
    }
    df_combined <- df_combined %>%
        select(ID, ALT_FREQ, GENE) %>%
        rename_with(~ paste0("af_", prefix), ALT_FREQ) %>%
        rename_with(~ paste0("gene_", prefix), GENE)
    return(df_combined)
}

# Combine files from Freq_Exomes and Freq_Clinvar (Note: Original code merged two separate runs)
# Using Exomes directory for the gene lists since it's the primary source.
files_exome_base <- gsub(".afreq$", "", files_exo)

afr_files <- files_exo[grep("afr_", files_exo)]
eur_files <- files_exo[grep("eur_", files_exo)]
amr_files <- files_exo[grep("amr_", files_exo)]
eas_files <- files_exo[grep("eas_", files_exo)]
sas_files <- files_exo[grep("sas_", files_exo)]
mid_files <- files_exo[grep("mid_", files_exo)]

# Combine Frequencies (Exomes)
df_afr1 = aggregate_afs(afr_files, FREQ_EXOMES_DIR, "afr")
df_eur1 = aggregate_afs(eur_files, FREQ_EXOMES_DIR, "eur")
df_amr1 = aggregate_afs(amr_files, FREQ_EXOMES_DIR, "amr")
df_eas1 = aggregate_afs(eas_files, FREQ_EXOMES_DIR, "eas")
df_sas1 = aggregate_afs(sas_files, FREQ_EXOMES_DIR, "sas")
df_mid1 = aggregate_afs(mid_files, FREQ_EXOMES_DIR, "mid")

# Merge into df_comb (Code used exactly as provided)
df_comb = merge(df_eur1,df_afr1,by="ID",all=T)
df_comb = merge(df_comb,df_amr1,by="ID",all=T)
df_comb = merge(df_comb,df_eas1,by="ID",all=T)
df_comb = merge(df_comb,df_mid1,by="ID",all=T)
df_comb = merge(df_comb,df_sas1,by="ID",all=T)

# Calculate PopMax AF (max_value) (Code used exactly as provided)
df_comb1 <- df_comb %>%
  mutate(
    max_value = do.call(pmax, c(df_comb[, c(2, 4, 6, 8, 10, 12)], na.rm = TRUE))
  )

# Apply PopMax filter <= 0.001 (Code used exactly as provided)
df_comb2 = df_comb1[which(df_comb1$max_value <= 0.001),]
print("AoU PopMax frequencies calculated and filtered.")


# --- 4. Load and Combine Annotated Files (Script 5 Output) ---
# NOTE: These sections were done in R in your original notebook (using fread inside loops)

# Define the file lists for aggregation
lof_wgs = list.files(LOF_DIR, pattern = "wgs.*tsv", full.names = TRUE)
lof_exome = list.files(LOF_DIR, pattern = "exomes.*tsv", full.names = TRUE)
missense_wgs = list.files(MISSENSE_DIR, pattern = "wgs.*tsv", full.names = TRUE)
missense_exome = list.files(MISSENSE_DIR, pattern = "exomes.*tsv", full.names = TRUE)
clinvar_wgs = list.files(CLINVAR_DIR, pattern = "wgs.*tsv", full.names = TRUE)
clinvar_exome = list.files(CLINVAR_DIR, pattern = "exomes.*tsv", full.names = TRUE)

# Read and combine LoF WGS and Exome separately (Code logic used exactly as provided)
df_lof_wgs = data.frame()
df_lof_exome = data.frame()
df_missense_wgs = data.frame()
df_missense_exome = data.frame()
df_clinvar_wgs = data.frame()
df_clinvar_exome = data.frame()

for (f in lof_wgs) { tmp = fread(f); df_lof_wgs = rbind(df_lof_wgs, tmp) }
for (f in lof_exome) { tmp = fread(f); df_lof_exome = rbind(df_lof_exome, tmp) }
for (f in missense_wgs) { tmp = fread(f); df_missense_wgs = rbind(df_missense_wgs, tmp) }
for (f in missense_exome) { tmp = fread(f); df_missense_exome = rbind(df_missense_exome, tmp) }
for (f in clinvar_wgs) { tmp = fread(f); df_clinvar_wgs = rbind(df_clinvar_wgs, tmp) }
for (f in clinvar_exome) { tmp = fread(f); df_clinvar_exome = rbind(df_clinvar_exome, tmp) }

# Final filtering by gene list (Code used exactly as provided)
genes_all = ['OBSCN','PLEKHM2','PRDM16','PSEN2','TNNI3K','LMNA','TNNT2','ACTN2','RIT1','NEXN','RYR2','TRIM63','ANKRD1','MYPN','NEBL','COX15','LDB3','RBM20','VCL','BAG3','CTNNA3','ILK','KCNQ1','MYBPC3','CRYAB','CSRP3','MYL2','ABCC9','CACNA1C','PTPN11','PKP2','MYH6','MYH7','TGFB3','TPM1','ALPK3','ACTC1','TJP1','CTF1','TCAP','GAA','JUP','DTNA','TTR','DSC2','DSG2','CDH2','MYOM1','TNNI3','CALR3','TTN','DES','JPH2','MYLK2','SCN5A','MYL3','CAV3','RAF1','TMEM43','TNNC1','SLC25A4','MYOZ2','PDLIM3','NKX2_5','SGCD','EYA4','LAMA4','DSP','MYO6','PLN','KCNH2','GATAD1','TBX20','FLNC','PRKAG2','KLF10','FXN','GLA','FHL1','LAMP2', 'NKX2-5'] # Added NKX2-5 from your filtering logic

df_lof_exomes1 = df_lof_exome[which(df_lof_exome$GENE %in% genes_all),]
df_missense_exomes1 = df_missense_exome[which(df_missense_exome$GENE %in% genes_all),]
df_clinvar_exomes1 = df_clinvar_exome[which(df_clinvar_exome$GENE %in% genes_all),]

df_lof_wgs1 =df_lof_wgs[which(df_lof_wgs$GENE %in% genes_all),]
df_missense_wgs1 =  df_missense_wgs[which(df_missense_wgs$GENE %in% genes_all),]
df_clinvar_wgs1 = df_clinvar_wgs[which(df_clinvar_wgs$GENE %in% genes_all),]

# Create ID field (Code used exactly as provided)
df_lof_exomes1$ID = paste0(df_lof_exomes1$CHROM,":",df_lof_exomes1$POS,":",df_lof_exomes1$REF,":",df_lof_exomes1$ALT)
df_lof_wgs1$ID = paste0(df_lof_wgs1$CHROM,":",df_lof_wgs1$POS,":",df_lof_wgs1$REF,":",df_lof_wgs1$ALT)
df_missense_exomes1$ID = paste0(df_missense_exomes1$CHROM,":",df_missense_exomes1$POS,":",df_missense_exomes1$REF,":",df_missense_exomes1$ALT)
df_missense_wgs1$ID = paste0(df_missense_wgs1$CHROM,":",df_missense_wgs1$POS,":",df_missense_wgs1$REF,":",df_missense_wgs1$ALT)
df_clinvar_exomes1$ID = paste0(df_clinvar_exomes1$CHROM,":",df_clinvar_exomes1$POS,":",df_clinvar_exomes1$REF,":",df_clinvar_exomes1$ALT)
df_clinvar_wgs1$ID = paste0(df_clinvar_wgs1$CHROM,":",df_clinvar_wgs1$POS,":",df_clinvar_wgs1$REF,":",df_clinvar_wgs1$ALT)
print("Annotated dataframes loaded and combined.")


# --- 5. Save Intermediate Filtered Data ---
save(df_lof_exomes1, df_missense_exomes1, df_clinvar_exomes1, 
     df_lof_wgs1, df_missense_wgs1, df_clinvar_wgs1, 
     df_comb2, gnomad1, 
     file=file.path(OUTPUT_DIR, "variants_categorized_12032025.rda"))

print("R Analysis Part 1 (Loading and PopMax) complete.")