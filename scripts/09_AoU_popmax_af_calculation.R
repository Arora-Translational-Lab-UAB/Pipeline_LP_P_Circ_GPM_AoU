# --- AoU PopMax Allele Frequency Calculation ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")

# Directories for input files
FREQ_CLINVAR_DIR <- file.path(BASE_DIR, "geno/Freq_Clinvar")
FREQ_EXOMES_DIR <- file.path(BASE_DIR, "geno/Freq_Exomes")

require(data.table)
library(dplyr)
library(splitstackshape)

# --- 2. Aggregate PLINK AF Files (ClinVar Subsets) ---

files_cln <- list.files(FREQ_CLINVAR_DIR, pattern = ".afreq", full.names = FALSE)
afr_files_cln <- files_cln[grep("afr_", files_cln)]
eur_files_cln <- files_cln[grep("eur_", files_cln)]
amr_files_cln <- files_cln[grep("amr_", files_cln)]
eas_files_cln <- files_cln[grep("eas_", files_cln)]
sas_files_cln <- files_cln[grep("sas_", files_cln)]
mid_files_cln <- files_cln[grep("mid_", files_cln)]

df_eur = data.frame()
df_afr = data.frame()
df_amr = data.frame()
df_eas = data.frame()
df_sas = data.frame()
df_mid= data.frame()

# Loop 1: Aggregate ClinVar Subsets (Code used exactly as provided)
for(i in 1:length(eur_files_cln)){
    tmp = fread(paste0("../../geno/Freq_Clinvar/",afr_files_cln[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(afr_files_cln[i]); df_afr = rbind(df_afr,tmp)
    tmp = fread(paste0("../../geno/Freq_Clinvar/",eas_files_cln[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(eas_files_cln[i]); df_eas = rbind(df_eas,tmp)
    tmp = fread(paste0("../../geno/Freq_Clinvar/",eur_files_cln[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(eur_files_cln[i]); df_eur = rbind(df_eur,tmp)
    tmp = fread(paste0("../../geno/Freq_Clinvar/",amr_files_cln[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(amr_files_cln[i]); df_amr = rbind(df_amr,tmp)
    tmp = fread(paste0("../../geno/Freq_Clinvar/",sas_files_cln[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(sas_files_cln[i]); df_sas = rbind(df_sas,tmp)
    tmp = fread(paste0("../../geno/Freq_Clinvar/",mid_files_cln[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(mid_files_cln[i]); df_mid = rbind(df_mid,tmp)
}

# --- 3. Aggregate PLINK AF Files (Exome Subsets - Appending to above DFs) ---

files_exo <- list.files(FREQ_EXOMES_DIR, pattern = ".afreq", full.names = FALSE)
afr_files_exo <- files_exo[grep("afr_", files_exo)]
eur_files_exo <- files_exo[grep("eur_", files_exo)]
amr_files_exo <- files_exo[grep("amr_", files_exo)]
eas_files_exo <- files_exo[grep("eas_", files_exo)]
sas_files_exo <- files_exo[grep("sas_", files_exo)]
mid_files_exo <- files_exo[grep("mid_", files_exo)]

# Loop 2: Aggregate Exome Subsets (Code used exactly as provided - Appends to existing DFs)
for(i in 1:length(eur_files_exo)){
    tmp = fread(paste0("../../geno/Freq_Exomes/",afr_files_exo[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(afr_files_exo[i]); df_afr = rbind(df_afr,tmp)
    tmp = fread(paste0("../../geno/Freq_Exomes/",eas_files_exo[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(eas_files_exo[i]); df_eas = rbind(df_eas,tmp)
    tmp = fread(paste0("../../geno/Freq_Exomes/",eur_files_exo[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(eur_files_exo[i]); df_eur = rbind(df_eur,tmp)
    tmp = fread(paste0("../../geno/Freq_Exomes/",amr_files_exo[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(amr_files_exo[i]); df_amr = rbind(df_amr,tmp)
    tmp = fread(paste0("../../geno/Freq_Exomes/",sas_files_exo[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(sas_files_exo[i]); df_sas = rbind(df_sas,tmp)
    tmp = fread(paste0("../../geno/Freq_Exomes/",mid_files_exo[i])); tmp = tmp[which(!(tmp$ALT_FREQ == 0)),]; tmp$GENE = paste0(mid_files_exo[i]); df_mid = rbind(df_mid,tmp)
}

# Convert to data.frame (Code used exactly as provided)
df_eur = as.data.frame(df_eur)
df_afr = as.data.frame(df_afr)
df_amr = as.data.frame(df_amr)
df_eas = as.data.frame(df_eas)
df_sas = as.data.frame(df_sas)
df_mid = as.data.frame(df_mid)


# --- 4. Rename Columns and Merge (Final AoU AF Data) ---

# Rename columns (Code used exactly as provided)
df_afr1 = df_afr[,c("ID","ALT_FREQS","GENE")]; colnames(df_afr1) = c("ID","af_afr","gene_afr")
df_amr1 = df_amr[,c("ID","ALT_FREQS","GENE")]; colnames(df_amr1) = c("ID","af_amr","gene_amr")
df_eas1 = df_eas[,c("ID","ALT_FREQS","GENE")]; colnames(df_eas1) = c("ID","af_eas","gene_eas")
df_sas1 = df_sas[,c("ID","ALT_FREQS","GENE")]; colnames(df_sas1) = c("ID","af_sas","gene_sas")
df_mid1 = df_mid[,c("ID","ALT_FREQS","GENE")]; colnames(df_mid1) = c("ID","af_mid","gene_mid")
df_eur1 = df_eur[,c("ID","ALT_FREQS","GENE")]; colnames(df_eur1) = c("ID","af_eur","gene_eur")

# Merge into df_comb (Code used exactly as provided)
df_comb = merge(df_eur1,df_afr1,by="ID",all=T)
df_comb = merge(df_comb,df_amr1,by="ID",all=T)
df_comb = merge(df_comb,df_eas1,by="ID",all=T)
df_comb = merge(df_comb,df_mid1,by="ID",all=T)
df_comb = merge(df_comb,df_sas1,by="ID",all=T)


# --- 5. Calculate AoU PopMax AF (max_value) ---

# Calculate PopMax AF (Code used exactly as provided)
df_comb1 <- df_comb %>%
  mutate(
    max_value = do.call(pmax, c(df_comb[, c(2, 4, 6, 8, 10, 12)], na.rm = TRUE))
  )

# Apply initial PopMax filter <= 0.001 (Code used exactly as provided)
df_comb2 = df_comb1[which(df_comb1$max_value <= 0.001),]


# --- 6. Final Combined PopMax Filter ---

# Load gnomAD filtered data (output from Script 5)
load(file.path(OUTPUT_DIR, "gnomad_popmax_filtered.rda"))

# Create the master list of variants passing PopMax filter from EITHER AoU or gnomAD
# NOTE: This uses df_comb2 (AoU PopMax filtered) and gnomad1 (gnomAD PopMax filtered)
af_popmax_aou_gnomad_filter <- bind_rows(
    df_comb2 %>% select(ID),
    gnomad1 %>% select(ID)
) %>%
    distinct(ID)

# Save results for next filtering step
save(df_comb1, df_comb2, gnomad1, af_popmax_aou_gnomad_filter, file = file.path(OUTPUT_DIR, "aou_gnomad_popmax_filters_v8_12032025.rda"))
print("AoU PopMax calculation and final filter list generation complete.")