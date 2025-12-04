# --- Gene-Specific and Functional Filtering ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
OUTPUT_DIR <- file.path(BASE_DIR, "analysis_data")
# NOTE: SCN5A annotation requires separate tool installation and execution (InterVar/ANNOVAR).
# These files must be manually generated before running this script.
SCN5A_ANNOVAR_DIR <- file.path(BASE_DIR, "annovar_intervar_out") # Placeholder directory for SCN5A files

require(data.table)
library(splitstackshape)
library(dplyr)

# Load previously cleaned and categorized variants
load(file.path(OUTPUT_DIR, "variants_categorized_and_cleaned_12032025.rda"))
print("Loaded cleaned dataframes (df_lof_exomes1, df_clinvar_exomes1, etc.).")


# --- 2. ClinVar Filtration (2-Star Evidence) ---

## Clinvar P/LP filtration based on evidence status (Code used exactly as provided)
clinvar_exomes_based_on_evidence = df_clinvar_exomes1[which(df_clinvar_exomes1$CLNSIG %in% c("Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic","Benign","Benign/Likely_benign","Likely_benign")),]
clinvar_wgs_based_on_evidence = df_clinvar_wgs1[which(df_clinvar_wgs1$CLNSIG %in% c("Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic","Benign","Benign/Likely_benign","Likely_benign")),]

# Filter for 2-star status (Code used exactly as provided)
clinvar_2star_exomes = clinvar_exomes_based_on_evidence[which(clinvar_exomes_based_on_evidence$CLNREVSTAT %in%  c("criteria_provided,_multiple_submitters,_no_conflicts","reviewed_by_expert_panel")),] 
clinvar_2star_wgs = clinvar_wgs_based_on_evidence[which(clinvar_wgs_based_on_evidence$CLNREVSTAT %in%  c("criteria_provided,_multiple_submitters,_no_conflicts","reviewed_by_expert_panel")),] 
print(paste("ClinVar filtered to 2-star+ status. Exomes:", nrow(clinvar_2star_exomes), "WGS:", nrow(clinvar_2star_wgs)))


# --- 3. TTN LoF Filtration (High-PSI Exons) ---

## Separate TTN LoF variants (Code used exactly as provided)
TTN_exomes = df_lof_exomes1[which(df_lof_exomes1$GENE %in% c("TTN")),]
TTN_wgs = df_lof_wgs1[which(df_lof_wgs1$GENE %in% c("TTN")),]

## Load Hi-PSI regions (NOTE: Collaborator must ensure this file is present locally)
hiPSI = fread("Exonic_Regions_PSI_90_percent.txt")

# Loop to filter TTN variants by Hi-PSI exon coordinates (Code used exactly as provided)
final_exomes = data.frame()
final_wgs = data.frame()
for(i in 1:nrow(hiPSI)){
  start = hiPSI$'Hg38 start'[i] - 10
  end = hiPSI$'Hg38 end'[i] + 10
  
  temp1 = TTN_exomes[which(TTN_exomes$POS >= start & TTN_exomes$POS <= end),]
  final_exomes = rbind(temp1,final_exomes)
    
  temp1 = TTN_wgs[which(TTN_wgs$POS >= start & TTN_wgs$POS <= end),]
  final_wgs = rbind(temp1,final_wgs)
}

# Final TTN variants inclusion (Code used exactly as provided)
TTN_hipsi_exomes_genes = TTN_exomes[which(TTN_exomes$ID %in% final_exomes$ID),]
TTN_hipsi_wgs_genes = TTN_wgs[which(TTN_wgs$ID %in% final_wgs$ID),]

## Recombine LoF variants (Code used exactly as provided)
df_lof_exomes2 = df_lof_exomes1[which(!(df_lof_exomes1$GENE %in% c("TTN"))),]
df_lof_wgs2 = df_lof_wgs1[which(!(df_lof_wgs1$GENE %in% c("TTN"))),]

df_lof_exomes_ttn_hipsi2 = rbind(df_lof_exomes2,TTN_hipsi_exomes_genes)
df_lof_wgs_ttn_hipsi2 = rbind(df_lof_wgs2,TTN_hipsi_wgs_genes)
print(paste("TTN LoF filtered by Hi-PSI regions and recombined."))


# --- 4. SCN5A Missense Filtration (REVEL >= 0.65 AND InterVar Pathogenic) ---

# NOTE: This section relies on external ANNOVAR/InterVar tool output files.
# Censored GCS paths are placeholders.

# Load SCN5A specific annotation files (Code used exactly as provided)
# ANNOVAR Output: Contains REVEL scores
SCN5A_annovar_exome = fread(file.path(SCN5A_ANNOVAR_DIR, "SCN5A_exome.hg38_multianno.txt"))
SCN5A_annovar_wgs = fread(file.path(SCN5A_ANNOVAR_DIR, "SCN5A_wgs.hg38_multianno.txt"))
# InterVar Output: Contains Pathogenicity calls
SCN5A_Intervar_exome = fread(file.path(SCN5A_ANNOVAR_DIR, "SCN5A_exome.hg38_multianno.txt.intervar"))
SCN5A_Intervar_wgs = fread(file.path(SCN5A_ANNOVAR_DIR, "SCN5A_wgs.hg38_multianno.txt.intervar"))


# Filter ANNOVAR output by REVEL score >= 0.65 and create ID (Code used exactly as provided)
SCN5A_annovar_exome = SCN5A_annovar_exome[which(SCN5A_annovar_exome$REVEL_score >= 0.65),]
SCN5A_annovar_exome$ID = paste0(SCN5A_annovar_exome$Chr,":",SCN5A_annovar_exome$Start,":",SCN5A_annovar_exome$Ref,":",SCN5A_annovar_exome$Alt)
SCN5A_annovar_wgs = SCN5A_annovar_wgs[which(SCN5A_annovar_wgs$REVEL_score >= 0.65),]
SCN5A_annovar_wgs$ID = paste0(SCN5A_annovar_exome$Chr,":",SCN5A_annovar_wgs$Start,":",SCN5A_annovar_wgs$Ref,":",SCN5A_annovar_wgs$Alt) # NOTE: Used SCN5A_annovar_exome$Chr in WGS ID calculation, maintaining original error/structure.

# Filter InterVar output by Pathogenic call and create ID (Code used exactly as provided)
SCN5A_Intervar_exome1 = SCN5A_Intervar_exome[grep("pathogenic|Pathogenic",SCN5A_Intervar_exome$`InterVar: InterVar and Evidence`),]
SCN5A_Intervar_wgs1 = SCN5A_Intervar_wgs[grep("pathogenic|Pathogenic",SCN5A_Intervar_wgs$`InterVar: InterVar and Evidence`),]
SCN5A_Intervar_exome1$ID = paste0("chr",SCN5A_Intervar_exome1$'#Chr',":",SCN5A_Intervar_exome1$Start,":",SCN5A_Intervar_exome1$Ref,":",SCN5A_Intervar_exome1$Alt)
SCN5A_Intervar_wgs1$ID = paste0("chr",SCN5A_Intervar_wgs1$'#Chr',":",SCN5A_Intervar_wgs1$Start,":",SCN5A_Intervar_wgs1$Ref,":",SCN5A_Intervar_wgs1$Alt)

# Combine REVEL and InterVar criteria to get final SCN5A variants to include (Code used exactly as provided)
SCN5A_exome_var_to_include = SCN5A_annovar_exome[which(SCN5A_annovar_exome$ID %in% SCN5A_Intervar_exome1$ID),]
SCN5A_wgs_var_to_include = SCN5A_annovar_wgs[which(SCN5A_annovar_wgs$ID %in% SCN5A_Intervar_wgs1$ID),]

# Intersect SCN5A filtered lists with original Missense lists (Code used exactly as provided)
SCN5A_revel_score_intervar_exome =  df_missense_exomes1[which(df_missense_exomes1$ID %in% SCN5A_exome_var_to_include$ID),]
SCN5A_revel_score_intervar_wgs = df_missense_wgs1[which(df_missense_wgs1$ID %in% SCN5A_wgs_var_to_include$ID),]
print(paste("SCN5A missense variants filtered by REVEL and InterVar consensus."))


# --- 5. General Missense Filtration (REVEL >= 0.65) ---

# Load general gnomAD annotation data containing REVEL scores (NOTE: Collaborator must ensure this file is present locally)
system("mv gnomad_annotations_all_01282024.tsv.bgz gnomad_annotations_all_01282024.tsv.gz") # Renaming if downloaded as .bgz

revel_scores = fread("gnomad_annotations_all_01282024.tsv.gz")

# Split multi-valued REVEL score column and filter by score >= 0.65 (Code used exactly as provided)
revel_scores = revel_scores[which(!(is.na(revel_scores$revel_score))),]
revel_scores$revel_score = gsub("\\]","",revel_scores$revel_score)
revel_scores$revel_score = gsub("\\[","",revel_scores$revel_score)
revel_scores = cSplit(revel_scores,"revel_score",",")

revel_scores1 = data.frame() # Initialize combined REVEL passes

# Loop through potential split REVEL columns (Code used exactly as provided)
revel_scores1 <- rbind(revel_scores1, revel_scores[which(revel_scores$revel_score_1 >= 0.65 & !is.na(revel_scores$revel_score_1)),])
revel_scores1 <- rbind(revel_scores1, revel_scores[which(revel_scores$revel_score_2 >= 0.65 & !is.na(revel_scores$revel_score_2)),])
revel_scores1 <- rbind(revel_scores1, revel_scores[which(revel_scores$revel_score_3 >= 0.65 & !is.na(revel_scores$revel_score_3)),])
revel_scores1 <- rbind(revel_scores1, revel_scores[which(revel_scores$revel_score_4 >= 0.65 & !is.na(revel_scores$revel_score_4)),])
revel_scores1 <- rbind(revel_scores1, revel_scores[which(revel_scores$revel_score_5 >= 0.65 & !is.na(revel_scores$revel_score_5)),])

# Create ID for merging (Code used exactly as provided)
revel_scores1$ID = paste0(revel_scores1$locus,":",revel_scores1$ref,":",revel_scores1$alt)

## Separate SCN5A from main missense list (Code used exactly as provided)
df_missense_exomes2 = df_missense_exomes1[which(!(df_missense_exomes1$GENE %in% c("SCN5A"))),]
df_missense_wgs2 = df_missense_wgs1[which(!(df_missense_wgs1$GENE %in% c("SCN5A"))),]

## Filter remaining missense variants by REVEL >= 0.65 (Code used exactly as provided)
df_missense_exomes3 = df_missense_exomes2[which(df_missense_exomes2$ID %in% revel_scores1$ID),]
df_missense_wgs3 = df_missense_wgs2[which(df_missense_wgs2$ID %in% revel_scores1$ID),]

## Recombine SCN5A and general missense (Code used exactly as provided)
df_missense_exomes4 = rbind(df_missense_exomes3,SCN5A_revel_score_intervar_exome)
df_missense_wgs4 = rbind(df_missense_wgs3,SCN5A_revel_score_intervar_wgs)
print(paste("General missense variants filtered by REVEL and recombined with SCN5A."))


# --- 6. Save Final Filtered Variants ---

save(clinvar_2star_exomes, clinvar_2star_wgs, 
     df_lof_exomes_ttn_hipsi2, df_lof_wgs_ttn_hipsi2, 
     df_missense_exomes4, df_missense_wgs4,
     file=file.path(OUTPUT_DIR, "verified_initial_filtered_v8_12032025.rda"))

print("Initial functional and gene-specific filtering (Script 9) complete.")