# --- Individual Carrier Status Generation and Summary ---

# --- 1. Configuration and Libraries ---
BASE_DIR <- getwd() 
FINAL_EXPORT_DIR <- file.path(BASE_DIR, "final_variant_lists")
PLINK_INPUT_DIR <- file.path(BASE_DIR, "geno/exome_subset")
RAW_OUTPUT_DIR <- file.path(BASE_DIR, "RAW") # Directory for PLINK .raw files
GCS_CARRIERS_BUCKET="gs://<YOUR_SECURE_BUCKET>/CMP_Gene_Data/Raw_Files/Processed_Carriers/" # Censored GCS path

dir.create(RAW_OUTPUT_DIR, showWarnings = FALSE)
dir.create(file.path(BASE_DIR, "final_summary_tables"), showWarnings = FALSE)

require(data.table)
library(dplyr) # Used for final carrier summary

# Define target gene list for PLINK subsetting (Censored/Consolidated list of 17 genes)
plink_target_genes = 'LMNA|TNNT2|ACTN2|RBM20|BAG3|MYBPC3|MYL2|PKP2|MYH7|ACTC1|TPM1|JUP|DSC2|DSG2|TNNI3|TTN|DES|TMEM43|SCN5A|MYL3|TNNC1|PLN|DSP|FLNC'


# --- 2. Prepare PLINK Extraction Files and Execute ---

# List Plink subset files and clean prefixes
exome_subset <- list.files(PLINK_INPUT_DIR, full.names = FALSE)
exome_subset <- exome_subset[grep(".bed",exome_subset)]
exome_subset1 <- exome_subset[grepl(plink_target_genes, exome_subset)]
exome_subset_prefixes <- gsub(".bed","",exome_subset1)

# List final variant CSVs
setwd(FINAL_EXPORT_DIR)
variants <- list.files()
variants <- variants[grep(".csv",variants)]
variants <- variants[grep("exomes",variants)]

# NOTE: The provided code relies on the exact order of files in the 'variants' list.
# Ensure the exported CSV file names match this order:
# 1. lof_arvc_...
# 2. lof_dcm_...
# 3. lof_hcm_...
# 4. lp_p_arvc_...
# 5. lp_p_dcm_...
# 6. lp_p_hcm_...
# 7. missense_arvc_...
# 8. missense_dcm_...
# 9. missense_hcm_...

# Execute PLINK2 --recode A for all gene/variant combinations (Code used exactly as provided)
for(i in 1:length(exome_subset_prefixes)){
    prefix = exome_subset_prefixes[i]
    gene = sub("([^_]+)_.*", "\\1", prefix)
    
    # NOTE: The index [1] to [9] must correspond to the files listed above.
    
    # LoF Extractions (1, 2, 3)
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[1]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_arvc_lof"))))
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[2]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_dcm_lof"))))
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[3]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_hcm_lof"))))
    
    # P/LP Extractions (4, 5, 6)
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[4]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_arvc_lp_p"))))
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[5]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_dcm_lp_p"))))
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[6]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_hcm_lp_p"))))
    
    # Missense Extractions (7, 8, 9)
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[7]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_arvc_missense"))))
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[8]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_dcm_missense"))))
    system(paste0("plink2 --bfile ", file.path(PLINK_INPUT_DIR, prefix)," --extract ",variants[9]," --recode A --out ",file.path(RAW_OUTPUT_DIR, paste0(gene,"_hcm_missense"))))
}
setwd(BASE_DIR) # Reset working directory
print("PLINK extraction to .raw files complete.")


# --- 3. Carriers Processing (Genotype Conversion and Aggregation) ---

list_raw_files = list.files(RAW_OUTPUT_DIR, pattern = ".raw", full.names = FALSE)

# Initialize the final merged dataframe using a guaranteed file's IID list
merged_init = fread(file.path(RAW_OUTPUT_DIR, "DSP_dcm_lp_p.raw"))
patho_final = data.frame(IID = merged_init$IID)

# Main aggregation loop: processes each .raw file (Code used exactly as provided)
for (i in 1:length(list_raw_files)){
    gene = list_raw_files[i]
    print(paste("Processing:", gene))
    merged = as.data.frame(fread(file.path(RAW_OUTPUT_DIR, gene)))
    
    var12 = merged  
    patho = data.frame(IID = merged$IID)
    
    # Loop 1: Invert Genotype Coding (2->0, 1->1, 0->2) - Assuming PLINK A/A=2, A/B=1, B/B=0
    for(j in 7:ncol(merged)){
        var11 = j
        # Invert: 2 (Homo Alt) -> 0 (Homo Ref), 1 (Het) -> 1 (Het), 0 (Homo Ref) -> 2 (Homo Alt)
        # NOTE: Your original inversion logic is: 2->0, 1->1, 0->2. This is a simple allele count flip.
        var12[,var11] <- ifelse(var12[,var11] == 2,0,ifelse(var12[,var11] == 1,1,ifelse(var12[,var11] == 0,2,NA)))
    }
    merged = var12
    
    # Loop 2: Summing up variants per gene/category (total carrier count)
    if(ncol(merged) > 7) {
        merged[,paste0(gene)] = rowSums(merged[,c(7:ncol(merged))],na.rm=T)
        merged = merged[,c("IID",paste0(gene))]
        patho = merge(patho,merged,by=c("IID"))
    } else {
        merged[,paste0(gene)] = merged[,7]
        merged = merged[,c("IID",paste0(gene))]
        patho = merge(patho,merged,by=c("IID"))
    }
    patho_final = merge(patho_final,patho,by=c("IID"),all.x=T)
}

# Final cleanup and rename
colnames(patho_final) = gsub(".raw","",colnames(patho_final))
carriers_final = as.data.frame(patho_final)
print("Genotype processing and aggregation complete.")


# --- 4. Calculate Final Carrier Flags and Prevalence ---

# Define column groups (Code used exactly as provided)
lp_p = colnames(carriers_final)[grep("lp_p",colnames(carriers_final))]
dcm_lp_p = lp_p[grep("dcm",lp_p)]
hcm_lp_p = lp_p[grep("hcm",lp_p)]
arvc_lp_p = lp_p[grep("arvc",lp_p)]

lof = colnames(carriers_final)[grep("lof",colnames(carriers_final))]
dcm_lof = lof[grep("dcm",lof)]
hcm_lof = lof[grep("hcm",lof)]
arvc_lof = lof[grep("arvc",lof)]

missense = colnames(carriers_final)[grep("missense",colnames(carriers_final))]
dcm_missense = missense[grep("dcm",missense)]
hcm_missense = missense[grep("hcm",missense)]
arvc_missense = missense[grep("arvc",missense)]

# Calculate total counts (Code used exactly as provided)
carriers_final$lp_p_carriers = rowSums(carriers_final[,lp_p] > 0,na.rm=T)
carriers_final$lof_carriers = rowSums(carriers_final[,lof] > 0,na.rm=T)
carriers_final$missense_carriers = rowSums(carriers_final[,missense] > 0,na.rm=T)

carriers_final$dcm_lp_p_carriers = rowSums(carriers_final[,dcm_lp_p] > 0,na.rm=T)
carriers_final$dcm_lof_carriers = rowSums(carriers_final[,dcm_lof] > 0,na.rm=T)
carriers_final$dcm_missense_carriers = rowSums(carriers_final[,dcm_missense] > 0,na.rm=T)

carriers_final$hcm_lp_p_carriers = rowSums(carriers_final[,hcm_lp_p] > 0,na.rm=T)
carriers_final$hcm_lof_carriers = ifelse(carriers_final[,hcm_lof] > 0,1,0) # Note: Special logic here
carriers_final$hcm_missense_carriers = rowSums(carriers_final[,hcm_missense] > 0,na.rm=T)

carriers_final$arvc_lp_p_carriers = rowSums(carriers_final[,arvc_lp_p] > 0,na.rm=T)
carriers_final$arvc_lof_carriers = rowSums(carriers_final[,arvc_lof] > 0,na.rm=T)
carriers_final$arvc_missense_carriers = ifelse(carriers_final[,arvc_missense] > 0,1,0) # Note: Special logic here

# Calculate overall carriers (P/LP/PuPV/LoF/Missense combined)
all_carriers = colnames(carriers_final)[!grepl("carriers",colnames(carriers_final))]
all_carriers = all_carriers[!all_carriers %in% c("IID")]
hcm = all_carriers[grep("hcm",all_carriers)]
dcm = all_carriers[grep("dcm",all_carriers)]
arvc = all_carriers[grep("arvc",all_carriers)]

carriers_final$lp_p_pupv_carriers = rowSums(carriers_final[,all_carriers] > 0,na.rm=T)
carriers_final$hcm_lp_p_pupv_carriers = rowSums(carriers_final[,hcm] > 0,na.rm=T)
carriers_final$dcm_lp_p_pupv_carriers = rowSums(carriers_final[,dcm] > 0,na.rm=T)
carriers_final$arvc_lp_p_pupv_carriers = rowSums(carriers_final[,arvc] > 0,na.rm=T)

# Save the final carriers file
save(carriers_final, file=file.path(BASE_DIR, "carriers_updated_v8_11042025.rda"))

# Upload final carriers file to GCS (Censored)
system(paste0("gsutil cp -r ", file.path(BASE_DIR, "carriers_updated_v8_11042025.rda"), " ", GCS_CARRIERS_BUCKET))

print("Final carrier generation (Script 13) complete.")