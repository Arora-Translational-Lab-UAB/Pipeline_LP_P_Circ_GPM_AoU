import hail as hl
import os

# --- 1. Configuration (CENSORED PLACEHOLDERS) ---
# NOTE: Collaborator MUST replace these placeholders with actual paths/IDs.
EXOME_MT_PATH = "gs://<CONTROLLED_DATASET_BUCKET>/exome/splitMT/hail.mt"
VDS_PATH = "gs://<CONTROLLED_DATASET_BUCKET>/vds/hail.vds"
VCF_OUTPUT_DIR = 'gs://<YOUR_SECURE_BUCKET>/CMP_Gene_Data/VCF_Format/'
PROJECT_ID = "<YOUR_PROJECT_ID>" 

hl.init(default_reference='GRCh38',min_block_size=100,gcs_requester_pays_configuration=PROJECT_ID)

# Full lists of 80 genes and loci for extraction
genes = ['OBSCN','PLEKHM2','PRDM16','PSEN2','TNNI3K','LMNA','TNNT2','ACTN2','RIT1','NEXN','RYR2','TRIM63','ANKRD1','MYPN','NEBL','COX15','LDB3','RBM20','VCL','BAG3','CTNNA3','ILK','KCNQ1','MYBPC3','CRYAB','CSRP3','MYL2','ABCC9','CACNA1C','PTPN11','PKP2','MYH6','MYH7','TGFB3','TPM1','ALPK3','ACTC1','TJP1','CTF1','TCAP','GAA','JUP','DTNA','TTR','DSC2','DSG2','CDH2','MYOM1','TNNI3','CALR3','TTN','DES','JPH2','MYLK2','SCN5A','MYL3','CAV3','RAF1','TMEM43','TNNC1','SLC25A4','MYOZ2','PDLIM3','NKX2_5','SGCD','EYA4','LAMA4','DSP','MYO6','PLN','KCNH2','GATAD1','TBX20','FLNC','PRKAG2','KLF10','FXN','GLA','FHL1','LAMP2']
locus = ['chr1:228208044-228378877','chr1:15681506-15734770','chr1:3069168-3438622','chr1:226870184-226927727','chr1:74235387-74544429','chr1:155682573-156540082','chr1:200959008-201777829','chr1:236264141-237164632','chr1:155497808-156311405','chr1:77488624-78343896','chr1:237042184-237833989','chr1:26051301-26068437','chr10:90912096-90921088','chr10:68087897-68212018','chr10:20779973-21293051','chr10:99294293-100132132','chr10:86268511-87136073','chr10:110244336-111239469','chr10:73595193-74521364','chr10:119251380-120077820','chr10:65912457-67763638','chr11:6603708-6610875','chr11:2444684-2849106','chr11:46952957-47774254','chr11:111508564-112323741','chr11:18782030-19610572','chr12:110510819-111321444','chr12:21397389-22342544','chr12:1569552-3097951','chr12:112018351-112909919','chr12:32390755-33296777','chr14:23381982-23408274','chr14:23012740-23835678','chr14:75958097-75983012','chr15:62642632-63471916','chr15:84417356-85273480','chr15:34390230-35195550','chr15:29699367-29969050','chr16:30895824-30903561','chr17:39665349-39666555','chr17:79701535-80519882','chr17:41354604-42186932','chr18:34493291-34891845','chr18:31157009-31998834','chr18:30658840-31502523','chr18:31098177-31949009','chr18:27932878-28177947','chr18:3066807-3220109','chr19:54751767-55557774','chr19:16479061-16496168','chr2:178125989-179230803','chr2:219018377-219826736','chr20:43706590-44587947','chr20:31819308-31834690','chr3:38148057-39049688','chr3:46435123-47282173','chr3:8333802-9241809','chr3:12183601-13064126','chr3:13725006-14543681','chr3:52051100-52854042','chr4:184743266-185550383','chr4:119135832-119187790','chr4:185500660-185535508','chr5:173232109-173235322','chr5:155728636-156767789','chr6:133240340-133532129','chr6:112107931-112254940','chr6:7141617-7986715','chr6:75349192-76319538','chr6:118148296-118961717','chr7:150944961-150978322','chr7:92447448-92494632','chr7:35199936-35254101','chr7:128430406-129259275','chr7:151156124-152277126','chr8:102648784-102655726','chr9:68635751-69479077','chrX:100993273-101808013','chrX:135746702-136611360','chrX:120026148-120869366']


# --- 2. Exome VCF Export (Sites Only) ---

def export_exome_vcf():
    """Reads Exome MT, filters by gene region, and exports sites-only VCFs."""
    print("--- Starting Exome VCF Export ---")
    mt = hl.read_matrix_table(EXOME_MT_PATH)

    if len(locus) > 0:
        for i in range(len(locus)):
            print(f"Processing Exome VCF for gene {i}: {genes[i]}")
            var = genes[i]
            loci = locus[i]
            
            # Filter the MatrixTable to the current locus
            aou_mt = hl.filter_intervals(mt,[hl.parse_locus_interval(loci)],keep=True)
            
            # Extract rows (sites only) and export as VCF
            row_table = aou_mt.rows()
            hl.export_vcf(row_table,f'{VCF_OUTPUT_DIR}{var}_Exome.vcf.bgz', overwrite=True)
    print("Exome VCF Export Complete.")


# --- 3. WGS VCF Export (Sites Only) ---

def export_wgs_vcf():
    """Reads WGS VDS, filters by gene region, and exports sites-only VCFs."""
    print("\n--- Starting WGS VCF Export ---")
    vds = hl.vds.read_vds(VDS_PATH)

    # Filter the VDS to contain only the 80 gene regions
    vds = hl.vds.filter_intervals(vds,[hl.parse_locus_interval(x) for x in locus])
    mt1 = vds.variant_data

    if len(locus) > 0:
        for i in range(len(locus)):
            print(f"Processing WGS VCF for gene {i}: {genes[i]}")
            var = genes[i]
            loci = locus[i]
            
            # Filter the variant data MT to the current locus
            aou_mt = hl.filter_intervals(mt1,[hl.parse_locus_interval(loci)],keep=True)
            
            # Extract rows (sites only) and export as VCF
            row_table = aou_mt.rows()
            hl.export_vcf(row_table,f'{VCF_OUTPUT_DIR}{var}wgs.vcf.bgz', overwrite=True)
    print("WGS VCF Export Complete.")


if __name__ == "__main__":
    export_exome_vcf()
    export_wgs_vcf()