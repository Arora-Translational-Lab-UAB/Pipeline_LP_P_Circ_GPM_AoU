import hail as hl
import os

# --- Configuration (UPDATE THESE PLACEHOLDERS) ---
# NOTE: Replace the project ID placeholder before running
PROJECT_ID = "<YOUR_PROJECT_ID>" 
GNOMAD_HT_PATH = "gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht"
EXPORT_PATH = "gs://<YOUR_SECURE_BUCKET>/CMP_Gene_Data/GnomAD_Annotations/gnomad_annotations_v4.1_export.tsv.bgz"
# You may need to copy the HT locally first depending on permissions/location:
# !gsutil -u $GOOGLE_PROJECT cp gs://gcp-public-data--gnomad/release/4.1/lof_curation/incomplete_penetrance_curation_results.csv .

# Initialize Hail
hl.init(default_reference='GRCh38',min_block_size=100,gcs_requester_pays_configuration=PROJECT_ID)
hl.cite_hail()

# --- 1. Load the gnomAD Hail Table ---
print(f"Reading gnomAD Table from: {GNOMAD_HT_PATH}")
ht = hl.read_table(GNOMAD_HT_PATH) 

# --- 2. Select and Flatten Top-Level Fields ---
# Select only the relevant fields, dropping complex structures not directly needed.
ht = ht.select(
    # Region flags
    region_flags = ht.region_flags,
    
    # Joint data (primary focus) - Select specific sub-fields
    joint_freq = ht.joint.freq.map(lambda x: x.select('AC', 'AF', 'AN', 'homozygote_count')),
    joint_faf = ht.joint.faf,
    joint_fafmax = ht.joint.fafmax.select('faf95_max', 'faf95_max_gen_anc', 'faf99_max', 'faf99_max_gen_anc'),
    joint_grpmax = ht.joint.grpmax.select('AC', 'AF', 'AN', 'homozygote_count', 'gen_anc'),
    
    # Exomes data 
    exomes_freq = ht.exomes.freq.map(lambda x: x.select('AC', 'AF', 'AN', 'homozygote_count')),
    exomes_faf = ht.exomes.faf,
    exomes_fafmax = ht.exomes.fafmax.select('faf95_max', 'faf95_max_gen_anc', 'faf99_max', 'faf99_max_gen_anc'),
    exomes_grpmax = ht.exomes.grpmax.select('AC', 'AF', 'AN', 'homozygote_count', 'gen_anc'),
    
    # Genomes data 
    genomes_freq = ht.genomes.freq.map(lambda x: x.select('AC', 'AF', 'AN', 'homozygote_count')),
    genomes_faf = ht.genomes.faf,
    genomes_fafmax = ht.genomes.fafmax.select('faf95_max', 'faf95_max_gen_anc', 'faf99_max', 'faf99_max_gen_anc'),
    genomes_grpmax = ht.genomes.grpmax.select('AC', 'AF', 'AN', 'homozygote_count', 'gen_anc'),
    
    # Frequency comparison stats
    freq_comparison_stats = ht.freq_comparison_stats
)


# --- 3. Flatten Annotations for Export ---
# Extract scalar and vector elements and create flat columns.
ht = ht.select(
    # Locus and Allele info
    ref = ht.alleles[0],
    alt = ht.alleles[1],
    is_snp = hl.is_snp(ht.alleles[0], ht.alleles[1]),

    # Joint data (Overall AF, AC, AN) - Index 0 is typically the aggregate/total frequency
    joint_AF = ht.joint_freq[0].AF,
    joint_AC = ht.joint_freq[0].AC,
    joint_AN = ht.joint_freq[0].AN,
    
    # Joint group max (Popmax) data
    joint_grpmax_AF = ht.joint_grpmax.AF,
    joint_grpmax_AC = ht.joint_grpmax.AC,
    joint_grpmax_AN = ht.joint_grpmax.AN,
    joint_grpmax_gen_anc = ht.joint_grpmax.gen_anc,

    # Filtering Allele Frequency (FAF) data
    joint_faf95 = ht.joint_faf[0].faf95,
    joint_faf99 = ht.joint_faf[0].faf99,
    joint_faf95_max = ht.joint_fafmax.faf95_max,
    joint_faf95_max_gen_anc = ht.joint_fafmax.faf95_max_gen_anc,
    joint_faf99_max = ht.joint_fafmax.faf99_max,
    joint_faf99_max_gen_anc = ht.joint_fafmax.faf99_max_gen_anc,

    # Region flags (Boolean filters)
    fail_interval_qc = ht.region_flags.fail_interval_qc,
    outside_broad_capture_region = ht.region_flags.outside_broad_capture_region,
    outside_ukb_capture_region = ht.region_flags.outside_ukb_capture_region,
    not_called_in_exomes = ht.region_flags.not_called_in_exomes,
    not_called_in_genomes = ht.region_flags.not_called_in_genomes,

    # Frequency comparison stats (For testing consistency between Exomes/Genomes)
    freq_cong_pval = ht.freq_comparison_stats.contingency_table_test[0].p_value,
    freq_cmh_pval = ht.freq_comparison_stats.cochran_mantel_haenszel_test.p_value,

    # Ancestry-specific allele frequencies (Arrays to be mapped)
    ancestry_AF = ht.joint_freq.map(lambda x: x.AF),
    ancestry_AC = ht.joint_freq.map(lambda x: x.AC),
    ancestry_AN = ht.joint_freq.map(lambda x: x.AN)
)

# --- 4. Add Ancestry-Specific Frequency Columns ---
# Label and extract elements from the frequency arrays.
# NOTE: Verify ancestry order matches the gnomAD HT structure!
ancestry_labels = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'sas', 'oth'] 

for i, label in enumerate(ancestry_labels):
    ht = ht.annotate(**{
        f'{label}_AF': ht.ancestry_AF[i],
        f'{label}_AC': ht.ancestry_AC[i],
        f'{label}_AN': ht.ancestry_AN[i]
    })

# Remove the temporary ancestry arrays and global fields that won't be exported
ht = ht.drop('ancestry_AF', 'ancestry_AC', 'ancestry_AN')
ht = ht.drop('freq_cong_odds_ratio', 'freq_cmh_chisq', 'freq_stat_union_pval', 
             'freq_stat_union_test_name', 'freq_stat_union_gen_ancs',
             'exomes_globals','genomes_globals','joint_globals', 'region_flags') # Dropping complex fields
ht = ht.drop('exomes_freq', 'exomes_faf', 'exomes_fafmax', 'exomes_grpmax', 'genomes_freq', 'genomes_faf', 'genomes_fafmax', 'genomes_grpmax', 'freq_comparison_stats')


# --- 5. Export Final Table ---

print(f"\nExporting final Hail Table to compressed TSV: {EXPORT_PATH}")
# Export to a gzipped TSV file, keeping the variant key (locus and alleles)
ht.export(EXPORT_PATH, header=True, delimiter='\t')

print("GnomAD annotation extraction and export complete.")