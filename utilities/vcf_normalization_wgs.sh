#!/bin/bash

# --- VCF Processing, Normalization, and Annotation Extraction ---

set -o pipefail
set -o errexit
set -u

echo "--- Starting VCF Processing and Annotation Extraction (Step 4 & 5) ---"

# --- 1. Configuration (Paths sourced from 03_software_setup.sh) ---
# NOTE: Collaborator MUST source the paths from 03_software_setup.sh before running this.
BASE_DIR="$(pwd)"

# Directories used in previous steps
ANNOT_VCF_DIR="${BASE_DIR}/annot"
OUTPUT_TSV_DIR="${BASE_DIR}/annot_snpeff"
NORMALIZED_VCF_DIR="${BASE_DIR}/normalized_vcf"

# Categorization directories
OUTPUT_CLINVAR_DIR="${OUTPUT_TSV_DIR}/clinvar"
OUTPUT_LOF_DIR="${OUTPUT_TSV_DIR}/lof"
OUTPUT_MISSENSE_DIR="${OUTPUT_TSV_DIR}/missense"

mkdir -p "$OUTPUT_TSV_DIR"
mkdir -p "$NORMALIZED_VCF_DIR"
mkdir -p "$OUTPUT_CLINVAR_DIR"
mkdir -p "$OUTPUT_LOF_DIR"
mkdir -p "$OUTPUT_MISSENSE_DIR"

# Tool Paths (Sourced/Expected from 03_software_setup.sh)
JAVA="${PATH_JAVA:-/usr/bin/java}" 
BCFTOOLS="${PATH_BCFTOOLS:-/usr/bin/bcftools}"
SNPSIFT_JAR="${PATH_SNPSIFT_JAR:-snpEff/SnpSift.jar}"
SNPEFF_VCF_SCRIPT="${PATH_SNPEFF_SCRIPT:-snpEff/scripts/vcfEffOnePerLine.pl}"
REFERENCE_FASTA="${PATH_FASTA:-GRCh38.d1.vd1.fa}"


# --- 2. WGS VCF Normalization (New Step) ---

echo "--- 2. Normalizing WGS VCFs ---"

# Iterate over raw WGS annotated VCFs
for vcf_file in "${ANNOT_VCF_DIR}"/*wgs.annot.vcf; do
    [ -e "$vcf_file" ] || continue
    
    base_name=$(basename "${vcf_file%.annot.vcf}")
    output_file="${NORMALIZED_VCF_DIR}/${base_name}_normalized.vcf.gz"
    
    echo "Normalizing and filtering PASS for: $vcf_file"
    
    # Pipeline: Filter PASS -> Normalize -> Output gzipped VCF
    # Note: Use bcftools view -f PASS | bcftools norm -m -both...
    
    # Corrected bcftools command based on your preferred logic:
    "$BCFTOOLS" view -f PASS "$vcf_file" | \
    "$BCFTOOLS" norm -m -both -f "$REFERENCE_FASTA" -Oz -o "$output_file"
done

# --- 3. VCF to Full TSV Conversion (Extraction) ---

echo "--- 3. Converting Annotated VCFs to Full TSV ---"

# Define VCFs to process: Normalized WGS VCFs and Annotated Exome VCFs
VCFS_TO_EXTRACT=("${NORMALIZED_VCF_DIR}"/*_normalized.vcf.gz "${ANNOT_VCF_DIR}"/*_Exome.annot.vcf)

for vcf_file in "${VCFS_TO_EXTRACT[@]}"; do
    [ -e "$vcf_file" ] || continue
    
    # Determine output file name based on file suffix
    base_name=$(basename "$vcf_file")
    
    if [[ "$base_name" == *"_Exome"* ]]; then
        output_file_name="${base_name%_Exome.annot.vcf}_exomes_full_annot.tsv"
        # Exome: contains ID field
        FIELDS="CHROM POS ID REF ALT \"ANN[*].ALLELE\" \"ANN[*].EFFECT\" \"ANN[*].GENE\" \"ANN[*].FEATURE\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" GENEINFO CLNSIG CLNREVSTAT"
    else
        output_file_name="${base_name%_normalized.vcf.gz}_wgs_full_annot.tsv"
        # WGS: positional fields only (ID removed during normalization)
        FIELDS="CHROM POS REF ALT \"ANN[*].ALLELE\" \"ANN[*].EFFECT\" \"ANN[*].GENE\" \"ANN[*].FEATURE\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" GENEINFO CLNSIG CLNREVSTAT"
    fi
    
    output_file="${OUTPUT_TSV_DIR}/${output_file_name}"
    
    echo "Extracting fields for: $base_name -> $output_file_name"

    # Pipeline: VCF -> One-consequence-per-line (SnpEff Perl script) -> Extract fields (SnpSift)
    "$BCFTOOLS" view "$vcf_file" | \
    "$SNPEFF_VCF_SCRIPT" | \
    "$JAVA" -jar "$SNPSIFT_JAR" extractFields - $FIELDS \
    > "$output_file"
done


# --- 4. Categorization by Clinical/Functional Impact ---

echo "--- 4. Categorizing Annotations (P/LP, LoF, Missense) ---"

# Define the standard output header for categorized files
HEADER="CHROM\tPOS\tID\tREF\tALT\tALLELE\tEFFECT\tGENE\tFEATURE\tHGVS_C\tHGVS_P\tGENEINFO\tCLNSIG\tCLNREVSTAT"

for f in "${OUTPUT_TSV_DIR}"/*_full_annot.tsv; do
    base_name=$(basename "$f")

    # 4.1 ClinVar Pathogenic/Likely Pathogenic (P/LP)
    outf_clinvar="${OUTPUT_CLINVAR_DIR}/${base_name/_full_annot.tsv/_clinvar_pathogenic.tsv}"
    echo -e "$HEADER" > "$outf_clinvar"
    grep -E 'Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic' "$f" | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' >> "$outf_clinvar"

    # 4.2 Loss of Function (LoF)
    outf_lof="${OUTPUT_LOF_DIR}/${base_name/_full_annot.tsv/_lof_variants.tsv}"
    echo -e "$HEADER" > "$outf_lof"
    grep -E 'stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|splice' "$f" | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' >> "$outf_lof"

    # 4.3 Missense
    outf_missense="${OUTPUT_MISSENSE_DIR}/${base_name/_full_annot.tsv/_missense_variants.tsv}"
    echo -e "$HEADER" > "$outf_missense"
    grep -E 'missense' "$f" | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' >> "$outf_missense"
done

echo "Annotation processing, normalization, and categorization complete. TSV files are in ${OUTPUT_TSV_DIR}."