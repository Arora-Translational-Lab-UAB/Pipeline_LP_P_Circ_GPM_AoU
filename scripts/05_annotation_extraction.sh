#!/bin/bash

# --- Annotation Extraction and Categorization Script ---
# This script converts annotated VCFs into tabular format and categorizes them.

set -o pipefail
set -o errexit
set -u

echo "--- Starting Annotation Extraction (Step 5) ---"

# --- 1. Configuration (Paths sourced from 03_software_setup.sh) ---
# NOTE: This script must be run after 03_software_setup.sh and 04_vcf_annotation.sh
BASE_DIR="$(pwd)"

# Directories used in previous steps
ANNOT_VCF_DIR="${BASE_DIR}/annot"
OUTPUT_TSV_DIR="${BASE_DIR}/annot_snpeff"

# Categorization directories
OUTPUT_CLINVAR_DIR="${OUTPUT_TSV_DIR}/clinvar"
OUTPUT_LOF_DIR="${OUTPUT_TSV_DIR}/lof"
OUTPUT_MISSENSE_DIR="${OUTPUT_TSV_DIR}/missense"

mkdir -p "$OUTPUT_TSV_DIR"
mkdir -p "$OUTPUT_CLINVAR_DIR"
mkdir -p "$OUTPUT_LOF_DIR"
mkdir -p "$OUTPUT_MISSENSE_DIR"

# Tool Paths (Sourced/Expected from 03_software_setup.sh)
JAVA="${PATH_JAVA:-/usr/bin/java}" 
SNPSIFT_JAR="${PATH_SNPSIFT_JAR:-snpEff/SnpSift.jar}"
SNPEFF_VCF_SCRIPT="${PATH_SNPEFF_SCRIPT:-snpEff/scripts/vcfEffOnePerLine.pl}"


# --- 2. VCF to Full TSV Conversion (All 80 Genes) ---

echo "--- 2. Converting Annotated VCFs to Full TSV (One Consequence Per Line) ---"

# Iterate over all annotated VCFs (Exome and WGS)
for vcf_file in "${ANNOT_VCF_DIR}"/*.annot.vcf; do
    [ -e "$vcf_file" ] || continue
    
    # Determine the output file name based on Exome/WGS and gene name
    base_name=$(basename "${vcf_file%.annot.vcf}")
    
    # Use standard logic to determine if it's WGS or Exome based on naming convention
    if [[ "$base_name" == *"_Exome"* ]]; then
        output_file_suffix="_exomes_full_annot.tsv"
        # Exome files contain the ID field
        FIELDS="CHROM POS ID REF ALT \"ANN[*].ALLELE\" \"ANN[*].EFFECT\" \"ANN[*].GENE\" \"ANN[*].FEATURE\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" GENEINFO CLNSIG CLNREVSTAT"
    else
        output_file_suffix="_wgs_full_annot.tsv"
        # WGS files from normalization may lack the original ID field, use positional fields
        FIELDS="CHROM POS REF ALT \"ANN[*].ALLELE\" \"ANN[*].EFFECT\" \"ANN[*].GENE\" \"ANN[*].FEATURE\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" GENEINFO CLNSIG CLNREVSTAT"
    fi
    
    output_file="${OUTPUT_TSV_DIR}/${base_name}${output_file_suffix}"
    
    echo "Extracting fields for: ${base_name} -> ${output_file}"

    # Pipeline: VCF -> One-consequence-per-line (SnpEff Perl script) -> Extract fields (SnpSift)
    cat "$vcf_file" | \
    "$SNPEFF_VCF_SCRIPT" | \
    "$JAVA" -jar "$SNPSIFT_JAR" extractFields - $FIELDS \
    > "$output_file"
done


# --- 3. Categorization by Clinical/Functional Impact ---
# Filters the full TSV output into specific functional categories.

echo "--- 3. Categorizing Annotations (P/LP, LoF, Missense) ---"

# Define the standard output header for categorized files
HEADER="CHROM\tPOS\tID\tREF\tALT\tALLELE\tEFFECT\tGENE\tFEATURE\tHGVS_C\tHGVS_P\tGENEINFO\tCLNSIG\tCLNREVSTAT"

for f in "${OUTPUT_TSV_DIR}"/*_full_annot.tsv; do
    base_name=$(basename "$f")

    # 3.1 ClinVar Pathogenic/Likely Pathogenic (P/LP)
    outf_clinvar="${OUTPUT_CLINVAR_DIR}/${base_name/_full_annot.tsv/_clinvar_pathogenic.tsv}"
    echo -e "$HEADER" > "$outf_clinvar"
    grep -E 'Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic' "$f" | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' >> "$outf_clinvar"
    echo "Categorized ClinVar P/LP: $(basename $outf_clinvar)"

    # 3.2 Loss of Function (LoF)
    outf_lof="${OUTPUT_LOF_DIR}/${base_name/_full_annot.tsv/_lof_variants.tsv}"
    echo -e "$HEADER" > "$outf_lof"
    grep -E 'stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|splice' "$f" | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' >> "$outf_lof"
    echo "Categorized LoF: $(basename $outf_lof)"

    # 3.3 Missense
    outf_missense="${OUTPUT_MISSENSE_DIR}/${base_name/_full_annot.tsv/_missense_variants.tsv}"
    echo -e "$HEADER" > "$outf_missense"
    grep -E 'missense' "$f" | \
    awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' >> "$outf_missense"
    echo "Categorized Missense: $(basename $outf_missense)"
done

echo "Annotation extraction and categorization complete. TSV files ready for R analysis."