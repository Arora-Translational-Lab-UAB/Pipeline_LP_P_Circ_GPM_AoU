#!/bin/bash

# --- PLINK Subsetting Script (Utility executed by dsub) ---

set -o pipefail
set -o errexit

chmod +x "${plink}"

bfile_prefix1="${input_bed1%.bed}"
bfile_prefix2="${input_bed2%.bed}"
bfile_prefix3="${input_bed3%.bed}"

region_start=$(echo "$region" | cut -d'-' -f1)
region_end=$(echo "$region" | cut -d'-' -f2)
region_string="${chrom}:${region}"

log_file="${gene}_${chrom}_${region}_processing.log"
exec > >(tee "${log_file}") 2>&1

echo "Processing ${gene} region ${region_string}..."

"${plink}" --bfile "${bfile_prefix1}" --chr "${chrom}" --from-bp "${region_start}" --to-bp "${region_end}" --make-bed --out "${gene}_${chrom}_${region}_1"
"${plink}" --bfile "${bfile_prefix2}" --chr "${chrom}" --from-bp "${region_start}" --to-bp "${region_end}" --make-bed --out "${gene}_${chrom}_${region}_2"
"${plink}" --bfile "${bfile_prefix3}" --chr "${chrom}" --from-bp "${region_start}" --to-bp "${region_end}" --make-bed --out "${gene}_${chrom}_${region}_3"

echo "Moving files to GCS outputs..."
mv "${gene}_${chrom}_${region}_1.bed" "${OUTPUT_PATH1}/"
mv "${gene}_${chrom}_${region}_1.bim" "${OUTPUT_PATH1}/"
mv "${gene}_${chrom}_${region}_1.fam" "${OUTPUT_PATH1}/"

mv "${gene}_${chrom}_${region}_2.bed" "${OUTPUT_PATH2}/"
mv "${gene}_${chrom}_${region}_2.bim" "${OUTPUT_PATH2}/"
mv "${gene}_${chrom}_${region}_2.fam" "${OUTPUT_PATH2}/"

mv "${gene}_${chrom}_${region}_3.bed" "${OUTPUT_PATH3}/"
mv "${gene}_${chrom}_${region}_3.bim" "${OUTPUT_PATH3}/"
mv "${gene}_${chrom}_${region}_3.fam" "${OUTPUT_PATH3}/"

echo "Pipeline completed for ${gene} at ${region_string}."
