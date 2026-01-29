# 🧬 Cardiomyopathy Variant Analysis Pipeline (AoU/gnomAD)

This repository contains the complete, modular pipeline for extracting, annotating, and filtering rare variants associated with cardiomyopathy (DCM, HCM, ARVC) from the **All of Us (AoU) Controlled Tier data**, integrated with gnomAD frequencies and functional prediction scores.

The workflow is strictly sequential, utilizing **Hail/Python**, **Bash/dsub**, and **R** to manage large-scale genomic data processing.

Please do cite this if you are planning to use this codes or results: Shetty NS, Pampana A, Gaonkar M, Patel N, Vekariya N, Smith JG, Kalra R, Chahal CAA, Semsarian C, Li P, Arora G, Arora P. Association of Pathogenic/Likely Pathogenic Genetic Variants for Cardiomyopathies With Clinical Outcomes: A Multiancestry Analysis in the All of Us Research Program. Circ Genom Precis Med. 2025 Jun;18(3):e005113. doi: 10.1161/CIRCGEN.124.005113. Epub 2025 May 28. PMID: 40433684.

---

## 📂 Repository Structure

The project is organized into three main directories to separate core logic from helper functions and data.

### 1. `scripts/`
Contains the **main execution pipeline** (Files 01-16). These are the drivers of the analysis and must be run in order.

### 2. `utilities/`
Contains auxiliary scripts and helper functions called by the main pipeline or used for pre-processing.

| Helper Script | Language | Purpose |
| :--- | :--- | :--- |
| **`plink_subsetter.sh`** | Bash | A utility script designed to slice large PLINK datasets into manageable subsets (by chromosome or region). It is essential for optimizing parallel processing steps. |
| **`vcf_normalization_wgs.sh`** | Bash | Standardizes WGS VCF files before annotation. It splits multiallelic sites and left-aligns indels (using `bcftools`) to ensure compatibility with SnpEff and other downstream tools. |

### 3. `data_files/`
Stores input datasets (e.g., raw VCFs, gene lists) and intermediate outputs. *Ensure this folder is configured according to security guidelines and contains no PII when pushing to version control.*

---

## ⚠️ ACTION REQUIRED: Security & Environment Setup

This repository contains **censored code** designed for use within controlled research environments (like the AoU Researcher Workbench).

1.  **Censored Paths:** All scripts contain placeholders (e.g., `<YOUR_SECURE_BUCKET>`, `<YOUR_PROJECT_ID>`). **You MUST replace these censored values** with your specific project credentials, GCS bucket paths, and resource links before execution.
2.  **No External Downloads:** The pipeline adheres to controlled environment policies. All required software and reference databases (Java, SnpEff, ClinVar, FASTA) must be pre-staged in your secure GCS bucket and are accessed via `gsutil cp`.

---

## ⚙️ Execution Workflow (Order of Operations)

The pipeline must be run sequentially. **Scripts are numbered 01 through 16.**

| Step | Script Filename | Language/Tool | Output/Purpose |
| :--- | :--- | :--- | :--- |
| **01** | `01_data_extraction.py` | Hail/Python | Extracts sites-only VCFs (WGS & Exome) for all 80 genes from the VDS/MT source files to GCS. |
| **02** | `02_dsub_plink_jobs.sh` | Bash/dsub/PLINK | Submits cloud jobs to filter large chromosomal PLINK files. **(Uses `plink_subsetter.sh`)** |
| **03** | `03_software_setup.sh` | Bash/gsutil | Downloads tool archives and reference files (FASTA, ClinVar) from GCS to local VM and sets environment paths. |
| **04** | `04_vcf_annotation.sh` | Bash/SnpEff | Annotates VCF files. **(Uses `vcf_normalization_wgs.sh` to prep data)** |
| **05** | `05_annotation_extraction.sh` | Bash/bcftools | Parses the annotated VCFs to extract specific annotation fields into TSV format for downstream processing. |
| **06** | `06_gnomad_annotations_extraction.py`| Python | Specifically targets and extracts data from the gnomAD database for population frequency comparisons. |
| **07** | `07_gnomad_annotations_processing.R` | R | Cleans the extracted gnomAD data and calculates the initial gnomAD PopMax AF lists. |
| **08** | `08_AoU_af_calculation_plink.sh` | Bash/PLINK | Calculates Allele Frequencies specifically for the All of Us (AoU) cohort using PLINK. |
| **09** | `09_AoU_popmax_af_calculation.R` | R | Calculates AoU PopMax AF and creates the master frequency filter list. |
| **10** | `10_af_aggregation_and_initial_filtration.R`| R | Aggregates calculated frequencies and performs a first-pass filtration based on frequency thresholds ($\le 0.001$). |
| **11** | `11_annotation_aggregation.R` | R | Merges the variant annotations with the frequency data into a unified dataset. |
| **12** | `12_gene_specific_filtering.R` | R | Applies advanced **functional filters** (e.g., TTN Hi-PSI, SCN5A REVEL/InterVar, ClinVar 2-star status). |
| **13** | `13_final_tiered_filtering.R` | R | Applies final disease-specific **FaF** and **Known Mechanism** filters, generating the final variant tiers. |
| **14** | `14_final_variant_export_and_cleanup.R` | R | Exports all final variant IDs into CSV files required for the final PLINK extraction steps. |
| **15** | `15_variant_counting_and_summary.R` | R | Generates summary tables of variant counts by gene and P/LP status. |
| **16** | `16_carrier_status_generation.R` | R/PLINK | Executes final PLINK extraction, processes genotypes, and determines individual carrier status. |

---

## 📜 Governance

### Code of Conduct

We are committed to providing a welcoming and inclusive environment for collaboration.

* **Be Respectful:** All communication should be constructive and respectful, regardless of differences in opinion or experience.
* **Be Professional:** Focus on the technical and research merits of the work. Harassment or discriminatory behavior will not be tolerated.
* **Be Clear:** Follow the project's documentation and processes when contributing code or reporting issues.

### LICENSE (MIT License Template)

The code in this repository is distributed under the MIT License.


