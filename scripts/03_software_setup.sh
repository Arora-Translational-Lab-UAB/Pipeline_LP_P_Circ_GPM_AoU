#!/bin/bash

# --- Script for Local VM Setup and Tool Installation ---
# This script copies required software and reference files from GCS and configures 
# environment variables locally, adhering to AOU policies against external downloads.

set -o pipefail
set -o errexit
set -u

# --- 1. Configuration (CENSORED PLACEHOLDERS) ---
# NOTE: Collaborator MUST set the Google Project ID and the path to the bucket 
# containing the tool archives and reference files.
GOOGLE_PROJECT="<YOUR_GOOGLE_PROJECT_ID>"
GCS_SOFTWARE_BUCKET="gs://<YOUR_SECURE_BUCKET>/Softwares"
GCS_CLINVAR_VCF="gs://<YOUR_SECURE_BUCKET>/Reference_Data/clinvar.vcf.gz"
GCS_CLINVAR_TBI="gs://<YOUR_SECURE_BUCKET>/Reference_Data/clinvar.vcf.gz.tbi"
GCS_FASTA="gs://<YOUR_SECURE_BUCKET>/Reference_Data/GRCh38.d1.vd1.fa"

# Local installation directory (assuming run from the root of the 'Annot' folder)
BASE_DIR="$(pwd)"
TOOL_INSTALL_DIR="${BASE_DIR}/tool_install"

mkdir -p "$TOOL_INSTALL_DIR"
cd "$TOOL_INSTALL_DIR"

echo "--- 2. Downloading Tools and References from GCS ---"

# 2.1 Download Tool Archives
echo "Copying software archives..."
gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_SOFTWARE_BUCKET}/openjdk-21+35_linux-x64_bin.tar.gz" .
gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_SOFTWARE_BUCKET}/bcftools.tar.gz" .
gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_SOFTWARE_BUCKET}/snpEff.tar.gz" .

# 2.2 Download Reference Files
echo "Copying ClinVar and FASTA reference files..."
gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_CLINVAR_VCF}" .
gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_CLINVAR_TBI}" .
gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_FASTA}" .


echo "--- 3. Tool Extraction and Configuration ---"

# 3.1 Java (OpenJDK) Setup
echo "Setting up OpenJDK..."
tar -xzf openjdk-21+35_linux-x64_bin.tar.gz
# Assuming extraction creates 'jdk-21' or similar; adjust path if necessary.
JAVA_HOME="${TOOL_INSTALL_DIR}/jdk-21"
export JAVA_HOME
export PATH="${JAVA_HOME}/bin:$PATH"
echo "JAVA_HOME set to: ${JAVA_HOME}"
java -version

# 3.2 SnpEff/SnpSift Setup (SnpSift is included in the SnpEff package)
echo "Setting up SnpEff/SnpSift..."
tar -xzf snpEff.tar.gz
SNPEFF_DIR="${TOOL_INSTALL_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_DIR}/SnpSift.jar"
# No explicit chmod is needed for JAR files, just ensure Java runs it.
java -jar "${SNPEFF_JAR}" -version

# 3.3 bcftools Setup
echo "Setting up bcftools..."
tar -xzf bcftools.tar.gz
# Assuming extraction creates a 'bcftools' folder and the binary is inside.
BCFTOOLS_EXE="${TOOL_INSTALL_DIR}/bcftools/bcftools" 
chmod +x "${BCFTOOLS_EXE}"
"${BCFTOOLS_EXE}" --version

echo "--- 4. Exporting Paths for Annotation Script Use ---"
# These variables can be sourced or reused by the next annotation script (e.g., 04_vcf_annotation.sh)
export PATH_JAVA="${JAVA_HOME}/bin/java"
export PATH_BCFTOOLS="${BCFTOOLS_EXE}"
export PATH_SNPEFF_JAR="${SNPEFF_JAR}"
export PATH_SNPSIFT_JAR="${SNPSIFT_JAR}"
export PATH_CLINVAR_GZ="${TOOL_INSTALL_DIR}/clinvar.vcf.gz"
export PATH_FASTA="${TOOL_INSTALL_DIR}/GRCh38.d1.vd1.fa"

echo "All annotation tools and references are installed and configured locally in: $TOOL_INSTALL_DIR"