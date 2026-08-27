### #!/usr/bin/env bash

# Exit immediately if a command fails (-e), treat unset variables as errors (-u),
# and propagate errors through pipes (-o pipefail).
set -euo pipefail


###############################################################################
# Configuration
###############################################################################

echo
echo "============================================================"
echo "        Genotype QC reference data preparation"
echo "============================================================"
echo

START_TIME=$(date +%s)

# Output/reference files
HIGH_LD_REGIONS="high-LD-regions-hg19-GRCh37.txt"
HRC_SITES="HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"

PHASE3_PGEN="all_phase3.pgen"
PHASE3_PVAR="all_phase3.pvar"
PHASE3_PSAM="all_phase3.psam"

# URLs
HIGH_LD_URL="https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt"
HRC_URL="ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
PGEN_URL="https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst"
PVAR_URL="https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst"
PSAM_URL="https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6"


###############################################################################
# 1. Download high-LD regions
###############################################################################

echo "[1/8] Downloading high-LD regions..."
echo "      File: ${HIGH_LD_REGIONS}"

curl -L --fail \
  -o "${HIGH_LD_REGIONS}" \
  "${HIGH_LD_URL}"

echo "      High-LD regions downloaded"
echo


###############################################################################
# 2. Download HRC reference sites
###############################################################################

echo "[2/8] Downloading HRC reference sites..."
echo "      File: ${HRC_SITES}"

curl -L --fail \
  -o "${HRC_SITES}" \
  "${HRC_URL}"

echo "      HRC sites downloaded"
echo


###############################################################################
# 3. Download 1000 Genomes Phase 3 PGEN file
###############################################################################

echo "[3/8] Downloading 1000 Genomes Phase 3 genotype data..."
echo "      Downloading: all_phase3.pgen.zst"

curl -L --fail \
  -o all_phase3.pgen.zst \
  "${PGEN_URL}"

echo "      Decompressing PGEN file..."

zstd -d --rm all_phase3.pgen.zst

echo "      PGEN file ready: ${PHASE3_PGEN}"
echo


###############################################################################
# 4. Download 1000 Genomes Phase 3 PVAR file
###############################################################################

echo "[4/8] Downloading 1000 Genomes Phase 3 variant information..."
echo "      Downloading: all_phase3.pvar.zst"

curl -L --fail \
  -o all_phase3.pvar.zst \
  "${PVAR_URL}"

echo "      Decompressing PVAR file..."

zstd -d --rm all_phase3.pvar.zst

echo "      PVAR file ready: ${PHASE3_PVAR}"
echo


###############################################################################
# 5. Download corrected sample information
###############################################################################

echo "[5/8] Downloading sample information..."
echo "      File: ${PHASE3_PSAM}"

curl -L --fail \
  -o all_phase3.psam \
  "${PSAM_URL}"

echo "      PSAM file downloaded"
echo


###############################################################################
# 6. Initial genotype filtering
###############################################################################

echo "[6/8] Performing initial genotype filtering..."
echo
echo "      Filters:"
echo "        - Missing genotype rate (--geno): 2%"
echo "        - Minor allele frequency (--maf): 5%"
echo "        - Biallelic variants only"
echo "        - SNPs only (A/C/G/T)"
echo

plink2 \
  --pfile all_phase3 \
  --geno 0.02 \
  --maf 0.05 \
  --max-alleles 2 \
  --snps-only just-acgt \
  --make-bed \
  --output-chr 26 \
  --out all_phase3_biallelic

echo
echo "      Initial filtering complete"
echo


###############################################################################
# 7. Sort variants and remove duplicate variant IDs
###############################################################################

echo "[7/8] Sorting variants and removing duplicate variant IDs..."
echo

# PLINK2 requires the split-chromosome correction to be performed
# as a standalone --make-pgen + --sort-vars operation.

echo "      Sorting variants..."

plink2 \
  --bfile all_phase3_biallelic \
  --make-pgen \
  --sort-vars \
  --out all_phase3_biallelic_sorted

echo "      Variants sorted"
echo

echo "      Removing duplicate variant IDs and assigning IDs to"
echo "      variants with missing IDs..."

# \$r and \$a are escaped so that Bash passes them literally to PLINK2.
# PLINK2 interprets them as reference/alternate alleles.

plink2 \
  --pfile all_phase3_biallelic_sorted \
  --make-bed \
  --new-id-max-allele-len 10 truncate \
  --rm-dup exclude-mismatch \
  --set-missing-var-ids @:#:\$r:\$a \
  --out all_phase3_nodup

echo
echo "      Duplicate handling complete"
echo


###############################################################################
# 8. LD pruning
###############################################################################

echo "[8/8] Performing LD pruning..."
echo
echo "      Parameters:"
echo "        - Exclude high-LD regions"
echo "        - Window: 200 variants"
echo "        - Step: 50 variants"
echo "        - r² threshold: 0.2"
echo

plink2 \
  --bfile all_phase3_nodup \
  --exclude range "${HIGH_LD_REGIONS}" \
  --indep-pairwise 200 50 0.2 \
  --out all_phase3_pruned

echo
echo "      LD pruning complete"
echo


###############################################################################
# Create final pruned dataset
###############################################################################

echo "Creating final pruned genotype dataset..."

plink2 \
  --bfile all_phase3_nodup \
  --extract all_phase3_pruned.prune.in \
  --make-bed \
  --out all_phase3_pruned_final \
  --output-chr 26

echo
echo "      Final pruned dataset created"
echo


###############################################################################
# Summary
###############################################################################

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "============================================================"
echo "                 Processing complete"
echo "============================================================"
echo
echo "Final genotype files:"
echo "  all_phase3_pruned_final.bed"
echo "  all_phase3_pruned_final.bim"
echo "  all_phase3_pruned_final.fam"
echo
echo "LD-pruning list:"
echo "  all_phase3_pruned.prune.in"
echo "  all_phase3_pruned.prune.out"
echo
echo
echo 
echo "NOTE: The following intermediate files are no longer required:"
echo "      Delete them if you no longer need them for"
echo "      additional analyses or to reproduce the pipeline."
echo
echo "  all_phase3_biallelic.bed"
echo "  all_phase3_biallelic.bim"
echo "  all_phase3_biallelic.fam"
echo "  all_phase3_biallelic_sorted.pgen"
echo "  all_phase3_biallelic_sorted.pvar"
echo "  all_phase3_biallelic_sorted.psam"
echo "  all_phase3_nodup.bed"
echo "  all_phase3_nodup.bim"
echo "  all_phase3_nodup.fam"
echo
echo
echo "Elapsed time: ${ELAPSED} seconds"
echo
echo "Data file processing complete"
echo
