#!/usr/bin/env bash
# =============================================================================
# collect_cromwell_results.sh
# 
# Collects and organises output files from a Cromwell genotype QC pipeline run.
# Results are written to a dated directory under the specified output directory.
#
# Usage:
#   bash collect_cromwell_results.sh <cromwell_run> <output_dir>
#
# Arguments:
#   cromwell_run   Full path to a specific Cromwell run directory (the UUID folder)
#   output_dir     Path where organised results will be written
#
# Output structure:
#   <output_dir>/genotype_qc_preimputation_YYYY-MM-DD/
#     qc_dataset/   Final QC'd PLINK files (bed/bim/fam)
#     plots/        PNG plots (MAC, heterozygosity, PCA)
#     reports/      QC reports (sex check, het, relatedness, freq, logs, ancestry)
#     pca/          PCA eigenvectors and eigenvalues
#     vcfs/         Imputation-ready VCFs per chromosome
#
# Example:
#   bash collect_cromwell_results.sh /data/cromwell-executions/myworkflow/3f2a1b4c-... ./results
#
# Notes:
#   - Sex check files are only present if X-chromosome SNPs were in the dataset
#   - Relatedness samples are flagged but NOT removed (user controls removal)
#   - Ancestry assignments file can be used with PLINK to filter by ancestry
# =============================================================================

set -euo pipefail

# -- Help message --------------------------------------------------------------
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 2 ]]; then
    sed -n '/^# ====/,/^# ====/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//';
    exit 0
fi


# -- Arguments -----------------------------------------------------------------
cromwell_run=$1
output_dir=$2

echo "Cromwell run directory : $cromwell_run"
echo "Output directory       : $output_dir"

BASE="$cromwell_run"

# -- Create output subdirectories ----------------------------------------------
OUTDIR="$output_dir"
mkdir -p "$OUTDIR"/{qc_dataset,plots,reports,vcfs,pca}
echo "Staging results in     : $OUTDIR"

# -- Pipeline log (WriteLog) ---------------------------------------------------
echo "Copying pipeline log..."
CALL="$BASE/call-WriteLog/execution"
cp "$CALL"/*_pipeline.log "$OUTDIR/"

# -- Final QC dataset (ChromosomeFilter; relatedness samples flagged but NOT removed) --
echo "Copying final QC dataset (PLINK bed/bim/fam)..."
CALL="$BASE/call-ChromosomeFilter/execution"
cp "$CALL"/*.bed "$CALL"/*.bim "$CALL"/*.fam "$OUTDIR/qc_dataset/"

# -- Threshold sweep plot (shows SNP/sample retention across thresholds) --------
echo "Copying threshold sweep plot..."
CALL="$BASE/call-ThresholdSweep/execution"
if [ -d "$CALL" ]; then
    cp "$CALL"/*.png "$OUTDIR/plots/" 2>/dev/null || true
fi

# -- Sex check report (step 3 - SexCheck / RemoveSexFails) --------------------
# Only present if the dataset contained X-chromosome SNPs
echo "Copying sex check reports (if present)..."
CALL="$BASE/call-SexCheck/execution"
[ -d "$CALL" ] && cp "$CALL"/*.sexcheck "$OUTDIR/reports/" 2>/dev/null || true
CALL="$BASE/call-RemoveSexFails/execution"
[ -d "$CALL" ] && cp "$CALL"/*.fam "$OUTDIR/reports/sex_fail_samples.fam" 2>/dev/null || true
CALL="$BASE/call-InitialVariantsPerChromosome/execution"
[ -d "$CALL" ] && cp "$CALL"/variants_per_chr_report.txt "$OUTDIR/reports/initial_variants_per_chr.txt" 2>/dev/null || true

# -- MAC plot & frequency files (step 2 - MacFilterSubset) --------------------------
echo "Copying MAC plots and frequency files..."
CALL="$BASE/call-MacFilter/execution"
cp "$CALL"/*.png "$OUTDIR/plots/"
cp "$CALL"/*.frq "$OUTDIR/reports/"

# -- HWE filter outputs --------------------------------------------------------
echo "Copying HWE filter reports (if present)..."
CALL="$BASE/call-HweFilterSubset/execution"
[ -d "$CALL" ] && cp "$CALL"/*.hwe "$OUTDIR/reports/" 2>/dev/null || true

# -- Heterozygosity report & fail list (step 7 - HetCheck / RemoveHetFails) ---
echo "Copying heterozygosity reports..."
CALL="$BASE/call-HeterozygosityCheck/execution"
cp "$CALL"/*.het "$OUTDIR/reports/"
cp "$CALL"/*.png "$OUTDIR/plots/" 2>/dev/null || true
cp "$CALL"/het_fail_ind.txt "$OUTDIR/reports/het_fail_samples.txt" 2>/dev/null || true

# -- Relatedness report (step 9 - RelatednessCheck) - skipped if run_relatedness_check=false
echo "Copying relatedness reports (if present)..."
CALL="$BASE/call-RelatednessCheck/execution"
if [ -d "$CALL" ]; then
    cp "$CALL"/*.genome "$OUTDIR/reports/" 2>/dev/null || true
    cp "$CALL"/*.king.cutoff.out.id "$OUTDIR/reports/relatedness_flagged_samples.txt" 2>/dev/null || true
fi

# -- LD-pruned SNP list used for ancestry PCA (step 3b - LdPruningPCA) --------
echo "Copying PCA LD-pruned SNP list..."
CALL="$BASE/call-LdPruningPCA/execution"
[ -d "$CALL" ] && cp "$CALL"/*.prune.in "$OUTDIR/pca/" 2>/dev/null || true

# -- Ancestry PCA against 1000 Genomes & ancestry assignments -----
echo "Copying ancestry PCA results and assignments..."
CALL="$BASE/call-AncestryPCA/execution"
if [ -d "$CALL" ]; then
    cp "$CALL"/*.eigenvec "$CALL"/*.eigenval "$OUTDIR/pca/" 2>/dev/null || true
    cp "$CALL"/*_ancestry_pca.png "$OUTDIR/plots/" 2>/dev/null || true
    for f in "$CALL"/*_ancestry_assignments.tsv; do
        [ -e "$f" ] || continue
        cp "$f" "$OUTDIR/reports/ancestry_assignments.tsv"
        cp "$f" "$OUTDIR/pca/"
    done
fi

# -- Imputation VCFs (step 11 - PrepareForImputation) -------------------------
echo "Copying imputation VCFs..."
CALL="$BASE/call-PrepareForImputation/execution"
cp "$CALL"/*_chr*.vcf.gz     "$OUTDIR/vcfs/"
cp "$CALL"/*_chr*.vcf.gz.tbi "$OUTDIR/vcfs/"
cp "$CALL"/check-bim.log     "$OUTDIR/reports/"

# -- Final variants-per-chromosome report (step 10b - FinalVariantsPerChromosome) --
CALL="$BASE/call-FinalVariantsPerChromosome/execution"
[ -d "$CALL" ] && cp "$CALL"/variants_per_chr_report.txt "$OUTDIR/reports/final_variants_per_chr.txt" 2>/dev/null || true

# -- Rename output directory with date ----------------------------------------
# outdir_renamed="${output_dir}_$(date +%Y-%m-%d)"
# echo "Renaming $OUTDIR -> $outdir_renamed"
# mv "$OUTDIR" "$outdir_renamed"

echo ""
echo "Done. Results written to: $output_dir"
echo ""
ls "$output_dir"/qc_dataset/ \
   "$output_dir"/plots/ \
   "$output_dir"/reports/ \
   "$output_dir"/pca/ \
   "$output_dir"/vcfs/
