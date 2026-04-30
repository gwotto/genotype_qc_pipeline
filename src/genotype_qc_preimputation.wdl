version 1.0

## =============================================================================
## genotype_qc_preimputation.wdl — Genotype QC Pre-imputation Pipeline
## =============================================================================
##
## PURPOSE
##   Performs standard pre-imputation QC on PLINK-format genotype data,
##   following the Anderson et al. (2010) protocol. Produces a clean binary
##   PLINK dataset suitable for association testing and imputation.
##
## PIPELINE STEPS (in order)
##   0.  [Optional] Convert text ped/map → binary bed/bim/fam
##   0b. Remove duplicate sample IDs
##   1.  SNP missingness filter        (--geno)
##   2.  Sample missingness filter     (--mind)
##   3.  Sex check & removal of fails  (skipped if no X-chromosome SNPs)
##   4.  MAF filter & monomorphic SNP removal
##   5.  Hardy-Weinberg equilibrium filter
##   6.  LD pruning (generates SNP list for steps 7 & 9)
##   7.  Heterozygosity outlier removal
##   8.  Restrict to autosomes (chr 1–22)
##   9.  Relatedness filtering (KING)
##   10. PCA (diagnostic, pre-imputation)
##   11. Prepare for imputation (TOPMed/HRC)
##
## EXECUTION
##   miniwdl run genotype_qc_preimputation.wdl -i inputs.json  # requires Docker
##   cromwell run genotype_qc_preimputation.wdl --inputs inputs.json
##
## DEPENDENCIES (must be on $PATH or supplied via *_bin inputs)
##   plink  >= 1.9
##   plink2 >= 2.0
##   Rscript >= 4.0
##
## REFERENCE
##   Anderson CA et al. (2010) Nat Protoc 5:1564–1573
##   https://doi.org/10.1038/nprot.2010.116
## =============================================================================


workflow genotype_qc_preimputation {

    input {
        # ── Input genotype files ───────────────────────────────────────────
        # Provide EITHER binary format (bed + bim + fam)
        # OR text format (ped + map). Binary is preferred.
        File? bed_file   # PLINK binary genotype matrix
        File? bim_file   # PLINK SNP information file
        File? fam_file   # PLINK sample information file

        File? ped_file   # PLINK text genotype file (alternative input)
        File? map_file   # PLINK SNP map file (alternative input)

        String output_prefix   # Prefix for all output files, e.g. "myStudy"

        # ── Reference and helper script files ─────────────────────────────
        File ld_regions              # High-LD genomic regions to exclude during pruning
                                     # Format: 4-column BED (chr start end name)
        File handle_duplicates_r     # R script: resolves duplicate sample IDs in .fam
        File check_heterozygosity_r  # R script: computes per-sample heterozygosity rates
        File heterozygosity_outliers_r  # R script: flags samples >3 SD from mean het rate
        File threshold_plot_r   # R script: plots SNP/sample counts across thresholds
        File pca_plot_r      # R script: plots PCA results
        File maf_plot_r      # R script: plots MAF distribution and barplot
    
        # ── Tool paths ─────────────────────────────────────────────────────
        # Change these if tools are not on $PATH
        String plink_bin    # PLINK 1.9 binary
        String plink2_bin   # PLINK 2.0 binary (used for KING relatedness)
        String rscript_bin  # R interpreter

        # ── QC thresholds ──────────────────────────────────────────────────
        # Typical values shown; adjust for your study.
        Float geno_threshold  # Max missing genotype rate per SNP    (e.g. 0.05 = 5%)
        Float mind_threshold  # Max missing genotype rate per sample (e.g. 0.05 = 5%)
        Float hwe_pvalue      # Min HWE p-value to retain SNP        (e.g. 1e-6)
        Float  maf_threshold   # Min MAF to retain a SNP (e.g. 0.01 = 1%)
        Float ld_r2           # Max r² for LD pruning                (e.g. 0.2)
        Int   ld_window_kb    # Sliding window size in kb            (e.g. 50)
        Int   ld_step         # Window step size in SNPs             (e.g. 5)
        Float pihat_min       # Min pi-hat to flag related pairs     (e.g. 0.2)
        Float king_cutoff     # KING kinship coefficient cutoff      (e.g. 0.0884 ≈ 3rd degree)

        # -- PCA ----------------------------------------------------------------
        Int    n_pcs     # Number of principal components to compute

        # -- Imputation preparation ----------------------------------------------
        File   check_bim_pl  # Will Rayner's HRC-1000G-check-bim.pl script
        File   hrc_ref_freq  # HRC r1.1 GRCh37 reference frequency file
        String perl_bin      # Perl interpreter
        String bcftools_bin  # bcftools binary
        String bgzip_bin   # bgzip binary
        String tabix_bin   # tabix binary

    }

    String pipeline_version = "1.0.0"

    # ── Step 0: Convert text format to binary (skipped if bed/bim/fam provided) ──
    # PLINK binary format is faster for all downstream steps.
    if (!defined(bed_file)) {
        call ConvertToBinary {
            input:
                ped_file      = select_first([ped_file]),
                map_file      = select_first([map_file]),
                output_prefix = output_prefix + "_QC0",
                plink_bin     = plink_bin
        }
    }

    # Resolve which files to carry forward (converted or original binary)
    File qc0_bed = select_first([ConvertToBinary.out_bed, bed_file])
    File qc0_bim = select_first([ConvertToBinary.out_bim, bim_file])
    File qc0_fam = select_first([ConvertToBinary.out_fam, fam_file])

    # ── Step 0b: Remove duplicate sample IDs ──────────────────────────────
    # Duplicate IDs cause downstream PLINK errors; the R script resolves them
    # by appending a suffix to make each ID unique.
    call HandleDuplicates {
        input:
            bed_file            = qc0_bed,
            bim_file            = qc0_bim,
            fam_file            = qc0_fam,
            handle_duplicates_r = handle_duplicates_r,
            output_prefix       = output_prefix + "_QC0",
            rscript_bin         = rscript_bin
    }

    call CountBimFam as LogStep0b {
        input:
            bim_file = HandleDuplicates.out_bim,
            fam_file = HandleDuplicates.out_fam,
            label    = "Step 0b  Remove duplicate IDs"
    }

    # ── Step 1: SNP missingness filter (--geno) ───────────────────────────
    # Removes SNPs with a call rate below (1 - geno_threshold).
    # Highly missing SNPs are likely to be poor-quality assays.
    call PlinkFilter as GenoFilter {
        input:
            bed_file      = HandleDuplicates.out_bed,
            bim_file      = HandleDuplicates.out_bim,
            fam_file      = HandleDuplicates.out_fam,
            plink_args    = "--geno " + geno_threshold,
            output_prefix = output_prefix + "_QC1",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep1 {
        input:
            bim_file = GenoFilter.out_bim,
            fam_file = GenoFilter.out_fam,
            label    = "Step 1   SNP missingness"
    }

    # ── Step 2: Sample missingness filter (--mind) ────────────────────────
    # Removes samples with a call rate below (1 - mind_threshold).
    # Highly missing samples indicate poor DNA quality or processing failure.
    call PlinkFilter as MindFilter {
        input:
            bed_file      = GenoFilter.out_bed,
            bim_file      = GenoFilter.out_bim,
            fam_file      = GenoFilter.out_fam,
            plink_args    = "--mind " + mind_threshold,
            output_prefix = output_prefix + "_QC2",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep2 {
        input:
            bim_file = MindFilter.out_bim,
            fam_file = MindFilter.out_fam,
            label    = "Step 2   Sample missingness"
    }

    # ── Threshold sweep: informational only, does not affect QC outputs ───
    # Runs PLINK across a range of --geno and --mind thresholds to show how
    # many SNPs/samples would be retained at each cutoff. Used to validate
    # the chosen thresholds. Results are plotted by an R script.
    call ThresholdSweep {
        input:
            bed_file         = HandleDuplicates.out_bed,
            bim_file         = HandleDuplicates.out_bim,
            fam_file         = HandleDuplicates.out_fam,
            geno_thresholds  = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2],
            mind_thresholds  = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2],
            geno_threshold   = geno_threshold,
            mind_threshold   = mind_threshold,
            output_prefix    = output_prefix + "_sweep",
            plink_bin        = plink_bin,
            rscript_bin      = rscript_bin,
            threshold_plot_r = threshold_plot_r
    }



    # ── Step 3: Sex check ─────────────────────────────────────────────────
    # Compares reported sex (from .fam) with X-chromosome heterozygosity.
    # Discordant samples may reflect sample swaps or genotyping errors.
    # Skipped entirely if the dataset contains no X-chromosome SNPs.

    call CountXSNPs {
        input:
            bim_file = MindFilter.out_bim
    }

    Boolean has_x_snps = CountXSNPs.n_x_snps > 0

    if (has_x_snps) {
        call SexCheck {
            input:
                bed_file      = MindFilter.out_bed,
                bim_file      = MindFilter.out_bim,
                fam_file      = MindFilter.out_fam,
                output_prefix = output_prefix + "_QC2",
                plink_bin     = plink_bin
        }

        # Remove samples that failed sex concordance check
        call RemoveSamples as RemoveSexFails {
            input:
                bed_file      = MindFilter.out_bed,
                bim_file      = MindFilter.out_bim,
                fam_file      = MindFilter.out_fam,
                remove_list   = SexCheck.problem_samples,
                output_prefix = output_prefix + "_QC3",
                plink_bin     = plink_bin
        }

        call CountBimFam as LogSexCheck {
            input:
                bim_file = RemoveSexFails.out_bim,
                fam_file = RemoveSexFails.out_fam,
                label    = "Step 3   Sex check"
        }
    }

    String sex_check_log_line = select_first([
        LogSexCheck.line,
        "Step 3   Sex check                            [skipped - no X SNPs]"
    ])


    # ── Step 4: MAF filter & monomorphic SNP removal ──────────────────────
    # Removes SNPs with MAF = 0 (monomorphic) and below maf_threshold.
    # Monomorphic SNPs carry no association signal and cause numerical issues.
    # Applied here — before HWE and LD pruning — so that:
    #   • HWE tests run only on common variants (reliable power)
    #   • LD pruning and heterozygosity checks use common variants only
    # Pre-imputation threshold: 0.01 recommended.
    call MafFilter {
        input:
            bed_file      = select_first([RemoveSexFails.out_bed, MindFilter.out_bed]),
            bim_file      = select_first([RemoveSexFails.out_bim, MindFilter.out_bim]),
            fam_file      = select_first([RemoveSexFails.out_fam, MindFilter.out_fam]),
            maf_threshold = maf_threshold,
            output_prefix = output_prefix + "_QC4",
            plink_bin     = plink_bin,
            rscript_bin   = rscript_bin,
            maf_plot_r    = maf_plot_r
    }

    call CountBimFam as LogStep4 {
        input:
            bim_file = MafFilter.out_bim,
            fam_file = MafFilter.out_fam,
            label    = "Step 4   MAF filter"
    }

    # ── Step 5: Hardy-Weinberg equilibrium filter ─────────────────────────
    # Removes SNPs that deviate significantly from HWE in controls.
    # Extreme HWE deviation often indicates genotyping error.
    # Applied after MAF filter (step 4): HWE tests are unreliable for rare
    # variants due to low expected counts at low frequencies.
    call PlinkFilter as HweFilter {
        input:
            bed_file      = MafFilter.out_bed,
            bim_file      = MafFilter.out_bim,
            fam_file      = MafFilter.out_fam,
            plink_args    = "--hwe " + hwe_pvalue,
            output_prefix = output_prefix + "_QC5",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep5 {
        input:
            bim_file = HweFilter.out_bim,
            fam_file = HweFilter.out_fam,
            label    = "Step 5   HWE filter"
    }


    # ── Step 6: LD pruning ────────────────────────────────────────────────
    # Generates a list of approximately independent SNPs by removing those
    # in high linkage disequilibrium. This pruned SNP set is used in the
    # heterozygosity check (step 7) and relatedness analysis (step 9) to
    # avoid inflated estimates from correlated markers.
    # High-LD genomic regions (e.g. MHC, inversions) are excluded first.
    call LdPruning {
        input:
            bed_file      = HweFilter.out_bed,
            bim_file      = HweFilter.out_bim,
            fam_file      = HweFilter.out_fam,
            ld_regions    = ld_regions,
            ld_window_kb  = ld_window_kb,
            ld_step       = ld_step,
            ld_r2         = ld_r2,
            output_prefix = "indepSNP",
            plink_bin     = plink_bin
    }

    # ── Step 7: Heterozygosity check ──────────────────────────────────────
    # Samples with unusually high or low heterozygosity are flagged.
    # High het → possible sample contamination.
    # Low het  → possible inbreeding or sample duplication.
    # Outliers are defined as >3 SD from the cohort mean.
    call HeterozygosityCheck {
        input:
            bed_file                  = HweFilter.out_bed,
            bim_file                  = HweFilter.out_bim,
            fam_file                  = HweFilter.out_fam,
            prune_in                  = LdPruning.prune_in,
            check_heterozygosity_r    = check_heterozygosity_r,
            heterozygosity_outliers_r = heterozygosity_outliers_r,
            plink_bin                 = plink_bin,
            rscript_bin               = rscript_bin
    }

    call RemoveSamples as RemoveHetFails {
        input:
            bed_file      = HweFilter.out_bed,
            bim_file      = HweFilter.out_bim,
            fam_file      = HweFilter.out_fam,
            remove_list   = HeterozygosityCheck.het_fail_ind,
            output_prefix = output_prefix + "_QC6",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep7 {
        input:
            bim_file = RemoveHetFails.out_bim,
            fam_file = RemoveHetFails.out_fam,
            label    = "Step 7   Heterozygosity filter"
    }

    # ── Step 8: Restrict to autosomes ─────────────────────────────────────
    # Retains only chromosomes 1–22 for downstream association analysis.
    # Sex chromosomes and mitochondrial SNPs require separate handling.
    call PlinkFilter as AutosomeFilter {
        input:
            bed_file      = RemoveHetFails.out_bed,
            bim_file      = RemoveHetFails.out_bim,
            fam_file      = RemoveHetFails.out_fam,
            plink_args    = "--chr 1-22",
            output_prefix = output_prefix + "_QC7",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep8 {
        input:
            bim_file = AutosomeFilter.out_bim,
            fam_file = AutosomeFilter.out_fam,
            label    = "Step 8   Autosome filter (chr 1-22)"
    }

    # ── Step 9: Relatedness filtering ─────────────────────────────────────
    # Identifies cryptically related samples using two methods:
    #   • pi-hat (PLINK --genome): proportion of IBD alleles shared
    #   • KING kinship coefficient (PLINK2 --king-cutoff): more robust in
    #     the presence of population stratification
    # Samples identified by KING as above the cutoff are removed.
    call RelatednessCheck {
        input:
            bed_file      = AutosomeFilter.out_bed,
            bim_file      = AutosomeFilter.out_bim,
            fam_file      = AutosomeFilter.out_fam,
            prune_in      = LdPruning.prune_in,
            pihat_min     = pihat_min,
            king_cutoff   = king_cutoff,
            output_prefix = output_prefix + "_QC8",
            plink_bin     = plink_bin,
            plink2_bin    = plink2_bin
    }

    call RemoveSamples as RemoveRelated {
        input:
            bed_file      = AutosomeFilter.out_bed,
            bim_file      = AutosomeFilter.out_bim,
            fam_file      = AutosomeFilter.out_fam,
            remove_list   = RelatednessCheck.king_cutoff_out,
            output_prefix = output_prefix + "_QC9",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep9 {
        input:
            bim_file = RemoveRelated.out_bim,
            fam_file = RemoveRelated.out_fam,
            label    = "Step 9   Relatedness filter (KING)"
    }

    # ── Step 10: PCA (pre-imputation, diagnostic only) ────────────────────
    # Computes principal components on the final QC-passed dataset using the
    # LD-pruned SNP list from step 6. Used to visualise population structure
    # and identify potential ancestry outliers BEFORE imputation.
    #
    # NOTE: These PCs are diagnostic only — no samples are removed here.
    #       PCA should be repeated after imputation to generate the PCs used
    #       as covariates in the GWAS model.
    call PCA {
        input:
            bed_file      = RemoveRelated.out_bed,
            bim_file      = RemoveRelated.out_bim,
            fam_file      = RemoveRelated.out_fam,
            prune_in      = LdPruning.prune_in,
            n_pcs         = n_pcs,
            output_prefix = output_prefix + "_QC9_pca",
            plink2_bin    = plink2_bin,
            rscript_bin   = rscript_bin,
            pca_plot_r    = pca_plot_r
    }

    # ── Step 11: Prepare for TOPMed imputation ────────────────────────────
    # Aligns strand orientation to the TOPMed reference panel using Will
    # Rayner's check-bim script. Removes SNPs not in the reference, ambiguous
    # A/T and C/G SNPs that cannot be strand-resolved, and SNPs with large
    # allele frequency differences vs the reference.
    # Outputs per-chromosome VCF files ready for upload to the TOPMed
    # imputation server (https://imputation.biodatacatalyst.nhlbi.nih.gov).
    call PrepareForImputation {
        input:
            bed_file        = RemoveRelated.out_bed,
            bim_file        = RemoveRelated.out_bim,
            fam_file        = RemoveRelated.out_fam,
            check_bim_pl    = check_bim_pl,
            hrc_ref_freq    = hrc_ref_freq,
            output_prefix   = output_prefix + "_imputation",
            plink_bin       = plink_bin,
            perl_bin        = perl_bin,
            bcftools_bin    = bcftools_bin,
            bgzip_bin       = bgzip_bin,
            tabix_bin       = tabix_bin
}


    # ── Pipeline log ──────────────────────────────────────────────────────
    # Collect per-step SNP/sample counts into a single human-readable log.
    Array[String] log_lines = [
        LogStep0b.line,
        LogStep1.line,
        LogStep2.line,
        sex_check_log_line,
        LogStep4.line,
        LogStep5.line,
        LdPruning.log_line,
        LogStep7.line,
        LogStep8.line,
        LogStep9.line,
        PrepareForImputation.log_line
    ]

    call WriteLog {
        input:
            pipeline_version = pipeline_version,
            output_prefix    = output_prefix,
            step_summaries   = log_lines
    }

    # ── Outputs ───────────────────────────────────────────────────────────
    output {
        # Pipeline run log — SNP/sample counts at each QC step
        File pipeline_log = WriteLog.log

        # Final QC-passed dataset — use these for association testing
        File final_bed = RemoveRelated.out_bed
        File final_bim = RemoveRelated.out_bim
        File final_fam = RemoveRelated.out_fam

        # Sex check outputs (only present if X-chromosome SNPs exist)
        File? sexcheck_report = SexCheck.sexcheck_report
        File? problem_samples = SexCheck.problem_samples

        # Heterozygosity check outputs
        File het_check_report = HeterozygosityCheck.r_check_het  # Per-sample het rates
        File het_fail_samples = HeterozygosityCheck.het_fail_ind  # Outlier sample list

        # MAF filter outputs (QC diagnostics; imputation freq computed separately)
        File maf_plot        = MafFilter.maf_plot
        File maf_freq_before = MafFilter.freq_before
        File maf_freq_after  = MafFilter.freq_after

        # Relatedness outputs
        File pihat_genome       = RelatednessCheck.pihat_genome       # All pairs above pihat_min
        File king_cutoff_out_id = RelatednessCheck.king_cutoff_out    # Samples removed by KING

        # LD-pruned SNP lists (used in steps 6 & 8; useful for PCA too)
        File prune_in  = LdPruning.prune_in
        File prune_out = LdPruning.prune_out

        # PCA outputs (diagnostic — inspect for population outliers)
        File pca_eigenvec  = PCA.eigenvec   # PC scores per sample
        File pca_eigenval  = PCA.eigenval   # Variance explained per PC
        File pca_plot      = PCA.pca_plot   # Scatter plots of PCs

        # Imputation-ready VCFs — upload these to the TOPMed server
        Array[File] imputation_vcfs     = PrepareForImputation.vcf_gz
        Array[File] imputation_vcf_tbis = PrepareForImputation.vcf_tbi
        File        check_bim_log       = PrepareForImputation.check_bim_log

    }
}


# =============================================================================
# TASKS
# =============================================================================
# Each task wraps a single command-line tool call. Tasks are reusable —
# PlinkFilter and RemoveSamples are called multiple times with different
# arguments. Inputs and outputs are explicit; no hidden state is shared
# between tasks.
# =============================================================================

## -----------------------------------------------------------------------------
## ConvertToBinary
## Convert PLINK text format (ped/map) to binary format (bed/bim/fam).
## Only called when bed/bim/fam files are not supplied as inputs.
## -----------------------------------------------------------------------------
task ConvertToBinary {
    input {
        File   ped_file       # Text genotype file
        File   map_file       # SNP position file
        String output_prefix
        String plink_bin
    }
    command <<<
        set -euo pipefail
        ~{plink_bin} \
            --ped ~{ped_file} \
            --map ~{map_file} \
            --make-bed \
            --out ~{output_prefix}
    >>>
    output {
        File out_bed = "~{output_prefix}.bed"
        File out_bim = "~{output_prefix}.bim"
        File out_fam = "~{output_prefix}.fam"
        File log     = "~{output_prefix}.log"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## HandleDuplicates
## Uses an R script to detect and resolve duplicate sample IDs in the .fam file.
## The bed/bim files are copied unchanged; only .fam IDs may be modified.
## -----------------------------------------------------------------------------
task HandleDuplicates {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   handle_duplicates_r   # R script that rewrites the .fam file
        String output_prefix
        String rscript_bin
    }

    command <<<
        set -euo pipefail
        cp ~{bed_file} ~{output_prefix}.bed
        cp ~{bim_file} ~{output_prefix}.bim
        cp ~{fam_file} ~{output_prefix}.fam
        ~{rscript_bin} --vanilla ~{handle_duplicates_r} ~{output_prefix}.fam
    >>>
    
    output {
        File out_bed = "~{output_prefix}.bed"
        File out_bim = "~{output_prefix}.bim"
        File out_fam = "~{output_prefix}.fam"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## PlinkFilter
## Generic PLINK 1.9 filter task. Accepts any valid PLINK flag(s) via
## plink_args and writes a new binary dataset.
## Used for: --geno, --mind, --hwe, --chr filtering.
## -----------------------------------------------------------------------------
task PlinkFilter {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        String plink_args    # Any PLINK filter flags, e.g. "--geno 0.05"
        String output_prefix
        String plink_bin
    }
    command <<<
        set -euo pipefail
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            ~{plink_args} \
            --make-bed \
            --out ~{output_prefix}
    >>>
    output {
        File out_bed = "~{output_prefix}.bed"
        File out_bim = "~{output_prefix}.bim"
        File out_fam = "~{output_prefix}.fam"
        File log     = "~{output_prefix}.log"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## ThresholdSweep
## Runs PLINK across a range of --geno and --mind thresholds and records how
## many SNPs and samples remain at each cutoff. An R script produces a ggplot
## showing retention curves, to help validate the chosen QC thresholds.
## This task is informational only — it does not affect any QC outputs.
## -----------------------------------------------------------------------------
task ThresholdSweep {
    input {
        File         bed_file
        File         bim_file
        File         fam_file
        Array[Float] geno_thresholds   # SNP missingness thresholds to test
        Array[Float] mind_thresholds   # Sample missingness thresholds to test
        Float        geno_threshold    # chosen threshold — passed to plot
        Float        mind_threshold    # chosen threshold — passed to plot
        String       output_prefix
        String       plink_bin
        String       rscript_bin
        File         threshold_plot_r
    }
    command <<<
        set -euo pipefail

        # Write header for results tables
        echo "threshold n_snps"   > geno_sweep.txt
        echo "threshold n_samples" > mind_sweep.txt

        # Sweep --geno thresholds: count remaining SNPs
        for T in ~{sep=' ' geno_thresholds}; do
            ~{plink_bin} \
                --bed ~{bed_file} \
                --bim ~{bim_file} \
                --fam ~{fam_file} \
                --geno $T \
                --make-bed \
                --out tmp_geno_$T > /dev/null 2>&1
            N=$(wc -l < tmp_geno_$T.bim)
            echo "$T $N" >> geno_sweep.txt
        done

        # Sweep --mind thresholds: count remaining samples
        for T in ~{sep=' ' mind_thresholds}; do
            ~{plink_bin} \
                --bed ~{bed_file} \
                --bim ~{bim_file} \
                --fam ~{fam_file} \
                --mind $T \
                --make-bed \
                --out tmp_mind_$T > /dev/null 2>&1
            N=$(wc -l < tmp_mind_$T.fam)
            echo "$T $N" >> mind_sweep.txt
        done
        ~{rscript_bin} --vanilla ~{threshold_plot_r} \
            ~{geno_threshold} \
            ~{mind_threshold}
    >>>
    output {
        File geno_sweep = "geno_sweep.txt"
        File mind_sweep = "mind_sweep.txt"
        File sweep_plot = "threshold_sweep.png"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## RemoveSamples
## Removes a list of samples from a PLINK dataset using --remove.
## The remove_list file must have two whitespace-separated columns: FID IID.
## Used for: sex check failures, heterozygosity outliers, related individuals.
## -----------------------------------------------------------------------------
task RemoveSamples {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   remove_list   # Two-column file: FID IID (one sample per line)
        String output_prefix
        String plink_bin
    }
    command <<<
        set -euo pipefail
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --remove ~{remove_list} \
            --make-bed \
            --out ~{output_prefix}
    >>>
    output {
        File out_bed = "~{output_prefix}.bed"
        File out_bim = "~{output_prefix}.bim"
        File out_fam = "~{output_prefix}.fam"
        File log     = "~{output_prefix}.log"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## CountXSNPs
## Counts the number of X-chromosome SNPs in a .bim file.
## Used as a guard: sex check is only meaningful when X SNPs are present.
## -----------------------------------------------------------------------------
task CountXSNPs {
    input {
        File bim_file   # PLINK .bim — chromosome is column 1
    }
    command <<<
        awk '$1=="X"{print}' ~{bim_file} | wc -l > x_count.txt
    >>>
    output {
        Int n_x_snps = read_int("x_count.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## SexCheck
## Runs PLINK --check-sex to compare reported sex with X-chromosome F-statistic.
## Writes a list of discordant samples (status != "OK") for removal.
## -----------------------------------------------------------------------------
task SexCheck {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        String output_prefix
        String plink_bin
    }
    command <<<
        set -euo pipefail
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --check-sex \
            --out ~{output_prefix}_sexcheck

        # Extract FID and IID of samples that failed (status column != "OK")
        awk '$5 != "OK" {print $1, $2}' ~{output_prefix}_sexcheck.sexcheck \
            > ~{output_prefix}_problem_samples.txt
    >>>
    output {
        File sexcheck_report = "~{output_prefix}_sexcheck.sexcheck"   # Full PLINK report
        File problem_samples = "~{output_prefix}_problem_samples.txt" # FID IID of failures
        File log             = "~{output_prefix}_sexcheck.log"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## LdPruning
## Generates two SNP lists using PLINK --indep-pairwise:
##   prune.in  — approximately independent SNPs (used downstream)
##   prune.out — correlated SNPs removed from consideration
## High-LD regions (e.g. MHC, chr8 inversion) are excluded before pruning
## to prevent these regions from dominating the independent SNP set.
## -----------------------------------------------------------------------------
task LdPruning {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   ld_regions     # Genomic regions to exclude (4-col BED format)
        Int    ld_window_kb   # Sliding window size in kilobases
        Int    ld_step        # Step size in number of SNPs
        Float  ld_r2          # r² threshold; SNP pairs above this are pruned
        String output_prefix
        String plink_bin
    }
    command <<<
        set -euo pipefail
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --exclude ~{ld_regions} --range \
            --indep-pairwise ~{ld_window_kb} kb ~{ld_step} ~{ld_r2} \
            --out ~{output_prefix}
        n_pruned=$(wc -l < ~{output_prefix}.prune.in)
        printf "%-42s  %7d independent SNPs\n" \
            "Step 6   LD pruning" "$n_pruned" > log_line.txt
    >>>
    output {
        File   prune_in  = "~{output_prefix}.prune.in"   # Independent SNPs
        File   prune_out = "~{output_prefix}.prune.out"  # Excluded SNPs
        File   log       = "~{output_prefix}.log"
        String log_line  = read_string("log_line.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## HeterozygosityCheck
## Computes per-sample heterozygosity rates on LD-pruned SNPs, then flags
## outliers using R scripts. Outlier threshold: mean ± 3 SD.
##   R_check.het     — PLINK heterozygosity output
##   het_fail_ind.txt — FID IID of outlier samples
## -----------------------------------------------------------------------------
task HeterozygosityCheck {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   prune_in                  # LD-pruned SNP list from LdPruning
        File   check_heterozygosity_r    # Computes observed heterozygosity rates
        File   heterozygosity_outliers_r # Flags samples outside mean ± 3 SD
        String plink_bin
        String rscript_bin
    }
    command <<<
        set -euo pipefail
        # Compute per-sample F-statistic and heterozygosity using pruned SNPs
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --extract ~{prune_in} \
            --het \
            --out R_check

        ~{rscript_bin} --vanilla ~{check_heterozygosity_r}
        ~{rscript_bin} --vanilla ~{heterozygosity_outliers_r}

        # Normalise output: strip quotes and keep FID + IID columns only
        sed 's/"//g' fail-het-qc.txt | awk '{print $1, $2}' > het_fail_ind.txt
    >>>
    output {
        File r_check_het  = "R_check.het"       # Full heterozygosity table
        File het_fail_ind = "het_fail_ind.txt"  # Outlier samples for removal
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## RelatednessCheck
## Identifies cryptically related sample pairs using two approaches:
##
##   1. PLINK 1.9 --genome  : computes pi-hat (proportion IBD) for all pairs
##      above pihat_min. Output is informational — not used for sample removal.
##
##   2. PLINK2 --king-cutoff : KING-robust kinship estimator. Samples above
##      the cutoff are written to king.cutoff.out.id and removed downstream.
##      Recommended cutoff: 0.0884 (3rd-degree relatives, ~12.5% IBD).
##
## Both analyses use the LD-pruned SNP list to reduce runtime and bias.
## -----------------------------------------------------------------------------
task RelatednessCheck {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   prune_in       # LD-pruned SNP list from LdPruning
        Float  pihat_min      # Min pi-hat to report a pair (e.g. 0.2)
        Float  king_cutoff    # Max kinship coefficient to retain a sample
        String output_prefix
        String plink_bin
        String plink2_bin
    }
    command <<<
        set -euo pipefail
        # Step 1: IBD estimation with PLINK 1.9
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --extract ~{prune_in} \
            --genome \
            --min ~{pihat_min} \
            --out ~{output_prefix}_pihat_min~{pihat_min}

        # Step 2: KING kinship with PLINK2 (robust to population stratification)
        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --extract ~{prune_in} \
            --king-cutoff ~{king_cutoff}
    >>>
    output {
        File pihat_genome    = "~{output_prefix}_pihat_min~{pihat_min}.genome"  # IBD pairs table
        File king_cutoff_out = "plink2.king.cutoff.out.id"  # Samples to REMOVE
        File king_cutoff_in  = "plink2.king.cutoff.in.id"   # Samples to KEEP
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## PCA
## Computes principal components using PLINK2 --pca on LD-pruned SNPs.
## The LD-pruned SNP list from LdPruning is reused here — no additional
## pruning is needed.
##
## Output:
##   .eigenvec — PC scores for each sample (used for plotting)
##   .eigenval — variance explained by each PC
##   pca_plot.png — scatter plots of PC1 vs PC2, PC1 vs PC3, PC2 vs PC3
##
## This task is diagnostic only. No samples are removed.
## Run PCA again after imputation to produce PCs for use as GWAS covariates.
## -----------------------------------------------------------------------------
task PCA {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   prune_in       # LD-pruned SNP list from LdPruning
        Int    n_pcs          # Number of PCs to compute (e.g. 10)
        String output_prefix
        String plink2_bin
        String rscript_bin
        File   pca_plot_r
    }
    command <<<
        set -euo pipefail

        # Extract LD-pruned SNPs and compute PCA
        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --extract ~{prune_in} \
            --pca ~{n_pcs} \
            --out ~{output_prefix}

        ~{rscript_bin} --vanilla ~{pca_plot_r} \
            ~{output_prefix}.eigenvec \
            ~{output_prefix}.eigenval \
            ~{n_pcs}
    >>>
    output {
        File eigenvec = "~{output_prefix}.eigenvec"
        File eigenval = "~{output_prefix}.eigenval"
        File pca_plot = "pca_plot.png"
        File log      = "~{output_prefix}.log"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## MafFilter
## Removes monomorphic SNPs (MAF = 0) and SNPs below maf_threshold.
## Steps:
##   1. Compute allele frequencies before filtering (--freq) → freq_before.frq
##   2. Remove monomorphic SNPs (--maf 0.0000001)
##   3. Remove SNPs below maf_threshold
##   4. Compute allele frequencies after filtering → freq_after.frq
##   5. R script plots MAF distribution and before/after barplot
## freq_before and freq_after are QC diagnostics only. The allele frequency
## file used for imputation prep is computed separately in PrepareForImputation
## on the final post-relatedness dataset.
## -----------------------------------------------------------------------------
task MafFilter {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        Float  maf_threshold
        String output_prefix
        String plink_bin
        String rscript_bin
        File   maf_plot_r
    }
    command <<<
        set -euo pipefail

        # Compute allele frequencies before any filtering
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --freq \
            --out freq_before

        # Remove monomorphic SNPs (MAF exactly 0)
        # is also implicitely done by the next step, so this seems redundant
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --maf 0.0000001 \
            --make-bed \
            --out tmp_no_monomorphic

        # Apply the chosen MAF threshold
        ~{plink_bin} \
            --bed tmp_no_monomorphic.bed \
            --bim tmp_no_monomorphic.bim \
            --fam tmp_no_monomorphic.fam \
            --maf ~{maf_threshold} \
            --make-bed \
            --out ~{output_prefix}

        # Compute allele frequencies after filtering (for plot)
        ~{plink_bin} \
            --bed ~{output_prefix}.bed \
            --bim ~{output_prefix}.bim \
            --fam ~{output_prefix}.fam \
            --freq \
            --out freq_after

        ~{rscript_bin} --vanilla ~{maf_plot_r} ~{maf_threshold}
    >>>
    output {
        File out_bed     = "~{output_prefix}.bed"
        File out_bim     = "~{output_prefix}.bim"
        File out_fam     = "~{output_prefix}.fam"
        File log         = "~{output_prefix}.log"
        File freq_before = "freq_before.frq"
        File freq_after  = "freq_after.frq"
        File maf_plot    = "maf_plot.png"
    }
    runtime { maxRetries: 1 }
}


## -----------------------------------------------------------------------------
## PrepareForImputation
## Aligns genotypes to the TOPMed reference panel and produces per-chromosome
## VCF files for upload to the TOPMed imputation server.
##
## Steps:
##   1. Compute allele frequencies on the final QC dataset (--freq)
##   2. Run HRC-1000G-check-bim.pl — produces Run-plink.sh with PLINK commands
##      that fix strand, remove problem SNPs, and split by chromosome
##   3. Execute the generated Run-plink.sh
##   4. Convert each per-chromosome PLINK dataset to bgzipped, tabix-indexed VCF
##
## The check-bim script removes:
##   • SNPs absent from the TOPMed reference panel
##   • Ambiguous SNPs (A/T and C/G) that cannot be strand-resolved
##   • SNPs with allele frequency difference > 0.2 vs reference
##   • SNPs with mismatched positions or alleles
##
## Reference:
##   Rayner W (2020) HRC or 1000G Imputation preparation and checking
##   https://www.well.ox.ac.uk/~wrayner/tools/
## -----------------------------------------------------------------------------
task PrepareForImputation {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   check_bim_pl      # Path to HRC-1000G-check-bim.pl
        File   hrc_ref_freq   # HRC reference frequency file (.tab.gz)
        String output_prefix
        String plink_bin
        String perl_bin
        String bcftools_bin
        String bgzip_bin
        String tabix_bin
    }
    command <<<
        set -euo pipefail

        # Decompress HRC reference file if needed
        # The check-bim script hardcodes /bin/gunzip which may not exist on macOS
        REF_FILE="~{hrc_ref_freq}"
        if [[ "$REF_FILE" == *.gz ]]; then
            gunzip -c "$REF_FILE" > hrc_ref.tab
            REF_FILE="hrc_ref.tab"
        fi

        # Step 1: Compute allele frequencies on the final QC-passed dataset
        # Computed here on the post-relatedness dataset so that frequencies
        # reflect the actual samples being submitted for imputation.
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --freq \
            --out final_freq

        # Step 2: Run check-bim script against TOPMed reference frequencies
        # Produces Run-plink.sh containing all fix/filter PLINK commands
        ~{perl_bin} ~{check_bim_pl} \
            -b ~{bim_file} \
            -f final_freq.frq \
            -r "$REF_FILE" \
            -h \
            -o ./

        ls -la

        # Rename the log to a predictable name for WDL output collection
        # check-bim names it based on the input bim stem
        # Rename log to predictable name — check-bim names it LOG-<stem>-HRC.txt
        mv $(ls LOG-*.txt | head -1) check-bim.log

        # the check-bim script's Run-plink.sh ends with `rm TEMP*`
        # which fails with set -e if no TEMP files exist.
        # Patch rm TEMP* to rm -f TEMP* — portable across macOS and Linux
        ~{perl_bin} -i -pe 's/^rm TEMP/rm -f TEMP/' Run-plink.sh

        # Step 3: Execute the generated PLINK commands
        # This performs strand flips, removes problem SNPs, and splits by chr
        bash Run-plink.sh

        # Step 4: Convert each per-chromosome PLINK file to VCF
        # TOPMed server requires one bgzipped, tabix-indexed VCF per chromosome
        for CHR in {1..22}; do
            PLINK_PREFIX=$(ls *-updated-chr${CHR}.bed 2>/dev/null \
                | sed 's/\.bed//' || true)

            # Skip chromosomes with no SNPs (ls returns nothing)
            if [ -z "$PLINK_PREFIX" ]; then
                continue
            fi

            # Sort by position (required by tabix) and convert to VCF
            # --real-ref-alleles preserves A1/A2 as ref/alt correctly
            ~{plink_bin} \
                --bfile "$PLINK_PREFIX" \
                --recode vcf \
                --real-ref-alleles \
                --out ~{output_prefix}_chr${CHR}

            # sort VCF by position (required by tabix)
            # bgzip the sorted VCF and index
            grep "^#" ~{output_prefix}_chr${CHR}.vcf > ~{output_prefix}_chr${CHR}_sorted.vcf
            grep -v "^#" ~{output_prefix}_chr${CHR}.vcf \
                | sort -k1,1V -k2,2n \
                >> ~{output_prefix}_chr${CHR}_sorted.vcf

            ~{bgzip_bin} ~{output_prefix}_chr${CHR}_sorted.vcf
            mv ~{output_prefix}_chr${CHR}_sorted.vcf.gz ~{output_prefix}_chr${CHR}.vcf.gz
            ~{tabix_bin} -p vcf ~{output_prefix}_chr${CHR}.vcf.gz

        done

        # Count total SNPs retained across all per-chromosome bim files
        # and format a log line matching the CountBimFam style
        n_snps=$(cat *-updated-chr*.bim 2>/dev/null | wc -l || echo 0)
        n_samples=$(wc -l < ~{fam_file})
        printf "%-42s  SNPs: %7d   Samples: %5d\n" \
            "Step 11  Imputation prep (check-bim)" "$n_snps" "$n_samples" \
            > step11_log.txt
    >>>
    output {
        Array[File] vcf_gz        = glob("~{output_prefix}_chr*.vcf.gz")
        Array[File] vcf_tbi       = glob("~{output_prefix}_chr*.vcf.gz.tbi")
        File        check_bim_log = "check-bim.log"
        String      log_line      = read_string("step11_log.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## CountBimFam
## Counts SNPs (.bim lines) and samples (.fam lines) in a PLINK binary dataset
## and formats a single-line summary for the pipeline log.
## -----------------------------------------------------------------------------
task CountBimFam {
    input {
        File   bim_file
        File   fam_file
        String label     # Step label, e.g. "Step 1   SNP missingness"
    }
    command <<<
        set -euo pipefail
        n_snps=$(wc -l < ~{bim_file})
        n_samples=$(wc -l < ~{fam_file})
        printf "%-42s  SNPs: %7d   Samples: %5d\n" \
            "~{label}" "$n_snps" "$n_samples" > line.txt
    >>>
    output {
        String line = read_string("line.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## WriteLog
## Assembles a human-readable pipeline run summary from per-step SNP/sample
## counts and writes it to a single log file.
## -----------------------------------------------------------------------------
task WriteLog {
    input {
        String        pipeline_version
        String        output_prefix
        Array[String] step_summaries
    }
    command <<<
        set -euo pipefail
        {
            printf "================================================================\n"
            printf " genotype_qc_preimputation  v%s\n" "~{pipeline_version}"
            printf " Run date : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
            printf " Prefix   : %s\n" "~{output_prefix}"
            printf "================================================================\n"
            printf "\n"
            cat ~{write_lines(step_summaries)}
            printf "\n"
            printf "================================================================\n"
        } > "~{output_prefix}_pipeline.log"
    >>>
    output {
        File log = "~{output_prefix}_pipeline.log"
    }
    runtime { maxRetries: 1 }
}
