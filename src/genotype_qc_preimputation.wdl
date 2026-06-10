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
##   2b. Report variants per chromosome
##   3.  Sex check & removal of fails  (skipped if no X-chromosome SNPs;
##       samples with unknown sex (sex = 0) are retained)
##   4.  Ancestry PCA against 1000 Genomes
##   5.  MAC filter & monomorphic SNP removal
##   6.  Hardy-Weinberg equilibrium filter
##   7.  LD pruning (generates SNP list for steps 7 & 9)
##   8.  Heterozygosity outlier removal
##   9.  Chromosome filter (configurable via chr_args, e.g. "--chr 1-22")
##   10. Relatedness filtering (KING)
##   10b.Report variants per chromosome (final, pre-imputation)
##   11. Prepare for imputation (TOPMed/HRC)
##
## CHROMOSOME ENCODING
##   Uses PLINK2 chromosome encoding throughout:
##     Autosomes:  1–22
##     X:          23
##     Y:          24
##     PAR (XY):   25
##     MT:         26
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
        # -- Input genotype files -------------------------------------------
        # Provide EITHER binary format (bed + bim + fam)
        # OR text format (ped + map). Binary is preferred.
        File? bed_file   # PLINK binary genotype matrix
        File? bim_file   # PLINK SNP information file
        File? fam_file   # PLINK sample information file

        File? ped_file   # PLINK text genotype file (alternative input)
        File? map_file   # PLINK SNP map file (alternative input)

        String output_prefix   # Prefix for all output files, e.g. "myStudy"

        # -- Reference and helper script files -----------------------------
        File ld_regions              # High-LD genomic regions to exclude during pruning
                                     # Format: 4-column BED (chr start end name)
        File handle_duplicates_r     # R script: resolves duplicate sample IDs in .fam
        File check_heterozygosity_r  # R script: computes per-sample heterozygosity rates
        File heterozygosity_outliers_r  # R script: flags samples >3 SD from mean het rate
        File threshold_plot_r   # R script: plots SNP/sample counts across thresholds
        File mac_plot_r      # R script: plots MAC distribution and barplot
    
        # -- Tool paths -----------------------------------------------------
        # Change these if tools are not on $PATH
        String plink_bin    # PLINK 1.9 binary
        String plink2_bin   # PLINK 2.0 binary (used for KING relatedness)
        String rscript_bin  # R interpreter

        # -- QC thresholds --------------------------------------------------
        # Typical values shown; adjust for your study.
        Float geno_threshold  # Max missing genotype rate per SNP    (e.g. 0.05 = 5%)
        Float mind_threshold  # Max missing genotype rate per sample (e.g. 0.05 = 5%)
        Float  hwe_pvalue      # Min HWE p-value to retain SNP        (e.g. 1e-6)
        Int    mac_threshold   # Min minor allele count to retain a SNP (e.g. 50)
        Float ld_r2           # Max r² for LD pruning                (e.g. 0.2)
        Int   ld_window_kb    # Sliding window size in kb            (e.g. 50)
        Int   ld_step         # Window step size in SNPs             (e.g. 5)
        Float pihat_min       # Min pi-hat to flag related pairs     (e.g. 0.2)
        Float king_cutoff     # KING kinship coefficient cutoff      (e.g. 0.0884 ≈ 3rd degree)

        # -- Chromosome filtering -----------------------------------------------
        String chr_args  # PLINK chromosome filter args (e.g. "--chr 1-22" or "--autosome")

        # -- PCA ----------------------------------------------------------------
        Int    n_pcs     # Number of principal components to compute

        # -- 1000 Genomes ancestry PCA --------------------------------------
        # Pre-processed 1000G Phase 3 (hg19) reference panel in PLINK binary
        # format. Should be LD-pruned and filtered to biallelic SNPs only
        # (see pipeline documentation for preparation steps).
        File   ref_1kg_bed       # 1000G reference panel .bed
        File   ref_1kg_bim       # 1000G reference panel .bim
        File   ref_1kg_fam       # 1000G reference panel .fam
        File   ref_1kg_psam      # 1000G sample metadata with SuperPop labels
        File   ancestry_pca_plot_r  # R script: plots study samples on 1000G PCA

        # -- Ancestry subset selection --------------------------------------
        # Comma-separated list of superpopulations to run variant/sample
        # detection on (e.g. "EUR,AFR"). Use "ALL" to run on entire cohort.
        String ancestry_populations = "ALL"

        # -- Imputation preparation ----------------------------------------------
        File   check_bim_pl  # Will Rayner's HRC-1000G-check-bim.pl script
        File   hrc_ref_freq  # HRC r1.1 GRCh37 reference frequency file
        String perl_bin      # Perl interpreter
        String bcftools_bin  # bcftools binary
        String bgzip_bin   # bgzip binary
        String tabix_bin   # tabix binary

    }

    String pipeline_version = "1.0.1"

    # -- Step 0: Convert text format to binary (skipped if bed/bim/fam provided) --
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

    # -- Step 0b: Remove duplicate sample IDs ------------------------------
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

    # -- Step 1: SNP missingness filter (--geno) ---------------------------
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

    # -- Step 2: Sample missingness filter (--mind) ------------------------
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

    # Report variants per chromosome before the sex-check step
    call VariantsPerChromosome as InitialVariantsPerChromosome {
        input:
            bim_file = MindFilter.out_bim,
            label    = "Step 2b  Variants per chromosome"
    }

    # -- Threshold sweep: informational only, does not affect QC outputs ---
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



    # -- Step 3: Sex check -------------------------------------------------
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

    # -- Step 4: Ancestry PCA against 1000 Genomes -----------------------
    # Merges the post-sex-check study data with the 1000G Phase 3 reference
    # panel and computes PCA jointly. Used to assign superpopulation ancestry
    # labels to study samples and identify ancestry outliers.
    # Input data must be on hg19/GRCh37 to match the 1000G reference.
    call AncestryPCA {
        input:
            bed_file           = select_first([RemoveSexFails.out_bed, MindFilter.out_bed]),
            bim_file           = select_first([RemoveSexFails.out_bim, MindFilter.out_bim]),
            fam_file           = select_first([RemoveSexFails.out_fam, MindFilter.out_fam]),
            ref_bed            = ref_1kg_bed,
            ref_bim            = ref_1kg_bim,
            ref_fam            = ref_1kg_fam,
            ref_psam           = ref_1kg_psam,
            n_pcs              = n_pcs,
            output_prefix      = output_prefix + "_ancestry_pca",
            plink_bin          = plink_bin,
            plink2_bin         = plink2_bin,
            rscript_bin        = rscript_bin,
            ancestry_pca_plot_r = ancestry_pca_plot_r
    }

    # -- Steps 5–8: For each QC step, detect problematic variants on the
    #               ancestry-defined subset and remove those variants from
    #               the full cohort before proceeding to the next step.
    # Note: LD pruning here is only used as a helper for heterozygosity
    #       detection (to generate a prune_in list). A separate pruning
    #       step is used for relatedness and optional final PCA.

    # 1) Build ancestry keep-list (used repeatedly)
    call MakeAncestryKeepList {
        input:
            assignments = AncestryPCA.assignments,
            fam_file    = select_first([RemoveSexFails.out_fam, MindFilter.out_fam]),
            populations = ancestry_populations,
            output_prefix = output_prefix
    }

    # --- Step 5: MAC detection on subset, removal on full cohort ----------
    call PlinkFilter as SubsetAncestry_MAC {
        input:
            bed_file      = select_first([RemoveSexFails.out_bed, MindFilter.out_bed]),
            bim_file      = select_first([RemoveSexFails.out_bim, MindFilter.out_bim]),
            fam_file      = select_first([RemoveSexFails.out_fam, MindFilter.out_fam]),
            plink_args    = "--keep " + MakeAncestryKeepList.keep_list,
            output_prefix = output_prefix + "_ancestry_subset",
            plink_bin     = plink_bin
    }

    call MacFilter as MacFilterSubset {
        input:
            bed_file      = SubsetAncestry_MAC.out_bed,
            bim_file      = SubsetAncestry_MAC.out_bim,
            fam_file      = SubsetAncestry_MAC.out_fam,
            mac_threshold = mac_threshold,
            output_prefix = output_prefix + "_QC5_subset",
            plink_bin     = plink_bin,
            rscript_bin   = rscript_bin,
            mac_plot_r    = mac_plot_r
    }

    call SnpListDiff as MacRemoved {
        input:
            before_bim = SubsetAncestry_MAC.out_bim,
            after_bim  = MacFilterSubset.out_bim,
            output_prefix = output_prefix + "_mac"
    }

    call PlinkFilter as RemoveVariantsAfterMAC {
        input:
            bed_file      = select_first([RemoveSexFails.out_bed, MindFilter.out_bed]),
            bim_file      = select_first([RemoveSexFails.out_bim, MindFilter.out_bim]),
            fam_file      = select_first([RemoveSexFails.out_fam, MindFilter.out_fam]),
            plink_args    = "--exclude " + MacRemoved.removed_snps,
            output_prefix = output_prefix + "_QC5",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep5 {
        input:
            bim_file = RemoveVariantsAfterMAC.out_bim,
            fam_file = RemoveVariantsAfterMAC.out_fam,
            label    = "Step 5   MAC filter (subset detection)"
    }

    # --- Step 6: HWE detection on subset (after MAC removals), removal on full
    call PlinkFilter as SubsetAncestry_HWE {
        input:
            bed_file      = RemoveVariantsAfterMAC.out_bed,
            bim_file      = RemoveVariantsAfterMAC.out_bim,
            fam_file      = RemoveVariantsAfterMAC.out_fam,
            plink_args    = "--keep " + MakeAncestryKeepList.keep_list,
            output_prefix = output_prefix + "_ancestry_subset_QC5",
            plink_bin     = plink_bin
    }

    call PlinkFilter as HweFilterSubset {
        input:
            bed_file      = SubsetAncestry_HWE.out_bed,
            bim_file      = SubsetAncestry_HWE.out_bim,
            fam_file      = SubsetAncestry_HWE.out_fam,
            plink_args    = "--hwe " + hwe_pvalue,
            output_prefix = output_prefix + "_QC6_subset",
            plink_bin     = plink_bin
    }

    call SnpListDiff as HweRemoved {
        input:
            before_bim = SubsetAncestry_HWE.out_bim,
            after_bim  = HweFilterSubset.out_bim,
            output_prefix = output_prefix + "_hwe"
    }

    call PlinkFilter as RemoveVariantsAfterHWE {
        input:
            bed_file      = RemoveVariantsAfterMAC.out_bed,
            bim_file      = RemoveVariantsAfterMAC.out_bim,
            fam_file      = RemoveVariantsAfterMAC.out_fam,
            plink_args    = "--exclude " + HweRemoved.removed_snps,
            output_prefix = output_prefix + "_QC6",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep6 {
        input:
            bim_file = RemoveVariantsAfterHWE.out_bim,
            fam_file = RemoveVariantsAfterHWE.out_fam,
            label    = "Step 6   HWE filter (subset detection)"
    }

    # --- Step 7: LD pruning only as helper for heterozygosity (do not remove LD SNPs)

    call LdPruning as LdPruningHet {
        input:
            bed_file      = RemoveVariantsAfterHWE.out_bed,
            bim_file      = RemoveVariantsAfterHWE.out_bim,
            fam_file      = RemoveVariantsAfterHWE.out_fam,
            ld_regions    = ld_regions,
            ld_window_kb  = ld_window_kb,
            ld_step       = ld_step,
            ld_r2         = ld_r2,
            output_prefix = output_prefix + "_indepSNP_het",
            plink_bin     = plink_bin
    }

    call HeterozygosityCheck as HeterozygosityCheck {
        input:
            bed_file                  = RemoveVariantsAfterHWE.out_bed,
            bim_file                  = RemoveVariantsAfterHWE.out_bim,
            fam_file                  = RemoveVariantsAfterHWE.out_fam,
            prune_in                  = LdPruningHet.prune_in,
            check_heterozygosity_r    = check_heterozygosity_r,
            heterozygosity_outliers_r = heterozygosity_outliers_r,
            plink_bin                 = plink_bin,
            rscript_bin               = rscript_bin
    }

    call RemoveSamples as RemoveHetFails {
        input:
            bed_file      = RemoveVariantsAfterHWE.out_bed,
            bim_file      = RemoveVariantsAfterHWE.out_bim,
            fam_file      = RemoveVariantsAfterHWE.out_fam,
            remove_list   = HeterozygosityCheck.het_fail_ind,
            output_prefix = output_prefix + "_QC8",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep8 {
        input:
            bim_file = RemoveHetFails.out_bim,
            fam_file = RemoveHetFails.out_fam,
            label    = "Step 8   Heterozygosity filter"
    }

    # -- Step 9: Restrict to selected chromosomes -------------------------------------
    # Retains only chromosomes 1–22 for downstream association analysis.
    # Sex chromosomes and mitochondrial SNPs require separate handling.
    # Chromosome filtering args are configurable (e.g. "--chr 1-23" or "--chr 1-10,12-22").
    call PlinkFilter as ChromosomeFilter {
        input:
            bed_file      = RemoveHetFails.out_bed,
            bim_file      = RemoveHetFails.out_bim,
            fam_file      = RemoveHetFails.out_fam,
            plink_args    = chr_args,
            output_prefix = output_prefix + "_QC9",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep9 {
        input:
            bim_file = ChromosomeFilter.out_bim,
            fam_file = ChromosomeFilter.out_fam,
            label    = "Step 9   Chromosome filter"
    }

    # -- Step 10: Relatedness filtering -------------------------------------
    # Identifies cryptically related samples using two methods:
    #   • pi-hat (PLINK --genome): proportion of IBD alleles shared
    #   • KING kinship coefficient (PLINK2 --king-cutoff): more robust in
    #     the presence of population stratification
    # Samples identified by KING as above the cutoff are flagged (not automatically removed).
    
    # Compute LD-pruned SNP list for relatedness check and optional final PCA.

    call LdPruning as LdPruningRelatedness {
        input:
            bed_file      = ChromosomeFilter.out_bed,
            bim_file      = ChromosomeFilter.out_bim,
            fam_file      = ChromosomeFilter.out_fam,
        ld_regions    = ld_regions,
        ld_window_kb  = ld_window_kb,
        ld_step       = ld_step,
        ld_r2         = ld_r2,
        output_prefix = output_prefix + "_indepSNP_related",
        plink_bin     = plink_bin
    }

    call RelatednessCheck {
    input:
        bed_file      = ChromosomeFilter.out_bed,
        bim_file      = ChromosomeFilter.out_bim,
        fam_file      = ChromosomeFilter.out_fam,
        prune_in      = LdPruningRelatedness.prune_in,
        pihat_min     = pihat_min,
        king_cutoff   = king_cutoff,
        output_prefix = output_prefix + "_QC10",
        plink_bin     = plink_bin,
        plink2_bin    = plink2_bin
    }


    call CountBimFam as LogStep10 {
    input:
        bim_file = ChromosomeFilter.out_bim,
        fam_file = ChromosomeFilter.out_fam,
        label    = "Step 10   Relatedness filter (KING)"
    }

String relatedness_log_line = LogStep10.line


    # -- Step 11: Prepare for TOPMed imputation ----------------------------
    # Aligns strand orientation to the TOPMed reference panel using Will
    # Rayner's check-bim script. Removes SNPs not in the reference, ambiguous
    # A/T and C/G SNPs that cannot be strand-resolved, and SNPs with large
    # allele frequency differences vs the reference.
    # Outputs per-chromosome VCF files ready for upload to the TOPMed
    # imputation server (https://imputation.biodatacatalyst.nhlbi.nih.gov).
    
    # Report final per-chromosome variant counts before imputation prep
    call VariantsPerChromosome as FinalVariantsPerChromosome {
        input:
            bim_file = select_first([ChromosomeFilter.out_bim, ChromosomeFilter.out_bim]),
            label    = "Step 10b Variants per chromosome (pre-imputation)"
    }

    call PrepareForImputation {
        input:
            bed_file        = select_first([ChromosomeFilter.out_bed, ChromosomeFilter.out_bed]),
            bim_file        = select_first([ChromosomeFilter.out_bim, ChromosomeFilter.out_bim]),
            fam_file        = select_first([ChromosomeFilter.out_fam, ChromosomeFilter.out_fam]),
            check_bim_pl    = check_bim_pl,
            hrc_ref_freq    = hrc_ref_freq,
            output_prefix   = output_prefix + "_imputation",
            plink_bin       = plink_bin,
            perl_bin        = perl_bin,
            bcftools_bin    = bcftools_bin,
            bgzip_bin       = bgzip_bin,
            tabix_bin       = tabix_bin
}


    # -- Pipeline log ------------------------------------------------------
    # Collect per-step SNP/sample counts into a single human-readable log.
    Array[String] log_lines = flatten([
        [
            LogStep0b.line,
            LogStep1.line,
            LogStep2.line,
            InitialVariantsPerChromosome.log_line
        ],
        InitialVariantsPerChromosome.log_lines,
        [
            sex_check_log_line,
            AncestryPCA.log_line,
            LogStep5.line,
            LogStep6.line,
            LdPruningHet.log_line,
            LogStep8.line,
            LogStep9.line,
            relatedness_log_line,
            FinalVariantsPerChromosome.log_line
        ],
        FinalVariantsPerChromosome.log_lines,
        [
            PrepareForImputation.log_line
        ]
    ])

    call WriteLog {
        input:
            pipeline_version = pipeline_version,
            output_prefix    = output_prefix,
            step_summaries   = log_lines
    }

    # -- Outputs -----------------------------------------------------------
    output {
        # Pipeline run log — SNP/sample counts at each QC step
        File pipeline_log = WriteLog.log

        # Final QC-passed dataset — use these for association testing
        File final_bed = select_first([ChromosomeFilter.out_bed, ChromosomeFilter.out_bed])
        File final_bim = select_first([ChromosomeFilter.out_bim, ChromosomeFilter.out_bim])
        File final_fam = select_first([ChromosomeFilter.out_fam, ChromosomeFilter.out_fam])

        # Sex check outputs (only present if X-chromosome SNPs exist)
        File? sexcheck_report = SexCheck.sexcheck_report
        File? problem_samples = SexCheck.problem_samples
        File? variants_per_chr_initial = InitialVariantsPerChromosome.report

        # Heterozygosity check outputs
        File het_check_report = HeterozygosityCheck.r_check_het  # Per-sample het rates
        File het_fail_samples = HeterozygosityCheck.het_fail_ind  # Outlier sample list

        # MAC filter outputs (QC diagnostics; imputation freq computed separately)
        File mac_plot        = MacFilterSubset.mac_plot
        File maf_freq_before = MacFilterSubset.freq_before
        File maf_freq_after  = MacFilterSubset.freq_after

        # Relatedness outputs (relatedness filter is optional)
        File? pihat_genome       = RelatednessCheck.pihat_genome       # All pairs above pihat_min
        File? king_cutoff_out_id = RelatednessCheck.king_cutoff_out    # Samples removed by KING

        # LD-pruned SNP lists (used in steps 8 & 10; useful for PCA too)
        File prune_in  = LdPruningHet.prune_in
        File prune_out = LdPruningHet.prune_out
        File? related_prune_in = LdPruningRelatedness.prune_in

        # Ancestry PCA against 1000 Genomes
        File ancestry_pca_eigenvec     = AncestryPCA.eigenvec       # PC scores (study + 1000G)
        File ancestry_pca_plot         = AncestryPCA.pca_plot        # Study samples overlaid on 1000G
        File ancestry_assignments      = AncestryPCA.assignments     # Per-sample superpopulation assignments

        # Per-chromosome variant counts (diagnostic — validates chromosome coding)
        File variants_per_chr_final = FinalVariantsPerChromosome.report

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
    ## MakeAncestryKeepList
    ## Creates a PLINK keep-list (FID IID) for study samples matching the user-
    ## supplied comma-separated superpopulation list, or all samples if "ALL".
    ## - assignments: TSV produced by AncestryPCA with at least an IID and a
    ##   predicted superpopulation column (header names: superpop|predicted_pop|predicted)
    ## - fam_file: study .fam to map IID → FID when assignments contain only IID
    ## -----------------------------------------------------------------------------
    task MakeAncestryKeepList {
        input {
            File assignments
            File fam_file
            String populations
            String output_prefix
        }
        command <<<
            set -euo pipefail

            POPS="~{populations}"
            OUT="~{output_prefix}_ancestry_keep.txt"

            if [ "$POPS" = "ALL" ]; then
                awk '{print $1, $2}' ~{fam_file} > "$OUT"
            else
                # Normalize comma/space separation
                PATT=$(echo "$POPS" | tr -d ' ')

                # Determine header column for population labels and extract matching IIDs
                awk -F"\t" -v patt="$PATT" '
                NR==1 { for(i=1;i<=NF;i++) hdr[$i]=i;
                          if ("superpop" in hdr) c=hdr["superpop"]; 
                          else if ("predicted_pop" in hdr) c=hdr["predicted_pop"]; 
                          else if ("predicted" in hdr) c=hdr["predicted"]; 
                          if ("IID" in hdr) iidc=hdr["IID"]; else iidc=2;
                          next }
                NR>1 {
                    n=split(patt, a, ",");
                    for (i=1;i<=n;i++) if ($c==a[i]) print $iidc
                }' ~{assignments} > selected_iids.txt

                # Map IIDs back to real FIDs using the provided .fam (IID -> FID)
                awk 'NR==FNR{iids[$0]=1; next} ($2 in iids){print $1, $2}' selected_iids.txt ~{fam_file} > "$OUT" || true
                rm -f selected_iids.txt || true
            fi

            wc -l < "$OUT" > keep_count.txt
        >>>
        output {
            File keep_list = "~{output_prefix}_ancestry_keep.txt"
            Int n_keep = read_int("keep_count.txt")
        }
        runtime { maxRetries: 1 }
    }

    ## -----------------------------------------------------------------------------
    ## SnpListDiff
    ## Given two .bim files (before and after filtering) write the SNP IDs present
    ## in the 'before' file but not in the 'after' file (one rsID per line).
    ## -----------------------------------------------------------------------------
    task SnpListDiff {
        input {
            File before_bim
            File after_bim
            String output_prefix
        }
        command <<<
            set -euo pipefail
            awk '{print $2}' ~{before_bim} | sort > before_snps.txt
            awk '{print $2}' ~{after_bim}  | sort > after_snps.txt
            comm -23 before_snps.txt after_snps.txt > ~{output_prefix}_removed_snps.txt
        >>>
        output {
            File removed_snps = "~{output_prefix}_removed_snps.txt"
        }
        runtime { maxRetries: 1 }
    }

    ## -----------------------------------------------------------------------------
    ## MergeSnpLists
    ## Concatenate an array of SNP-list files and output a single unique, sorted
    ## SNP exclude file suitable for PLINK --exclude.
    ## -----------------------------------------------------------------------------
    task MergeSnpLists {
        input {
            Array[File] snp_lists
            String output_prefix
        }
        command <<<
            set -euo pipefail
            > combined.tmp
            for f in ~{sep=' ' snp_lists}; do
                cat "$f" >> combined.tmp
            done
            sort -u combined.tmp > ~{output_prefix}_remove_snps_union.txt
        >>>
        output {
            File union_snps = "~{output_prefix}_remove_snps_union.txt"
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
## Compatible with both chrX and chr23 chromosome coding formats.
## -----------------------------------------------------------------------------
task CountXSNPs {
    input {
        File bim_file   # PLINK .bim — chromosome is column 1
    }
    command <<<
        # Count X-chromosome SNPs in both old (X) and new (23) notation
        x_old=$(awk '$1=="X"{print}' ~{bim_file} | wc -l)
        x_new=$(awk '$1=="23"{print}' ~{bim_file} | wc -l)
        echo $((x_old + x_new)) > x_count.txt
    >>>
    output {
        Int n_x_snps = read_int("x_count.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## VariantsPerChromosome
## Counts and reports the number of variants per chromosome.
## Recognizes both PLINK1 (1-22, X, Y, MT) and PLINK2 (1-22, 23, 24, 25, 26) 
## chromosome coding schemes.
##
## Output includes:
##   - A formatted report of variants per chromosome (chromosomes with 0
##     variants are omitted from the report)
##   - Separate counts for autosomes, sex chromosomes, and mitochondrial DNA
##   - A single-line summary for the pipeline log
##   - A multi-line breakdown for inclusion in the pipeline log
## -----------------------------------------------------------------------------
task VariantsPerChromosome {
    input {
        File   bim_file
        String label    # Step label for the log line
    }
    command <<<
        set -euo pipefail
        
        # Initialize output files
        > variants_per_chr_report.txt
        > log_lines.txt
        
        # Count and report variants per chromosome (1-22), skip zeros
        total_autosomes=0
        # Derive chromosome list from the bim file so we handle whatever
        # chromosomes were selected in the config (not just 1-22).
        # Exclude non-autosomal chromosomes here; they are handled separately below.
        CHROMS=$(awk '{print $1}' ~{bim_file} | sort -u -V | grep -vE '^(X|Y|MT|XY|23|24|25|26)$')

        for CHR in $CHROMS; do
            count=$(awk -v chr="$CHR" '$1==chr' ~{bim_file} | wc -l)
            if [ "$count" -gt 0 ]; then
                label=$(printf "Chromosome %2d:" "$CHR")
                printf "  %-20s %8d variants\n" "$label" "$count" >> variants_per_chr_report.txt
                printf "  %-20s %8d variants\n" "$label" "$count" >> log_lines.txt
            fi
            total_autosomes=$((total_autosomes + count))
        done
        
        # Count sex chromosomes (both PLINK1 and PLINK2 notations, combined)
        x_plink1=$(awk '$1=="X"' ~{bim_file} | wc -l)
        x_plink2=$(awk '$1=="23"' ~{bim_file} | wc -l)
        x_total=$((x_plink1 + x_plink2))
        
        y_plink1=$(awk '$1=="Y"' ~{bim_file} | wc -l)
        y_plink2=$(awk '$1=="24"' ~{bim_file} | wc -l)
        y_total=$((y_plink1 + y_plink2))
        
        # Pseudoautosomal region (PLINK2 only)
        xy_plink2=$(awk '$1=="25"' ~{bim_file} | wc -l)
        
        # Mitochondrial DNA
        mt_plink1=$(awk '$1=="MT"' ~{bim_file} | wc -l)
        mt_plink2=$(awk '$1=="26"' ~{bim_file} | wc -l)
        mt_total=$((mt_plink1 + mt_plink2))
        
        # Only report sex chromosomes and related categories if non-zero
        if [ "$x_total" -gt 0 ]; then
            printf "  %-20s %8d variants\n" "Chromosome 23 (X):" "$x_total" >> variants_per_chr_report.txt
            printf "  %-20s %8d variants\n" "Chromosome 23 (X):" "$x_total" >> log_lines.txt
        fi
        if [ "$y_total" -gt 0 ]; then
            printf "  %-20s %8d variants\n" "Chromosome 24 (Y):" "$y_total" >> variants_per_chr_report.txt
            printf "  %-20s %8d variants\n" "Chromosome 24 (Y):" "$y_total" >> log_lines.txt
        fi
        if [ "$xy_plink2" -gt 0 ]; then
            printf "  %-20s %8d variants\n" "Chromosome 25 (PAR):" "$xy_plink2" >> variants_per_chr_report.txt
            printf "  %-20s %8d variants\n" "Chromosome 25 (PAR):" "$xy_plink2" >> log_lines.txt
        fi
        if [ "$mt_total" -gt 0 ]; then
            printf "  %-20s %8d variants\n" "Chromosome 26 (MT):" "$mt_total" >> variants_per_chr_report.txt
            printf "  %-20s %8d variants\n" "Chromosome 26 (MT):" "$mt_total" >> log_lines.txt
        fi
        
        # Compute total variants and format log line
        total_variants=$((total_autosomes + x_total + y_total + xy_plink2 + mt_total))
        printf "%-42s  Total: %8d variants\n" \
            "~{label}" "$total_variants" >> variants_per_chr_report.txt
        
        # Create a concise single-line summary for the pipeline log
        printf "%-42s  Total: %8d variants (A:%d X:%d Y:%d XY:%d MT:%d)\n" \
            "~{label}" "$total_variants" "$total_autosomes" "$x_total" "$y_total" "$xy_plink2" "$mt_total" \
            > log_line.txt
    >>>
    output {
        File          report    = "variants_per_chr_report.txt"
        String        log_line  = read_string("log_line.txt")
        Array[String] log_lines = read_lines("log_lines.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## SexCheck
## Runs PLINK --check-sex to compare reported sex with X-chromosome F-statistic.
## Writes a list of discordant samples (status != "OK") for removal.
## Aware of chromosome coding: recognizes both X (PLINK1) and 23 (PLINK2) as 
## the sex chromosome for F-stat calculation.
## Outputs problem samples with annotated chromosome information.
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

        # Extract FID and IID of samples that failed sex check.
        # PLINK already handles sex=0 (unknown) samples appropriately in --check-sex.
        awk '
            BEGIN { status_col = 5 }

            # Skip header/comment lines in the PLINK sexcheck output.
            /^#/ || /^$/ { next }

            # If the output contains a header line, capture the STATUS column index.
            # Otherwise, default to using column 5 as the status column.
            $1 == "FID" && $2 == "IID" {
                for (i = 1; i <= NF; i++) {
                    header[$i] = i
                }
                if ("STATUS" in header) {
                    status_col = header["STATUS"]
                }
                next
            }

            # Flag all samples with a non-OK sex-check status.
            ($status_col != "OK") {
                print $1, $2
            }' ~{output_prefix}_sexcheck.sexcheck \
            > ~{output_prefix}_problem_samples.txt

        # Check which sex chromosome coding was used and annotate report
        x_count_old=$(awk '$1=="X"' ~{bim_file} | wc -l)
        x_count_new=$(awk '$1=="23"' ~{bim_file} | wc -l)
        
        if [ $x_count_new -gt 0 ]; then
            CHR_CODING="PLINK2 (chr 23=X)"
        else
            CHR_CODING="PLINK1 (chr X)"
        fi
        
        # Add chromosome coding annotation to the report
        {
            printf "## Sex check report (chromosome coding: %s)\n" "$CHR_CODING"
            printf "## Run on X-chromosome variants: %d SNPs\n" "$((x_count_old + x_count_new))"
            printf "## Samples with status != OK are flagged as discordant\n"
            printf "\n"
            cat ~{output_prefix}_sexcheck.sexcheck
        } > ~{output_prefix}_sexcheck_annotated.sexcheck
    >>>
    output {
        File sexcheck_report = "~{output_prefix}_sexcheck_annotated.sexcheck"   # Annotated report with chr coding
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
## MacFilter
## Removes monomorphic SNPs (MAC = 0) and SNPs below mac_threshold.
## Steps:
##   1. Compute allele frequencies before filtering (--freq) → freq_before.frq
##   2. Remove monomorphic SNPs (--mac 0)
##   3. Remove SNPs below mac_threshold
##   4. Compute allele frequencies after filtering → freq_after.frq
##   5. R script plots MAC distribution and before/after barplot
## freq_before and freq_after are QC diagnostics only. The allele frequency
## file used for imputation prep is computed separately in PrepareForImputation
## on the final post-relatedness dataset.
## -----------------------------------------------------------------------------
task MacFilter {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        Int    mac_threshold
        String output_prefix
        String plink_bin
        String rscript_bin
        File   mac_plot_r
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

        # Remove monomorphic SNPs (MAC exactly 0)
        # is also implicitely done by the next step, so this seems redundant
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --mac 1 \
            --make-bed \
            --out tmp_no_monomorphic

        # Apply the chosen MAC threshold
        ~{plink_bin} \
            --bed tmp_no_monomorphic.bed \
            --bim tmp_no_monomorphic.bim \
            --fam tmp_no_monomorphic.fam \
            --mac ~{mac_threshold} \
            --make-bed \
            --out ~{output_prefix}

        # Compute allele frequencies after filtering (for plot)
        ~{plink_bin} \
            --bed ~{output_prefix}.bed \
            --bim ~{output_prefix}.bim \
            --fam ~{output_prefix}.fam \
            --freq \
            --out freq_after

        ~{rscript_bin} --vanilla ~{mac_plot_r} ~{mac_threshold}
    >>>
    output {
        File out_bed     = "~{output_prefix}.bed"
        File out_bim     = "~{output_prefix}.bim"
        File out_fam     = "~{output_prefix}.fam"
        File log         = "~{output_prefix}.log"
        File freq_before = "freq_before.frq"
        File freq_after  = "freq_after.frq"
        File mac_plot    = "mac_plot.png"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## AncestryPCA
## Merges study data with 1000 Genomes Phase 3 reference panel and runs a
## joint PCA to assign superpopulation ancestry labels to study samples.
##
## Steps:
##   1. Extract SNPs present in both study and 1000G reference (.bim overlap)
##   2. Merge study data with 1000G using the overlapping SNP list
##      - Strand flips are attempted automatically; unresolvable SNPs removed
##   3. Run PCA on the merged dataset
##   4. R script assigns superpopulation labels using a random forest trained
##      on 1000G samples with known labels, and plots study samples overlaid
##      on 1000G reference clusters
##
## Outputs:
##   eigenvec        — PC scores for all samples (study + 1000G)
##   pca_plot        — PC1 vs PC2 scatter plot coloured by population
##   assignments     — TSV: study sample IID, predicted superpop, probability
## -----------------------------------------------------------------------------
task AncestryPCA {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   ref_bed          # 1000G reference .bed
        File   ref_bim          # 1000G reference .bim
        File   ref_fam          # 1000G reference .fam
        File   ref_psam         # 1000G sample metadata with SuperPop column
        Int    n_pcs
        String output_prefix
        String plink_bin
        String plink2_bin
        String rscript_bin
        File   ancestry_pca_plot_r
    }
    command <<<
        set -euo pipefail

        # Step 1: find SNPs present in both study and 1000G bim files
        awk '{print $2}' ~{bim_file} | sort > study_snps.txt
        awk '{print $2}' ~{ref_bim}  | sort > ref_snps.txt
        comm -12 study_snps.txt ref_snps.txt > overlap_snps.txt

        echo "Overlapping SNPs: $(wc -l < overlap_snps.txt)"

        # Step 2: extract overlapping SNPs from both datasets
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --extract overlap_snps.txt \
            --make-bed \
            --out study_overlap

        ~{plink_bin} \
            --bed ~{ref_bed} \
            --bim ~{ref_bim} \
            --fam ~{ref_fam} \
            --extract overlap_snps.txt \
            --make-bed \
            --out ref_overlap

        # Step 3: merge — first attempt
        ~{plink_bin} \
            --bed study_overlap.bed \
            --bim study_overlap.bim \
            --fam study_overlap.fam \
            --bmerge ref_overlap.bed ref_overlap.bim ref_overlap.fam \
            --make-bed \
            --out merged_study_1kg || true

        # If strand flip errors, attempt flip and retry
        if [ -f merged_study_1kg-merge.missnp ]; then
            echo "Strand flip issues detected — flipping and retrying merge"

            ~{plink_bin} \
                --bed study_overlap.bed \
                --bim study_overlap.bim \
                --fam study_overlap.fam \
                --flip merged_study_1kg-merge.missnp \
                --make-bed \
                --out study_overlap_flipped

            ~{plink_bin} \
                --bed study_overlap_flipped.bed \
                --bim study_overlap_flipped.bim \
                --fam study_overlap_flipped.fam \
                --bmerge ref_overlap.bed ref_overlap.bim ref_overlap.fam \
                --make-bed \
                --out merged_study_1kg || true

            # If still unresolvable SNPs remain after flip, exclude them
            if [ -f merged_study_1kg-merge.missnp ]; then
                echo "Excluding unresolvable SNPs after flip"
                ~{plink_bin} \
                    --bed study_overlap_flipped.bed \
                    --bim study_overlap_flipped.bim \
                    --fam study_overlap_flipped.fam \
                    --exclude merged_study_1kg-merge.missnp \
                    --make-bed \
                    --out study_overlap_clean

                ~{plink_bin} \
                    --bed ref_overlap.bed \
                    --bim ref_overlap.bim \
                    --fam ref_overlap.fam \
                    --exclude merged_study_1kg-merge.missnp \
                    --make-bed \
                    --out ref_overlap_clean

                ~{plink_bin} \
                    --bed study_overlap_clean.bed \
                    --bim study_overlap_clean.bim \
                    --fam study_overlap_clean.fam \
                    --bmerge ref_overlap_clean.bed ref_overlap_clean.bim ref_overlap_clean.fam \
                    --make-bed \
                    --out merged_study_1kg
            fi
        fi

        # Step 4: run PCA on merged dataset
        ~{plink2_bin} \
            --bed merged_study_1kg.bed \
            --bim merged_study_1kg.bim \
            --fam merged_study_1kg.fam \
            --pca ~{n_pcs} \
            --out ~{output_prefix}

        # Step 5: plot and assign ancestry labels
        ~{rscript_bin} --vanilla ~{ancestry_pca_plot_r} \
            ~{output_prefix}.eigenvec \
            ~{ref_psam} \
            ~{n_pcs} \
            ~{output_prefix}

        # Create a single-line summary of ancestry assignments for the pipeline log
        if [ -f ~{output_prefix}_ancestry_assignments.tsv ]; then
            counts=$(awk -F'\t' '
                NR==1 { for(i=1;i<=NF;i++) hdr[$i]=i; \
                         if ("superpop" in hdr) c=hdr["superpop"]; \
                         else if ("predicted_pop" in hdr) c=hdr["predicted_pop"]; \
                         else if ("predicted" in hdr) c=hdr["predicted"]; \
                         else c=(NF>=2?2:1); next }
                NR>1 && c { count[$c]++ }
                END { sep=""; for (k in count) { printf "%s%s:%d", sep, k, count[k]; sep=" " } }
            ' ~{output_prefix}_ancestry_assignments.tsv)

            printf "%-42s  %s\n" "Step 4   Ancestry PCA" "$counts" > log_line.txt
        else
            printf "%-42s  %s\n" "Step 4   Ancestry PCA" "[no assignments]" > log_line.txt
        fi
    >>>
    output {
        File eigenvec    = "~{output_prefix}.eigenvec"
        File eigenval    = "~{output_prefix}.eigenval"
        File pca_plot    = "~{output_prefix}_ancestry_pca.png"
        File assignments = "~{output_prefix}_ancestry_assignments.tsv"
        File log         = "~{output_prefix}.log"
        String log_line  = read_string("log_line.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## PrepareForImputation
## Checks and prepares the final QC dataset for upload to the TOPMed imputation
## server. Produces per-chromosome VCF files.
##
## Steps:
##   1. Compute allele frequencies on the final QC dataset (--freq)
##   2. Run HRC-1000G-check-bim.pl — produces Run-plink.sh with PLINK commands
##      that fix strand, remove problem SNPs, and split by chromosome
##   3. Execute the generated Run-plink.sh
##   4. Convert each per-chromosome PLINK dataset to bgzipped, tabix-indexed VCF
##
## The check-bim script removes:
##   - SNPs absent from the TOPMed reference panel
##   - Ambiguous SNPs (A/T and C/G) that cannot be strand-resolved
##   - SNPs with allele frequency difference > 0.2 vs reference
##   - SNPs with mismatched positions or alleles
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

        # Ensure the generated Run-plink.sh uses the configured PLINK binary
        # rather than relying on `plink` being on the PATH.
        PLINK_BIN="~{plink_bin}"
        # Use a different delimiter (#) so replacement paths with / don't
        # break the perl s/// expression. Keep \b word-boundaries intact.
        ~{perl_bin} -i -pe 's#\bplink\b#'"$PLINK_BIN"'#g' Run-plink.sh

        # Step 3: Execute the generated PLINK commands
        # This performs strand flips, removes problem SNPs, and splits by chr
        bash Run-plink.sh

        # Step 4: Convert each per-chromosome PLINK file to VCF
        # TOPMed server requires one bgzipped, tabix-indexed VCF per chromosome
        # Derive chromosome list from the bim file so we handle whatever
        # chromosomes were selected in the config (not just 1-22).
        CHROMS=$(awk '{print $1}' ~{bim_file} | sort -u -V)

        for CHR in $CHROMS; do
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
