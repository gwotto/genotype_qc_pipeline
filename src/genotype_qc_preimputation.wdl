version 1.0

## =============================================================================
## genotype_qc_preimputation.wdl — Genotype QC Pre-imputation Pipeline
## =============================================================================
##
## PURPOSE
##   Performs standard pre-imputation QC on PLINK-format genotype data,
##   following the Anderson et al. (2010) protocol. Produces a clean binary
##   PLINK dataset and per-chromosome VCFs ready for imputation servers
##   (HRC, TOPMed).
##
## PIPELINE STEPS (in order)
##   0.  [Optional] Convert text ped/map → binary bed/bim/fam
##   0b. Handle duplicate sample IDs (suffix-based deduplication)
##   1.  SNP missingness filter        (--geno)
##   2.  Sample missingness filter     (--mind)
##   2b. Report variants per chromosome (informational)
##   —   Threshold sweep visualisation (informational plots only)
##   3.  Sex check & removal of fails  (skipped if no X-chromosome SNPs;
##       samples with sex code 0 are retained as "unknown")
##   4.  Ancestry PCA against 1000 Genomes (reference-only PCA + projection)
##         a. SNPs overlapping study and reference are identified by rsID
##            and positional concordance (guards against build mismatches).
##         b. Optionally restrict reference to a subset of superpopulations
##            via pca_reference_populations.
##         c. PCA is computed on the reference panel only; variant weights
##            (biallelic-var-wts) are saved and used to project ALL study
##            samples (including related) into the reference PC space.
##         d. A random forest classifier trained on the reference eigenvec
##            assigns each study sample a superpopulation label and
##            confidence probability. Samples below ancestry_prob_threshold
##            are marked "unassigned".
##   5.  MAC filter & monomorphic SNP removal (full cohort)
##   6.  HWE filter (detection on test_populations subset, removal from full cohort)
##   7a. LD pruning (helper for heterozygosity check only)
##   7b. Heterozygosity outlier removal (detection and removal within subset)
##   8.  Chromosome filter (configurable via chr_args) — FINAL QC DATASET
##   9.  Relatedness check (KING kinship; samples flagged, not removed)
##   10. Prepare for imputation (HRC-1000G-check-bim.pl; per-chr VCFs)
##
## KEY DESIGN DECISIONS
##   • Reference-only PCA: PC axes are defined by 1000G structure alone,
##     so they are stable regardless of study size or relatedness.
##   • MAC: applied to full cohort (rarity is a global property).
##   • HWE: detected within test_populations, removed from full cohort
##     (prevents false violations in admixed samples).
##   • Heterozygosity: detected and removed within subset only (preserves
##     samples from untested populations).
##   • Relatedness: flagged on full cohort; user decides which pairs to break.
##
## CHROMOSOME ENCODING
##   Uses PLINK chromosome encoding throughout:
##     Autosomes 1–22 | X = 23 | Y = 24 | PAR (XY) = 25 | MT = 26
##
## EXECUTION
##   java -jar cromwell-92.jar run src/genotype_qc_preimputation.wdl \
##        -i config/genotype_qc_preimputation_inputs.json
##
## DEPENDENCIES (must be on $PATH or supplied via *_bin inputs)
##   plink 1.9, plink2 2.0, Rscript ≥ 4.0, perl, bcftools, bgzip, tabix
##
## REFERENCE
##   Anderson CA et al. (2010) Nat Protoc 5:1564–1573
##   https://doi.org/10.1038/nprot.2010.116
## =============================================================================



workflow genotype_qc_preimputation {

     String pipeline_version = "v2026-07.1"

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
        File create_qc_table_r # R sript: creates table of ancestry assignments and relatedness checks
    
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
        Int   ld_step         # Window step size in SNPs             (e.g. 1)
        Float pihat_min       # Min pi-hat to flag related pairs     (e.g. 0.2)
        Float king_cutoff     # KING kinship coefficient cutoff      (e.g. 0.0884 ≈ 3rd degree)

        # -- Chromosome filtering -----------------------------------------------
        String chr_args  # PLINK chromosome filter args (e.g. "--chr 1-22" or "--autosome")

        # -- PCA ----------------------------------------------------------------
        ## PCA to assign ancestry
        Int    n_pcs     # Number of principal components to compute
        Float  ancestry_prob_threshold = 0.5  # Probability threshold for ancestry assignment (0.0-1.0)

        # Number of within-cohort principal components emitted as association
        # covariates for the final (imputation-prepared) dataset. Computed on
        # unrelated samples and projected onto all samples.
        Int    n_covariate_pcs

        # Additional deliverable: an unrelated single-ancestry subset of the
        # imputation-ready dataset, exported in PLINK 1 and PLINK 2 formats with
        # its own within-subset covariate PCs. Superpopulation label to select
        # (must match AncestryPCA labels: EUR/AFR/EAS/SAS/AMR).
        String subset_population

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
        String test_populations

        # Comma-separated list of 1000G superpopulations to include in the
        # ancestry PCA reference panel (e.g. "EUR,AFR,EAS,SAS,AMR").
        # Use "ALL" to include all 2504 phase-3 samples.
        String pca_reference_populations

        # -- Imputation preparation ----------------------------------------------
        File   check_bim_pl  # Will Rayner's HRC-1000G-check-bim.pl script
        File   hrc_ref_freq  # HRC r1.1 GRCh37 reference frequency file
        String perl_bin      # Perl interpreter
        String bcftools_bin  # bcftools binary
        String bgzip_bin   # bgzip binary
        String tabix_bin   # tabix binary

    }

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
    # Principal components are computed on the LD-pruned 1000G Phase 3
    # reference panel alone (saving PLINK2 variant weights), then ALL study
    # samples are projected onto those reference axes with --score. Because the
    # study cohort never enters the eigendecomposition, related study samples
    # cannot distort the PCs — so no relatedness filtering is needed here (the
    # relatedness check runs later, at step 9, on the QC'd variant set).
    # Only autosomal SNPs shared by study and reference are used.
    call AncestryPCA {
        input:
            bed_file            = select_first([RemoveSexFails.out_bed, MindFilter.out_bed]),
            bim_file            = select_first([RemoveSexFails.out_bim, MindFilter.out_bim]),
            fam_file            = select_first([RemoveSexFails.out_fam, MindFilter.out_fam]),
            ref_bed             = ref_1kg_bed,
            ref_bim             = ref_1kg_bim,
            ref_fam             = ref_1kg_fam,
            ref_psam            = ref_1kg_psam,
            ld_regions          = ld_regions,
            ld_window_kb        = ld_window_kb,
            ld_step             = ld_step,
            ld_r2               = ld_r2,
            n_pcs               = n_pcs,
            ancestry_prob_threshold    = ancestry_prob_threshold,
            output_prefix              = output_prefix + "_ancestry_pca",
            pca_reference_populations  = pca_reference_populations,
            plink_bin                  = plink_bin,
            plink2_bin                 = plink2_bin,
            rscript_bin                = rscript_bin,
            ancestry_pca_plot_r        = ancestry_pca_plot_r
    }

    # --- Step 5: MAC filter (full cohort) ----------------------------------------

    call MacFilter {
        input:
            bed_file      = select_first([RemoveSexFails.out_bed, MindFilter.out_bed]),
            bim_file      = select_first([RemoveSexFails.out_bim, MindFilter.out_bim]),
            fam_file      = select_first([RemoveSexFails.out_fam, MindFilter.out_fam]),
            mac_threshold = mac_threshold,
            output_prefix = output_prefix + "_QC5",
            plink_bin     = plink_bin,
            rscript_bin   = rscript_bin,
            mac_plot_r    = mac_plot_r
    }

    call CountBimFam as LogStep5 {
        input:
            bim_file = MacFilter.out_bim,
            fam_file = MacFilter.out_fam,
            label    = "Step 5   MAC filter"
    }

    # --- Step 6: HWE filter (subset detection, full cohort removal) -----------
    # Population structure can cause deviation from HWE after pooling (Wahlund effect).
    # Here, the ancestry group(s) are chosen by the user in the config file, e.g. 
    # "EUR" or "ALL". The selected population is used to detect variants deviating from 
    # HWE. These SNPs are then removed from the entire cohort.
    
    # Build ancestry keep-list for HWE and heterozygosity subsetting
    call MakeAncestryKeepList {
        input:
            assignments = AncestryPCA.assignments,
            fam_file    = MacFilter.out_fam,
            populations = test_populations,
            output_prefix = output_prefix
    }

    call PlinkFilter as SubsetAncestry_HWE {
        input:
            bed_file      = MacFilter.out_bed,
            bim_file      = MacFilter.out_bim,
            fam_file      = MacFilter.out_fam,
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
            bed_file      = MacFilter.out_bed,
            bim_file      = MacFilter.out_bim,
            fam_file      = MacFilter.out_fam,
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

    # --- Step 7a: LD pruning as a helper for the heterozygosity (step 7b) and
    # relatedness (step 9) checks only — the LD SNPs are NOT removed from the
    # dataset. High-LD regions are excluded, which suits both uses.

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
            label         = "Step 7a  LD pruning (heterozygosity + relatedness)",
            plink_bin     = plink_bin
    }

    # --- Step 7b: Heterozygosity filter (subset detection, subset removal) ---
    # Detect heterozygosity outliers only in user-defined subpopulation(s). Different 
    # and admixed ancestries could have heterozygosity rates that might count as 
    # outlyers in the selected population, but should not be removed.
    
    call PlinkFilter as SubsetAncestry_HET {
        input:
            bed_file      = RemoveVariantsAfterHWE.out_bed,
            bim_file      = RemoveVariantsAfterHWE.out_bim,
            fam_file      = RemoveVariantsAfterHWE.out_fam,
            plink_args    = "--keep " + MakeAncestryKeepList.keep_list,
            output_prefix = output_prefix + "_ancestry_subset_HET",
            plink_bin     = plink_bin
    }

    call HeterozygosityCheck as HeterozygosityCheck {
        input:
            bed_file                  = SubsetAncestry_HET.out_bed,
            bim_file                  = SubsetAncestry_HET.out_bim,
            fam_file                  = SubsetAncestry_HET.out_fam,
            prune_in                  = LdPruningHet.prune_in,
            check_heterozygosity_r    = check_heterozygosity_r,
            heterozygosity_outliers_r = heterozygosity_outliers_r,
            plink_bin                 = plink_bin,
            rscript_bin               = rscript_bin
    }

    call RemoveSamples as RemoveHetFailsFromSubset {
        input:
            bed_file      = SubsetAncestry_HET.out_bed,
            bim_file      = SubsetAncestry_HET.out_bim,
            fam_file      = SubsetAncestry_HET.out_fam,
            remove_list   = HeterozygosityCheck.het_fail_ind,
            output_prefix = output_prefix + "_QC7_het_subset",
            plink_bin     = plink_bin
    }

    # Merge het-filtered subset back with other populations (if test_populations != "ALL")
    # When test_populations = "ALL": use the het-filtered subset directly (all het failures removed)
    # When testing a specific population: combine het-filtered subset + non-tested populations
    call MergeKeepLists as MergeHetKeepLists {
        input:
            tested_fam      = RemoveHetFailsFromSubset.out_fam,
            all_fam         = RemoveVariantsAfterHWE.out_fam,
            tested_keep     = MakeAncestryKeepList.keep_list,
            is_all_pops     = (test_populations == "ALL"),
            output_prefix   = output_prefix + "_het_merged_keep"
    }

    call PlinkFilter as MergeHetFiltered {
        input:
            bed_file      = RemoveVariantsAfterHWE.out_bed,
            bim_file      = RemoveVariantsAfterHWE.out_bim,
            fam_file      = RemoveVariantsAfterHWE.out_fam,
            plink_args    = "--keep " + MergeHetKeepLists.merged_keep,
            output_prefix = output_prefix + "_QC7",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep7 {
        input:
            bim_file = MergeHetFiltered.out_bim,
            fam_file = MergeHetFiltered.out_fam,
            label    = "Step 7b  Heterozygosity filter (subset detection, subset removal)"
    }

    # -- Step 8: Restrict to selected chromosomes -------------------------------------
    # Retains only chromosomes 1–22 for downstream association analysis.
    # Sex chromosomes and mitochondrial SNPs require separate handling.
    # Chromosome filtering args are configurable (e.g. "--chr 1-23" or "--chr 1-10,12-22").
    call PlinkFilter as ChromosomeFilter {
        input:
            bed_file      = MergeHetFiltered.out_bed,
            bim_file      = MergeHetFiltered.out_bim,
            fam_file      = MergeHetFiltered.out_fam,
            plink_args    = chr_args,
            output_prefix = output_prefix + "_QC8",
            plink_bin     = plink_bin
    }

    call CountBimFam as LogStep8 {
        input:
            bim_file = ChromosomeFilter.out_bim,
            fam_file = ChromosomeFilter.out_fam,
            label    = "Step 8   Chromosome filter"
    }

    # -- Step 9: Relatedness check -----------------------------------------
    # KING kinship on the final QC dataset, so that kinship is estimated from
    # MAC- and HWE-filtered variants rather than the raw call set. Reuses the
    # LD-pruned SNP list from step 7a (high-LD regions already excluded);
    # --extract simply intersects it with the variants remaining in QC8.
    # Related samples are FLAGGED, never removed — users apply --remove in
    # downstream association tests. The flags feed steps 12 and 13 and the
    # sample QC status table.
    call RelatednessCheck {
        input:
            bed_file      = ChromosomeFilter.out_bed,
            bim_file      = ChromosomeFilter.out_bim,
            fam_file      = ChromosomeFilter.out_fam,
            prune_in      = LdPruningHet.prune_in,
            pihat_min     = pihat_min,
            king_cutoff   = king_cutoff,
            output_prefix = output_prefix + "_relatedness",
            plink_bin     = plink_bin,
            plink2_bin    = plink2_bin
    }

    call CountBimFam as LogStep9 {
        input:
            bim_file = ChromosomeFilter.out_bim,
            fam_file = ChromosomeFilter.out_fam,
            label    = "Step 9   Relatedness check              (related flagged, not removed)"
    }

    # -- Step 10: Prepare for TOPMed imputation ----------------------------
    # Aligns strand orientation to the HRC reference panel using Will
    # Rayner's check-bim script. Removes SNPs not in the reference, ambiguous
    # A/T and C/G SNPs that cannot be strand-resolved, and SNPs with large
    # allele frequency differences vs the reference.
    # Outputs per-chromosome VCF files ready for upload to the TOPMed
    # imputation server (https://imputation.biodatacatalyst.nhlbi.nih.gov).
    
    call PrepareForImputation {
        input:
            bed_file        = ChromosomeFilter.out_bed,
            bim_file        = ChromosomeFilter.out_bim,
            fam_file        = ChromosomeFilter.out_fam,
            check_bim_pl    = check_bim_pl,
            hrc_ref_freq    = hrc_ref_freq,
            output_prefix   = output_prefix + "_imputation",
            plink_bin       = plink_bin,
            perl_bin        = perl_bin,
            bcftools_bin    = bcftools_bin,
            bgzip_bin       = bgzip_bin,
            tabix_bin       = tabix_bin
    }

    # Report per-chromosome variant counts AFTER imputation prep (check-bim),
    # i.e. on the harmonised variant set actually written to the VCFs.
    call VariantsPerChromosome as PostImputationVariantsPerChromosome {
        input:
            bim_file = PrepareForImputation.combined_bim,
            label    = "Step 10a Variants per chromosome (harmonised)"
    }

    # -- Step 11: PLINK 2 (pgen/pvar/psam) copy of the combined dataset ----
    call ConvertToPlink2 as CombinedToPlink2 {
        input:
            bed_file      = PrepareForImputation.combined_bed,
            bim_file      = PrepareForImputation.combined_bim,
            fam_file      = PrepareForImputation.combined_fam,
            output_prefix = output_prefix + "_imputation_combined",
            label         = "Step 11  PLINK 2 conversion (combined)",
            plink2_bin    = plink2_bin
    }

    # -- Step 12: Within-cohort covariate PCA on the combined dataset ------
    # PCs are computed on unrelated samples (KING-flagged relatives excluded)
    # and projected onto all samples, giving association covariates for the
    # whole cohort on a common scale.
    call CovariatePCA {
        input:
            bed_file      = PrepareForImputation.combined_bed,
            bim_file      = PrepareForImputation.combined_bim,
            fam_file      = PrepareForImputation.combined_fam,
            related_list  = RelatednessCheck.king_cutoff_out,
            ld_regions    = ld_regions,
            ld_window_kb  = ld_window_kb,
            ld_step       = ld_step,
            ld_r2         = ld_r2,
            n_pcs         = n_covariate_pcs,
            output_prefix = output_prefix + "_covariate_pca",
            plink_bin     = plink_bin,
            plink2_bin    = plink2_bin
    }

    # -- Create comprehensive QC status table ------------------------------
    # Combines ancestry, relatedness, and the covariate PCs into one per-sample
    # table so users can subset the cohort (by ancestry / relatedness) and pull
    # covariates from a single file.
    call CreateQCStatusTable {
        input:
            ancestry_assignments = AncestryPCA.assignments,
            relatedness_flagged  = RelatednessCheck.king_cutoff_out,
            pca_covariates       = CovariatePCA.covariates,
            output_prefix        = output_prefix,
            rscript_bin          = rscript_bin,
            create_qc_table_r    = create_qc_table_r
    }

    # -- Step 13: Unrelated single-ancestry subset ---
    # Select samples assigned to `subset_population` that are NOT KING-flagged
    # as related, emit bed/bim/fam (in SubsetUnrelated) plus a PLINK 2 copy, and
    # compute covariate PCs *within* the subset (correct for within-population
    # association; the whole-cohort PCs mostly capture between-ancestry structure).
    call SubsetUnrelatedPopulation as SubsetUnrelated {
        input:
            bed_file      = PrepareForImputation.combined_bed,
            bim_file      = PrepareForImputation.combined_bim,
            fam_file      = PrepareForImputation.combined_fam,
            assignments   = AncestryPCA.assignments,
            related_list  = RelatednessCheck.king_cutoff_out,
            population    = subset_population,
            output_prefix = output_prefix + "_unrelated_" + subset_population,
            plink2_bin    = plink2_bin
    }

    call ConvertToPlink2 as SubsetToPlink2 {
        input:
            bed_file      = SubsetUnrelated.out_bed,
            bim_file      = SubsetUnrelated.out_bim,
            fam_file      = SubsetUnrelated.out_fam,
            output_prefix = output_prefix + "_unrelated_" + subset_population,
            plink2_bin    = plink2_bin
    }

    call CovariatePCA as SubsetCovariatePCA {
        input:
            bed_file      = SubsetUnrelated.out_bed,
            bim_file      = SubsetUnrelated.out_bim,
            fam_file      = SubsetUnrelated.out_fam,
            related_list  = RelatednessCheck.king_cutoff_out,
            ld_regions    = ld_regions,
            ld_window_kb  = ld_window_kb,
            ld_step       = ld_step,
            ld_r2         = ld_r2,
            n_pcs         = n_covariate_pcs,
            output_prefix = output_prefix + "_unrelated_" + subset_population + "_pca",
            plink_bin     = plink_bin,
            plink2_bin    = plink2_bin
    }

    call CountBimFam as LogSubset {
        input:
            bim_file = SubsetUnrelated.out_bim,
            fam_file = SubsetUnrelated.out_fam,
            label    = "Step 13  Unrelated " + subset_population + " subset"
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
            LogStep7.line,
            LogStep8.line,
            LogStep9.line,
            PrepareForImputation.log_line,
            PostImputationVariantsPerChromosome.log_line
        ],
        PostImputationVariantsPerChromosome.log_lines,
        [
            CombinedToPlink2.log_line,
            CovariatePCA.log_line,
            LogSubset.line
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
        File final_bed = ChromosomeFilter.out_bed
        File final_bim = ChromosomeFilter.out_bim
        File final_fam = ChromosomeFilter.out_fam

        # Sex check outputs (only present if X-chromosome SNPs exist)
        File? sexcheck_report = SexCheck.sexcheck_report
        File? problem_samples = SexCheck.problem_samples
        File? variants_per_chr_initial = InitialVariantsPerChromosome.report

        # Heterozygosity check outputs
        File het_check_report = HeterozygosityCheck.r_check_het  # Per-sample het rates
        File het_fail_samples = HeterozygosityCheck.het_fail_ind  # Outlier sample list

        # MAC filter outputs (QC diagnostics; imputation freq computed separately)
        File mac_plot        = MacFilter.mac_plot
        File maf_freq_before = MacFilter.freq_before
        File maf_freq_after  = MacFilter.freq_after

        # Relatedness outputs (step 9, on the final QC dataset; flagged not removed)
        File pihat_genome       = RelatednessCheck.pihat_genome    # All pairs above pihat_min
        File king_cutoff_out_id = RelatednessCheck.king_cutoff_out # Flagged related samples
        File king_cutoff_in_id  = RelatednessCheck.king_cutoff_in  # Unrelated samples (for downstream use)

        # LD-pruned SNP lists
        File pca_prune_in = AncestryPCA.prune_in   # LD-pruned reference SNPs used for ancestry PCA
        File prune_in     = LdPruningHet.prune_in  # used for heterozygosity + relatedness checks
        File prune_out    = LdPruningHet.prune_out

        # Threshold sweep plot — shows SNP/sample retention across filtering thresholds
        File threshold_plot = ThresholdSweep.sweep_plot

        # Sample QC status table — comprehensive reference for filtering decisions
        File sample_qc_status_table = CreateQCStatusTable.qc_status_table  # FID, IID, ancestry, ancestry_prob, related

        # Ancestry PCA against 1000 Genomes
        File ancestry_pca_eigenvec     = AncestryPCA.eigenvec       # PC scores (study + 1000G)
        File ancestry_pca_plot         = AncestryPCA.pca_plot        # Study samples overlaid on 1000G
        File ancestry_assignments      = AncestryPCA.assignments     # Per-sample superpopulation assignments

        # Per-chromosome variant counts (post check-bim — matches the VCFs)
        File variants_per_chr_postprep = PostImputationVariantsPerChromosome.report

        # Imputation-ready VCFs — upload these to the TOPMed server
        Array[File] imputation_vcfs     = PrepareForImputation.vcf_gz
        Array[File] imputation_vcf_tbis = PrepareForImputation.vcf_tbi
        File        check_bim_log       = PrepareForImputation.check_bim_log

        # Combined (all-chromosome) imputation-ready dataset — PLINK equivalent
        # of the per-chromosome VCFs, in both PLINK 1 and PLINK 2 formats
        File combined_bed  = PrepareForImputation.combined_bed
        File combined_bim  = PrepareForImputation.combined_bim
        File combined_fam  = PrepareForImputation.combined_fam
        File combined_pgen = CombinedToPlink2.out_pgen
        File combined_pvar = CombinedToPlink2.out_pvar
        File combined_psam = CombinedToPlink2.out_psam

        # Within-cohort covariate PCA (unrelated axes, all samples projected)
        File covariate_pcs      = CovariatePCA.covariates    # FID IID PC1..PCn
        File covariate_eigenvec = CovariatePCA.eigenvec
        File covariate_eigenval = CovariatePCA.eigenval

        # Unrelated single-ancestry subset (subset_population, default EUR):
        # PLINK 1 (bed/bim/fam) + PLINK 2 (pgen/pvar/psam) and within-subset covariate PCs
        File subset_keep       = SubsetUnrelated.keep_list      # FID IID of the subset
        File subset_bed        = SubsetUnrelated.out_bed
        File subset_bim        = SubsetUnrelated.out_bim
        File subset_fam        = SubsetUnrelated.out_fam
        File subset_pgen       = SubsetToPlink2.out_pgen
        File subset_pvar       = SubsetToPlink2.out_pvar
        File subset_psam       = SubsetToPlink2.out_psam
        File subset_covariates = SubsetCovariatePCA.covariates  # FID IID PC1..PCn (within-subset)
        File subset_eigenval   = SubsetCovariatePCA.eigenval
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
                # 
                awk 'NR==FNR{iids[$0]=1; next} ($2 in iids){print $1, $2}' selected_iids.txt ~{fam_file} > "$OUT" || true
                rm -f selected_iids.txt || true
            fi

            # Validate that keep-list is not empty
            n_lines=$(wc -l < "$OUT")
            if [ "$n_lines" -eq 0 ]; then
                echo "ERROR: No samples found for populations: $POPS" >&2
                echo "Check that population names match those in the ancestry assignments file" >&2
                exit 1
            fi
            
            echo $n_lines > keep_count.txt
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
                --out tmp_geno_$T --silent
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
                --out tmp_mind_$T --silent
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


    # Helper task: Merge het-filtered tested population with non-tested populations
task MergeKeepLists {
    input {
        File tested_fam      # FAM file with het-filtered tested population
        File all_fam         # Original FAM file with all populations
        File tested_keep     # Keep-list of tested populations
        Boolean is_all_pops  # Whether testing all populations
        String output_prefix
    }
        
    command <<<
        set -euo pipefail
        
        if [ "~{is_all_pops}" = "true" ]; then
            # When testing all populations, just use the het-filtered FAM
            awk '{print $1, $2}' ~{tested_fam} > ~{output_prefix}_merged_keep.txt
        else
            # When testing a subset:
            # 1. Extract IIDs from tested population keep-list
            awk '{print $2}' ~{tested_keep} > tested_iids.txt
            
            # 2. Get het-filtered tested samples from tested_fam
            awk '{print $1, $2}' ~{tested_fam} > tested_filtered.txt
            
            # 3. Get non-tested populations (those NOT in the tested_keep list)
            awk 'NR==FNR{iids[$0]=1; next} !($2 in iids) {print $1, $2}' tested_iids.txt ~{all_fam} > other_pops.txt
            
            # 4. Combine: het-filtered tested population + non-tested populations
            cat tested_filtered.txt other_pops.txt > ~{output_prefix}_merged_keep.txt
            
            rm -f tested_iids.txt tested_filtered.txt other_pops.txt
        fi
    >>>
        
    output {
        File merged_keep = "~{output_prefix}_merged_keep.txt"
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
        set -euo pipefail
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
                chr_label=$(printf "Chromosome %2d:" "$CHR")
                printf "  %-20s %8d variants\n" "$chr_label" "$count" >> variants_per_chr_report.txt
                printf "  %-20s %8d variants\n" "$chr_label" "$count" >> log_lines.txt
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
        String label = "LD pruning"   # Label for the pipeline log line
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
            "~{label}" "$n_pruned" > log_line.txt
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
            --king-cutoff ~{king_cutoff} \
            --out plink2
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

        # Apply MAC threshold (also removes monomorphic SNPs since mac_threshold >= 1)
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
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
##   1. Restrict to autosomal SNPs (chr 1-22); intersect by rsID, filter to
##      positionally concordant SNPs, then intersect with LD-pruned list
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
        File   ref_bed              # 1000G reference .bed
        File   ref_bim              # 1000G reference .bim
        File   ref_fam              # 1000G reference .fam
        File   ref_psam             # 1000G sample metadata with SuperPop column
        File   ld_regions           # High-LD regions to exclude during pruning
        Int    ld_window_kb         # Sliding window size in kb
        Int    ld_step              # Step size in SNPs
        Float  ld_r2                # r² pruning threshold
        Int    n_pcs
        Float  ancestry_prob_threshold # ancestry assignment probability threshold
        String output_prefix
        String pca_reference_populations  # "ALL" or comma-separated e.g. "EUR,AFR,EAS,SAS,AMR"
        String plink_bin
        String plink2_bin
        String rscript_bin
        File   ancestry_pca_plot_r
    }
    command <<<
        set -euo pipefail

        # Step 1: find autosomal SNPs (chr 1-22 only) present in both datasets.
        # Non-autosomal SNPs excluded: sex-chr heterozygosity and PAR effects
        # distort principal components and are uninformative for ancestry inference.

        # (a) rsID intersection restricted to autosomes
        awk '$1+0 >= 1 && $1+0 <= 22 {print $2}' ~{bim_file} | sort > study_snps.txt
        awk '$1+0 >= 1 && $1+0 <= 22 {print $2}' ~{ref_bim}  | sort > ref_snps.txt
        comm -12 study_snps.txt ref_snps.txt > rsid_overlap.txt
        n_rsid=$(wc -l < rsid_overlap.txt)
        echo "Autosomal SNPs overlapping by rsID: $n_rsid"

        # (b) positional concordance check — guards against rsID collisions or
        #     genome build mismatches between study and reference.
        awk '$1+0 >= 1 && $1+0 <= 22 {print $2, $1"_"$4}' ~{bim_file} | sort -k1,1 > study_id_pos.txt
        awk '$1+0 >= 1 && $1+0 <= 22 {print $2, $1"_"$4}' ~{ref_bim}  | sort -k1,1 > ref_id_pos.txt
        join -1 1 -2 1 study_id_pos.txt ref_id_pos.txt \
            | awk '$2 == $3 {print $1}' | sort > overlap_snps.txt
        n_concordant=$(wc -l < overlap_snps.txt)
        n_discord=$(( n_rsid - n_concordant ))
        if [ "$n_discord" -gt 0 ]; then
            echo "WARNING: Excluded $n_discord rsID-matched SNPs with discordant chr:pos (build mismatch?)" >&2
        fi

        # Step 2: extract overlapping SNPs from the 1000G reference.
        # PCA is computed on the reference only — PC axes are defined by 1000G
        # population structure, independent of study composition or relatedness.
        ~{plink_bin} \
            --bed ~{ref_bed} \
            --bim ~{ref_bim} \
            --fam ~{ref_fam} \
            --extract overlap_snps.txt \
            --chr 1-22 \
            --make-bed \
            --out ref_overlap

        # Step 2b: optionally restrict reference to specified superpopulations.
        # psam column layout: #IID(1) PAT(2) MAT(3) SEX(4) SuperPop(5) Population(6)
        # plink1 --keep expects a two-column FID/IID file; FID is taken from the .fam.
        ref_pca=ref_overlap
        if [ "~{pca_reference_populations}" != "ALL" ]; then
            awk -v pops="~{pca_reference_populations}" '
                BEGIN { n=split(pops, a, ","); for (i=1;i<=n;i++) keep[a[i]]=1 }
                NR>1 && $5 in keep { print $1 }
            ' ~{ref_psam} | sort > pop_iids.txt

            # Resolve FID from the fam file (1000G fams often have FID=0)
            awk 'NR==FNR { want[$1]=1; next } $2 in want { print $1, $2 }' \
                pop_iids.txt ~{ref_fam} > ref_keep.txt

            n_ref_keep=$(wc -l < ref_keep.txt)
            echo "Reference samples after population filter (~{pca_reference_populations}): $n_ref_keep"

            ~{plink_bin} --bfile ref_overlap --keep ref_keep.txt --make-bed --out ref_overlap_pop
            ref_pca=ref_overlap_pop
        fi

        # Step 3: LD prune the reference overlap.
        ~{plink_bin} \
            --bfile ${ref_pca} \
            --exclude ~{ld_regions} --range \
            --indep-pairwise ~{ld_window_kb} kb ~{ld_step} ~{ld_r2} \
            --out ref_prune

        n_pruned=$(wc -l < ref_prune.prune.in)
        echo "Reference SNPs after LD pruning: $n_pruned"

        # Step 4: PCA on the 1000G reference, saving variant weights for projection.
        ~{plink2_bin} \
            --bfile ${ref_pca} \
            --extract ref_prune.prune.in \
            --pca ~{n_pcs} biallelic-var-wts \
            --out ~{output_prefix}

        # Step 4a: reference allele frequencies for projection standardisation.
        # These must match the frequencies used during PCA so that projected
        # study scores land on the same scale as the reference eigenvec scores.
        ~{plink2_bin} \
            --bfile ${ref_pca} \
            --extract ref_prune.prune.in \
            --freq \
            --out ref_prune

        # Step 4b: project ALL study samples (related and unrelated) onto the
        # reference PC space using the variant weights from step 4.
        # eigenvec.var column layout: #CHROM(1) ID(2) MAJ(3) NONMAJ(4) PC1(5)...PCN(N+4)
        # Column 2 = variant ID (rsID), column 4 = scored allele (non-major).
        # PLINK2 resolves strand orientation automatically from the allele column.
        # No chr:pos renaming is needed: the reference BIM uses rsIDs throughout.
        score_end_col=$(( ~{n_pcs} + 4 ))
        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --extract ref_prune.prune.in \
            --chr 1-22 \
            --read-freq ref_prune.afreq \
            --score ~{output_prefix}.eigenvec.var 2 4 header no-mean-imputation variance-standardize \
            --score-col-nums 5-${score_end_col} \
            --out all_study_projected

        # Step 5: plot and assign ancestry labels.
        # The R script uses the 1000G eigenvec to build reference clusters and
        # the all-sample sscore for classification (includes related samples).
        ~{rscript_bin} --vanilla ~{ancestry_pca_plot_r} \
            ~{output_prefix}.eigenvec \
            all_study_projected.sscore \
            ~{ref_psam} \
            ~{n_pcs} \
            ~{output_prefix} \
            ~{output_prefix}.eigenval \
            ~{ancestry_prob_threshold}

        # Write pipeline log lines: SNP filter chain then ancestry assignment counts
        printf "%-42s  rsID: %d  concordant: %d  LD-pruned: %d\n" \
            "Step 4   Ancestry PCA (SNPs)" "$n_rsid" "$n_concordant" "$n_pruned" > log_line.txt

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

            printf "%-42s  %s\n" "Step 4   Ancestry PCA (ancestry)" "$counts" >> log_line.txt
        else
            printf "%-42s  %s\n" "Step 4   Ancestry PCA (ancestry)" "[no assignments]" >> log_line.txt
        fi
    >>>
    output {
        File eigenvec     = "~{output_prefix}.eigenvec"
        File eigenvec_var = "~{output_prefix}.eigenvec.var"
        File eigenval     = "~{output_prefix}.eigenval"
        File sscore       = "all_study_projected.sscore"
        File pca_plot     = "~{output_prefix}_ancestry_pca.png"
        File assignments  = "~{output_prefix}_ancestry_assignments.tsv"
        File prune_in     = "ref_prune.prune.in"
        File prune_out    = "ref_prune.prune.out"
        File log          = "~{output_prefix}.log"
        String log_line   = read_string("log_line.txt")
    }
    runtime { maxRetries: 1 }
}

# Task: Create comprehensive sample QC status table
# Combines ancestry assignments and relatedness flags into a single reference table
# Output: FID, IID, ancestry, ancestry_prob, related
task CreateQCStatusTable {
    input {
        File ancestry_assignments      # From AncestryPCA.assignments
        File relatedness_flagged       # From RelatednessCheck.king_cutoff_out
        File? pca_covariates           # From CovariatePCA.covariates (FID IID PC1..PCn)
        String output_prefix
        String rscript_bin
        File create_qc_table_r
    }

    command <<<
        set -euo pipefail
        ~{rscript_bin} ~{create_qc_table_r} \
            ~{ancestry_assignments} \
            ~{relatedness_flagged} \
            ~{output_prefix} \
            ~{default="NA" pca_covariates}
    >>>
        
    output {
        File qc_status_table = "~{output_prefix}_sample_qc_status.tsv"
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

        # Step 5: merge the per-chromosome harmonized filesets back into a
        # single dataset spanning all chromosomes. This carries exactly the
        # variant set written to the per-chromosome imputation VCFs (post
        # strand-fix, post HRC-panel filtering), so the combined bed/bim/fam
        # is the PLINK-format equivalent of the uploaded VCFs.
        ls *-updated-chr*.bed 2>/dev/null | sed 's/\.bed$//' > combined_stems.txt || true
        n_stems=$(wc -l < combined_stems.txt)
        if [ "$n_stems" -eq 0 ]; then
            echo "ERROR: no per-chromosome harmonized filesets found to combine" >&2
            exit 1
        elif [ "$n_stems" -eq 1 ]; then
            stem=$(head -1 combined_stems.txt)
            ~{plink_bin} --bfile "$stem" --make-bed --out ~{output_prefix}_combined
        else
            awk '{print $1".bed "$1".bim "$1".fam"}' combined_stems.txt > merge_list.txt
            ~{plink_bin} --merge-list merge_list.txt --make-bed --out ~{output_prefix}_combined
        fi

        # Count total SNPs retained across all per-chromosome bim files
        # and format a log line matching the CountBimFam style
        n_snps=$(cat *-updated-chr*.bim 2>/dev/null | wc -l || echo 0)
        n_samples=$(wc -l < ~{fam_file})
        printf "%-42s  SNPs: %7d   Samples: %5d\n" \
            "Step 10  Imputation prep (check-bim)" "$n_snps" "$n_samples" \
            > step11_log.txt
    >>>
    output {
        Array[File] vcf_gz        = glob("~{output_prefix}_chr*.vcf.gz")
        Array[File] vcf_tbi       = glob("~{output_prefix}_chr*.vcf.gz.tbi")
        File        check_bim_log = "check-bim.log"
        String      log_line      = read_string("step11_log.txt")
        # Combined (all-chromosome) harmonized dataset — PLINK equivalent of the VCFs
        File        combined_bed  = "~{output_prefix}_combined.bed"
        File        combined_bim  = "~{output_prefix}_combined.bim"
        File        combined_fam  = "~{output_prefix}_combined.fam"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## ConvertToPlink2
## Converts a PLINK 1.x binary dataset (bed/bim/fam) to PLINK 2 format
## (pgen/pvar/psam). Used to emit the combined imputation-ready dataset in
## PLINK 2 format alongside the PLINK 1 bed/bim/fam.
## -----------------------------------------------------------------------------
task ConvertToPlink2 {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        String output_prefix
        String label = "PLINK 2 conversion"   # Label for the pipeline log line
        String plink2_bin
    }
    command <<<
        set -euo pipefail
        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --make-pgen \
            --out ~{output_prefix}

        # pvar/psam carry ## and #CHROM header lines, so skip them when counting.
        n_snps=$(awk '!/^#/' ~{output_prefix}.pvar | wc -l)
        n_samples=$(awk '!/^#/' ~{output_prefix}.psam | wc -l)
        printf "%-42s  SNPs: %7d   Samples: %5d\n" \
            "~{label}" "$n_snps" "$n_samples" > log_line.txt
    >>>
    output {
        File   out_pgen = "~{output_prefix}.pgen"
        File   out_pvar = "~{output_prefix}.pvar"
        File   out_psam = "~{output_prefix}.psam"
        File   log      = "~{output_prefix}.log"
        String log_line = read_string("log_line.txt")
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## SubsetUnrelatedPopulation
## Builds a keep-list of samples assigned to a given superpopulation that are
## NOT flagged as related (KING), then extracts that subset from a PLINK binary
## dataset. Sample IDs are unique (guaranteed by HandleDuplicates), so related
## individuals are excluded by IID.
##   assignments  : AncestryPCA TSV (FID IID superpop probability ...)
##   related_list : KING .id file of related samples (header lines start with #)
## -----------------------------------------------------------------------------
task SubsetUnrelatedPopulation {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   assignments
        File   related_list
        String population       # superpopulation label, e.g. "EUR"
        String output_prefix
        String plink2_bin
    }
    command <<<
        set -euo pipefail
        POP="~{population}"

        # IIDs assigned to the requested superpopulation.
        awk -F'\t' -v pop="$POP" '
            NR==1 { for (i=1;i<=NF;i++) h[$i]=i;
                    c    = ("superpop" in h) ? h["superpop"] : \
                           (("predicted_pop" in h) ? h["predicted_pop"] : h["predicted"]);
                    iidc = ("IID" in h) ? h["IID"] : 2; next }
            $c == pop { print $iidc }' ~{assignments} > pop_iids.txt

        # Related IIDs to exclude (skip header lines; IID is the last column).
        grep -v '^#' ~{related_list} 2>/dev/null | awk '{print $NF}' > related_iids.txt || true

        # Unrelated population IIDs = population IIDs minus related IIDs.
        awk 'NR==FNR{rel[$1]=1; next} !($1 in rel){print $1}' \
            related_iids.txt pop_iids.txt > keep_iids.txt

        # Resolve FID/IID from the genotype .fam.
        awk 'NR==FNR{keep[$1]=1; next} ($2 in keep){print $1, $2}' \
            keep_iids.txt ~{fam_file} > ~{output_prefix}_keep.txt

        n=$(wc -l < ~{output_prefix}_keep.txt)
        if [ "$n" -eq 0 ]; then
            echo "ERROR: no unrelated $POP samples found (check ancestry labels / relatedness)" >&2
            exit 1
        fi
        echo "Unrelated $POP samples retained: $n"

        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --keep ~{output_prefix}_keep.txt \
            --make-bed \
            --out ~{output_prefix}
    >>>
    output {
        File out_bed   = "~{output_prefix}.bed"
        File out_bim   = "~{output_prefix}.bim"
        File out_fam   = "~{output_prefix}.fam"
        File keep_list = "~{output_prefix}_keep.txt"
        File log       = "~{output_prefix}.log"
    }
    runtime { maxRetries: 1 }
}

## -----------------------------------------------------------------------------
## CovariatePCA
## Computes within-cohort principal components on the final (imputation-prepared)
## dataset for use as association-analysis covariates.
##
## Design:
##   • PC axes are estimated on LD-pruned autosomal SNPs using UNRELATED samples
##     only (related individuals — flagged by the KING check — are removed so
##     they do not distort the axes). Variant weights are saved.
##   • ALL samples (related included) are then projected onto those axes, so
##     every participant receives covariates on a common scale.
##
## This differs from AncestryPCA: those PCs are projected onto 1000 Genomes to
## infer ancestry, whereas these are within-cohort PCs capturing residual
## structure to adjust for in association testing.
##
## Outputs:
##   eigenvec    — PC scores of the unrelated samples used to define the axes
##   eigenval    — eigenvalues (variance explained)
##   covariates  — FID IID PC1..PCn for ALL samples (plink2 --covar compatible)
## -----------------------------------------------------------------------------
task CovariatePCA {
    input {
        File   bed_file
        File   bim_file
        File   fam_file
        File   related_list    # Samples to exclude from axis estimation (KING flagged)
        File   ld_regions      # High-LD regions to exclude during pruning
        Int    ld_window_kb
        Int    ld_step
        Float  ld_r2
        Int    n_pcs
        String output_prefix
        String plink_bin
        String plink2_bin
    }
    command <<<
        set -euo pipefail

        # (1) LD prune autosomal SNPs, excluding high-LD regions.
        ~{plink_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --chr 1-22 \
            --exclude ~{ld_regions} --range \
            --indep-pairwise ~{ld_window_kb} kb ~{ld_step} ~{ld_r2} \
            --out cov_prune

        # (2) Decide whether any related samples need removing from axis estimation.
        # The KING .id file carries a header (lines starting with '#'); count the
        # data rows to know whether --remove is needed.
        n_related=$(grep -vc '^#' ~{related_list} || true)
        if [ "${n_related:-0}" -gt 0 ]; then
            REMOVE_ARG="--remove ~{related_list}"
        else
            REMOVE_ARG=""
        fi

        # (3) Reference allele frequencies from the unrelated samples, used both
        # for the PCA and to standardise the projection consistently.
        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            ${REMOVE_ARG} \
            --extract cov_prune.prune.in \
            --freq \
            --out cov_ref_freq

        # (4) PCA on the unrelated samples; save variant weights for projection.
        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            ${REMOVE_ARG} \
            --extract cov_prune.prune.in \
            --read-freq cov_ref_freq.afreq \
            --pca ~{n_pcs} biallelic-var-wts \
            --out ~{output_prefix}

        # (5) Project ALL samples (related included) onto the unrelated axes.
        # eigenvec.var: col 2 = variant ID, col 4 = scored allele, cols 5.. = PCs.
        score_end=$(( ~{n_pcs} + 4 ))
        ~{plink2_bin} \
            --bed ~{bed_file} \
            --bim ~{bim_file} \
            --fam ~{fam_file} \
            --extract cov_prune.prune.in \
            --read-freq cov_ref_freq.afreq \
            --score ~{output_prefix}.eigenvec.var 2 4 header no-mean-imputation variance-standardize \
            --score-col-nums 5-${score_end} \
            --out ~{output_prefix}_projected

        # (6) Build the covariate table: FID, IID, PC1..PCn for all samples.
        # The projected .sscore ends with the n PC score columns; take the first
        # two columns (FID, IID) and the last n columns.
        awk -v n=~{n_pcs} '
            NR==1 {
                printf "FID\tIID";
                for (i=1;i<=n;i++) printf "\tPC%d", i;
                printf "\n";
                next
            }
            {
                printf "%s\t%s", $1, $2;
                for (i=NF-n+1; i<=NF; i++) printf "\t%s", $i;
                printf "\n"
            }' ~{output_prefix}_projected.sscore > ~{output_prefix}_covariates.tsv

        n_prune=$(wc -l < cov_prune.prune.in)
        printf "%-42s  %d PCs on %d LD-pruned SNPs (unrelated), projected to all\n" \
            "Step 12  Covariate PCA" "~{n_pcs}" "$n_prune" > cov_log.txt
    >>>
    output {
        File   eigenvec     = "~{output_prefix}.eigenvec"
        File   eigenval     = "~{output_prefix}.eigenval"
        File   eigenvec_var = "~{output_prefix}.eigenvec.var"
        File   covariates   = "~{output_prefix}_covariates.tsv"
        String log_line      = read_string("cov_log.txt")
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
