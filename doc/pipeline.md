# Genotype QC Pre-Imputation Pipeline

## Overview

This WDL workflow performs comprehensive quality control (QC) on genotype data before imputation. It implements a multi-stage filtering approach that balances data retention with stringency, including ancestry-aware filtering and relatedness detection.

**Current Version:** 2.0 (Updated June 2026)

## Pipeline Steps

### Step 0: Data Preparation
- **ConvertToBinary** (optional): Converts PLINK text format (ped/map) to binary format (bed/bim/fam)
- **HandleDuplicates**: Resolves duplicate sample IDs and prepares data for analysis

### Step 1: SNP Missingness Filter (--geno)
- Removes SNPs with ≥X% missing genotypes across samples
- Threshold: `geno_threshold` (default: 0.03 = 3%)
- **Rationale**: Unreliable genotype calls

### Step 2: Sample Missingness Filter (--mind)
- Removes samples with ≥X% missing genotypes
- Threshold: `mind_threshold` (default: 0.05 = 5%)
- **Rationale**: Low-quality samples

### **Threshold Sweep Visualization** (Informational)
- Runs PLINK with multiple thresholds (0.01, 0.02, 0.05, 0.10, 0.20) for both geno and mind
- Produces plots showing SNP and sample retention across thresholds
- **Purpose**: Helps validate that chosen thresholds represent reasonable trade-offs
- **Output**: `threshold_plot.png`

### Step 3: Sex Check
- Validates reported vs. inferred sex using X-chromosome heterozygosity
- **Conditional**: Only runs if X-chromosome SNPs are present
- Removes samples with inconsistent sex assignments
- **Output**: Sex check report and list of failed samples (if any)

### Step 4: Ancestry PCA Against 1000 Genomes
- Projects study samples onto 1000 Genomes Phase 3 principal components
- Assigns superpopulation labels (EUR, AFR, EAS, SAS, AMR) to each sample
- **Inputs Required**:
  - 1000G reference panel (bed/bim/fam format)
  - 1000G sample metadata with superpopulation labels
- **Outputs**:
  - Per-sample superpopulation assignments: `ancestry_assignments.tsv`
  - PCA eigenvectors and eigenvalues
  - PCA plot with study samples overlaid on 1000G

### Step 5: MAC Filter (Full Cohort)
- Removes monomorphic SNPs and variants with minor allele count < X
- Threshold: `mac_threshold` (default: 50)
- **Rationale**: Rare variants tested on full cohort (not ancestry subsets) because rarity is a cohort-level property
- **Output**: MAC distribution plot with allele frequencies (not counts)

### Step 6: Hardy-Weinberg Equilibrium Filter (Ancestry Subset Detection, Full Cohort Removal)
- **Detection**: Tests HWE separately within each user-specified ancestry group
- **Rationale**: Allele frequencies differ by ancestry; testing on subset avoids false HWE violations in admixed samples
- **Removal**: Removes failed SNPs from ENTIRE cohort (prevents artifacts in any ancestry group)
- Threshold: `hwe_pvalue` (default: 1e-6)
- **Inputs**: 
  - `ancestry_populations`: Comma-separated list (e.g., "EUR,AFR") or "ALL" for entire cohort
  - Ancestry assignments from Step 4

### Step 7a: LD Pruning (Helper for Heterozygosity)
- Generates independent SNP list using specified LD threshold
- Used only as input for heterozygosity detection (Step 7b)
- Parameters: `ld_r2` (default: 0.1), `ld_window_kb` (default: 50), `ld_step` (default: 5)

### Step 7b: Heterozygosity Outlier Removal (Ancestry Subset Detection, Subset Removal)
- **Detection**: Identifies samples with heterozygosity >3 SD from mean within ancestry group
- **Removal**: Removes outlier SAMPLES (not variants) only from tested ancestry group
- **Merge**: Combines het-filtered ancestry group with untested populations
- **Rationale**: Heterozygosity is sample-level; removing from subset avoids discarding samples from untested groups
- **Outputs**: 
  - Per-sample heterozygosity rates
  - Heterozygosity plots and outlier list

### Step 8: Chromosome Filter
- Restricts analysis to specified chromosomes (default: 1-23, autosomes + X)
- Parameter: `chr_args` (e.g., "--chr 1-22" for autosomes only)

### Step 9: Relatedness Check (KING)
- Detects cryptic relatedness using KING algorithm
- **Important**: Samples are FLAGGED, not automatically removed
- Threshold: `king_cutoff` (default: 0.0884 ≈ 3rd-degree relative)
- Threshold: `pihat_min` (default: 0.2 for IBD pairs)
- **Outputs**:
  - Relatedness report: `king_cutoff_related_samples.txt` (samples to optionally remove)
  - Detailed IBD pairs: `pihat.genome`
- **User Decision**: Users can manually remove related samples using this file in downstream PLINK association tests

### Step 10: Prepare for Imputation
- Harmonizes SNPs with imputation reference (HRC or TOPMed)
- Checks allele strand, ref/alt assignments
- Creates VCFs per chromosome (GRCh37/GRCh38 compatible)
- Outputs ready for imputation server submission

## Configuration File

### Location
`config/genotype_qc_preimputation_inputs.json`

### Input Variables (35 total)

#### Data Input Paths
- **ped_file**: Path to PLINK PED file (text genotypes) - use if data is in text format
- **map_file**: Path to PLINK MAP file - use with ped_file
- **bed_file**: Path to PLINK BED file (binary) - use if data is already binary
- **bim_file**: Path to PLINK BIM file - use with bed_file
- **fam_file**: Path to PLINK FAM file - use with bed_file
- **Note**: Provide EITHER (ped_file + map_file) OR (bed_file + bim_file + fam_file)

#### Output Settings
- **output_prefix**: Prefix for all output files (default: "array_qc")

#### QC Thresholds
- **geno_threshold**: Max SNP missingness (0.0-1.0, default: 0.03 = 3%)
- **mind_threshold**: Max sample missingness (0.0-1.0, default: 0.05 = 5%)
- **hwe_pvalue**: Min HWE p-value to retain SNP (default: 1e-6)
- **mac_threshold**: Min minor allele count (default: 50)
- **pihat_min**: Min IBD coefficient for relatedness (default: 0.2)
- **king_cutoff**: KING kinship cutoff (default: 0.0884 for 3rd degree)

#### LD Pruning Parameters
- **ld_r2**: Max r² for LD pruning (default: 0.1)
- **ld_window_kb**: Sliding window size in kb (default: 50)
- **ld_step**: Window step size in SNPs (default: 5)
- **ld_regions**: Path to file with high-LD regions to exclude (BED format)

#### Ancestry & Population Selection
- **ancestry_populations**: Population(s) for ancestry-aware filtering
  - Format: Comma-separated list (e.g., "EUR,AFR") or "ALL"
  - Default: "ALL"
  - **Important**: Values must match superpopulation labels in 1000G assignments (EUR, AFR, EAS, SAS, AMR)
  - **Effect**: MAC and HWE detection will only use specified populations

#### Reference Files (1000 Genomes)
- **ref_1kg_bed/bim/fam**: Path to 1000G Phase 3 reference panel (binary PLINK format)
  - Should be LD-pruned and filtered to biallelic SNPs only
- **ref_1kg_psam**: 1000G sample metadata file with SuperPop labels (EUR, AFR, EAS, SAS, AMR)
- **n_pcs**: Number of principal components to compute (default: 10)

#### Imputation Reference
- **hrc_ref_freq**: Path to HRC or 1000G allele frequency file (typically .gz)
  - Example: HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
- **check_bim_pl**: Path to Will Rayner's HRC-1000G-check-bim.pl script

#### Chromosome Selection
- **chr_args**: PLINK chromosome filter arguments (default: "--chr 1-23")
  - Examples: "--chr 1-22" (autosomes only), "--chr 1-22 X" (auto + X)

#### Script Paths
- **handle_duplicates_r**: Path to R script for duplicate handling
- **check_heterozygosity_r**: Path to R script for heterozygosity computation
- **heterozygosity_outliers_r**: Path to R script for outlier detection
- **ancestry_pca_plot_r**: Path to R script for PCA visualization
- **mac_plot_r**: Path to R script for MAC distribution plot
- **threshold_plot_r**: Path to R script for threshold sweep visualization

#### Tool Paths
- **plink_bin**: Path to PLINK 1.9 binary (or command name if in PATH)
- **plink2_bin**: Path to PLINK 2.0 binary (for KING relatedness)
- **rscript_bin**: Path to Rscript binary (or command name if in PATH)
- **perl_bin**: Path to Perl interpreter (or command name if in PATH)
- **bcftools_bin**: Path to bcftools binary (or command name if in PATH)
- **bgzip_bin**: Path to bgzip binary (or command name if in PATH)
- **tabix_bin**: Path to tabix binary (or command name if in PATH)

## Running the Pipeline

### Step 1: Prepare Input Data
```bash
# Ensure you have PLINK binary files (bed/bim/fam) or text files (ped/map)
ls -lh your_data.*
```

### Step 2: Update Configuration
```bash
# Edit config/genotype_qc_preimputation_inputs.json
# - Set input file paths
# - Update tool paths if needed
# - Adjust QC thresholds if desired
# - Specify ancestry_populations (e.g., "EUR" or "EUR,AFR")
```

### Step 3: Prepare Reference Files
```bash
# Download 1000 Genomes reference panel
# Ensure HRC/TOPMed frequency file is available
# Verify script paths point to correct R and Perl scripts
```

### Step 4: Run Workflow with Cromwell
```bash
java -jar ~/apps/cromwell-92.jar run \
  src/genotype_qc_preimputation.wdl \
  -i config/genotype_qc_preimputation_inputs.json
```

### Step 5: Collect Results
```bash
# After workflow completes:
./src/collect_cromwell_results.bash
# Results will be in results/ directory
```

## Output Files

### QC Dataset
- **final_bed/bim/fam**: Final QC'd dataset ready for imputation
- **final.vcf.gz**: Final dataset in VCF format (bgzipped)

### QC Reports
- **threshold_plot.png**: SNP and sample retention across thresholds
- **mac_plot.png**: Allele frequency distribution before/after filtering
- **ancestry_assignments.tsv**: Per-sample superpopulation assignments
- **relatedness_flagged_samples.txt**: Samples flagged as related (for user decision)
- **sex_check_report**: Sex consistency report
- **heterozygosity plots and outlier list**

### Imputation-Ready VCFs
- **per-chromosome VCFs**: `chr{1-23}.vcf.gz` (ready for imputation server)

## Important Notes

### Ancestry-Aware Filtering Strategy
- **MAC filter**: Full cohort (rarity is global property)
- **HWE filter**: Detected on ancestry subset, removed from full cohort
  - Rationale: Avoids false HWE violations in admixed samples
- **Heterozygosity filter**: Detected on ancestry subset, removed from subset only
  - Rationale: Preserves samples from untested populations
- **Relatedness**: Detected on full cohort, samples flagged (not removed)
  - Rationale: User controls removal decisions in downstream analysis

### Handling Related Samples
The pipeline **flags** related samples but does NOT automatically remove them. To use the relatedness information:

1. View flagged samples in `relatedness_flagged_samples.txt`
2. Decide which pairs to break based on your study design
3. Use PLINK's `--remove` flag in association tests

### Using Ancestry Assignments Downstream
After QC, use the ancestry assignments file to:
1. Filter samples in PLINK association tests: `plink --keep ancestry_assignments.tsv --filter-cases`
2. Perform ancestry-stratified analyses
3. Adjust for ancestry in statistical models

## Troubleshooting

### Empty Ancestry Keep-List
**Error**: "No samples found for populations: EUR"
**Cause**: Population names in config don't match 1000G assignments
**Fix**: Check `ancestry_assignments.tsv` for correct population labels; edit config

### Missing X Chromosome
**Info**: Sex check step is skipped
**Cause**: Input data has no X chromosome SNPs
**Fix**: This is normal; proceed with analysis (sex check optional)

### Relatedness Check Takes Too Long
**Cause**: KING kinship computation is expensive on large datasets
**Solution**: Consider increasing `pihat_min` threshold or running on subset

## References

- PLINK 1.9: https://www.cog-genomics.org/plink/1.9/
- PLINK 2.0: https://www.cog-genomics.org/plink/2.0/
- KING relatedness: https://www.kingrelatedness.com/
- Will Rayner's HRC-1000G-check-bim: http://www.well.ox.ac.uk/~wrayner/
- 1000 Genomes Project: https://www.internationalgenome.org/
