# GWAS Genotype QC Pipeline

A WDL workflow that cleans raw genotype array data before genome-wide association analysis. It combines published workflows and follows best practices in the field combining PLINK commands with custom R scripts for diagnostics and sample filtering. The pipeline is designed to be modular and adaptable to different datasets and QC thresholds. It follows a standard set of QC steps combining methods from Anderson et al. (2010) and Marees et al. (2018), as well from the pipeline developed for the DIVERGE study by the Kuchenbaaecker lab at UCL (https://github.com/DIVERGEstudy). 

---

## What this pipeline does

Raw SNP array data contains errors from genotyping noise, sample mixups,
and population structure. This pipeline removes those problems in a fixed
order, producing a dataset suitable for GWAS.

```
Input genotypes (bed/bim/fam or ped/map)
    │
    ├─ Step 0  : Convert to binary format (if needed)
    ├─ Step 0b : Remove duplicate sample IDs
    ├─ Step 1  : Remove low-quality SNPs      --geno (missingness per SNP)
    ├─ Step 2  : Remove low-quality samples   --mind (missingness per sample)
    ├─ Step 3  : Sex check                    flag sample swaps / genotyping errors
    ├─ Step 4  : MAF filter                   remove monomorphic & rare SNPs
    ├─ Step 5  : Hardy-Weinberg filter        remove SNPs that fail equilibrium test
    ├─ Step 6  : LD pruning                   select independent SNPs (used in steps 7 & 9)
    ├─ Step 7  : Heterozygosity filter        remove contaminated / inbred samples
    ├─ Step 8  : Restrict to autosomes        chromosomes 1–22 only
    ├─ Step 9  : Relatedness filter           remove one of each related pair
    ├─ Step 10 : PCA                          diagnostic population structure check
    └─ Step 11 : Prepare for imputation       strand-align and export per-chr VCFs
                                                    │
                                              Final QC dataset
```

The pipeline graph:

![genotype_qc_preimputation pipeline graph](genotype_qc_preimputation.png "genotype_qc_preimputation pipeline graph")


A flow diagram of the pipeline:

![genotype_qc_preimputation pipeline flow diagram](genotype_qc_preimputation_mm.png "genotype_qc_preimputation pipeline flow diagram")

---

## Step-by-step explanation

### Step 0 — Format conversion *(optional)*
If you supply text-format `.ped`/`.map` files, PLINK converts them to binary
`.bed`/`.bim`/`.fam`. Binary format is faster for every subsequent step.
Skipped if you supply binary files directly.

### Step 0b — Duplicate sample IDs
PLINK requires unique sample identifiers. An R script scans the `.fam` file
and appends a numeric suffix to any duplicated IDs.

### Step 1 — SNP missingness (`--geno`)
SNPs with too many missing genotype calls across samples are removed.
A typical threshold is 5% missing (`geno_threshold = 0.05`).
Highly missing SNPs usually reflect a failed assay probe.

### Step 2 — Sample missingness (`--mind`)
Samples with too many missing genotype calls across SNPs are removed.
A typical threshold is also 5% (`mind_threshold = 0.05`).
Highly missing samples reflect poor DNA quality or plate failures.

### Step 3 — Sex check
PLINK uses X-chromosome heterozygosity to infer biological sex and compares
it to the reported sex in the `.fam` file. Discordant samples are removed —
they likely represent sample swaps or labelling errors.  
Skipped automatically if the dataset has no X-chromosome SNPs.

### Step 4 — MAF filter
Removes monomorphic SNPs (MAF = 0) and SNPs below `maf_threshold`.
Monomorphic SNPs carry no association signal and can cause numerical issues
in downstream analyses. A typical pre-imputation threshold is 1%
(`maf_threshold = 0.01`). MAF filtering should be repeated after imputation
on the imputed dataset.

This step is placed before HWE testing and LD pruning for two reasons:
- HWE tests are unreliable for rare variants due to low expected counts;
  applying the MAF filter first means HWE is tested only where it has power.
- LD pruning and heterozygosity checks are more stable when based on common
  variants only.

An R script produces a MAF distribution plot and a before/after barplot for
quality inspection. The frequency files (`maf_freq_before`, `maf_freq_after`)
are QC diagnostics; the allele frequency file used for imputation strand
alignment (step 11) is computed separately on the final post-relatedness
dataset.

### Step 5 — Hardy-Weinberg equilibrium (HWE)
SNPs that deviate significantly from HWE are removed. Extreme deviation
(typical threshold p < 1×10⁻⁶, `hwe_pvalue = 1e-6`) usually indicates
genotyping error rather than true selection. Applied after the MAF filter
(step 4) so that HWE is only tested on common variants, where the test has
reliable power.

### Step 6 — LD pruning
Generates a list of approximately independent SNPs by removing pairs in high
linkage disequilibrium (correlated because they are physically close on the
chromosome). This pruned SNP set is not used to filter the dataset here —
it is only used as input to steps 7 and 9 to ensure those analyses are not
biased by correlated markers. The same list is useful for PCA.

| Parameter | Meaning | Typical value |
|---|---|---|
| `ld_window_kb` | Sliding window size | 50 kb |
| `ld_step` | Step size (SNPs) | 5 |
| `ld_r2` | r² threshold | 0.2 |

High-LD genomic regions (e.g. the MHC locus on chr 6, the chr 8 inversion)
are excluded before pruning so they do not dominate the independent SNP set.
The exclusion list is supplied via the `ld_regions` input file. The file
bundled with this pipeline (`ref/high-LD-regions.txt`) is the hg19/GRCh37
region list from the supplementary data of Anderson et al. (2010) at https://static-content.springer.com/esm/art%3A10.1038%2Fnprot.2010.116/MediaObjects/41596_2010_BFnprot2010116_MOESM396_ESM.zip. If your data are
aligned to a different genome build, supply the appropriate coordinates
instead.

### Step 7 — Heterozygosity check
Using the LD-pruned SNPs, per-sample heterozygosity rates are computed.
Outliers beyond ±3 SD of the cohort mean are flagged and removed.

- High heterozygosity: possible sample contamination (DNA from two individuals)
- Low heterozygosity: possible inbreeding or an accidental sample duplicate

The R scripts for this step (`check_heterozygosity_rate.R` and
`heterozygosity_outliers_list.R`) are taken from the practical GWAS tutorial
by Marees et al. (2018).

### Step 8 — Autosome filter
Retains only chromosomes 1–22. Sex chromosomes and mitochondrial variants
require specialised handling not included in this standard GWAS QC workflow.

### Step 9 — Relatedness filtering *(optional)*
Cryptically related individuals (e.g. undisclosed family members, sample
duplicates) violate the independence assumption in standard GWAS. Two methods are run on the LD-pruned SNPs:

| Method | Tool | Metric | Use |
|---|---|---|---|
| pi-hat | PLINK 1.9 `--genome` | Proportion of alleles identical by descent | Informational report |
| KING kinship | PLINK2 `--king-cutoff` | Kinship coefficient | Decides which samples to remove |

KING is preferred for removal decisions because it is robust to population
stratification. A cutoff of 0.0884 corresponds approximately to 3rd-degree
relatives (cousins).

Set `run_relatedness_check = false` in the configuration json file to skip this step entirely — for example, when working with a cohort containing family trios. When skipped, the autosome-filtered dataset from step 8 passes directly to PCA.

### Step 10 — PCA *(diagnostic only)*
Principal components are computed on the final QC-passed dataset using the
LD-pruned SNP list from step 6 (no additional pruning needed). Scatter plots
of PC1–PC3 are produced to visualise population structure and identify
potential ancestry outliers before imputation.

No samples are removed at this step. PCA should be repeated after imputation
to generate the PCs used as covariates in the GWAS model.

### Step 11 — Prepare for imputation
Aligns strand orientation to the TOPMed/HRC reference panel using Will
Rayner's `HRC-1000G-check-bim.pl` script, then exports one bgzipped,
tabix-indexed VCF per autosome for upload to the TOPMed imputation server
(https://imputation.biodatacatalyst.nhlbi.nih.gov).

The check-bim script removes:
- SNPs absent from the reference panel
- Ambiguous A/T and C/G SNPs that cannot be strand-resolved
- SNPs with allele frequency difference > 0.2 vs the reference
- SNPs with mismatched positions or alleles

---

## Inputs

### Genome build

It is assumed that the input dataset is already aligned to the same genome build as the HRC reference panel (GRCh37). If this is not the case, the pipeline does not perform any liftover itself. Instead it corrects variant positions and alleles to match the HRC reference panel in step 11, which is the critical step for imputation. The `HRC-1000G-check-bim.pl` script will remove any SNPs that cannot be aligned to the reference, so it is important to ensure the input dataset is as close as possible to the reference build to avoid excessive SNP loss at this step. The ld_regions file used in step 6 (LD pruning) also needs to be on the same build as the input data, to remove the correct high-LD regions before pruning. 

If your data are on GRCh38, you can use the HRC reference files lifted over to GRCh38 by the Michigan Imputation Server (https://imputationserver.sph.umich.edu/index.html#!pages/download/reference). 



### Required files

| Input | Description |
|---|---|
| `bed_file`, `bim_file`, `fam_file` | PLINK binary genotype files **(or)** |
| `ped_file`, `map_file` | PLINK text genotype files |
| `ld_regions` | High-LD genomic regions to exclude from LD pruning. The file bundled with this pipeline (`ref/high-LD-regions.txt`) lists hg19/GRCh37 coordinates from Anderson et al. (2010). Replace with hg38 coordinates if your data use GRCh38. |
| `handle_duplicates_r` | R script to resolve duplicate sample IDs |
| `check_heterozygosity_r` | R script to compute heterozygosity rates |
| `heterozygosity_outliers_r` | R script to flag heterozygosity outliers |
| `threshold_plot_r` | R script to plot SNP/sample counts across missingness thresholds |
| `maf_plot_r` | R script to plot MAF distribution |
| `pca_plot_r` | R script to plot PCA results |
| `check_bim_pl` | Will Rayner's `HRC-1000G-check-bim.pl` script |
| `hrc_ref_freq` | HRC r1.1 GRCh37 reference frequency file |

### Tool paths

By default the pipeline expects `plink`, `plink2`, `Rscript`, `perl`, `bcftools`, `bgzip`, and `tabix` to be on `$PATH`. If they are not, override these inputs in the JSON:

| Input | Default | Description |
|---|---|---|
| `plink_bin` | `"plink"` | Path to PLINK 1.9 binary |
| `plink2_bin` | `"plink2"` | Path to PLINK2 binary |
| `rscript_bin` | `"Rscript"` | Path to Rscript binary |
| `perl_bin` | `"perl"` | Path to Perl binary |
| `bcftools_bin` | `"bcftools"` | Path to bcftools binary |
| `bgzip_bin` | `"bgzip"` | Path to bgzip binary |
| `tabix_bin` | `"tabix"` | Path to tabix binary |


### QC thresholds

QC thresholds are defined in the json configuration file. The table lists values that are commonly used in similar pipelines.

| Parameter | Meaning | Typical value | Remarks |
|---|---|---|
| `geno_threshold` | Max SNP missingness rate | `0.05` | |
| `mind_threshold` | Max sample missingness rate | `0.05` | |
| `hwe_pvalue` | Min HWE p-value | `1e-6` | |
| `maf_threshold` | Min MAF to retain a SNP | `0.01` | |
| `ld_r2` | Max r² for LD pruning | `0.2` | |
| `ld_window_kb` | LD pruning window size (kb) | `50` | |
| `ld_step` | LD pruning step size (SNPs) | `5` | |
| `pihat_min` | Min pi-hat to report a related pair | `0.2` | Only used if `run_relatedness_check = true` |
| `king_cutoff` | KING kinship cutoff for removal | `0.0884` | Only used if `run_relatedness_check = true` |
| `n_pcs` | Number of principal components to compute | `10` | |

---

## Outputs

Final QC-passed files are copied to `results/` after the run (see
[How to run](#how-to-run) below). Intermediate per-step files remain in the
executor's working directory (`_miniwdl_run/` or `cromwell-executions/`) and
can be deleted once the run is verified.

| Output | Description | Remarks |
|---|---|
| `pipeline_log` | Human-readable run summary: version, date, SNP/sample counts at each QC step | |
| `final_bed/bim/fam` | QC-passed dataset — use this for downstream analyses | |
| `sexcheck_report` | Full PLINK sex check results (if X SNPs present) | |
| `problem_samples` | Samples removed at sex check step | |
| `het_check_report` | Per-sample heterozygosity rates | |
| `het_fail_samples` | Samples removed at heterozygosity step | |
| `maf_plot` | MAF distribution and before/after barplot | |
| `maf_freq_before` | Allele frequencies before MAF filter (QC diagnostic) | |
| `maf_freq_after` | Allele frequencies after MAF filter (QC diagnostic) | |
| `pihat_genome` | All related pairs above `pihat_min` (informational) | Only present if `run_relatedness_check = true`|
| `king_cutoff_out_id` | Samples removed at relatedness step | Only present if `run_relatedness_check = true` |
| `prune_in` | LD-pruned SNP list (useful for PCA) | |
| `prune_out` | SNPs excluded by LD pruning | |
| `pca_eigenvec` | PC scores per sample | |
| `pca_eigenval` | Variance explained per PC | |
| `pca_plot` | Scatter plots of PC1–PC3 | |
| `imputation_vcfs` | Per-chromosome bgzipped VCFs for imputation server upload | |
| `imputation_vcf_tbis` | Tabix indices for imputation VCFs | |
| `check_bim_log` | Log from `HRC-1000G-check-bim.pl` | |

---


## How to run

### 1. Install dependencies

Option 1: conda (recommended):

Create the environment from the provided `environment.yml` file in the
repository root. This installs all required tools (PLINK 1.9, PLINK2,
R ≥ 4.0, Perl, bcftools, bgzip, tabix) and Java in a single step:

```bash
conda env create -f environment.yml
conda activate genotype-qc
```

To update an existing environment after the file changes:

```bash
conda env update -f environment.yml --prune
```

R packages (ggplot2, dplyr, tidyr, readr, patchwork, scales) are included
in `environment.yml` via `r-*` conda packages. If you prefer to install them
from CRAN instead:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales"))
```

The Perl module `Term::ReadKey` (required by Will Rayner's check-bim script)
is not on conda; install it via cpan:

```bash
cpan App::cpanminus
cpanm Term::ReadKey
```

Cromwell is not included in the conda environment (it is a standalone JAR).
Download it from https://github.com/broadinstitute/cromwell/releases and
place it somewhere on your `$PATH`, or reference it by full path. Java ≥ 17
is required (provided by the `openjdk` package in the environment).

Option 2: UCL Myriad HPC modules:
```bash
module load plink/1.9
module load plink2/2.0
module load r/4.2.1
```

Option 3: Manual install:

- PLINK 1.9 and PLINK2: download binaries from
  https://www.cog-genomics.org/plink/ and place them on `$PATH`, or supply
  their paths via `plink_bin` / `plink2_bin` in the inputs JSON.
- R ≥ 4.0: https://cran.r-project.org/. Install required packages from CRAN:
  ```r
  install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales"))
  ```
- bcftools / bgzip / tabix: install via your system package manager, e.g.
  ```bash
  # Debian/Ubuntu
  sudo apt install bcftools tabix
  # macOS
  brew install bcftools htslib
  ```
  or download from https://www.htslib.org/download/.
- Perl and `Term::ReadKey` (required by the check-bim script):
  ```bash
  cpan App::cpanminus
  cpanm Term::ReadKey
  ```
- Java ≥ 17 (for Cromwell): https://adoptium.net/ or via your package manager.
- Cromwell: download the standalone JAR from
  https://github.com/broadinstitute/cromwell/releases and place it on `$PATH`
  or reference it by full path when running the pipeline.



### 2. Create an inputs JSON file

```json
{
  "genotype_qc_preimputation.bed_file":                  "data/myStudy.bed",
  "genotype_qc_preimputation.bim_file":                  "data/myStudy.bim",
  "genotype_qc_preimputation.fam_file":                  "data/myStudy.fam",
  "genotype_qc_preimputation.output_prefix":             "myStudy",
  "genotype_qc_preimputation.ld_regions":                "ref/high-LD-regions.txt",
  "genotype_qc_preimputation.handle_duplicates_r":       "scripts/handle_duplicates.R",
  "genotype_qc_preimputation.check_heterozygosity_r":    "scripts/check_heterozygosity_rate.R",
  "genotype_qc_preimputation.heterozygosity_outliers_r": "scripts/heterozygosity_outliers_list.R",
  "genotype_qc_preimputation.threshold_plot_r":          "scripts/threshold_plot.R",
  "genotype_qc_preimputation.maf_plot_r":                "scripts/maf_plot.R",
  "genotype_qc_preimputation.pca_plot_r":                "scripts/pca_plot.R",
  "genotype_qc_preimputation.check_bim_pl":              "ref/HRC-1000G-check-bim.pl",
  "genotype_qc_preimputation.hrc_ref_freq":              "ref/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz",
  "genotype_qc_preimputation.geno_threshold":            0.05,
  "genotype_qc_preimputation.mind_threshold":            0.05,
  "genotype_qc_preimputation.hwe_pvalue":                1e-6,
  "genotype_qc_preimputation.maf_threshold":             0.01,
  "genotype_qc_preimputation.ld_r2":                     0.2,
  "genotype_qc_preimputation.ld_window_kb":              50,
  "genotype_qc_preimputation.ld_step":                   5,
  "genotype_qc_preimputation.run_relatedness_check":     true,
  "genotype_qc_preimputation.pihat_min":                 0.2,
  "genotype_qc_preimputation.king_cutoff":               0.0884,
  "genotype_qc_preimputation.n_pcs":                     10
}
```

### 3. Run with miniwdl

miniwdl requires Docker to be running (on macOS, Docker Desktop must be
started first). If Docker is not available, use Cromwell instead — it runs
tasks directly on the host when no `docker:` runtime attribute is set, which
is the case for all tasks in this pipeline.

```bash
miniwdl run src/genotype_qc_preimputation.wdl -i inputs.json

# Copy final outputs to results/
RUNDIR=$(ls -td _miniwdl_run/*/genotype_qc_preimputation | head -1)
mkdir -p results
cp "$RUNDIR/out/final_bed"/*.bed \
   "$RUNDIR/out/final_bim"/*.bim \
   "$RUNDIR/out/final_fam"/*.fam \
   "$RUNDIR/out/pipeline_log"/*_pipeline.log \
   results/
```
### 4. Run with Cromwell

```bash
cromwell run src/genotype_qc_preimputation.wdl --inputs inputs.json
```

Once the run completes, collect all outputs into a timestamped subdirectory
of `results/`. The script below uses the Cromwell run UUID as the directory
name so each run's outputs are kept separate and can be traced back to the
exact execution.

```bash
# Root of the most recent Cromwell run
RUNID=$(ls -t cromwell-executions/genotype_qc_preimputation/ | head -1)
BASE="cromwell-executions/genotype_qc_preimputation/$RUNID"

# Destination: results/<uuid>/  (unique per run, traceable to Cromwell logs)
OUTDIR="results/$RUNID"
mkdir -p "$OUTDIR"/{qc_dataset,plots,reports,vcfs,pca,ld_pruning}

# ── Pipeline log (WriteLog) ───────────────────────────────────────────────────
CALL="$BASE/call-WriteLog/execution"
cp "$CALL"/*_pipeline.log "$OUTDIR/"

# ── Final QC dataset (step 9 — RemoveRelated) ────────────────────────────────
CALL="$BASE/call-RemoveRelated/execution"
cp "$CALL"/*.bed "$CALL"/*.bim "$CALL"/*.fam "$OUTDIR/qc_dataset/"

# ── Sex check report (step 3 — SexCheck / RemoveSexFails) ────────────────────
# Only present if the dataset contained X-chromosome SNPs
CALL="$BASE/call-SexCheck/execution"
[ -d "$CALL" ] && cp "$CALL"/*.sexcheck "$OUTDIR/reports/" 2>/dev/null || true

CALL="$BASE/call-RemoveSexFails/execution"
[ -d "$CALL" ] && cp "$CALL"/*.fam "$OUTDIR/reports/sex_fail_samples.fam" 2>/dev/null || true

# ── Heterozygosity report & fail list (step 7 — HetCheck / RemoveHetFails) ───
CALL="$BASE/call-HetCheck/execution"
cp "$CALL"/*.het "$OUTDIR/reports/"
cp "$CALL"/*.png "$OUTDIR/plots/" 2>/dev/null || true

CALL="$BASE/call-RemoveHetFails/execution"
cp "$CALL"/*.txt "$OUTDIR/reports/het_fail_samples.txt" 2>/dev/null || true

# ── MAF plot & frequency files (step 4 — MafFilter) ─────────────────────────
CALL="$BASE/call-MafFilter/execution"
cp "$CALL"/*.png "$OUTDIR/plots/"
cp "$CALL"/*.frq "$OUTDIR/reports/"

# ── Relatedness report (step 9 — RelatednessCheck) ──────────────────────────
CALL="$BASE/call-RelatednessCheck/execution"
cp "$CALL"/*.genome "$OUTDIR/reports/"
cp "$CALL"/*.king.cutoff.out.id "$OUTDIR/reports/" 2>/dev/null || true

# ── LD pruning SNP lists (step 6 — LdPruning) ────────────────────────────────
CALL="$BASE/call-LdPruning/execution"
cp "$CALL"/*.prune.in "$CALL"/*.prune.out "$OUTDIR/ld_pruning/"

# ── PCA (step 10 — PCA) ───────────────────────────────────────────────────────
CALL="$BASE/call-PCA/execution"
cp "$CALL"/*.eigenvec "$CALL"/*.eigenval "$OUTDIR/pca/"
cp "$CALL"/*.png "$OUTDIR/plots/" 2>/dev/null || true

# ── Imputation VCFs (step 11 — PrepareForImputation) ─────────────────────────
CALL="$BASE/call-PrepareForImputation/execution"
cp "$CALL"/*_chr*.vcf.gz     "$OUTDIR/vcfs/"
cp "$CALL"/*_chr*.vcf.gz.tbi "$OUTDIR/vcfs/"
cp "$CALL"/check-bim.log     "$OUTDIR/reports/"

echo "Done. Results written to: $OUTDIR"
ls "$OUTDIR"/qc_dataset/ "$OUTDIR"/plots/ "$OUTDIR"/reports/ \
   "$OUTDIR"/pca/ "$OUTDIR"/vcfs/ "$OUTDIR"/ld_pruning/
```

If you prefer a human-readable name over the UUID, you can rename the
directory afterwards:

```bash
mv "results/$RUNID" "results/genotype_qc_preimputation_$(date +%Y%m%d)"
```

After collecting outputs, verify sample and SNP counts decrease
monotonically through the pipeline:

```bash
# Samples at each step (.fam lines = number of samples)
wc -l "$BASE"/call-*/execution/*.fam | sort -k2

# SNPs at each step (.bim lines = number of SNPs)
wc -l "$BASE"/call-*/execution/*.bim | sort -k2
```

---

## Troubleshooting

Each PLINK step produces a `.log` file in the executor working directory. Check these first if a step fails:

```bash
# miniwdl — logs are alongside each task's outputs
find _miniwdl_run -name "*.log" | sort

# Cromwell
find cromwell-executions -name "*.log" | sort
```

Common failure causes:

| Error | Likely cause |
|---|---|
| `Error: No SNPs pass --geno` | Threshold too strict, or input data quality is very poor |
| `Error: No samples pass --mind` | As above for samples |
| Sex check produces no output | No X-chromosome SNPs in dataset (pipeline skips this automatically) |
| KING produces empty output | Too few SNPs in pruned set; try relaxing `ld_r2` |
| `check-bim` produces no VCFs | No SNPs passed reference alignment; check genome build and `hrc_ref_freq` |

---

## Expected resource usage

For a typical dataset (~500k SNPs, ~2,000 samples):

| Resource | Requirement |
|---|---|
| CPUs | 1 (pipeline is single-threaded) |
| RAM | 4–8 GB |
| Runtime | 15–30 minutes |
| Disk | ~3× the size of input files |

Larger datasets (e.g. whole-genome imputed data, >100k samples) will require
proportionally more RAM and runtime.

---

## Dependencies

All tools must be on `$PATH` or their paths supplied via the `*_bin` inputs
in the JSON (see [Tool paths](#tool-paths-optional) above).

| Tool | Version | Purpose |
|---|---|---|
| PLINK 1.9 | ≥ 1.9 | Most QC steps |
| PLINK2 | ≥ 2.0 | KING relatedness, PCA |
| R / Rscript | ≥ 4.0 | Heterozygosity, duplicate, MAF, PCA scripts |
| Perl | any recent | HRC check-bim script |
| bcftools / bgzip / tabix | any recent | VCF conversion for imputation |
| miniwdl or Cromwell | any recent | WDL executor |

---

## Reference

Anderson CA et al. (2010). Data quality control in genetic case-control
association studies. Nature Protocols 5, 1564–1573.  
https://doi.org/10.1038/nprot.2010.116

Marees AT et al. (2018). A tutorial on conducting genome-wide association
studies: Quality control and statistical analysis. International Journal of
Methods in Psychiatric Research 27, e1608.  
https://doi.org/10.1002/mpr.1608
