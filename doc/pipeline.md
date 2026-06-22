# Genotype QC Pre-Imputation Pipeline

**Version 0.2.0 — June 2026**

## Overview

`src/genotype_qc_preimputation.wdl` is a WDL 1.0 workflow that performs
standard pre-imputation quality control on genotyping array data. It follows
the Anderson et al. (2010) protocol with ancestry-aware extensions, and
produces a clean PLINK binary dataset and per-chromosome VCFs ready for
submission to imputation servers (HRC, TOPMed).

The pipeline is run with Cromwell:

```bash
java -jar cromwell-92.jar run src/genotype_qc_preimputation.wdl \
     -i config/genotype_qc_preimputation_inputs.json
```

Results are collected with:

```bash
bash src/collect_cromwell_results.bash <cromwell-run-uuid-dir> ./results
```
<!--
This is commented out

![Pipeline flowchart](genotype_qc_preimputation_mm.png)

![Pipeline graph](genotype_qc_preimputation.png)
-->

---

## Pipeline Steps

### Step 0 — Format conversion (optional)

If the input is in PLINK text format (ped/map), it is converted to binary
(bed/bim/fam). Skip this by supplying bed/bim/fam directly.

### Step 0b — Duplicate sample IDs

Duplicate IDs in the FAM file cause downstream PLINK errors. The R script
`handle_duplicates.R` appends a numeric suffix to make each ID unique and
writes a log of affected samples.

### Step 1 — SNP missingness filter (`--geno`)

Removes variants with a missing genotype rate above `geno_threshold`.
Run before sample filtering so that poorly genotyped SNPs do not inflate
per-sample missingness.

### Step 2 — Sample missingness filter (`--mind`)

Removes samples with a missing genotype rate above `mind_threshold`.

### Threshold sweep (informational)

Runs PLINK across a grid of geno/mind values (0.01, 0.02, 0.05, 0.10,
0.20) and plots variant and sample retention. Used to validate that the
chosen thresholds represent sensible trade-offs. No data are removed.

Output: `threshold_plot.png`

### Step 3 — Sex check

Infers sex from X-chromosome heterozygosity and compares it to the FAM
file. Samples where the reported and inferred sex disagree are removed.
Skipped automatically if the dataset has no X-chromosome SNPs. Samples
with sex code 0 ("unknown") are retained.

### Step 4 — Ancestry PCA

PCA is computed on the 1000 Genomes Phase 3 reference panel only, then
all study samples (including related pairs) are projected into that
reference PC space. This keeps the PC axes stable and independent of
study composition.

**Procedure:**

1. SNPs are matched between study and reference by rsID. A positional
   concordance check discards rsID-matched variants that sit at different
   chr:pos coordinates (catches genome-build mismatches).
2. If `pca_reference_populations` is not `"ALL"`, the reference panel is
   filtered to the listed superpopulations before PCA.
3. The filtered reference is LD-pruned (same parameters as step 7a).
4. PLINK2 `--pca biallelic-var-wts` computes `n_pcs` PCs on the reference
   and saves per-variant weights in `.eigenvec.var`.
5. Reference allele frequencies are saved (`--freq`) for use in step 6.
6. All study samples are projected with `--score … variance-standardize`.
   The raw `SCORE_AVG` values are rescaled to the eigenvec scale before
   classification: `eigenvec_k = SCORE_AVG_k × 2 / sqrt(eigenvalue_k)`.
7. A random forest (500 trees, `randomForest` R package) is trained on the
   reference eigenvec scores and their known superpopulation labels, then
   applied to all study samples. Each sample receives:
   - `superpop`: predicted superpopulation (`EUR`, `AFR`, `EAS`, `SAS`,
     `AMR`), or `"unassigned"` if the prediction probability is below
     `ancestry_prob_threshold`.
   - `probability`: maximum class probability from the random forest.

**Outputs:**

| File | Description |
|------|-------------|
| `*_ancestry_assignments.tsv` | FID, IID, superpop, probability |
| `*_ancestry_pca.png` | Grid of PC pair plots (PC1/2 … PC9/10), study samples coloured by predicted superpopulation |
| `*.eigenvec` | Reference-panel PC scores |
| `*.eigenvec.var` | Per-variant PC weights |
| `*.eigenval` | PC eigenvalues |
| `all_study_projected.sscore` | Raw projected scores for all study samples |

### Step 5 — MAC filter

Removes monomorphic variants and variants with minor allele count below
`mac_threshold`. Applied to the full cohort because rarity is a global
property of the dataset.

Output: `mac_plot.png`

### Step 6 — HWE filter

Hardy-Weinberg equilibrium is tested separately within each superpopulation
listed in `test_populations`. Variants failing HWE in any tested population
are removed from the entire cohort. Testing on a subset prevents false
violations caused by population stratification in admixed samples.

### Step 7a — LD pruning (helper)

Produces an independent SNP list using `ld_window_kb`, `ld_step`, and
`ld_r2`. Used only as input to the heterozygosity check; not applied to
the final dataset.

### Step 7b — Heterozygosity outlier removal

Per-sample heterozygosity (F-statistic) is computed on LD-pruned variants.
Samples more than 3 SD from the mean within the tested ancestry subset are
flagged and removed — but only from that subset. The removed samples are
not propagated to the rest of the cohort, preserving samples from
populations not covered by `test_populations`.

### Step 8 — Chromosome filter

Restricts the dataset to the chromosomes listed in `chr_args`. The output
of this step is the **final QC dataset** used for all downstream analysis.

### Step 9 — Relatedness check (KING)

KING kinship coefficients are computed on the full post-QC dataset. Pairs
with kinship above `king_cutoff` are flagged. Samples are **not removed**
automatically. The output file lists one member of each related pair and
can be supplied to PLINK `--remove` in downstream association tests.

Output: `relatedness_flagged_samples.txt`

### Step 10 — Prepare for imputation

Will Rayner's `HRC-1000G-check-bim.pl` harmonises the dataset against the
HRC or TOPMed frequency reference: checks strand, ref/alt assignment, and
frequency concordance. The script generates per-chromosome VCFs (bgzipped
and tabix-indexed) ready for submission.

---

## Configuration

All parameters live in `config/genotype_qc_preimputation_inputs.json`.
Every key is prefixed with `genotype_qc_preimputation.`.

### Input files

| Parameter | Description |
|-----------|-------------|
| `bed_file`, `bim_file`, `fam_file` | PLINK binary input. Supply these **or** ped/map, not both. |
| `ped_file`, `map_file` | PLINK text input (converted to binary automatically). |
| `output_prefix` | String prefix for all output files (e.g. `"array_qc"`). |

### Reference files

| Parameter | Description |
|-----------|-------------|
| `ld_regions` | High-LD genomic regions excluded from LD pruning. Tab-separated, four columns: chr start end name (hg19). Download with `bash src/download_resources.bash`. |
| `ref_1kg_bed`, `ref_1kg_bim`, `ref_1kg_fam` | 1000 Genomes Phase 3 reference panel in PLINK binary format. Must be autosomal, biallelic SNPs only. |
| `ref_1kg_psam` | 1000G sample metadata file. Must contain a `SuperPop` column with values `EUR`, `AFR`, `EAS`, `SAS`, `AMR`. |
| `hrc_ref_freq` | HRC r1.1 GRCh37 allele frequency file (`.tab.gz`). Download with `bash src/download_resources.bash`. |

### Script paths

All scripts in `src/`. Set these to absolute paths.

| Parameter | Script |
|-----------|--------|
| `handle_duplicates_r` | `src/handle_duplicates.R` |
| `check_heterozygosity_r` | `src/check_heterozygosity_rate.R` |
| `heterozygosity_outliers_r` | `src/heterozygosity_outliers_list.R` |
| `ancestry_pca_plot_r` | `src/ancestry_pca_plot.R` |
| `mac_plot_r` | `src/mac_plot.R` |
| `threshold_plot_r` | `src/threshold_plot.R` |
| `create_qc_table_r` | `src/create_qc_status_table.R` |
| `check_bim_pl` | `src/HRC-1000G-check-bim.pl` |

### Tool paths

Default to bare command names (assumes tools are on `$PATH`). Override
with absolute paths if needed.

| Parameter | Tool | Minimum version |
|-----------|------|-----------------|
| `plink_bin` | PLINK 1.9 | 1.9 |
| `plink2_bin` | PLINK 2.0 | 2.0 |
| `rscript_bin` | R | 4.0 |
| `perl_bin` | Perl | 5.x |
| `bcftools_bin` | bcftools | — |
| `bgzip_bin` | bgzip | — |
| `tabix_bin` | tabix | — |

### QC thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `geno_threshold` | `0.03` | Maximum per-SNP missing rate (0–1). |
| `mind_threshold` | `0.05` | Maximum per-sample missing rate (0–1). |
| `hwe_pvalue` | `1e-6` | Minimum HWE p-value to retain a variant. |
| `mac_threshold` | `50` | Minimum minor allele count across the full cohort. |
| `pihat_min` | `0.2` | Minimum IBD coefficient for flagging related pairs. |
| `king_cutoff` | `0.0884` | KING kinship cutoff (~3rd-degree relative). |

### LD pruning parameters

Used in two places: ancestry PCA (step 4) and heterozygosity check
(step 7a). Same parameters apply to both.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ld_window_kb` | `200` | Sliding window size in kilobases. |
| `ld_step` | `1` | Window step size in SNPs. |
| `ld_r2` | `0.1` | Maximum r² for a variant to be retained. |

### Ancestry PCA parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_pcs` | `10` | Number of principal components to compute and use for ancestry classification. |
| `pca_reference_populations` | `"ALL"` | Comma-separated list of 1000G superpopulations whose samples are included in the reference PCA. Possible values: any subset of `EUR,AFR,EAS,SAS,AMR`, or `"ALL"` to use all 2504 Phase 3 samples. Restricting to relevant populations can sharpen cluster separation for studies with limited ancestry diversity. |
| `ancestry_prob_threshold` | `0.5` | Random forest prediction probability below which a study sample is labelled `"unassigned"` rather than assigned to a superpopulation. Range 0–1; increase for stricter assignment. |

### Ancestry-aware filtering parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `test_populations` | `"ALL"` | Comma-separated list of superpopulations on which HWE and heterozygosity detection are run (e.g. `"EUR"`, `"EUR,AFR"`). Use `"ALL"` to test the full cohort without stratification. Values must match the `superpop` column of the ancestry assignments file. |

### Imputation and chromosome parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `chr_args` | `"--chr 1-23"` | PLINK chromosome filter applied at step 8. Use `"--chr 1-22"` for autosomes only or `"--autosome"` as an alternative. |

---

## Output structure

After running `bash src/collect_cromwell_results.bash <run-dir> ./results`:

```
results/
  qc_dataset/          Final PLINK files (bed/bim/fam) — output of step 8
  vcfs/                chr1–23.vcf.gz + .tbi — imputation-ready VCFs
  pca/                 eigenvec, eigenval, eigenvec.var, sscore,
                       ancestry_assignments.tsv, ancestry_pca.png
  plots/               threshold_plot.png, mac_plot.png, het plots
  reports/             sex check, het outlier list, relatedness flags,
                       HWE fail list, per-chr freq files, pipeline log
```

**Key files:**

| File | Description |
|------|-------------|
| `pca/ancestry_assignments.tsv` | Per-sample: FID, IID, superpop, probability. Use to select ancestry strata for downstream analysis. |
| `reports/relatedness_flagged_samples.txt` | One member of each related pair. Pass to `plink --remove` in association tests if needed. |
| `qc_dataset/*.bed/bim/fam` | Final QC'd genotype data. |
| `vcfs/chr*.vcf.gz` | Per-chromosome VCFs for imputation servers. |

---

## Reference files — how to obtain

```bash
# High-LD regions (hg19) and HRC frequency file
bash src/download_resources.bash

# 1000 Genomes Phase 3 reference panel
# Obtain a biallelic-SNP-only, hg19, PLINK binary version of the 1000G
# Phase 3 panel. The companion .psam file must contain a SuperPop column
# with values: EUR, AFR, EAS, SAS, AMR. This can be obtained from 
# https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
```

---

## Design notes

### Why reference-only PCA?

Computing PCA on the 1000G reference and projecting study samples keeps
the PC axes stable regardless of study size, population composition, or
the presence of related individuals. Related pairs would otherwise deflate
PC variance and distort ancestry clusters.

### Why SCORE_AVG needs rescaling

PLINK2 `--score` with `variance-standardize` outputs `SCORE_AVG =
SCORE_SUM / ALLELE_CT`, where `ALLELE_CT = 2 × (non-missing variants)`.
Because `SCORE_SUM` scales with the singular value of the genotype matrix
for each PC, the per-PC conversion is:

```
eigenvec_k = SCORE_AVG_k × 2 / sqrt(eigenvalue_k)
```

The `M_scored / M_pca` missingness factor cancels exactly in the division
by `ALLELE_CT`, so the formula is correct for any genotyping coverage.

### Ancestry-aware filtering strategy

| Filter | Applied to | Why |
|--------|-----------|-----|
| MAC | Full cohort | Rarity is a global property |
| HWE | Detected within `test_populations`, removed from full cohort | Avoids false violations from stratification in admixed samples |
| Heterozygosity | Detected and removed within `test_populations` only | Preserves samples from untested populations |
| Relatedness | Flagged on full cohort | User decides which samples to exclude |

---

## References

- Anderson CA et al. (2010) *Nat Protoc* 5:1564–1573. doi:10.1038/nprot.2010.116
- Chang CC et al. (2015) *GigaScience* 4:7. doi:10.1186/s13742-015-0047-8 (PLINK)
- Manichaikul A et al. (2010) *Bioinformatics* 26:2867–2873. doi:10.1093/bioinformatics/btq559 (KING)
- The 1000 Genomes Project Consortium (2015) *Nature* 526:68–74. doi:10.1038/nature15393
- Marees AT et al. (2018) *Int J Methods Psychiatr Res* 27:e1608. doi:10.1002/mpr.1608
