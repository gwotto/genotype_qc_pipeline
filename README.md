# Genotype QC Pre-Imputation Pipeline

WDL 1.0 workflow for quality control of genotyping array data before
imputation, following the Anderson et al. (2010) protocol with
ancestry-aware extensions.

**Full documentation:** `doc/pipeline.md`

---

## Quick start

### 1. Install dependencies

```bash
conda env create -f environment.yml
conda activate genotype_qc
```

Required tools (all provided by the conda environment): `plink` 1.9,
`plink2`, `Rscript` ≥ 4.0, `perl`, `bcftools`, `bgzip`, `tabix`, and
`openjdk` ≥ 11 (for Cromwell).

### 2. Download reference files

```bash
bash src/download_resources.bash
```

This downloads:
- `high-LD-regions-hg19-GRCh37.txt` — regions excluded from LD pruning
- `HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz` — HRC imputation reference frequencies

You also need a 1000 Genomes Phase 3 reference panel in PLINK binary
format (hg19, biallelic SNPs only) with an accompanying `.psam` file
containing a `SuperPop` column (`EUR`, `AFR`, `EAS`, `SAS`, `AMR`). These can be 
obtained from https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3

### 3. Configure

Edit `config/genotype_qc_preimputation_inputs.json`:

- Set input file paths (`bed_file`/`bim_file`/`fam_file` or
  `ped_file`/`map_file`).
- Set absolute paths for all scripts in `src/` and all reference files.
- Adjust QC thresholds if needed (defaults are suitable for most arrays).
- Set `pca_reference_populations` to control which 1000G superpopulations
  define the ancestry PCA space (e.g. `"EUR,AFR,EAS,SAS"` or `"ALL"`).
- Set `test_populations` to the superpopulation(s) on which HWE and
  heterozygosity outlier detection should run (e.g. `"EUR"` or `"ALL"`).

See `doc/pipeline.md` for the complete parameter reference.

### 4. Run

```bash
java -jar cromwell-92.jar run src/genotype_qc_preimputation.wdl \
     -i config/genotype_qc_preimputation_inputs.json
```

### 5. Collect results

```bash
bash src/collect_cromwell_results.bash <cromwell-run-uuid-dir> ./results
```

Outputs land in `results/qc_dataset/` (final PLINK files),
`results/vcfs/` (imputation-ready per-chromosome VCFs),
`results/pca/` (ancestry assignments and PCA plots), 
`results/subset/` (subset of unrelated EUR samples), and
`results/reports/` (sex check, het outliers, relatedness flags, log).

---

## Pipeline summary

| Step | What it does |
|------|-------------|
| 0 | Convert ped/map → bed/bim/fam (optional) |
| 0b | Deduplicate sample IDs |
| 1 | SNP missingness filter (`--geno`) |
| 2 | Sample missingness filter (`--mind`) |
| 3 | Sex check (X-chr heterozygosity) |
| 4 | Ancestry PCA: reference-only PCA on 1000G -> project all study samples -> random forest ancestry assignment |
| 5 | MAC filter and monomorphic SNP removal (full cohort) |
| 6 | HWE filter (detected within ancestry subset, removed from full cohort) |
| 7 | LD pruning + heterozygosity outlier removal (within ancestry subset) |
| 8 | Restrict to selected chromosomes |
| 9 | Relatedness check (KING; samples flagged, not removed) |
| 10 | Harmonise with HRC/TOPMed; produce per-chromosome VCFs |
| 11 | BED and PLINK 2 copy of the combined dataset |
| 12 | Within-cohort covariate PCA on the combined dataset |
| 13 | Unrelated single-ancestry subset |

---

## Key outputs

| File | Description |
|------|-------------|
| `results/qc_dataset/` | Final QC'd genotype data in bed and Plink 2 format |
| `results/reports/*_sample_qc_status.tsv` | reference for filtering e.g. on ancestry and relatednes |
| `results/pca/*_covariate_pca_covariates.tsv` | Plink format PCA covariates of the QC'ed dataset |
| `results/subset/` | Final QC'd genotype data, subset of unrelated EUR samples | 
| `results/subset/*_unrelated_EUR_pca_covariates.tsv`| PCA covariates for subset of unrelated EUR samples |
| `results/vcfs/chr*.vcf.gz` | Per-chromosome VCFs for imputation servers |

---

## Tools and references

- PLINK 1.9/2.0: https://www.cog-genomics.org/plink/
- KING relatedness: https://www.kingrelatedness.com/
- HRC-1000G-check-bim (Will Rayner): https://www.well.ox.ac.uk/~wrayner/tools/
- 1000 Genomes Phase 3: https://www.internationalgenome.org/
- Anderson CA et al. (2010) *Nat Protoc* doi:10.1038/nprot.2010.116
- Chang CC et al. (2015) *GigaScience* doi:10.1186/s13742-015-0047-8
- Manichaikul A et al. (2010) *Bioinformatics* doi:10.1093/bioinformatics/btq559

## License

BSD 3-Clause. See [LICENSE](LICENSE).
