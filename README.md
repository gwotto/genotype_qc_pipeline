## Genotyping Array Data QC and Imputation Pipelines

This repository contains WDL workflows for quality control and imputation of genotyping array data.

## Main Pipeline: Genotype QC Pre-Imputation

A comprehensive workflow for QC of genotyping array data before imputation (fully implemented in `src/`).

The pipeline follows best practices in the field and utilizes established tools and software to ensure high-quality results.

**Documentation**: See `doc/pipeline.md` for detailed step-by-step documentation and configuration options.

### Key Features

- **Ancestry-aware filtering**: HWE detection on ancestry subsets, removal from full cohort
- **Heterozygosity outlier removal**: Detected on ancestry subsets, removed from tested population only
- **Relatedness detection**: KING-based relatedness check with user-controlled removal
- **Ancestry assignments**: Per-sample superpopulation labels (EUR, AFR, EAS, SAS, AMR)
- **Threshold visualization**: Plots showing SNP/sample retention across filtering thresholds
- **Imputation-ready output**: Per-chromosome VCFs ready for imputation servers

### Quick Start

```bash
# 1. Update configuration file
nano config/genotype_qc_preimputation_inputs.json

# 2. Run the workflow
java -jar cromwell-92.jar run src/genotype_qc_preimputation.wdl \
  -i config/genotype_qc_preimputation_inputs.json

# 3. Collect results
./src/collect_cromwell_results.bash
```

### Configuration

All pipeline parameters are defined in `config/genotype_qc_preimputation_inputs.json` (35 configuration variables).

Key parameters:
- `ancestry_populations`: Specify which populations to test (e.g., "EUR" or "EUR,AFR" or "ALL")
- `geno_threshold`: SNP missingness threshold (default: 0.03)
- `mind_threshold`: Sample missingness threshold (default: 0.05)
- `mac_threshold`: Minor allele count threshold (default: 50)
- `hwe_pvalue`: Hardy-Weinberg equilibrium p-value threshold (default: 1e-6)

See `doc/pipeline.md` for complete parameter documentation.

### Output Files

- **final.bed/bim/fam**: QC'd PLINK binary files
- **chr1-23.vcf.gz**: Per-chromosome VCFs ready for imputation
- **ancestry_assignments.tsv**: Per-sample superpopulation assignments
- **relatedness_flagged_samples.txt**: Samples flagged as related (for user decision)
- **threshold_plot.png**: SNP and sample retention across filtering thresholds
- **mac_plot.png**: Allele frequency distribution

See `doc/pipeline.md` for complete output file descriptions.

## Resources and Tools Used

1. **PLINK 1.9 and 2.0**: https://www.cog-genomics.org/plink/
2. **KING relatedness**: https://www.kingrelatedness.com/
3. **Will Rayner's HRC-1000G-check-bim**: https://www.well.ox.ac.uk/~wrayner/tools/
4. **1000 Genomes Project Phase 3**: https://www.internationalgenome.org/
5. **High-LD regions**:
   - Anderson et al (2010), doi: 10.1038/nprot.2010.116
   - Syed et al (2025), doi: 10.1101/2025.11.25.690541
   - https://dougspeed.com/high-ld-regions/
6. **Heterozygosity scripts**: Marees et al (2018), doi: 10.1002/mpr.1608
7. **HRC reference frequencies**: ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/

## Helper Script

A download helper script `src/download_resources.sh` downloads the high-LD-regions and HRC reference frequency files.


## File Structure

- `src/` — WDL workflow and R/Perl scripts
- `doc/pipeline.md` — Detailed pipeline documentation
- `config/` — Configuration template and example inputs
- `ref/` — Reference files (1000 Genomes, HRC frequencies, high-LD regions)

## Citation

If you use this pipeline, please cite the relevant tools:
- PLINK: Chang et al (2015), doi: 10.1186/s13742-015-0047-8
- KING: Manichaikul et al (2010), doi: 10.1093/bioinformatics/btq559
- 1000 Genomes: The 1000 Genomes Project Consortium (2015), doi: 10.1038/nature15393

## License

This repository is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.
