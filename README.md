## Genotyping Array Data QC and Imputation Pipelines

A list of pipelines for quality control and imputation of genotyping array data in `doc/genotype-qc-pipeline.md`

## Genotyping Array Data QC and Imputation Workflow

A draft of a pipeline for quality control and imputation of genotyping array data (work in progess) in `src/`.

The pipeline follows best practices in the field and utilizes established tools and software to ensure high-quality results.

### Resources used

1. Liftover tool: https://hgdownload.soe.ucsc.edu/admin/exe/
2. Liftover chain files: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
3. Pre-imputation QC tools: https://www.well.ox.ac.uk/~wrayner/tools/
4. HRC reference sites for quality control of variants: ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
5. PLINK: https://www.cog-genomics.org/plink/2.0/
6. Michigan Imputation Server: https://imputationserver.sph.umich.edu/index.html


# Data and scripts used

1. high-LD-regions.txt
  - Regions of high LD, from Anderson et al (2010),doi: 10.1038/nprot.2010.116. Link to dataset: https://static-content.springer.com/esm/art%3A10.1038%2Fnprot.2010.116/MediaObjects/41596_2010_BFnprot2010116_MOESM396_ESM.zip
  - Regions of high LD, from Syed et al (2025), doi: 10.1101/2025.11.25.690541. Link to dataset:
    https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt
  - Regions of high LD,   https://dougspeed.com/high-ld-regions/
2. Hetrozygosity check R scripts: `check_heterozygosity_rate.R` and `heterozygosity_outliers_list.R` from GWAS tutorial Marees et al (2018), doi: 10.1002/mpr.1608
3. Script by Will Rayner for pre-imputation checks (strand alignment, allele frequency checks, etc.), version 4.3: https://www.chg.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
4. HRC r1.1 GRCh37 reference frequency file. Download: ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz


There is a download helper script `src/download_resources.sh` that downloads the high-LD-regions and HRC reference frequency files to the current working directory:

```bash
cd ref/  # or wherever you want the files
bash ../src/download_resources.sh
```

 
# Run pipeline

```{bash}

# Run the genotyping QC and imputation pipeline using Cromwell
java -jar cromwell-85.jar run src/genotype_qc_imputation.wdl -i inputs.json

```

# TODO

- pruning and pcr before and after relationship check.

# License
This repository is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.
