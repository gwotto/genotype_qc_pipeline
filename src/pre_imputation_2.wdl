version 1.0

workflow SNP_QC_Pipeline {
    input {
        File bim_file
        File bed_file
        File fam_file
        File liftOver_chain
        File hrc_reference_tab  # Required for HRC-1000G-check-bim.pl
        String sample_id
    }

    # 1. Convert BIM to BED format for LiftOver
    call CreateLiftOverBed {
        input: bim = bim_file, name = sample_id
    }

    # 2. Perform LiftOver (hg18 -> hg19)
    call RunLiftOver {
        input: 
            input_bed = CreateLiftOverBed.bed_for_lift,
            chain_file = liftOver_chain,
            name = sample_id
    }

    # 3. Update BIM coordinates and filter out unmapped variants
    call UpdateBimCoordinates {
        input:
            old_bim = bim_file,
            lifted_bed = RunLiftOver.lifted_bed,
            name = sample_id
    }

    # 4. Pre-Imputation QC (Missingness, Sex Check, HWE, MAF, Heterozygosity)
    call QualityControl {
        input:
            bed = bed_file,
            fam = fam_file,
            bim = UpdateBimCoordinates.new_bim,
            name = sample_id
    }

    output {
        File cleaned_bed = QualityControl.qc_bed
        File cleaned_bim = QualityControl.qc_bim
        File cleaned_fam = QualityControl.qc_fam
        File unmapped_snps = RunLiftOver.unmapped
    }
}

task CreateLiftOverBed {
    input { 
        File bim 
        String name
        }
    command <<<
        # BIM is: Chr, ID, GD, Pos. LiftOver BED needs: Chr, Start, End, ID
        # Note: Handled 1-based to 0-based conversion (Pos-1)
        awk '{print "chr"$1, $4-1, $4, $2}' ~{bim} > ~{name}.for_lift.bed
    >>>
    output { File bed_for_lift = "for_lift.bed" }
    # runtime { docker: "ubuntu:latest" }
    runtime {cpu: 4
    memory: "8G"
  }
}

task RunLiftOver {
    input {
        File input_bed
        File chain_file
        String name
    }
    command <<<
        liftOver ~{input_bed} ~{chain_file} ~{name}_hg19.bed unmapped.txt
    >>>
    output {
        File lifted_bed = "~{name}_hg19.bed"
        File unmapped = "unmapped.txt"
    }
    runtime { docker: "biocontainers/ucsc-liftover:v357_cv1" }
}

task UpdateBimCoordinates {
    input {
        File old_bim
        File lifted_bed
        String name
    }
    command <<<
        # Join lifted coordinates back to the BIM file
        # This replaces the missing liftmap.py script
        awk 'FNR==NR{a[$4]=$3; next} ($2 in a){$4=a[$2]; print $0}' ~{lifted_bed} ~{old_bim} > ~{name}_hg19.bim
    >>>
    output { File new_bim = "~{name}_hg19.bim" }
    runtime { docker: "ubuntu:latest" }
}

task QualityControl {
    input {
        File bed
        File bim
        File fam
        String name
    }
    command <<<
        # 1. Missingness (Sample > 2%, SNP > 3%)
        plink --bfile ~{name}_hg19 --mind 0.02 --geno 0.03 --make-bed --out step1
        
        # 2. Sex Check (F < 0.2 female, F > 0.8 male)
        plink --bfile step1 --check-sex 0.2 0.8 --out sex_info
        
        # 3. HWE (P < 1x10-6) and MAF (< 1%)
        plink --bfile step1 --hwe 1e-6 --maf 0.01 --make-bed --out step3
        
        # 4. Heterozygosity (> 3 SD)
        plink --bfile step3 --het --out het_info
        awk '{if ($6 <= -3 || $6 >= 3) print $1, $2}' het_info.het > outliers.txt
        plink --bfile step3 --remove outliers.txt --make-bed --out ~{name}_cleaned
    >>>
    output {
        File qc_bed = "~{name}_cleaned.bed"
        File qc_bim = "~{name}_cleaned.bim"
        File qc_fam = "~{name}_cleaned.fam"
    }
    runtime { docker: "biocontainers/plink1.9:v1.90b6.21-210416_cv1" }
}
