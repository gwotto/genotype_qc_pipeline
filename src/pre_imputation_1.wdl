version 1.0

workflow PreImputationTOPMed {

  input {
    File raw_vcf
    String input_build                # e.g. "GRCh37" or "GRCh38"
    File reference_fasta
    File reference_fasta_index
    File chain_file                   # required if liftover needed
    Boolean do_liftover = true
  }

  call ValidateVCF {
    input:
      vcf = raw_vcf
  }

  if (do_liftover) {
    call LiftoverVCF {
      input:
        vcf = ValidateVCF.validated_vcf,
        chain = chain_file,
        reference = reference_fasta
    }
  }

  File qc_input_vcf = select_first([
    LiftoverVCF.lifted_vcf,
    ValidateVCF.validated_vcf
  ])

  call VariantQC {
    input:
      vcf = qc_input_vcf
  }

  call SampleQC {
    input:
      vcf = VariantQC.filtered_vcf
  }

  call PrepareTOPMed {
    input:
      vcf = SampleQC.qc_vcf
  }

  output {
    File topmed_ready_vcf = PrepareTOPMed.final_vcf
    File topmed_index     = PrepareTOPMed.final_index
    File qc_report        = SampleQC.qc_report
  }
}



task ValidateVCF {
  input {
    File vcf
  }

  command {
    bcftools view ${vcf} -Oz -o validated.vcf.gz
    tabix -p vcf validated.vcf.gz
  }

  output {
    File validated_vcf = "validated.vcf.gz"
    File validated_index = "validated.vcf.gz.tbi"
  }

  runtime {
    docker: "biocontainers/bcftools:v1.17"
  }
}


task LiftoverVCF {
  input {
    File vcf
    File chain
    File reference
  }

  command {
    picard LiftoverVcf \
      I=${vcf} \
      O=lifted.vcf.gz \
      CHAIN=${chain} \
      REFERENCE_SEQUENCE=${reference} \
      REJECT=rejected.vcf.gz

    tabix -p vcf lifted.vcf.gz
  }

  output {
    File lifted_vcf = "lifted.vcf.gz"
    File rejected_vcf = "rejected.vcf.gz"
  }

  runtime {
    docker: "broadinstitute/picard:latest"
  }
}


task VariantQC {
  input {
    File vcf
  }

  command {
    bcftools filter \
      -e 'F_MISSING > 0.05 || INFO/MAF < 0.01' \
      ${vcf} -Oz -o variants_qc.vcf.gz

    tabix -p vcf variants_qc.vcf.gz
  }

  output {
    File filtered_vcf = "variants_qc.vcf.gz"
  }

  runtime {
    docker: "biocontainers/bcftools:v1.17"
  }
}


task SampleQC {
  input {
    File vcf
  }

  command {
    plink2 \
      --vcf ${vcf} \
      --mind 0.05 \
      --make-bed \
      --out qc

    plink2 \
      --bfile qc \
      --recode vcf bgz \
      --out qc_final

    tabix -p vcf qc_final.vcf.gz
  }

  output {
    File qc_vcf = "qc_final.vcf.gz"
    File qc_report = "qc.log"
  }

  runtime {
    docker: "biocontainers/plink2:v2.00a3"
  }
}


task PrepareTOPMed {
  input {
    File vcf
  }

  command {
    bcftools norm \
      -m -any \
      ${vcf} -Oz -o topmed_ready.vcf.gz

    tabix -p vcf topmed_ready.vcf.gz
  }

  output {
    File final_vcf = "topmed_ready.vcf.gz"
    File final_index = "topmed_ready.vcf.gz.tbi"
  }

  runtime {
    docker: "biocontainers/bcftools:v1.17"
  }
}
