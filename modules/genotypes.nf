#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  Genotype processing: sample list, PCA, PLINK beds
// --------------------------------------------------

// List samples in the VCF
process VCFSamples {

    label 'small'
    container "${params.container_wgs}"
    publishDir "${params.out}/genotypes", mode: 'copy'

    input:
        path vcf
        path vcf_index

    output:
        path "vcf_samples.txt", emit: samples

    script:
        """
        bcftools query -l ${vcf} > vcf_samples.txt
        """
}

// Genome-wide genotype PCA on LD-pruned autosomal SNPs
process GenotypePCA {

    label 'plink'
    container "${params.container_wgs}"
    publishDir "${params.out}/genotypes", mode: 'copy'

    input:
        path vcf
        path vcf_index
        path metadata

    output:
        path "${params.prefix}.geno_pca.eigenvec", emit: eigenvec
        path "${params.prefix}.geno_pca.eigenval"
        path "${params.prefix}.geno_pca.log"

    script:
        """
        # Samples with metadata (first column)
        tail -n +2 ${metadata} | cut -f1 > keep.txt

        # LD prune
        plink2 \\
            --vcf ${vcf} \\
            --vcf-half-call missing \\
            --allow-extra-chr \\
            --chr ${params.pca_chroms} \\
            --keep keep.txt \\
            --snps-only just-acgt \\
            --max-alleles 2 \\
            --set-all-var-ids @:# \\
            --rm-dup force-first \\
            --geno ${params.pca_geno} \\
            --maf ${params.pca_maf} \\
            --indep-pairwise 50 5 0.2 \\
            --threads ${task.cpus} \\
            --memory ${task.memory.toMega() - 1000} \\
            --out prune

        # PCA on pruned SNPs
        plink2 \\
            --vcf ${vcf} \\
            --vcf-half-call missing \\
            --allow-extra-chr \\
            --chr ${params.pca_chroms} \\
            --keep keep.txt \\
            --snps-only just-acgt \\
            --max-alleles 2 \\
            --set-all-var-ids @:# \\
            --rm-dup force-first \\
            --extract prune.prune.in \\
            --pca ${params.n_pca} \\
            --threads ${task.cpus} \\
            --memory ${task.memory.toMega() - 1000} \\
            --out ${params.prefix}.geno_pca
        """
}

// PLINK bed per chromosome for QC-passing samples
process MakeGenotypeBed {

    label 'plink'
    container "${params.container_wgs}"
    publishDir "${params.out}/genotypes/bed", mode: 'copy'

    input:
        val chromosome
        path vcf
        path vcf_index
        path samples

    output:
        tuple val(chromosome),
              path("${params.prefix}.${chromosome}.bed"),
              path("${params.prefix}.${chromosome}.bim"),
              path("${params.prefix}.${chromosome}.fam"), emit: plink

    script:
        """
        # Treat all samples as diploid on X (avoids haploid handling of unknown sex)
        printf "#IID\\tSEX\\n" > sex.txt
        awk -v OFS='\\t' '{print \$1, 2}' ${samples} >> sex.txt

        plink2 \\
            --vcf ${vcf} \\
            --vcf-half-call missing \\
            --allow-extra-chr \\
            --chr ${chromosome} \\
            --keep ${samples} \\
            --update-sex sex.txt \\
            --snps-only just-acgt \\
            --max-alleles 2 \\
            --set-all-var-ids @:# \\
            --rm-dup force-first \\
            --geno ${params.geno_missing} \\
            --maf ${params.maf} \\
            --make-bed \\
            --threads ${task.cpus} \\
            --memory ${task.memory.toMega() - 1000} \\
            --out ${params.prefix}.${chromosome}
        """
}
