#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  RNA QC, normalization, covariates
// --------------------------------------------------

process PrepExpression {

    label 'small'
    container "${params.container_r}"
    publishDir "${params.out}/expression", mode: 'copy'

    input:
        path counts
        path metadata
        path gtf
        path chrom_file
        path vcf_samples
        path geno_pcs
        path script

    output:
        path "${params.prefix}.expression.bed", emit: bed
        path "${params.prefix}.covariates.K*.txt", emit: covariates
        path "${params.prefix}.samples.txt", emit: sample_list
        path "${params.prefix}.interaction.txt", emit: interaction, optional: true
        path "${params.prefix}.covariates_int.K*.txt", emit: covariates_int, optional: true
        path "${params.prefix}.sample_qc.tsv"
        path "${params.prefix}.expr_pc_variance.tsv"
        path "${params.prefix}.gene_r2_known_covariates.tsv"
        path "${params.prefix}.expression_qc.pdf"
        path "${params.prefix}.prep_objects.rds", emit: rds
        path "prep_expression.log"

    script:
        def iv = params.interaction
        def interaction_var = (iv == null || iv instanceof Boolean || iv.toString().trim() in ['', 'false', 'null']) ? '' : iv.toString()
        """
        Rscript ${script} \\
            --counts ${counts} \\
            --metadata ${metadata} \\
            --gtf ${gtf} \\
            --chromosomes ${chrom_file} \\
            --genotyped_samples ${vcf_samples} \\
            --geno_pcs ${geno_pcs} \\
            --n_geno_pcs ${params.n_geno_pcs} \\
            --expr_pcs "${params.expr_pcs}" \\
            --sample_col "${params.sample_col}" \\
            --formula "${params.covariate_formula}" \\
            --factor_covariates "${params.factor_covariates}" \\
            --cpm_threshold ${params.cpm_threshold} \\
            --count_threshold ${params.count_threshold} \\
            --sample_frac ${params.sample_frac} \\
            --min_lib_size ${params.min_lib_size} \\
            --lib_mad ${params.lib_mad} \\
            --pca_outlier_p ${params.pca_outlier_p} \\
            --norm_method ${params.norm_method} \\
            --interaction "${interaction_var}" \\
            --out_prefix ${params.prefix} 2>&1 | tee prep_expression.log
        """
}
