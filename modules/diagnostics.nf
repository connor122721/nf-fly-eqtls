#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  Gene annotation, expression diagnostics, Manhattan plots
// --------------------------------------------------

// FBgn -> gene symbol / CG ID from FlyBase fbgn_annotation_ID table
process PrepGeneMap {

    label 'small'
    container "${params.container_r}"
    publishDir "${params.out}/expression", mode: 'copy'

    input:
        path annotation
        path expression_bed
        path script

    output:
        path "${params.prefix}.fbgn_to_symbol.tsv", emit: map

    script:
        """
        Rscript ${script} \\
            --annotation ${annotation} \\
            --genes ${expression_bed} \\
            --out ${params.prefix}.fbgn_to_symbol.tsv
        """
}

// Expression QC / normalization / covariate diagnostics
process ExpressionDiagnostics {

    label 'small'
    container "${params.container_r}"
    publishDir "${params.out}/expression/diagnostics", mode: 'copy'

    input:
        path prep_rds
        path eigenvec
        path eigenval
        path gene_map
        path script

    output:
        path "${params.prefix}.*"

    script:
        """
        Rscript ${script} \\
            --prep ${prep_rds} \\
            --geno_pcs ${eigenvec} \\
            --gene_map ${gene_map} \\
            --group "${params.plot_group}" \\
            --out_prefix ${params.prefix}
        """
}

// Collapse nominal parquet files to one row per SNP for genome-wide Manhattan
process ReduceNominal {

    label 'small'
    container "${params.container_tensorqtl}"
    publishDir "${params.out}/summary", mode: 'copy'

    input:
        path parquet_files
        path script

    output:
        path "${params.prefix}.nominal_snp_min_p.txt.gz", emit: snp

    script:
        """
        python3 ${script} ${parquet_files} \\
            --keep_p ${params.manhattan_keep_p} \\
            --thin_frac ${params.manhattan_thin_frac} \\
            --out ${params.prefix}.nominal_snp_min_p.txt.gz
        """
}

// Manhattan + QQ plots for main and interaction models
process PlotManhattan {

    label 'small'
    container "${params.container_r}"
    publishDir "${params.out}/plots", mode: 'copy'

    input:
        path main_table
        path interaction_table
        path snp_table
        path chrom_file
        path script

    output:
        path "${params.prefix}.manhattan_qq.pdf"

    script:
        def int_arg = interaction_table.name != 'NO_INTERACTION' ? "--interaction ${interaction_table}" : ''
        """
        Rscript ${script} \\
            --main ${main_table} \\
            ${int_arg} \\
            --snp ${snp_table} \\
            --chromosomes ${chrom_file} \\
            --fdr ${params.fdr} \\
            --n_label ${params.n_label} \\
            --out_prefix ${params.prefix}
        """
}
