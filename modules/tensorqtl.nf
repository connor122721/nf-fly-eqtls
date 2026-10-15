#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  tensorQTL cis mapping (per chromosome x covariate model)
// --------------------------------------------------

// cis-eQTL with permutations (eGene discovery)
process TensorQTLCis {

    label 'tensorqtl'
    container "${params.container_tensorqtl}"
    publishDir "${params.out}/tensorqtl_cis", mode: 'copy'

    input:
        tuple val(chromosome), path(bed), path(bim), path(fam), val(k), path(covariates)
        path expression_bed

    output:
        tuple val(k), path("${params.prefix}.chr${chromosome}.K${k}.cis_qtl.txt.gz"), emit: cis

    script:
        """
        # Phenotypes on this chromosome only
        awk -v c="${chromosome}" 'NR==1 || \$1==c' ${expression_bed} > pheno.${chromosome}.bed

        python3 -m tensorqtl \\
            ${bed.baseName} \\
            pheno.${chromosome}.bed \\
            ${params.prefix}.chr${chromosome}.K${k} \\
            --covariates ${covariates} \\
            --mode cis \\
            --window ${params.cis_window} \\
            --maf_threshold ${params.maf} \\
            --permutations ${params.permutations} \\
            --seed ${params.seed}
        """
}

// cis-eQTL nominal associations (all variant-gene pairs) for the best model
process TensorQTLNominal {

    label 'tensorqtl'
    container "${params.container_tensorqtl}"
    publishDir "${params.out}/tensorqtl_nominal", mode: 'copy'

    input:
        tuple val(chromosome), path(bed), path(bim), path(fam), val(k), path(covariates)
        path expression_bed

    output:
        path "${params.prefix}.chr${chromosome}.K${k}.cis_qtl_pairs.*.parquet"

    script:
        """
        awk -v c="${chromosome}" 'NR==1 || \$1==c' ${expression_bed} > pheno.${chromosome}.bed

        python3 -m tensorqtl \\
            ${bed.baseName} \\
            pheno.${chromosome}.bed \\
            ${params.prefix}.chr${chromosome}.K${k} \\
            --covariates ${covariates} \\
            --mode cis_nominal \\
            --window ${params.cis_window} \\
            --maf_threshold ${params.maf}
        """
}

// SNP x interaction (e.g. condition) model, lead association per gene with eigenMT
// expr ~ g + i + g*i + covariates (covariates exclude the interaction variable)
process TensorQTLInteraction {

    label 'tensorqtl'
    container "${params.container_tensorqtl}"
    publishDir "${params.out}/tensorqtl_interaction", mode: 'copy'

    input:
        tuple val(chromosome), path(bed), path(bim), path(fam), val(k), path(covariates)
        path expression_bed
        path interaction

    output:
        path "${params.prefix}.chr${chromosome}.K${k}.interaction.cis_qtl_top_assoc.txt.gz", emit: top
        path "${params.prefix}.chr${chromosome}.K${k}.interaction.cis_qtl_pairs.*.parquet", optional: true

    script:
        def all_pairs = params.interaction_all_pairs ? '' : '--best_only'
        """
        awk -v c="${chromosome}" 'NR==1 || \$1==c' ${expression_bed} > pheno.${chromosome}.bed

        python3 -m tensorqtl \\
            ${bed.baseName} \\
            pheno.${chromosome}.bed \\
            ${params.prefix}.chr${chromosome}.K${k}.interaction \\
            --covariates ${covariates} \\
            --interaction ${interaction} \\
            --mode cis_nominal \\
            --window ${params.cis_window} \\
            --maf_threshold ${params.maf} \\
            --maf_threshold_interaction ${params.maf_threshold_interaction} \\
            ${all_pairs}
        """
}
