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
