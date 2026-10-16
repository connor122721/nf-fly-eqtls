#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pool cis results across chromosomes, FDR, and pick best K
process SummarizeCis {

    label 'small'
    container "${params.container_r}"
    publishDir "${params.out}/summary", mode: 'copy'

    input:
        path cis_files
        path gene_map
        path script

    output:
        path "best_k", emit: best_k
        path "${params.prefix}.egenes_by_K.tsv"
        path "${params.prefix}.egenes_by_K.pdf"
        path "${params.prefix}.egenes_by_chrom.K*.tsv"
        path "${params.prefix}.egenes.K*.tsv"
        path "${params.prefix}.cis_eqtl.K*.all_chroms.txt.gz", emit: table

    script:
        """
        ls *.cis_qtl.txt.gz > cis.list

        Rscript ${script} \\
            --cis_list cis.list \\
            --fdr ${params.fdr} \\
            --gene_map ${gene_map} \\
            --out_prefix ${params.prefix}
        """
}
