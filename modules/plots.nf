#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  eQTL plots for top gene-SNP pairs
// --------------------------------------------------

// Pool interaction results and pick top pairs from both models
process SelectTopPairs {

    label 'small'
    container "${params.container_r}"
    publishDir "${params.out}/summary", mode: 'copy'

    input:
        path main_table
        path interaction_files
        path script

    output:
        path "${params.prefix}.plot_pairs.tsv", emit: pairs
        path "plot_variants.txt", emit: variants
        path "${params.prefix}.interaction.*", optional: true

    script:
        """
        ls *.interaction.cis_qtl_top_assoc.txt.gz > interaction.list 2>/dev/null || true

        Rscript ${script} \\
            --main ${main_table} \\
            --interaction_list interaction.list \\
            --n_plot ${params.n_plot} \\
            --fdr ${params.fdr} \\
            --out_prefix ${params.prefix}
        """
}

// Genotype dosages for the plotted variants (all chromosomes)
process ExtractPlotGenotypes {

    label 'plink'
    container "${params.container_wgs}"

    input:
        path variants
        path plink_files

    output:
        path "plot_genotypes.raw", emit: raw
        path "plot_variants.bim", emit: bim

    script:
        """
        for b in *.bed; do
            p=\$(basename \$b .bed)
            # Skip chromosomes with none of the plotted variants
            if ! awk 'NR==FNR {v[\$1]; next} \$2 in v {f=1; exit} END {exit !f}' ${variants} \$p.bim; then
                continue
            fi
            plink2 \\
                --bfile \$p \\
                --chr-set ${params.plink_chr_set} \\
                --allow-extra-chr \\
                --extract ${variants} \\
                --export A \\
                --out geno.\$p
            awk 'NR==FNR {v[\$1]; next} \$2 in v' ${variants} \$p.bim >> plot_variants.bim
        done

        # Join per-chromosome .raw files by sample (same sample order in all beds)
        files=(geno.*.raw)
        cp \${files[0]} plot_genotypes.raw
        for f in "\${files[@]:1}"; do
            cut -f7- \$f | paste plot_genotypes.raw - > tmp.raw && mv tmp.raw plot_genotypes.raw
        done
        """
}

// Genotype vs normalized expression plots
process PlotEQTLs {

    label 'small'
    container "${params.container_r}"
    publishDir "${params.out}/plots", mode: 'copy'

    input:
        path pairs
        path raw
        path bim
        path expression_bed
        path metadata
        path script

    output:
        path "${params.prefix}.eqtl_plots.*.pdf"

    script:
        """
        Rscript ${script} \\
            --pairs ${pairs} \\
            --genotypes ${raw} \\
            --bim ${bim} \\
            --expression ${expression_bed} \\
            --metadata ${metadata} \\
            --sample_col "${params.sample_col}" \\
            --group "${params.plot_group}" \\
            --out_prefix ${params.prefix}
        """
}
