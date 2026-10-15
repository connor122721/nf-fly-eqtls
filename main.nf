#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 *  nf-fly-eqtls: cis-eQTL mapping in Drosophila (dm3) with tensorQTL
 *  Derived from nf-eqtls (Connor Murray), simplified for fly cage experiments.
 *
 *  1) List VCF samples and run genome-wide genotype PCA (LD-pruned autosomes)
 *  2) RNA QC + TMM/INT normalization + covariates (~ cage_name + condition [+ geno PCs] + K expression PCs)
 *  3) PLINK bed per chromosome for QC-passing samples
 *  4) tensorQTL cis (permutations) per chromosome x K
 *  5) Pool chromosomes, FDR, pick K with most eGenes
 *  6) tensorQTL cis_nominal per chromosome for best K
 *  7) SNP x condition interaction model (tensorQTL interaction, eigenMT) for best K
 *  8) Genotype vs normalized expression plots for top eGenes / interaction pairs
 */

include { VCFSamples; GenotypePCA; MakeGenotypeBed } from './modules/genotypes.nf'
include { PrepExpression } from './modules/expression.nf'
include { TensorQTLCis; TensorQTLNominal; TensorQTLInteraction } from './modules/tensorqtl.nf'
include { SummarizeCis } from './modules/summarize.nf'
include { SelectTopPairs; ExtractPlotGenotypes; PlotEQTLs } from './modules/plots.nf'

workflow {

    // Main chromosomes
    chroms = Channel
        .fromPath(params.chrom_file)
        .splitText()
        .map { it.trim() }
        .filter { it }

    vcf = file(params.vcf, checkIfExists: true)
    vcf_index = file("${params.vcf}.tbi", checkIfExists: true)

    // 1) Genotype sample list + PCA
    vcf_samples = VCFSamples(vcf, vcf_index)
    geno_pca = GenotypePCA(vcf, vcf_index, file(params.metadata, checkIfExists: true))

    // 2) Expression QC, normalization, covariates
    prep = PrepExpression(
        file(params.counts, checkIfExists: true),
        file(params.metadata, checkIfExists: true),
        file(params.gtf, checkIfExists: true),
        file(params.chrom_file, checkIfExists: true),
        vcf_samples.samples,
        geno_pca.eigenvec,
        file("${projectDir}/scripts/prep_expression.R"))

    // 3) PLINK beds for QC-passing samples
    plink = MakeGenotypeBed(chroms, vcf, vcf_index, prep.sample_list)

    // One covariate file per K expression PCs: [K, file]
    covs = prep.covariates
        .flatten()
        .map { f -> tuple((f.name =~ /\.K(\d+)\.txt$/)[0][1] as Integer, f) }

    // [chrom, bed, bim, fam, K, covariates]
    qtl_input = plink.combine(covs)

    // 4) cis-eQTL permutation pass for every chromosome x K
    TensorQTLCis(qtl_input, prep.bed)

    // 5) Summarize and choose best K
    summary = SummarizeCis(
        TensorQTLCis.out.cis.map { it[1] }.collect(),
        file("${projectDir}/scripts/summarize_cis.R"))

    best_k = params.best_k != null ?
        Channel.value(params.best_k as Integer) :
        summary.best_k.map { it.text.trim() as Integer }

    // 6) Nominal pass for best K
    nominal_input = qtl_input
        .combine(best_k)
        .filter { it[4] == it[6] }
        .map { it[0..5] }

    TensorQTLNominal(nominal_input, prep.bed)

    // 7) SNP x interaction model for best K (covariates without the interaction variable)
    if (params.interaction) {
        int_input = plink
            .combine(prep.covariates_int.flatten()
                .map { f -> tuple((f.name =~ /\.K(\d+)\.txt$/)[0][1] as Integer, f) })
            .combine(best_k)
            .filter { it[4] == it[6] }
            .map { it[0..5] }

        int_files = TensorQTLInteraction(int_input, prep.bed, prep.interaction).top.collect()
    } else {
        int_files = Channel.value([])
    }

    // 8) eQTL plots for top pairs
    top = SelectTopPairs(
        summary.table,
        int_files,
        file("${projectDir}/scripts/select_top_pairs.R"))

    plot_geno = ExtractPlotGenotypes(
        top.variants,
        plink.map { [it[1], it[2], it[3]] }.flatten().collect())

    PlotEQTLs(
        top.pairs,
        plot_geno.raw,
        plot_geno.bim,
        prep.bed,
        file(params.metadata),
        file("${projectDir}/scripts/plot_eqtls.R"))
}
