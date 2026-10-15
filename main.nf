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
 */

include { VCFSamples; GenotypePCA; MakeGenotypeBed } from './modules/genotypes.nf'
include { PrepExpression } from './modules/expression.nf'
include { TensorQTLCis; TensorQTLNominal } from './modules/tensorqtl.nf'
include { SummarizeCis } from './modules/summarize.nf'

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
}
