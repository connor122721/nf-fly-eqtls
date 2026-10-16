#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 *  nf-fly-eqtls: cis-eQTL mapping in Drosophila (dm3) with tensorQTL
 *  Derived from nf-eqtls (Connor Murray), simplified for fly cage experiments.
 *
 *  1) List VCF samples and run genome-wide genotype PCA (LD-pruned autosomes)
 *  2) RNA QC + TMM/INT normalization + covariates (~ cage_name + condition [+ geno PCs] + K expression PCs)
 *  3) FlyBase FBgn -> gene symbol map; expression diagnostics
 *  4) PLINK bed per chromosome for QC-passing samples
 *  5) tensorQTL cis (permutations) per chromosome x K
 *  6) Pool chromosomes, FDR, pick K with most eGenes
 *  7) tensorQTL cis_nominal per chromosome for best K
 *  8) SNP x condition interaction model (tensorQTL interaction, eigenMT) for best K
 *  9) eQTL plots (genotype vs expression), Manhattan and QQ plots
 */

include { VCFSamples; GenotypePCA; MakeGenotypeBed } from './modules/genotypes.nf'
include { PrepExpression } from './modules/expression.nf'
include { TensorQTLCis; TensorQTLNominal; TensorQTLInteraction } from './modules/tensorqtl.nf'
include { SummarizeCis } from './modules/summarize.nf'
include { SelectTopPairs; ExtractPlotGenotypes; PlotEQTLs } from './modules/plots.nf'
include { PrepGeneMap; ExpressionDiagnostics; ReduceNominal; PlotManhattan } from './modules/diagnostics.nf'

// Optional string params: "", false, or a bare --flag (parsed as true) all mean "off"
def optParam(v) {
    (v == null || v instanceof Boolean || v.toString().trim() in ['', 'false', 'null']) ? '' : v.toString()
}

// [K, file] from <prefix>.covariates[_int].K<k>.txt
def withK(ch) {
    ch.flatten().map { f -> tuple((f.name =~ /\.K(\d+)\.txt$/)[0][1] as Integer, f) }
}

workflow {

    // Main chromosomes
    chroms = Channel
        .fromPath(params.chrom_file)
        .splitText()
        .map { it.trim() }
        .filter { it }

    vcf = file(params.vcf, checkIfExists: true)
    vcf_index = file("${params.vcf}.tbi", checkIfExists: true)
    no_int = file("${projectDir}/assets/NO_INTERACTION")

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

    // 3) Gene symbols + expression diagnostics
    gene_map = optParam(params.gene_annotation) ?
        PrepGeneMap(file(optParam(params.gene_annotation)), prep.bed, file("${projectDir}/scripts/prep_gene_map.R")).map :
        Channel.value(file("${projectDir}/assets/NO_GENE_MAP"))

    ExpressionDiagnostics(
        prep.rds,
        geno_pca.eigenvec,
        geno_pca.eigenval,
        gene_map,
        file("${projectDir}/scripts/expression_diagnostics.R"))

    // 4) PLINK beds for QC-passing samples
    plink = MakeGenotypeBed(chroms, vcf, vcf_index, prep.sample_list)

    // [chrom, bed, bim, fam, K, covariates]
    covs = withK(prep.covariates)
    qtl_input = plink.combine(covs)

    // 5) cis-eQTL permutation pass for every chromosome x K
    TensorQTLCis(qtl_input, prep.bed)

    // 6) Summarize and choose best K
    summary = SummarizeCis(
        TensorQTLCis.out.cis.map { it[1] }.collect(),
        gene_map,
        file("${projectDir}/scripts/summarize_cis.R"))

    best_k = params.best_k != null ?
        Channel.value(params.best_k as Integer) :
        summary.best_k.map { it.text.trim() as Integer }

    // 7) Nominal pass for best K
    nominal_input = qtl_input
        .combine(best_k)
        .filter { it[4] == it[6] }
        .map { it[0..5] }

    nominal = TensorQTLNominal(nominal_input, prep.bed)

    best_cov = covs.combine(best_k).filter { it[0] == it[2] }.map { it[1] }.first()

    // 8) SNP x interaction model for best K (covariates without the interaction variable)
    if (optParam(params.interaction)) {
        int_covs = withK(prep.covariates_int)
        int_input = plink
            .combine(int_covs)
            .combine(best_k)
            .filter { it[4] == it[6] }
            .map { it[0..5] }

        int_files = TensorQTLInteraction(int_input, prep.bed, prep.interaction).top.collect()
        best_cov_int = int_covs.combine(best_k).filter { it[0] == it[2] }.map { it[1] }.first()
    } else {
        int_files = Channel.value([])
        best_cov_int = Channel.value(no_int)
    }

    // 9a) eQTL plots for top pairs
    top = SelectTopPairs(
        summary.table,
        int_files,
        gene_map,
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
        best_cov,
        best_cov_int,
        file("${projectDir}/scripts/plot_eqtls.R"))

    // 9b) Manhattan + QQ
    snp = ReduceNominal(nominal.collect(), file("${projectDir}/scripts/reduce_nominal.py"))

    PlotManhattan(
        summary.table,
        top.int_table.ifEmpty(no_int),
        snp.snp,
        file(params.chrom_file),
        file("${projectDir}/scripts/plot_manhattan.R"))
}
