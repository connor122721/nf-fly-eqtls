#!/usr/bin/env Rscript
# Connor Murray
# Pool interaction (SNP x condition) results across chromosomes and pick the
# top gene-SNP pairs from the main and interaction models for plotting.

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--main", required=TRUE, help="Pooled cis eGene table (summarize_cis.R output)")
parser$add_argument("--interaction_list", default="", help="File listing *.cis_qtl_top_assoc.txt.gz interaction files")
parser$add_argument("--n_plot", type="integer", default=20, help="Top pairs per model to plot")
parser$add_argument("--fdr", type="double", default=0.05, help="FDR threshold")
parser$add_argument("--gene_map", default="", help="FBgn -> symbol map (prep_gene_map.R)")
parser$add_argument("--out_prefix", default="fly", help="Output prefix")
args <- parser$parse_args()

read_gz <- function(f) fread(cmd=paste("zcat", shQuote(f)))

gm <- if (args$gene_map != "" && file.exists(args$gene_map) && file.size(args$gene_map) > 0) {
  fread(args$gene_map)[, .(phenotype_id=query_fbgn, gene_symbol, annotation_ID)]
} else NULL
annotate <- function(dt) {
  if (!is.null(gm) && !"gene_symbol" %in% colnames(dt)) dt <- merge(dt, gm, by="phenotype_id", all.x=TRUE, sort=FALSE)
  if (!"gene_symbol" %in% colnames(dt)) dt[, gene_symbol := phenotype_id]
  dt[is.na(gene_symbol), gene_symbol := phenotype_id]
  dt
}

# Main model: lead SNP per gene, most significant genes first
main <- annotate(read_gz(args$main))
main <- main[order(pval_beta, pval_nominal)]
top.main <- main[seq_len(min(args$n_plot, .N)),
                 .(model="main", phenotype_id, gene_symbol, variant_id, pval=pval_nominal, fdr, beta=slope)]
cat("Main eGenes (FDR <", args$fdr, "):", sum(main$fdr < args$fdr, na.rm=TRUE), "\n")

top.int <- NULL
if (args$interaction_list != "" && file.exists(args$interaction_list)) {
  files <- readLines(args$interaction_list)
  files <- files[files != ""]
  if (length(files) > 0) {
    int <- rbindlist(lapply(files, function(f) {
      dt <- read_gz(f)
      dt[, chrom := sub(".*\\.chr([^.]+)\\.K[0-9]+\\.interaction\\.cis_qtl_top_assoc\\.txt\\.gz$", "\\1", basename(f))]
      dt
    }), fill=TRUE)

    # eigenMT-corrected per-gene p-values (pval_emt), BH across all chromosomes
    int[, fdr_gi := p.adjust(pval_emt, method="BH")]
    int <- annotate(int)
    int <- int[order(pval_emt, pval_gi)]
    fwrite(int, paste0(args$out_prefix, ".interaction.top_assoc.all_chroms.txt.gz"), sep="\t")
    fwrite(int[fdr_gi < args$fdr], paste0(args$out_prefix, ".interaction.sig_genes.fdr", args$fdr, ".tsv"), sep="\t")
    fwrite(int[, .(n_genes=.N, n_sig=sum(fdr_gi < args$fdr, na.rm=TRUE)), by=chrom],
           paste0(args$out_prefix, ".interaction.sig_by_chrom.tsv"), sep="\t")
    cat("SNP x interaction genes (FDR <", args$fdr, "):", sum(int$fdr_gi < args$fdr, na.rm=TRUE), "\n")

    top.int <- int[seq_len(min(args$n_plot, .N)),
                   .(model="interaction", phenotype_id, gene_symbol, variant_id, pval=pval_gi, fdr=fdr_gi, beta=b_gi)]
  }
}

pairs <- rbind(top.main, top.int)
fwrite(pairs, paste0(args$out_prefix, ".plot_pairs.tsv"), sep="\t")
writeLines(unique(pairs$variant_id), "plot_variants.txt")
