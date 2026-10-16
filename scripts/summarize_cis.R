#!/usr/bin/env Rscript
# Connor Murray
# Summarize tensorQTL cis (permutation) output across chromosomes and
# expression-PC models (K); pick the K that maximizes eGenes; write the
# pooled, FDR-corrected eGene table for that K.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--cis_list", required=TRUE, help="File listing *.cis_qtl.txt.gz files")
parser$add_argument("--fdr", type="double", default=0.05, help="eGene FDR threshold")
parser$add_argument("--gene_map", default="", help="FBgn -> symbol map (prep_gene_map.R)")
parser$add_argument("--out_prefix", default="fly", help="Output prefix")
args <- parser$parse_args()

files <- readLines(args$cis_list)
files <- files[files != ""]

# Files are named <prefix>.chr<CHROM>.K<k>.cis_qtl.txt.gz
cis <- rbindlist(lapply(files, function(f) {
  dt <- fread(cmd=paste("zcat", shQuote(f)))
  m <- regmatches(basename(f), regexec("\\.chr([^.]+)\\.K([0-9]+)\\.cis_qtl\\.txt\\.gz$", basename(f)))[[1]]
  dt[, chrom := m[2]]
  dt[, K := as.integer(m[3])]
  dt
}), fill=TRUE)

# Pool chromosomes within each model before FDR correction
cis[, fdr := p.adjust(pval_beta, method="BH"), by=K]
summ <- cis[, .(n_genes=.N, n_egenes=sum(fdr < args$fdr, na.rm=TRUE)), by=K][order(K)]
print(summ)

# Best K: most eGenes (ties -> fewest PCs)
best_k <- summ[order(-n_egenes, K)]$K[1]
cat("Best K:", best_k, "\n")
writeLines(as.character(best_k), "best_k")
fwrite(summ, paste0(args$out_prefix, ".egenes_by_K.tsv"), sep="\t")

# Per-chromosome summary for the best model
fwrite(cis[K == best_k, .(n_genes=.N, n_egenes=sum(fdr < args$fdr, na.rm=TRUE)), by=chrom],
       paste0(args$out_prefix, ".egenes_by_chrom.K", best_k, ".tsv"), sep="\t")
best <- cis[K == best_k][order(fdr, pval_beta)]
if (args$gene_map != "" && file.exists(args$gene_map) && file.size(args$gene_map) > 0) {
  gm <- fread(args$gene_map)[, .(phenotype_id=query_fbgn, gene_symbol, annotation_ID)]
  best <- merge(best, gm, by="phenotype_id", all.x=TRUE, sort=FALSE)[order(fdr, pval_beta)]
}
fwrite(best, paste0(args$out_prefix, ".cis_eqtl.K", best_k, ".all_chroms.txt.gz"), sep="\t")
fwrite(best[fdr < args$fdr], paste0(args$out_prefix, ".egenes.K", best_k, ".fdr", args$fdr, ".tsv"), sep="\t")

pdf(paste0(args$out_prefix, ".egenes_by_K.pdf"), width=6, height=4)
print(ggplot(summ, aes(x=K, y=n_egenes)) + geom_point() + geom_line() +
        geom_vline(xintercept=best_k, linetype="dashed", color="red") + theme_bw() +
        labs(x="Expression PCs in model", y=paste0("eGenes (FDR < ", args$fdr, ")"),
             title="cis-eQTL discovery by number of expression PCs"))
dev.off()
