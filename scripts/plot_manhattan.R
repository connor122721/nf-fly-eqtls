#!/usr/bin/env Rscript
# Connor Murray
# Manhattan, QQ, and summary plots for the main cis-eQTL and SNP x condition
# interaction models.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--main", required=TRUE, help="Pooled main-model eGene table (with gene_symbol)")
parser$add_argument("--interaction", default="", help="Pooled interaction table (with gene_symbol)")
parser$add_argument("--snp", default="", help="Reduced nominal SNP table (reduce_nominal.py)")
parser$add_argument("--chromosomes", required=TRUE, help="Chromosome order file")
parser$add_argument("--fdr", type="double", default=0.05)
parser$add_argument("--n_label", type="integer", default=15, help="Top genes to label")
parser$add_argument("--out_prefix", default="fly")
args <- parser$parse_args()

read_gz <- function(f) fread(cmd=paste("zcat", shQuote(f)))
chroms <- trimws(readLines(args$chromosomes)); chroms <- chroms[chroms != ""]
split_var <- function(dt) {
  dt[, chrom := sub(":[0-9]+$", "", variant_id)]
  dt[, pos := as.numeric(sub("^.*:", "", variant_id))]
  dt[chrom %in% chroms]
}

main <- split_var(read_gz(args$main))
if (!"gene_symbol" %in% colnames(main)) main[, gene_symbol := phenotype_id]
int <- NULL
if (args$interaction != "" && file.exists(args$interaction)) {
  int <- split_var(read_gz(args$interaction))
  if (!"gene_symbol" %in% colnames(int)) int[, gene_symbol := phenotype_id]
}
snp <- if (args$snp != "" && file.exists(args$snp)) split_var(read_gz(args$snp)) else NULL

# Cumulative genome coordinates
all.pos <- rbind(main[, .(chrom, pos)], if (!is.null(snp)) snp[, .(chrom, pos)], if (!is.null(int)) int[, .(chrom, pos)])
clen <- all.pos[, .(len=max(pos)), by=chrom][match(chroms, chrom)][!is.na(chrom)]
clen[, offset := cumsum(c(0, head(len, -1)))]
clen[, mid := offset + len / 2]
add_cum <- function(dt) {
  dt <- merge(dt, clen[, .(chrom, offset)], by="chrom")
  dt[, cum := pos + offset]
  dt[, chrom_i := match(chrom, chroms) %% 2]
  dt
}
xaxis <- scale_x_continuous(breaks=clen$mid, labels=clen$chrom, expand=c(0.01, 0))

manhattan <- function(dt, p, sig, title, ylab) {
  dt <- add_cum(copy(dt))
  dt[, lp := -log10(get(p))]
  dt[, sig := sig]
  lab <- dt[sig == TRUE][order(-lp)][!duplicated(gene_symbol)][seq_len(min(args$n_label, sum(dt$sig)))]
  ggplot(dt, aes(x=cum, y=lp)) +
    geom_point(aes(color=factor(chrom_i)), size=0.6, alpha=0.6, show.legend=FALSE) +
    scale_color_manual(values=c("0"="grey30", "1"="grey65")) +
    geom_point(data=dt[sig == TRUE], color="firebrick", size=0.9) +
    geom_text(data=lab, aes(label=gene_symbol), size=2.6, vjust=-0.6, check_overlap=TRUE, fontface="italic") +
    xaxis + theme_bw() + theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
    labs(x="Chromosome", y=ylab, title=title)
}

qq <- function(p, title) {
  p <- sort(p[!is.na(p) & p > 0])
  d <- data.table(exp=-log10(ppoints(length(p))), obs=-log10(p))
  lambda <- median(qchisq(p, 1, lower.tail=FALSE)) / qchisq(0.5, 1)
  ggplot(d, aes(x=exp, y=obs)) + geom_point(size=0.6) + geom_abline(color="red") + theme_bw() +
    labs(x="Expected -log10 p", y="Observed -log10 p", title=title, subtitle=sprintf("lambda = %.2f", lambda))
}

pages <- list()
add <- function(p) pages[[length(pages) + 1]] <<- p

n.main <- sum(main$fdr < args$fdr, na.rm=TRUE)
if (!is.null(snp)) {
  # Label the SNP-level plot with eGene lead SNPs
  eg <- main[fdr < args$fdr, .(variant_id, gene_symbol)]
  s2 <- merge(snp, eg[!duplicated(variant_id)], by="variant_id", all.x=TRUE)
  s2[is.na(gene_symbol), gene_symbol := phenotype_id]
  add(manhattan(s2, "pval_nominal", s2$variant_id %in% eg$variant_id,
                "cis-eQTL Manhattan (all tested SNPs; min p across genes; non-significant SNPs thinned)",
                "-log10 nominal p") + labs(subtitle="Red: lead SNPs of eGenes"))
}
add(manhattan(main, "pval_nominal", main$fdr < args$fdr,
              sprintf("Main model: lead SNP per gene (%d eGenes at FDR < %.2g, red)", n.main, args$fdr),
              "-log10 nominal p (lead SNP)"))
add(qq(main$pval_beta, "Main model QQ: permutation (beta) p-value per gene"))

if (!is.null(int)) {
  n.int <- sum(int$fdr_gi < args$fdr, na.rm=TRUE)
  add(manhattan(int, "pval_gi", int$fdr_gi < args$fdr,
                sprintf("SNP x condition: lead interaction SNP per gene (%d genes at FDR < %.2g, red)", n.int, args$fdr),
                "-log10 interaction p (lead SNP)"))
  add(qq(int$pval_emt[int$pval_emt < 1], "Interaction QQ: eigenMT-corrected p per gene (p < 1)"))
}

# eGenes per chromosome
cnt <- rbind(main[, .(model="main", n=sum(fdr < args$fdr, na.rm=TRUE), tested=.N), by=chrom],
             if (!is.null(int)) int[, .(model="interaction", n=sum(fdr_gi < args$fdr, na.rm=TRUE), tested=.N), by=chrom])
cnt[, chrom := factor(chrom, levels=chroms)]
add(ggplot(cnt, aes(x=chrom, y=n, fill=model)) + geom_col(position="dodge") + theme_bw() +
      geom_text(aes(label=sprintf("%d/%d", n, tested)), position=position_dodge(0.9), vjust=-0.3, size=2.6) +
      labs(x="Chromosome", y=sprintf("Significant genes (FDR < %.2g)", args$fdr), title="Significant genes per chromosome"))

# Lead SNP position relative to TSS, and effect sizes
if ("start_distance" %in% colnames(main)) {
  add(ggplot(main[fdr < args$fdr], aes(x=start_distance / 1000)) + geom_histogram(bins=80) + theme_bw() +
        labs(x="Lead SNP distance to TSS (kb)", y="eGenes", title="eGene lead SNP position relative to TSS"))
}
if ("slope" %in% colnames(main)) {
  add(ggplot(main[fdr < args$fdr], aes(x=slope)) + geom_histogram(bins=80) + theme_bw() +
        labs(x="Effect size (slope, normalized expression per allele)", y="eGenes", title="eGene effect sizes"))
}
if ("af" %in% colnames(main)) {
  add(ggplot(main, aes(x=pmin(af, 1 - af), fill=fdr < args$fdr)) + geom_histogram(bins=50, position="identity", alpha=0.6) +
        theme_bw() + labs(x="Lead SNP minor allele frequency", y="Genes", fill="eGene", title="Lead SNP MAF"))
}

pdf(paste0(args$out_prefix, ".manhattan_qq.pdf"), width=11, height=5)
for (p in pages) print(p)
invisible(dev.off())
