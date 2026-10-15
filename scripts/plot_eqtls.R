#!/usr/bin/env Rscript
# Connor Murray
# eQTL plots: genotype (x) vs normalized expression (y) for the top
# gene-SNP pairs from the main and SNP x condition interaction models.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--pairs", required=TRUE, help="plot_pairs.tsv from select_top_pairs.R")
parser$add_argument("--genotypes", required=TRUE, help="plink2 --export A .raw file")
parser$add_argument("--bim", required=TRUE, help="Concatenated .bim lines for plotted variants")
parser$add_argument("--expression", required=TRUE, help="tensorQTL phenotype BED")
parser$add_argument("--metadata", required=TRUE, help="Sample metadata")
parser$add_argument("--sample_col", default="sample", help="Sample ID column in metadata")
parser$add_argument("--group", default="condition", help="Metadata column for coloring / facets")
parser$add_argument("--out_prefix", default="fly", help="Output prefix")
args <- parser$parse_args()

pairs <- fread(args$pairs)
meta <- fread(args$metadata, colClasses="character")
setnames(meta, args$sample_col, "sample")
has.group <- args$group %in% colnames(meta)

# Expression: genes x samples -> long
expr <- fread(args$expression)
setnames(expr, 1:4, c("chr", "start", "end", "phenotype_id"))
expr <- expr[phenotype_id %in% pairs$phenotype_id]
expr.long <- melt(expr[, -(1:3)], id.vars="phenotype_id", variable.name="sample", value.name="expression")
expr.long[, sample := as.character(sample)]

# Genotypes: columns <variant>_<counted allele>, values = counted-allele dosage
geno <- fread(args$genotypes)
geno.cols <- setdiff(colnames(geno), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
geno.long <- melt(geno[, c("IID", geno.cols), with=FALSE], id.vars="IID",
                  variable.name="col", value.name="dosage")
geno.long[, col := as.character(col)]
geno.long[, variant_id := sub("_[^_]+$", "", col)]
geno.long[, counted := sub("^.*_", "", col)]
setnames(geno.long, "IID", "sample")

# Other allele from the bim (A1/A2 columns)
bim <- fread(args$bim, header=FALSE, col.names=c("chr", "variant_id", "cm", "pos", "a1", "a2"))
bim <- unique(bim, by="variant_id")
geno.long <- merge(geno.long, bim[, .(variant_id, a1, a2)], by="variant_id")
geno.long[, other := ifelse(counted == a1, a2, a1)]
geno.long[, genotype := factor(
  ifelse(dosage == 0, paste0(other, other), ifelse(dosage == 1, paste0(other, counted), paste0(counted, counted))))]

plot_pair <- function(i) {
  pr <- pairs[i]
  d <- merge(expr.long[phenotype_id == pr$phenotype_id],
             geno.long[variant_id == pr$variant_id & !is.na(dosage)], by="sample")
  if (nrow(d) == 0) return(NULL)
  if (has.group) d <- merge(d, meta[, c("sample", args$group), with=FALSE], by="sample")
  lv <- unique(d[order(dosage), as.character(genotype)])
  d[, genotype := factor(as.character(genotype), levels=lv)]
  ttl <- sprintf("%s  |  %s", pr$phenotype_id, pr$variant_id)
  sub <- sprintf("%s model: p = %.2g, FDR = %.2g, beta = %.2f", pr$model, pr$pval, pr$fdr, pr$beta)

  if (pr$model == "interaction" && has.group) {
    d[, grp := as.character(get(args$group))]
    p <- ggplot(d, aes(x=genotype, y=expression, fill=grp)) +
      geom_boxplot(outlier.shape=NA, alpha=0.5, position=position_dodge(0.8)) +
      geom_point(aes(color=grp), position=position_jitterdodge(jitter.width=0.2, dodge.width=0.8),
                 size=0.8, alpha=0.6) +
      stat_summary(aes(group=grp, color=grp), fun=median, geom="line",
                   position=position_dodge(0.8), linewidth=0.8) +
      labs(fill=args$group, color=args$group)
  } else {
    p <- ggplot(d, aes(x=genotype, y=expression)) +
      geom_boxplot(outlier.shape=NA, alpha=0.4)
    if (has.group) {
      d[, grp := as.character(get(args$group))]
      p <- p + geom_jitter(data=d, aes(color=grp), width=0.2, size=0.8, alpha=0.6) + labs(color=args$group)
    } else {
      p <- p + geom_jitter(width=0.2, size=0.8, alpha=0.6)
    }
  }
  n.tab <- d[, .N, by=genotype]
  p + theme_bw() +
    scale_x_discrete(labels=setNames(paste0(n.tab$genotype, "\n(n=", n.tab$N, ")"), n.tab$genotype)) +
    labs(x="Genotype", y="Normalized expression", title=ttl, subtitle=sub)
}

for (m in unique(pairs$model)) {
  pdf(paste0(args$out_prefix, ".eqtl_plots.", m, ".pdf"), width=7, height=5)
  for (i in which(pairs$model == m)) {
    p <- plot_pair(i)
    if (!is.null(p)) print(p)
  }
  dev.off()
}
