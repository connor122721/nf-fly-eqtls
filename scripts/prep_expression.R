#!/usr/bin/env Rscript
# Connor Murray
# QC, normalization, and covariate construction for fly cis-eQTL mapping
#
# Steps:
#   1) Read featureCounts matrix + metadata, keep samples with RNA, metadata, and genotypes
#   2) Sample QC: library size / detected genes (MAD) and expression PCA outliers (Mahalanobis)
#   3) Gene QC (GTEx-style): CPM >= threshold AND counts >= threshold in >= frac of samples
#   4) TMM normalization (edgeR) -> log2 CPM -> inverse normal transform per gene
#   5) Gene TSS from GTF -> tensorQTL phenotype BED (main chromosomes only)
#   6) Covariates: model.matrix(formula) (default ~ cage_name + condition)
#      + optional genotype PCs + K expression PCs (PCs of residuals after known covariates)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(argparse)
})

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--counts", required=TRUE, help="Gene x sample count matrix (featureCounts style, first column gene ID)")
parser$add_argument("--metadata", required=TRUE, help="Sample metadata TSV")
parser$add_argument("--gtf", required=TRUE, help="GTF used for counting (gene_id must match count matrix)")
parser$add_argument("--chromosomes", required=TRUE, help="File with chromosomes to keep, one per line")
parser$add_argument("--genotyped_samples", default=NULL, help="Samples present in the VCF, one per line")
parser$add_argument("--geno_pcs", default=NULL, help="plink2 .eigenvec file")
parser$add_argument("--n_geno_pcs", type="integer", default=0, help="Number of genotype PCs to include")
parser$add_argument("--expr_pcs", default="0,5,10,20", help="Comma list of expression PC counts (one covariate file each)")
parser$add_argument("--sample_col", default="sample", help="Sample ID column in metadata")
parser$add_argument("--formula", default="~ cage_name + condition", help="Known covariate model")
parser$add_argument("--factor_covariates", default="cage_name,condition", help="Metadata columns treated as factors")
parser$add_argument("--cpm_threshold", type="double", default=1, help="CPM threshold for gene filter")
parser$add_argument("--count_threshold", type="integer", default=6, help="Read count threshold for gene filter")
parser$add_argument("--sample_frac", type="double", default=0.2, help="Fraction of samples passing thresholds")
parser$add_argument("--min_lib_size", type="double", default=0, help="Absolute minimum library size")
parser$add_argument("--lib_mad", type="double", default=3, help="Drop samples this many MADs below median log10 library size / detected genes")
parser$add_argument("--pca_outlier_p", type="double", default=0.05, help="Bonferroni-corrected Mahalanobis p-value cutoff for expression PCA outliers")
parser$add_argument("--pca_outlier_npc", type="integer", default=5, help="Number of expression PCs used for outlier detection")
parser$add_argument("--norm_method", default="int", choices=c("int", "logcpm"), help="int: TMM logCPM + inverse normal transform; logcpm: TMM logCPM only")
parser$add_argument("--out_prefix", default="fly", help="Output prefix")
args <- parser$parse_args()

prefix <- args$out_prefix
msg <- function(...) cat(sprintf(...), "\n", sep="")

### 1) Read data ###

# Count matrix: drop featureCounts annotation columns if present
counts.dt <- fread(args$counts, header=TRUE)
setnames(counts.dt, 1, "gene_id")
drop.cols <- intersect(c("Chr", "Start", "End", "Strand", "Length"), colnames(counts.dt))
if (length(drop.cols) > 0) counts.dt[, (drop.cols) := NULL]
counts.dt <- counts.dt[!duplicated(gene_id)]
counts <- as.matrix(counts.dt[, -1, with=FALSE])
rownames(counts) <- counts.dt$gene_id
storage.mode(counts) <- "double"
msg("Count matrix: %d genes x %d samples", nrow(counts), ncol(counts))

# Metadata
meta <- fread(args$metadata, header=TRUE, colClasses="character")
if (!args$sample_col %in% colnames(meta)) stop("Metadata lacks sample column: ", args$sample_col)
setnames(meta, args$sample_col, "sample")
meta <- meta[!duplicated(sample)]

# Covariates used in the model must be non-missing
form <- as.formula(args$formula)
model.vars <- all.vars(form)
missing.vars <- setdiff(model.vars, colnames(meta))
if (length(missing.vars) > 0) stop("Formula variables not in metadata: ", paste(missing.vars, collapse=", "))
for (v in model.vars) meta <- meta[!is.na(get(v)) & get(v) != "" & get(v) != "NA"]

# Factor vs numeric covariates
factor.vars <- trimws(strsplit(args$factor_covariates, ",")[[1]])
factor.vars <- factor.vars[factor.vars != ""]
for (v in model.vars) {
  if (v %in% factor.vars) meta[, (v) := factor(get(v))]
  else meta[, (v) := as.numeric(get(v))]
}

# Sample intersection: RNA ∩ metadata ∩ genotypes
samps <- intersect(colnames(counts), meta$sample)
msg("Samples with RNA + metadata: %d (RNA only: %d; metadata only: %d)",
    length(samps), length(setdiff(colnames(counts), meta$sample)),
    length(setdiff(meta$sample, colnames(counts))))
if (!is.null(args$genotyped_samples)) {
  geno.samps <- fread(args$genotyped_samples, header=FALSE)$V1
  msg("Samples with RNA + metadata but no genotypes (dropped): %d", length(setdiff(samps, geno.samps)))
  samps <- intersect(samps, geno.samps)
}

# Genotype PCs
geno.pc <- NULL
if (args$n_geno_pcs > 0) {
  if (is.null(args$geno_pcs)) stop("--n_geno_pcs > 0 requires --geno_pcs")
  geno.pc <- fread(args$geno_pcs, header=TRUE)
  setnames(geno.pc, sub("^#", "", colnames(geno.pc)))
  if ("FID" %in% colnames(geno.pc) && "IID" %in% colnames(geno.pc)) geno.pc[, FID := NULL]
  setnames(geno.pc, "IID", "sample")
  pc.cols <- paste0("PC", seq_len(args$n_geno_pcs))
  if (!all(pc.cols %in% colnames(geno.pc))) stop("Requested more genotype PCs than in ", args$geno_pcs)
  geno.pc <- geno.pc[, c("sample", pc.cols), with=FALSE]
  setnames(geno.pc, pc.cols, paste0("geno_", pc.cols))
  msg("Samples without genotype PCs (dropped): %d", length(setdiff(samps, geno.pc$sample)))
  samps <- intersect(samps, geno.pc$sample)
}
msg("Samples entering QC: %d", length(samps))
if (length(samps) < 20) stop("Too few samples to continue")
counts <- counts[, samps, drop=FALSE]

### 2-4) Sample QC, gene QC, normalization ###

# Gene filter + TMM + logCPM (+ INT)
normalize_counts <- function(cnt) {
  y <- DGEList(cnt)
  cpm.raw <- cpm(y)
  keep <- rowMeans(cpm.raw >= args$cpm_threshold) >= args$sample_frac &
          rowMeans(cnt >= args$count_threshold) >= args$sample_frac
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method="TMM")
  logcpm <- cpm(y, log=TRUE, prior.count=1)
  list(y=y, logcpm=logcpm, n_genes=sum(keep))
}

# Rank-based inverse normal transform (per gene)
inverse_normal <- function(mat) {
  out <- t(apply(mat, 1, function(x) qnorm((rank(x, ties.method="average") - 0.5) / length(x))))
  dimnames(out) <- dimnames(mat)
  out
}

# Library size / detected gene outliers
qc <- data.table(sample=colnames(counts),
                 lib_size=colSums(counts),
                 detected=colSums(counts > 0))
qc[, log_lib := log10(lib_size + 1)]
qc[, low_lib := log_lib < median(log_lib) - args$lib_mad * mad(log_lib) | lib_size < args$min_lib_size]
qc[, low_detected := detected < median(detected) - args$lib_mad * mad(detected)]
msg("Low library size: %d; low detected genes: %d", sum(qc$low_lib), sum(qc$low_detected))

# Expression PCA outliers (Mahalanobis distance on top PCs)
pass1 <- qc[!low_lib & !low_detected]$sample
norm1 <- normalize_counts(counts[, pass1, drop=FALSE])
npc <- min(args$pca_outlier_npc, length(pass1) - 1)
pca1 <- prcomp(t(norm1$logcpm), center=TRUE, scale.=FALSE)
pcx <- pca1$x[, seq_len(npc), drop=FALSE]
mal <- mahalanobis(pcx, colMeans(pcx), cov(pcx))
mal.p <- pchisq(mal, df=npc, lower.tail=FALSE)
qc[, mahal_p := mal.p[match(sample, names(mal.p))]]
qc[, pca_outlier := !is.na(mahal_p) & mahal_p < args$pca_outlier_p / length(pass1)]
msg("Expression PCA outliers: %d", sum(qc$pca_outlier))

qc[, keep := !low_lib & !low_detected & !pca_outlier]
samps <- qc[keep == TRUE]$sample
msg("Samples passing QC: %d", length(samps))

# Final normalization on passing samples
norm <- normalize_counts(counts[, samps, drop=FALSE])
msg("Genes passing expression filter: %d of %d", norm$n_genes, nrow(counts))
expr <- if (args$norm_method == "int") inverse_normal(norm$logcpm) else norm$logcpm

### 5) Gene TSS from GTF ###

gtf <- fread(args$gtf, sep="\t", header=FALSE, quote="", select=c(1, 3, 4, 5, 7, 9),
             col.names=c("chr", "feature", "start", "end", "strand", "attr"))
gtf <- gtf[!grepl("^#", chr)]
gtf[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attr)]
gtf <- gtf[gene_id != attr]

# Match GTF chromosome naming (e.g. chr2L) to genotype naming (e.g. 2L)
chroms <- trimws(readLines(args$chromosomes))
chroms <- chroms[chroms != ""]
if (sum(sub("^chr", "", gtf$chr) %in% chroms) > sum(gtf$chr %in% chroms)) gtf[, chr := sub("^chr", "", chr)]

# One TSS per gene: most upstream transcript start; drop genes on >1 chromosome/strand
genes <- gtf[, .(n_loc=uniqueN(paste(chr, strand)), chr=chr[1], strand=strand[1],
                 gstart=min(start), gend=max(end)), by=gene_id]
genes <- genes[n_loc == 1 & chr %in% chroms]
genes[, tss := ifelse(strand == "+", gstart, gend)]
msg("Genes with TSS on %s: %d", paste(chroms, collapse=","), nrow(genes))

expr <- expr[rownames(expr) %in% genes$gene_id, , drop=FALSE]
msg("Expressed genes on main chromosomes: %d", nrow(expr))

# tensorQTL BED: 0-based start = TSS - 1, end = TSS
bed <- genes[match(rownames(expr), gene_id), .(`#chr`=chr, start=tss - 1L, end=tss, phenotype_id=gene_id)]
bed <- cbind(bed, as.data.table(round(expr, 6)))
bed[, chr.o := match(`#chr`, chroms)]
setorder(bed, chr.o, start)
bed[, chr.o := NULL]
fwrite(bed, paste0(prefix, ".expression.bed"), sep="\t", quote=FALSE)

### 6) Covariates ###

meta.s <- meta[match(samps, sample)]
mm <- model.matrix(form, data=meta.s)

# Drop aliased (collinear) columns, e.g. cage nested within condition
qr.mm <- qr(mm)
if (qr.mm$rank < ncol(mm)) {
  dropped <- colnames(mm)[qr.mm$pivot[(qr.mm$rank + 1):ncol(mm)]]
  msg("WARNING: dropping aliased covariate columns: %s", paste(dropped, collapse=", "))
  mm <- mm[, qr.mm$pivot[seq_len(qr.mm$rank)], drop=FALSE]
}
known <- mm[, colnames(mm) != "(Intercept)", drop=FALSE]
if (!is.null(geno.pc)) {
  gp <- as.matrix(geno.pc[match(samps, sample), -1, with=FALSE])
  known <- cbind(known, gp)
}
colnames(known) <- make.names(colnames(known))
rownames(known) <- samps
msg("Known covariates (%d): %s", ncol(known), paste(colnames(known), collapse=", "))

# Expression PCs from residuals after removing known covariates
X <- cbind(1, known)
resid <- lm.fit(X, t(expr))$residuals
pca.res <- prcomp(resid, center=TRUE, scale.=FALSE)
expr.pcs <- pca.res$x
colnames(expr.pcs) <- paste0("expr_PC", seq_len(ncol(expr.pcs)))

ks <- sort(unique(as.integer(strsplit(args$expr_pcs, ",")[[1]])))
max.k <- length(samps) - ncol(known) - 10
if (any(ks > max.k)) {
  msg("WARNING: dropping K > %d (too many covariates for %d samples)", max.k, length(samps))
  ks <- ks[ks <= max.k]
}
for (k in ks) {
  covk <- if (k > 0) cbind(known, expr.pcs[, seq_len(k), drop=FALSE]) else known
  if (ncol(covk) == 0) {
    msg("WARNING: K=%d has no covariates; skipping", k)
    next
  }
  out <- data.table(id=colnames(covk), t(covk))
  setnames(out, c("id", samps))
  fwrite(out, paste0(prefix, ".covariates.K", k, ".txt"), sep="\t", quote=FALSE)
}

### Outputs: sample list, QC tables, plots ###

writeLines(samps, paste0(prefix, ".samples.txt"))
qc.out <- merge(qc[, .(sample, lib_size, detected, mahal_p, low_lib, low_detected, pca_outlier, keep)],
                meta, by="sample", all.x=TRUE)
fwrite(qc.out, paste0(prefix, ".sample_qc.tsv"), sep="\t")
fwrite(data.table(PC=seq_along(pca.res$sdev), var_explained=pca.res$sdev^2 / sum(pca.res$sdev^2)),
       paste0(prefix, ".expr_pc_variance.tsv"), sep="\t")

# Variance explained per gene by known covariates
r2 <- 1 - colSums(resid^2) / colSums(sweep(t(expr), 2, colMeans(t(expr)))^2)
fwrite(data.table(gene_id=rownames(expr), r2_known=r2), paste0(prefix, ".gene_r2_known_covariates.tsv"), sep="\t")

col.var <- intersect(c("condition", model.vars), colnames(meta))[1]
shape.var <- setdiff(intersect(c("cage_name", model.vars), colnames(meta)), col.var)[1]

pdf(paste0(prefix, ".expression_qc.pdf"), width=8, height=6)

print(ggplot(qc.out, aes(x=log10(lib_size), fill=keep)) +
        geom_histogram(bins=50) + theme_bw() +
        labs(x="log10 library size", y="Samples", title="Library size"))

print(ggplot(qc.out, aes(x=log10(lib_size), y=detected, color=keep)) +
        geom_point(alpha=0.7) + theme_bw() +
        labs(x="log10 library size", y="Genes detected (count > 0)", title="Sample QC"))

pdat <- data.table(sample=rownames(pca1$x), PC1=pca1$x[, 1], PC2=pca1$x[, 2])
pdat <- merge(pdat, qc.out, by="sample")
v1 <- round(100 * pca1$sdev[1:2]^2 / sum(pca1$sdev^2), 1)
p <- ggplot(pdat, aes(x=PC1, y=PC2)) + theme_bw() +
  labs(x=paste0("PC1 (", v1[1], "%)"), y=paste0("PC2 (", v1[2], "%)"),
       title="Expression PCA before covariate removal (logCPM)")
if (!is.na(col.var)) {
  pdat[, col := as.character(get(col.var))]
  pdat[, shp := if (!is.na(shape.var)) as.character(get(shape.var)) else "all"]
  p <- p + geom_point(data=pdat, aes(color=col, shape=shp), alpha=0.7) +
    labs(color=col.var, shape=ifelse(is.na(shape.var), "", shape.var))
} else {
  p <- p + geom_point(alpha=0.7)
}
p <- p + geom_point(data=pdat[pca_outlier == TRUE], shape=4, size=3, color="black")
print(p)

rdat <- data.table(sample=samps, PC1=pca.res$x[, 1], PC2=pca.res$x[, 2])
rdat <- merge(rdat, meta.s, by="sample")
v2 <- round(100 * pca.res$sdev[1:2]^2 / sum(pca.res$sdev^2), 1)
pr <- ggplot(rdat, aes(x=PC1, y=PC2)) + theme_bw() +
  labs(x=paste0("PC1 (", v2[1], "%)"), y=paste0("PC2 (", v2[2], "%)"),
       title="Expression PCA of residuals after known covariates")
if (!is.na(col.var)) {
  pr <- pr + geom_point(aes(color=as.character(.data[[col.var]])), alpha=0.7) + labs(color=col.var)
} else {
  pr <- pr + geom_point(alpha=0.7)
}
print(pr)

sdat <- data.table(PC=1:min(30, length(pca.res$sdev)))
sdat[, var := 100 * pca.res$sdev[PC]^2 / sum(pca.res$sdev^2)]
print(ggplot(sdat, aes(x=PC, y=var)) + geom_point() + geom_line() + theme_bw() +
        labs(x="Residual expression PC", y="% variance explained", title="Scree plot"))

print(ggplot(data.table(r2=r2), aes(x=r2)) + geom_histogram(bins=50) + theme_bw() +
        labs(x="Per-gene R^2", y="Genes", title=paste("Variance explained by", args$formula)))

dev.off()
msg("Done.")
