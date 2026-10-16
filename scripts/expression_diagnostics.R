#!/usr/bin/env Rscript
# Connor Murray
# Expression diagnostics for eQTL readiness: normalization, gene filtering,
# sample structure, covariate associations, residual PCs, genotype PCs.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--prep", required=TRUE, help="prep_objects.rds from prep_expression.R")
parser$add_argument("--geno_pcs", default=NULL, help="plink2 .eigenvec file")
parser$add_argument("--gene_map", default=NULL, help="FBgn -> symbol map (prep_gene_map.R)")
parser$add_argument("--group", default="condition", help="Main grouping column for colors")
parser$add_argument("--out_prefix", default="fly", help="Output prefix")
args <- parser$parse_args()

obj <- readRDS(args$prep)
prefix <- args$out_prefix
meta <- copy(obj$meta)
qc <- copy(obj$qc)
grp <- if (args$group %in% colnames(meta)) args$group else NULL
other.vars <- setdiff(obj$model_vars, grp)

# Group label (e.g. condition x cage) for ordering samples
meta[, group_label := if (length(obj$model_vars) > 0)
  do.call(paste, c(lapply(obj$model_vars, function(v) as.character(get(v))), sep=" | ")) else "all"]
# Group labels for all samples entering QC (including those that failed)
meta.all <- if (!is.null(obj$meta_all)) copy(obj$meta_all) else copy(meta)
meta.all[, group_label := if (length(obj$model_vars) > 0)
  do.call(paste, c(lapply(obj$model_vars, function(v) as.character(get(v))), sep=" | ")) else "all"]
qc <- merge(qc, meta.all[, .(sample, group_label)], by="sample", all.x=TRUE)

pages <- list()
add <- function(p) pages[[length(pages) + 1]] <<- p

### 1) Library size and detection per group ###
add(ggplot(qc[!is.na(group_label)], aes(x=group_label, y=log10(lib_size), fill=group_label)) +
      geom_violin(alpha=0.4) + geom_jitter(width=0.15, size=0.6, alpha=0.6) + theme_bw() +
      theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="none") +
      labs(x=NULL, y="log10 library size", title="Library size by group (samples entering QC)"))

add(ggplot(qc[!is.na(group_label)], aes(x=group_label, y=detected, fill=group_label)) +
      geom_violin(alpha=0.4) + geom_jitter(width=0.15, size=0.6, alpha=0.6) + theme_bw() +
      theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="none") +
      labs(x=NULL, y="Genes with count > 0", title="Detected genes by group"))

add(ggplot(qc, aes(x=-log10(mahal_p), fill=pca_outlier)) + geom_histogram(bins=60) + theme_bw() +
      labs(x="-log10 Mahalanobis p (expression PCs)", y="Samples", title="Expression PCA outlier test"))

# QC failure counts by group
fail <- qc[, .(pass=sum(keep), low_lib=sum(low_lib), low_detected=sum(low_detected),
               pca_outlier=sum(pca_outlier, na.rm=TRUE)), by=group_label]
fwrite(fail, paste0(prefix, ".qc_by_group.tsv"), sep="\t")
add(ggplot(melt(fail, id.vars="group_label"), aes(x=group_label, y=value, fill=variable)) +
      geom_col(position="dodge") + theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1)) +
      labs(x=NULL, y="Samples", fill=NULL, title="Sample QC outcome by group"))

### 2) Normalization ###
ns <- as.data.table(obj$norm_samples, keep.rownames="sample")
ns <- merge(ns, meta, by="sample")
add(ggplot(ns, aes(x=group_label, y=norm.factors, fill=group_label)) +
      geom_boxplot(outlier.shape=NA, alpha=0.4) + geom_jitter(width=0.15, size=0.6) + theme_bw() +
      theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="none") +
      geom_hline(yintercept=1, linetype="dashed") +
      labs(x=NULL, y="TMM normalization factor", title="TMM factors by group"))

# Per-sample density before (raw logCPM, no TMM) and after TMM
cnt <- obj$counts[obj$genes_kept, obj$samples]
raw.lcpm <- cpm(cnt, log=TRUE, prior.count=1)
set.seed(1)
show <- if (ncol(cnt) > 120) sample(colnames(cnt), 120) else colnames(cnt)
dens <- rbind(
  data.table(stage="1. logCPM (library size only)", melt(as.data.table(raw.lcpm[, show], keep.rownames="g"), id.vars="g")),
  data.table(stage="2. logCPM after TMM", melt(as.data.table(obj$logcpm[, show], keep.rownames="g"), id.vars="g")))
dens <- merge(dens, meta[, c("sample", grp), with=FALSE], by.x="variable", by.y="sample")
add(ggplot(dens, aes(x=value, group=variable)) +
      geom_density(do.call(aes, if (!is.null(grp)) list(color=as.name(grp)) else list()), linewidth=0.2, alpha=0.5) +
      facet_wrap(~stage, ncol=1) + theme_bw() +
      labs(x="log2 CPM", y="Density", title="Per-sample expression distributions (filtered genes, up to 120 samples)"))

# Relative log expression (RLE) after TMM, samples ordered by group
rle <- obj$logcpm - apply(obj$logcpm, 1, median)
rle.s <- data.table(sample=colnames(rle), t(apply(rle, 2, quantile, c(0.1, 0.25, 0.5, 0.75, 0.9))))
setnames(rle.s, c("sample", "q10", "q25", "q50", "q75", "q90"))
rle.s <- merge(rle.s, meta[, .(sample, group_label)], by="sample")
rle.s <- rle.s[order(group_label, q50)]
rle.s[, x := .I]
add(ggplot(rle.s, aes(x=x)) +
      geom_linerange(aes(ymin=q10, ymax=q90, color=group_label), alpha=0.4) +
      geom_linerange(aes(ymin=q25, ymax=q75, color=group_label), linewidth=0.8) +
      geom_point(aes(y=q50), size=0.4) + geom_hline(yintercept=0) + theme_bw() +
      labs(x="Samples (ordered by group)", y="RLE (logCPM - gene median)", color=NULL,
           title="Relative log expression after TMM (should center on 0)"))

### 3) Gene filter ###
all.lcpm <- cpm(obj$counts[, obj$samples], log=TRUE, prior.count=1)
gstat <- data.table(gene_id=rownames(all.lcpm), mean=rowMeans(all.lcpm),
                    sd=apply(all.lcpm, 1, sd), kept=rownames(all.lcpm) %in% obj$genes_kept)
add(ggplot(gstat, aes(x=mean, fill=kept)) + geom_histogram(bins=80, alpha=0.7, position="identity") +
      theme_bw() + labs(x="Mean log2 CPM", y="Genes", title=sprintf("Expression filter: %d of %d genes kept",
                                                                   sum(gstat$kept), nrow(gstat))))
add(ggplot(gstat, aes(x=mean, y=sd)) + geom_point(aes(color=kept), size=0.4, alpha=0.4) +
      geom_smooth(method="loess", se=FALSE, color="black", formula=y ~ x) + theme_bw() +
      labs(x="Mean log2 CPM", y="SD log2 CPM", title="Mean-variance relationship"))

tss <- obj$tss[gene_id %in% rownames(obj$expr)]
add(ggplot(tss[, .N, by=chr], aes(x=chr, y=N)) + geom_col() + theme_bw() +
      geom_text(aes(label=N), vjust=-0.3) +
      labs(x="Chromosome", y="Genes tested", title="Expressed genes per chromosome (phenotypes for tensorQTL)"))

### 4) Sample structure ###
# Sample-sample correlation on normalized expression, ordered by group
ord <- meta[order(group_label)]$sample
cc <- cor(obj$expr[, ord])
cdt <- melt(as.data.table(cc, keep.rownames="s1"), id.vars="s1", variable.name="s2", value.name="r")
cdt[, s1 := factor(s1, levels=ord)][, s2 := factor(s2, levels=ord)]
brk <- meta[order(group_label), .(pos=.N), by=group_label][, pos := cumsum(pos)]
add(ggplot(cdt, aes(x=s1, y=s2, fill=r)) + geom_raster() +
      scale_fill_viridis_c() +
      geom_vline(xintercept=brk$pos + 0.5, color="white", linewidth=0.3) +
      geom_hline(yintercept=brk$pos + 0.5, color="white", linewidth=0.3) +
      theme_minimal() + theme(axis.text=element_blank(), axis.ticks=element_blank()) +
      labs(x=NULL, y=NULL, fill="Pearson r",
           title="Sample-sample correlation (normalized expression)",
           subtitle=paste("Blocks:", paste(brk$group_label, collapse="; "))))

# Mean correlation of each sample to all others (low = possible outlier)
mc <- data.table(sample=rownames(cc), mean_r=(rowSums(cc) - 1) / (ncol(cc) - 1))
mc <- merge(mc, meta[, .(sample, group_label)], by="sample")
add(ggplot(mc, aes(x=group_label, y=mean_r)) + geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.15, size=0.6) + theme_bw() +
      theme(axis.text.x=element_text(angle=30, hjust=1)) +
      labs(x=NULL, y="Mean correlation to other samples", title="Sample connectivity"))

# PCA before covariate removal (logCPM of QC-passing samples)
pca <- prcomp(t(obj$logcpm), center=TRUE)
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pd <- data.table(sample=rownames(pca$x), pca$x[, 1:min(10, ncol(pca$x))])
pd <- merge(pd, meta, by="sample")
pd <- merge(pd, qc[, .(sample, lib_size, detected)], by="sample")
# Scatter of two columns colored by a covariate (columns copied so the loop does not
# depend on ggplot's lazy evaluation)
scatter <- function(dt, xa, ya, cv, xlab, ylab, title) {
  d <- data.table(x=dt[[xa]], y=dt[[ya]], col=dt[[cv]])
  num <- is.numeric(d$col)
  if (!num) d[, col := as.character(col)]
  p <- ggplot(d, aes(x=x, y=y, color=col)) + geom_point(size=1, alpha=0.7) + theme_bw() +
    labs(x=xlab, y=ylab, color=cv, title=title)
  if (num) p <- p + scale_color_viridis_c(trans=if (cv == "lib_size") "log10" else "identity")
  p
}
for (pp in list(c(1, 2), c(3, 4), c(5, 6))) {
  xa <- paste0("PC", pp[1]); ya <- paste0("PC", pp[2])
  if (!all(c(xa, ya) %in% colnames(pd))) next
  for (cv in unique(c(obj$model_vars, "lib_size"))) {
    add(scatter(pd, xa, ya, cv, sprintf("%s (%.1f%%)", xa, ve[pp[1]]), sprintf("%s (%.1f%%)", ya, ve[pp[2]]),
                sprintf("Expression PCA before covariate removal, colored by %s", cv)))
  }
}

### 5) Covariate associations with PCs (R^2 heatmaps) ###
assoc_r2 <- function(pcs, dt, vars) {
  rbindlist(lapply(vars, function(v) {
    x <- dt[[v]]
    rbindlist(lapply(colnames(pcs), function(pc) {
      f <- if (is.numeric(x)) lm(pcs[, pc] ~ x) else lm(pcs[, pc] ~ factor(x))
      data.table(covariate=v, PC=pc, r2=summary(f)$r.squared)
    }))
  }))
}
dt.a <- merge(data.table(sample=obj$samples), meta, by="sample", sort=FALSE)
dt.a <- merge(dt.a, qc[, .(sample, lib_size, detected)], by="sample", sort=FALSE)
dt.a[, log_lib := log10(lib_size)]
vars <- unique(c(obj$model_vars, "log_lib", "detected"))

geno <- NULL
if (!is.null(args$geno_pcs) && file.exists(args$geno_pcs)) {
  geno <- fread(args$geno_pcs)
  setnames(geno, sub("^#", "", colnames(geno)))
  if ("FID" %in% colnames(geno) && "IID" %in% colnames(geno)) geno[, FID := NULL]
  setnames(geno, "IID", "sample")
  gcols <- grep("^PC", colnames(geno), value=TRUE)[1:min(5, sum(grepl("^PC", colnames(geno))))]
  setnames(geno, gcols, paste0("geno_", gcols))
  dt.a <- merge(dt.a, geno[, c("sample", paste0("geno_", gcols)), with=FALSE], by="sample", all.x=TRUE, sort=FALSE)
  vars <- c(vars, paste0("geno_", gcols))
}
setkey(dt.a, NULL)
dt.a <- dt.a[match(obj$samples, sample)]

npc <- min(15, ncol(pca$x))
raw.pcs <- pca$x[obj$samples, 1:npc, drop=FALSE]
res.pcs <- obj$expr_pcs[obj$samples, 1:min(15, ncol(obj$expr_pcs)), drop=FALSE]
heat <- rbind(data.table(stage="Before covariate removal", assoc_r2(raw.pcs, dt.a, vars)),
              data.table(stage="Residual PCs (covariates in model)", assoc_r2(res.pcs, dt.a, vars)))
heat[, PC := factor(sub("expr_", "", PC), levels=paste0("PC", 1:15))]
fwrite(heat, paste0(prefix, ".pc_covariate_r2.tsv"), sep="\t")
add(ggplot(heat, aes(x=PC, y=covariate, fill=r2)) + geom_tile() +
      geom_text(aes(label=ifelse(r2 >= 0.05, sprintf("%.2f", r2), "")), size=2.5, color="white") +
      scale_fill_viridis_c(limits=c(0, 1)) + facet_wrap(~stage, ncol=1) + theme_bw() +
      labs(x=NULL, y=NULL, fill="R^2",
           title="Variance in expression PCs explained by each covariate",
           subtitle="Residual PCs should show ~0 for modeled covariates; strong hits elsewhere = hidden structure"))

### 6) Residual PCs ###
rv <- 100 * obj$resid_sdev^2 / sum(obj$resid_sdev^2)
rvd <- data.table(PC=seq_len(min(100, length(rv))))
rvd[, var := rv[PC]][, cum := cumsum(var)]
add(ggplot(rvd, aes(x=PC)) + geom_col(aes(y=var), fill="grey60") +
      geom_line(aes(y=cum / max(cum) * max(var)), color="red") + theme_bw() +
      scale_y_continuous(sec.axis=sec_axis(~ . / max(rvd$var) * max(rvd$cum), name="Cumulative % (red)")) +
      labs(x="Residual expression PC", y="% variance", title="Residual expression PCs 1-100 (candidate covariates)"))

rp <- data.table(sample=obj$samples, PC1=res.pcs[, 1], PC2=res.pcs[, min(2, ncol(res.pcs))])
rp <- merge(rp, meta, by="sample")
for (cv in obj$model_vars) {
  add(scatter(rp, "PC1", "PC2", cv, "Residual PC1", "Residual PC2",
              sprintf("Residual PCs colored by %s (should be mixed)", cv)))
}

### 7) Variance explained per gene by each modeled covariate ###
vp <- rbindlist(lapply(obj$model_vars, function(v) {
  x <- dt.a[[v]]
  X <- if (is.numeric(x)) cbind(1, x) else model.matrix(~ factor(x))
  fit <- lm.fit(X, t(obj$expr))
  r2 <- 1 - colSums(fit$residuals^2) / colSums(sweep(t(obj$expr), 2, colMeans(t(obj$expr)))^2)
  data.table(covariate=v, gene_id=rownames(obj$expr), r2=r2)
}))
vp <- rbind(vp, data.table(covariate=paste("Full:", obj$formula), gene_id=names(obj$r2_known), r2=obj$r2_known))
add(ggplot(vp, aes(x=r2)) + geom_histogram(bins=60) + facet_wrap(~covariate, scales="free_y", ncol=1) +
      theme_bw() + labs(x="Per-gene R^2", y="Genes", title="Expression variance explained per gene"))

# Top covariate-responsive genes (annotated)
gm <- if (!is.null(args$gene_map) && file.exists(args$gene_map) && file.size(args$gene_map) > 0) fread(args$gene_map) else NULL
top.vp <- vp[order(-r2), head(.SD, 25), by=covariate]
if (!is.null(gm)) top.vp <- merge(top.vp, gm[, .(gene_id=query_fbgn, gene_symbol)], by="gene_id", all.x=TRUE)
fwrite(top.vp[order(covariate, -r2)], paste0(prefix, ".top_genes_by_covariate.tsv"), sep="\t")

### 8) Genotype PCs ###
if (!is.null(geno)) {
  gp <- merge(geno, meta, by="sample")
  ge <- NULL
  ev <- sub("\\.eigenvec$", ".eigenval", args$geno_pcs)
  if (file.exists(ev)) ge <- scan(ev, quiet=TRUE)
  lab <- function(i) if (is.null(ge)) paste0("PC", i) else sprintf("PC%d (eigenvalue %.1f)", i, ge[i])
  for (cv in obj$model_vars) {
    add(scatter(gp, "geno_PC1", "geno_PC2", cv, lab(1), lab(2), sprintf("Genotype PCA colored by %s", cv)))
  }
}

pdf(paste0(prefix, ".expression_diagnostics.pdf"), width=9, height=6.5)
for (p in pages) print(p)
invisible(dev.off())
cat("Wrote", length(pages), "diagnostic pages\n")
