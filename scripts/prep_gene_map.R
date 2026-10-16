#!/usr/bin/env Rscript
# Connor Murray
# FBgn -> gene symbol / CG annotation ID map from FlyBase
# fbgn_annotation_ID_fb_*.tsv.gz. dm3-era FBgns that were later merged/renamed
# are resolved through the secondary FBgn column.

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--annotation", required=TRUE, help="FlyBase fbgn_annotation_ID_*.tsv(.gz)")
parser$add_argument("--genes", required=TRUE, help="tensorQTL expression BED (gene IDs to map)")
parser$add_argument("--out", default="fbgn_to_symbol.tsv", help="Output map")
args <- parser$parse_args()

f <- args$annotation
lines <- if (grepl("\\.gz$", f)) readLines(gzcon(file(f, "rb"))) else readLines(f)
hdr <- grep("^##?gene_symbol", lines, value=TRUE)
body <- lines[!grepl("^#", lines) & lines != ""]
fb <- fread(text=paste(body, collapse="\n"), sep="\t", header=FALSE, fill=TRUE, colClasses="character")

# Column layout: symbol, [organism], primary FBgn, secondary FBgn(s), annotation ID, secondary annotation IDs
if (length(hdr) > 0) {
  setnames(fb, seq_along(strsplit(sub("^#+", "", hdr[1]), "\t")[[1]]),
           make.names(strsplit(sub("^#+", "", hdr[1]), "\t")[[1]]))
}
cn <- colnames(fb)
col_sym <- 1
col_org <- if (any(grepl("organism", cn))) grep("organism", cn)[1] else if (all(fb[[2]][1:min(50, nrow(fb))] %like% "^D[a-z]{3}$")) 2 else NA
col_pri <- which(sapply(fb, function(x) mean(grepl("^FBgn[0-9]+$", head(x, 500))) > 0.9))[1]
col_sec <- col_pri + 1
col_ann <- col_pri + 2
if (!is.na(col_org)) fb <- fb[fb[[col_org]] == "Dmel"]

map <- data.table(gene_symbol=fb[[col_sym]], primary_fbgn=fb[[col_pri]],
                  secondary=fb[[col_sec]], annotation_ID=fb[[col_ann]])

# Long table: every primary and secondary FBgn -> current gene
sec <- map[secondary != "", .(query_fbgn=trimws(unlist(strsplit(secondary, ",")))),
           by=.(gene_symbol, primary_fbgn, annotation_ID)]
long <- rbind(map[, .(query_fbgn=primary_fbgn, gene_symbol, primary_fbgn, annotation_ID, match="primary")],
              sec[, .(query_fbgn, gene_symbol, primary_fbgn, annotation_ID, match="secondary")])
long <- long[order(query_fbgn, match)][!duplicated(query_fbgn)]

genes <- fread(args$genes, select=4)[[1]]
out <- merge(data.table(query_fbgn=genes), long, by="query_fbgn", all.x=TRUE)
out[is.na(match), `:=`(match="unmapped", gene_symbol=query_fbgn)]
fwrite(out, args$out, sep="\t")
cat(sprintf("Mapped %d genes: %d primary, %d secondary (renamed/merged), %d unmapped\n",
            nrow(out), sum(out$match == "primary"), sum(out$match == "secondary"), sum(out$match == "unmapped")))
