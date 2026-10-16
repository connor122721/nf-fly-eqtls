#!/usr/bin/env python3
"""
Connor Murray
Reduce tensorQTL cis_nominal parquet files to one row per SNP (best gene,
minimum nominal p) for genome-wide Manhattan plots. SNPs with p < --keep_p
are all kept; the rest are randomly thinned to --thin_frac.
"""
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("parquet", nargs="+", help="*.cis_qtl_pairs.*.parquet")
parser.add_argument("--keep_p", type=float, default=1e-3)
parser.add_argument("--thin_frac", type=float, default=0.05)
parser.add_argument("--seed", type=int, default=1)
parser.add_argument("--out", default="nominal_snp_min_p.txt.gz")
args = parser.parse_args()

rng = np.random.default_rng(args.seed)
out = []
for f in args.parquet:
    df = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "pval_nominal", "slope"])
    df = df.sort_values("pval_nominal").drop_duplicates("variant_id")
    keep = (df["pval_nominal"] < args.keep_p).values | (rng.random(len(df)) < args.thin_frac)
    out.append(df[keep])
    print(f"{f}: {len(df)} SNPs -> {keep.sum()} kept")

res = pd.concat(out)
res[["chrom", "pos"]] = res["variant_id"].str.rsplit(":", n=1, expand=True)
res["pos"] = res["pos"].astype(int)
res.to_csv(args.out, sep="\t", index=False, compression="gzip")
