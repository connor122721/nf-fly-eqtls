# nf-fly-eqtls: cis-eQTL mapping in *Drosophila* with tensorQTL

![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-brightgreen?style=for-the-badge)
![License: MIT](https://img.shields.io/badge/License-MIT-blue?style=for-the-badge)
![Profile: Slurm](https://img.shields.io/badge/Profile-Slurm-orange?style=for-the-badge)

A Nextflow DSL2 pipeline for cis-eQTL mapping of single-fly RNA-seq and WGS from a desiccation-tolerance cage experiment (dm3). It is a simplified version of [nf-eqtls](https://github.com/connor122721/nf-eqtls): the sQTL/LeafCutter, kinship, coloc, and TOPMed-specific modules are gone. What's left goes from genotypes and counts to tensorQTL cis-eQTLs on the main chromosomes (2L, 2R, 3L, 3R, X).

## Pipeline

```mermaid
flowchart TD
    V["filtered.all.vcf.gz"] --> S["List VCF samples"]
    V --> P["Genotype PCA<br/>(LD-pruned autosomes)"]
    C["dm3.featureCountMatrix"] --> R["R: sample & gene QC<br/>TMM + inverse-normal<br/>covariates"]
    M["dm3.metadata.tsv"] --> R
    G["dm3.ensGene.gtf"] --> R
    S --> R
    P --> R
    R -- "QC-passing samples" --> B["PLINK bed per chromosome"]
    V --> B
    R -- "expression BED + covariates (K = 0,5,10,...)" --> T["tensorQTL cis (permutations)<br/>chromosome x K"]
    B --> T
    T --> Z["Pool chromosomes, FDR,<br/>pick K with most eGenes"]
    Z --> N["tensorQTL cis_nominal<br/>for best K"]
    B --> N
```

### Expression QC and normalization (`scripts/prep_expression.R`)

1. Keeps the samples that appear in the count matrix, the metadata, and the VCF. Every drop is logged.
2. **Sample QC:** drops samples more than `lib_mad` MADs below the median log10 library size or median number of detected genes. It then drops expression-PCA outliers: Mahalanobis distance on PCs 1–5, Bonferroni p < `pca_outlier_p`.
3. **Gene QC (GTEx-style):** keeps genes with CPM ≥ `cpm_threshold` and counts ≥ `count_threshold` in ≥ `sample_frac` of samples.
4. **Normalization:** edgeR TMM, then log2 CPM, then a rank-based inverse normal transform per gene (`norm_method = "int"`). Use `"logcpm"` to skip the transform.
5. **Phenotype BED:** one TSS per gene (the most upstream transcript start, strand-aware) from the GTF. The `chr` prefix is stripped automatically to match the VCF (`chr2L` → `2L`). Only genes on the chromosomes listed in `chromosomes` are kept.
6. **Covariates:** `model.matrix(covariate_formula)`, by default `~ cage_name + condition`, with `cage_name` and `condition` treated as factors. Optional extras are `n_geno_pcs` genotype PCs and K expression PCs. The expression PCs are computed on the residuals left after regressing out the known covariates, so they don't just rediscover condition/cage. Aliased columns (e.g. cage nested within condition) are detected and dropped with a warning. One covariate file is written per K in `expr_pcs`.

QC outputs: `fly.sample_qc.tsv`, `fly.expression_qc.pdf` (library size, PCA by condition/cage, residual PCA, scree, per-gene R² of the known covariates), and `prep_expression.log`.

### Genotypes (`modules/genotypes.nf`)

- plink2 reads the VCF with `--vcf-half-call missing`, because the VCF contains half calls such as `./0` and `1/.`. It keeps biallelic ACGT SNPs and applies `--geno` (`geno_missing`) and `--maf` filters. Variant IDs are set to `CHR:POS`.
- All samples are set to diploid on X (`--update-sex` as female) so plink2 does not apply haploid handling to unknown-sex flies.
- Genotype PCA runs on LD-pruned autosomal SNPs (`--indep-pairwise 50 5 0.2`), with `pca_geno` and `pca_maf` filters.

### QTL mapping (`modules/tensorqtl.nf`, `scripts/summarize_cis.R`)

- tensorQTL `cis` runs per chromosome × K, using `cis_window` around the TSS (default ±100 kb) and `permutations`.
- The summary step pools all chromosomes for each K, applies BH FDR to `pval_beta`, and counts eGenes. It picks the K with the most eGenes (ties go to fewer PCs) and writes `fly.cis_eqtl.K<best>.all_chroms.txt.gz`.
- tensorQTL `cis_nominal` then runs for that K and writes parquet files per chromosome. You can force a K with `--best_k 10`; it must be one of the values in `expr_pcs`.

## Usage

```bash
cd /project/berglandlab/connor/fly_QTLs_2026
git clone https://github.com/connor122721/nf-fly-eqtls.git

module load nextflow apptainer
nextflow run nf-fly-eqtls/main.nf -profile slurm -resume \
    --out /scratch/$USER/fly_eqtl/output

# Example overrides
nextflow run nf-fly-eqtls/main.nf -profile slurm -resume \
    --covariate_formula "~ cage_name + condition" \
    --n_geno_pcs 3 \
    --expr_pcs "0,5,10,15,20,30" \
    --cis_window 1000000
```

All parameters and their defaults are in `nextflow.config`.

### Inputs

| param | file | notes |
|---|---|---|
| `counts` | `dm3.featureCountMatrix.fixedNames` | first column gene ID (FBgn); featureCounts `Chr/Start/End/Strand/Length` columns are dropped if present |
| `metadata` | `dm3.metadata.tsv` | `sample`, `cage_name`, `fly_number`, `condition`; sample ID must be the first column |
| `vcf` | `filtered.all.vcf.gz` (+ `.tbi`) | chromosome names `2L`, `2R`, ... |
| `gtf` | `dm3.ensGene.gtf` | `gene_id` must match the count matrix |
| `chrom_file` | `chromosomes` | `2L 2R 3L 3R X` |

### Containers

| param | contents |
|---|---|
| `container_wgs` | bcftools, plink2 (`library://connmurr243/wgs/wgs_common`) |
| `container_tensorqtl` | tensorQTL, pandas_plink (`definition/fly_tensorqtl_python.def`, same as the TOPCHeF image) |
| `container_r` | R + edgeR, data.table, ggplot2, argparse (`definition/fly_eqtl_r.def`) |

The R image is new. Build and push it once before the first run:

```bash
apptainer build fly_eqtl_r.sif definition/fly_eqtl_r.def
apptainer push fly_eqtl_r.sif library://connmurr243/wgs/fly_eqtl_r.sif:latest
# or point params.container_r at the local .sif path
```

### Outputs (`params.out`)

```
genotypes/        vcf_samples.txt, fly.geno_pca.eigenvec, bed/fly.<chr>.{bed,bim,fam}
expression/       fly.expression.bed, fly.covariates.K*.txt, fly.samples.txt, QC tables + PDF
tensorqtl_cis/    fly.chr<chr>.K<k>.cis_qtl.txt.gz
summary/          fly.egenes_by_K.{tsv,pdf}, best_k, fly.cis_eqtl.K<best>.all_chroms.txt.gz
tensorqtl_nominal/ fly.chr<chr>.K<best>.cis_qtl_pairs.<chr>.parquet
```

## Notes

- **Cage vs condition:** if each cage appears in only one condition, `condition` is aliased with `cage_name` and is dropped automatically, with a warning in `prep_expression.log`. In that case cage absorbs the condition effect. Check the log to see which columns were kept.
- **Missingness:** tensorQTL mean-imputes missing genotypes. Tighten `geno_missing` if the low-coverage calls are noisy.
- **Cis window:** the default is ±100 kb, chosen because fly genes are compact and LD decays quickly. Increase it with `--cis_window` if you expect longer-range LD in the cage population.

## Citation

Derived from [connor122721/nf-eqtls](https://github.com/connor122721/nf-eqtls). Uses [tensorQTL](https://github.com/broadinstitute/tensorqtl) (Taylor-Weiner et al. 2019, *Genome Biology*), [edgeR](https://bioconductor.org/packages/edgeR), and [PLINK 2](https://www.cog-genomics.org/plink/2.0/).
