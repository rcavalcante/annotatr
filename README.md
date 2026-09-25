[![check](https://github.com/rcavalcante/annotatr/actions/workflows/check.yml/badge.svg?branch=devel)](https://github.com/rcavalcante/annotatr/actions/workflows/check.yml)

# annotatr

`annotatr` annotates genomic regions (e.g. differentially methylated regions, ChIP-seq peaks, or SNPs) to genomic features: genes, CpG islands, enhancers, ENCODE cCREs, lncRNAs, chromatin states, and custom annotations. It then summarizes and plots the annotated regions.

## Installation

`annotatr` is a [Bioconductor package](https://bioconductor.org/packages/annotatr):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("annotatr")
```

## Built-in annotations

Which annotations are available for which genome builds, from `builtin_annotations_table()`:

<!-- Generated with builtin_annotations_table(). Update when annotations change. -->
|genome      | genes | mane | canonical | cpgs | enhancers | chromatin | lncrna | ccres |
|:-----------|:-----:|:----:|:---------:|:----:|:---------:|:---------:|:------:|:-----:|
|dm3         |   ✓   |      |           |      |           |           |        |       |
|dm6         |   ✓   |      |     ✓     |      |           |           |        |       |
|danRer10    |   ✓   |      |           |  ✓   |           |           |        |       |
|danRer11    |   ✓   |      |     ✓     |  ✓   |           |           |        |       |
|galGal5     |   ✓   |      |           |  ✓   |           |           |        |       |
|hg19        |   ✓   |      |           |  ✓   |     ✓     |     ✓     |   ✓    |       |
|hg38        |   ✓   |  ✓   |     ✓     |  ✓   |     ✓     |           |   ✓    |   ✓   |
|mm9         |   ✓   |      |           |  ✓   |     ✓     |           |        |       |
|mm10        |   ✓   |      |           |  ✓   |     ✓     |           |   ✓    |   ✓   |
|mm39        |   ✓   |      |     ✓     |  ✓   |           |           |        |       |
|oviariramb2 |   ✓   |      |     ✓     |  ✓   |           |           |        |       |
|rn4         |   ✓   |      |           |  ✓   |           |           |        |       |
|rn5         |   ✓   |      |           |  ✓   |           |           |        |       |
|rn6         |   ✓   |      |           |  ✓   |           |           |        |       |
|rn7         |   ✓   |      |     ✓     |  ✓   |           |           |        |       |

Custom annotations can be read from BED files with `read_annotations()`, and gene annotations can be built from any `TxDb` or `EnsDb` with `build_txdb_annotations()`.

## Documentation

See the [package vignette](https://bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html) for a fully worked through use case, and `news(package = "annotatr")` for the changes in each version.

## Citation

Cavalcante RG, Sartor MA. annotatr: genomic regions in context. *Bioinformatics* (2017) 33(15):2381-2383. [doi:10.1093/bioinformatics/btx183](https://doi.org/10.1093/bioinformatics/btx183)
