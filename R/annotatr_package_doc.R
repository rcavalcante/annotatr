#' annotatr: Annotation of Genomic Regions to Functional Annotations
#'
#' Given a set of genomic sites/regions (e.g. ChIP-seq peaks, CpGs, differentially methylated CpGs or regions, SNPs, etc.) it is often of interest to investigate the intersecting functional annotations. Such annotations include those relating to gene models (promoters, 5'UTRs, exons, introns, and 3'UTRs), CpGs (CpG islands, CpG shores, CpG shelves), the non-coding genome, and enhancers. The annotatr package provides an easy way to summarize and visualize the intersection of genomic sites/regions with the above functional annotations.
#'
#' @docType package
#' @name annotatr
#'
#' @import AnnotationDbi
#' @import AnnotationHub
#' @import dplyr
#' @import ggplot2
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors endoapply
#' @import methods
#' @importFrom readr read_tsv
#' @importFrom reshape2 melt
#' @importFrom regioneR randomizeRegions
#' @import rtracklayer
#' @importClassesFrom S4Vectors Hits Rle
#' @importFrom stats as.formula
#' @importFrom utils combn data
NULL
