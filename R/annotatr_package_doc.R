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
#' @importFrom IRanges IRanges endoapply
#' @import methods
#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @importFrom readr read_tsv
#' @importFrom reshape2 melt
#' @importFrom regioneR randomizeRegions
#' @import rtracklayer
#' @importClassesFrom S4Vectors Hits Rle
#' @importFrom stats as.formula
#' @import TxDb.Dmelanogaster.UCSC.dm3.ensGene
#' @import TxDb.Dmelanogaster.UCSC.dm6.ensGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import TxDb.Mmusculus.UCSC.mm9.knownGene
#' @import TxDb.Mmusculus.UCSC.mm10.knownGene
#' @import TxDb.Rnorvegicus.UCSC.rn4.ensGene
#' @import TxDb.Rnorvegicus.UCSC.rn5.refGene
#' @import TxDb.Rnorvegicus.UCSC.rn6.refGene
#' @importFrom utils combn data
NULL
