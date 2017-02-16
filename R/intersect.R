#' A function to intersect user region data with annotation data
#'
#' Annotate genomic regions to selected genomic annotations while preserving the data associated with the genomic regions.
#'
#' @param regions The GRanges object returned by \code{read_regions()}.
#' @param annotations A character vector of annotations to build. Valid annotation codes are listed with \code{builtin_annotations()}. The "basicgenes" shortcut builds the following regions: 1-5Kb upstream of TSSs, promoters, 5UTRs, exons, introns, and 3UTRs. The "cpgs" shortcut builds the following regions: CpG islands, shores, shelves, and interCGI regions. NOTE: Shortcuts need to be appended by the genome, e.g. \code{hg19_basicgenes}.
#' Custom annotations whose names are of the form \code{[genome]_custom_[name]} should also be included. Custom annotations should be read in and converted to \code{GRanges} with \code{read_annotations()}. They can be for a \code{supported_genome()}, or for an unsupported genome.
#' @param minoverlap A scalar, positive integer, indicating the minimum required overlap of regions with annotations.
#' @param ignore.strand Logical indicating whether strandedness should be respected in findOverlaps(). Default FALSE.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{GRanges} where the \code{granges} are from the regions, and the \code{mcols} include the \code{mcols} from the regions and a column with the annotation \code{GRanges}.
#'
#' @examples
#'    ### An example using Gm12878 Stat3 ChIP-seq from ENCODE and an annotation shortcut
#'    file = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')
#'    r = read_regions(con = file, genome = 'hg19')
#'
#'    # Select and build annotations
#'    annots = 'hg19_cpg_islands'
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
#'
#'    a = annotate_regions(
#'        regions = r,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#'    ### An example with built-in and custom annotations
#'    a = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
#'    read_annotations(con = a, name = 'TFBS', genome = 'hg19')
#'
#'    r_file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
#'    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
#'    r = read_regions(con = r_file, genome = 'hg19', extraCols = extraCols, rename_score = 'coverage')
#'
#'    # Select and build annotations
#'    annots = c('hg19_cpg_islands','hg19_custom_TFBS')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
#'
#'    a = annotate_regions(
#'        regions = r,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#' @export
annotate_regions = function(regions, annotations, minoverlap = 1L, ignore.strand = TRUE, quiet = FALSE) {
    # Checks before moving forward
    if(class(regions)[1] != "GRanges") {
        stop('Error in annotate_regions(...): regions object is not GRanges.')
    }

    if(class(annotations)[1] != "GRanges") {
        stop('Error in annotate_regions(...): annotations object is not GRanges. Use build_annotations(...) to construct the annotations before calling annotate_regions(...).')
    }

    # Perform the intersections
    if(!quiet) {
        message('Annotating...')
    }

    intersections = GenomicRanges::findOverlaps(regions, annotations, minoverlap = minoverlap, ignore.strand = ignore.strand)

    if(length(intersections) > 0) {
        gr = regions[S4Vectors::queryHits(intersections)]
        GenomicRanges::mcols(gr)$annot = annotations[S4Vectors::subjectHits(intersections)]
        return(gr)
    } else {
        stop('No annotations intersect the regions.')
    }
}
