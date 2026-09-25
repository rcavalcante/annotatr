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
#'    r_file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
#'    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
#'    r = read_regions(con = r_file, extraCols = extraCols, rename_score = 'coverage')
#'
#'    # Get premade CpG annotations
#'    data('annotations', package = 'annotatr')
#'
#'    a = annotate_regions(
#'        regions = r,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#' @export
annotate_regions = function(regions, annotations, minoverlap = 1L, ignore.strand = TRUE, quiet = FALSE) {
    # Checks before moving forward
    if(!methods::is(regions, "GRanges")) {
        stop('Error in annotate_regions(...): regions object is not GRanges.')
    }

    if(!methods::is(annotations, "GRanges")) {
        stop('Error in annotate_regions(...): annotations object is not GRanges. Use build_annotations(...) to construct the annotations before calling annotate_regions(...).')
    }

    # Check the genomes and chromosome names match, with clearer errors than
    # findOverlaps() gives
    check_regions_genome(regions, annotations, quiet = quiet)

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

#' Function to check regions and annotations are from the same genome
#'
#' Gives an error if the regions and annotations have different genomes, or no chromosome names in common (e.g. \code{2} and \code{chr2}). If the regions have no genome, it can't check the genomes match, so it gives a message once per session suggesting \code{read_regions(genome = ...)}.
#'
#' @param regions A \code{GRanges} object of regions.
#' @param annotations A \code{GRanges} object of annotations.
#' @param quiet Print the message about regions without a genome (FALSE) or not (TRUE).
#'
#' @return \code{NULL}, invisibly, if the checks pass.
check_regions_genome = function(regions, annotations, quiet = FALSE) {
    regions_seqlevels = Seqinfo::seqlevelsInUse(regions)
    common = intersect(regions_seqlevels, Seqinfo::seqlevels(annotations))

    if(length(regions_seqlevels) > 0 && length(common) == 0) {
        stop(sprintf(paste(
            'Error: The regions and annotations have no chromosome names in common, e.g. %s in the regions and %s in the annotations.',
            'If only the naming style differs, change the regions to match, e.g. with GenomeInfoDb::seqlevelsStyle(regions) = \'UCSC\' for chr1, chr2, etc.'),
            paste(utils::head(regions_seqlevels, 3), collapse = ', '),
            paste(utils::head(Seqinfo::seqlevels(annotations), 3), collapse = ', ')))
    }

    regions_genome = Seqinfo::genome(regions)[common]
    annotations_genome = Seqinfo::genome(annotations)[common]
    conflicts = !is.na(regions_genome) & !is.na(annotations_genome) & regions_genome != annotations_genome
    if(any(conflicts)) {
        stop(sprintf('Error: The regions are from genome %s but the annotations are from genome %s. Use regions and annotations from the same genome.',
            paste(unique(regions_genome[conflicts]), collapse = ', '),
            paste(unique(annotations_genome[conflicts]), collapse = ', ')))
    }

    annotations_genome = unique(stats::na.omit(annotations_genome))
    if(!quiet && all(is.na(regions_genome)) && length(annotations_genome) > 0) {
        rlang::inform(
            sprintf(paste(
                'The regions have no genome, so annotate_regions() can\'t check they match the annotations (%s).',
                'Use read_regions(genome = ...) to set it.'),
                paste(annotations_genome, collapse = ', ')),
            .frequency = 'once', .frequency_id = 'annotatr_regions_without_genome')
    }

    return(invisible(NULL))
}
