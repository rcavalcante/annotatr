#' Randomize Regions
#'
#' \code{randomize_regions} is deprecated and will be removed in a future release. Regions placed uniformly at random across the genome are a poor null model for most data, because real regions come from a much narrower set of possible locations (e.g. mappable sequence, or the CpGs an assay measures). Instead, annotate a background of the regions your data could have come from (e.g. all tested CpGs) with \code{annotate_regions()}, and pass it as \code{annotated_random} to \code{summarize_annotations()}, \code{plot_annotation()}, or \code{plot_categorical()}. For a statistical test, see \code{regioneR::permTest()} or \code{regioneR::resampleRegions()}.
#'
#' \code{randomize_regions} is a wrapper function for \code{regioneR::randomizeRegions()} that simplifies the creation of randomized regions for an input set of regions read with \code{read_regions()}. It relies on the \code{seqlengths} of \code{regions} in order to build the appropriate \code{genome} object for \code{regioneR::randomizeRegions()}.
#'
#' NOTE: The data associated with the input \code{regions} are not passed on to the random regions.
#'
#' @param regions A \code{GRanges} object from \code{read_regions}.
#' @param allow.overlaps A logical stating whether random regions can overlap each other (TRUE) or not (FALSE). Default TRUE.
#' @param per.chromosome A logical stating whether the random regions should remain on the same chromosome (TRUE) or not (FALSE). Default TRUE.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{GRanges} object of randomized regions based on \code{regions} from \code{read_regions()}. NOTE: Data associated with the original regions is not attached to the randomized regions.
#'
#' @seealso \code{\link{annotatr-deprecated}}
#'
#' @examples
#'    # Create random region set based on ENCODE ChIP-seq data
#'    file = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')
#'    r = read_regions(con = file, genome = 'hg19')
#'
#'    random_r = suppressWarnings(randomize_regions(regions = r))
#'
#' @export
randomize_regions = function(regions, allow.overlaps = TRUE, per.chromosome = TRUE, quiet = FALSE) {

    .Deprecated(msg = paste(
        'randomize_regions() is deprecated and will be removed in a future release.',
        'Random regions are a poor null model for most data. Instead, annotate a background',
        'of the regions your data could have come from (e.g. all tested CpGs) and pass it as',
        'annotated_random, or use regioneR::permTest() or regioneR::resampleRegions().',
        'See the "Comparing to a Background" section of the annotatr vignette.'))

    ########################################################################
    # Argument parsing and error handling
    if(!methods::is(regions, "GRanges")) {
        stop('Error: regions must have class GRanges. The best way to ensure this is to pass the result of read_regions() into this function.')
    }

    # Get the genome from the regions
    genome = unique(Seqinfo::genome(regions))

    if(is.na(genome)) {
        stop('Error: regions GRanges object must have a valid genome to randomize its regions.')
    } else {
        chr_lengths = Seqinfo::Seqinfo(genome = genome)
        chr_lengths = Seqinfo::seqlengths(chr_lengths)

        df_genome = data.frame(
            'chr' = names(chr_lengths),
            'start' = rep.int(1, length(chr_lengths)),
            'end' = as.numeric(chr_lengths),
            stringsAsFactors = FALSE)
    }

    if(!quiet) {
        message('Randomizing regions...')
    }

    # Randomize the regions
    randomized = regioneR::randomizeRegions(A = regions, genome = df_genome,
        per.chromosome = per.chromosome, allow.overlaps = allow.overlaps)

    # Sort the randomized
    randomized = sort(randomized)

    return(randomized)
}
