#' Summarize annotation counts
#'
#' Given a \code{GRanges} of annotated regions, count the number of regions in each annotation type. If \code{annotated_random} is not \code{NULL}, then the same is computed for the random regions.
#'
#' If a region is annotated to multiple annotations of the same \code{annot.type}, the region will only be counted once. For example, if a region were annotated to multiple exons, it would only count once toward the exons, but if it were annotated to an exon and an intron, it would count towards both.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}.
#' @param annotated_random The \code{GRanges} result of \code{annotate_regions()} on the randomized regions created from \code{randomize_regions()}.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{tbl_df} of the number of regions per annotation type.
#'
#' @examples
#'    ### An example of ChIP-seq peaks with signalValue
#'
#'    # Select and build annotations
#'    annots = c('hg19_cpg_islands','hg19_cpg_shores')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
#'
#'    file = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')
#'    r = read_regions(con = file, genome = 'hg19')
#'
#'    a = annotate_regions(
#'        regions = r,
#'        annotations = annotations,
#'        ignore.strand = TRUE,
#'        quiet = FALSE)
#'
#'    rnd = randomize_regions(regions = r)
#'
#'    rnd_annots = annotate_regions(
#'        regions = rnd,
#'        annotations = annotations,
#'        ignore.strand = TRUE,
#'        quiet = FALSE)
#'
#'    # Summarize the annotated regions without randomized regions
#'    s = summarize_annotations(annotated_regions = a)
#'
#'    # Summarize the annotated regions with randomized regions
#'    s_rnd = summarize_annotations(
#'        annotated_regions = a,
#'        annotated_random = rnd_annots)
#'
#' @export
summarize_annotations = function(annotated_regions, annotated_random, quiet = FALSE) {
    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_regions = as.data.frame(annotated_regions)

    ########################################################################
    # If a region has multiple annotation types that are the same, count only one
    # from each type of annotation
    annotated_regions = dplyr::distinct_(
        dplyr::ungroup(annotated_regions),
        .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)

    # Tally over data and random regions if annotated_random isn't null,
    # otherwise tally over data only
    if(!missing(annotated_random)) {
        # Tidy the GRanges into a tbl_df for use with dplyr functions
        annotated_random = as.data.frame(annotated_random)

        # If a region has multiple annotation types that are the same, count only one
        # from each type of annotation
        annotated_random = dplyr::distinct_(
            dplyr::ungroup(annotated_random),
            .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)

        if(!quiet) {
            message('Counting annotation types in data and random regions')
        }

        combined_annots = dplyr::bind_rows('Data' = annotated_regions, 'Random Regions' = annotated_random, .id = 'data_type')

        agg = dplyr::tally(
            dplyr::group_by_(combined_annots, .dots=c('data_type', 'annot.type'))
        )
    } else {
        if(!quiet) {
            message('Counting annotation types')
        }

        # Tally over the normal data
        agg = dplyr::tally(
            dplyr::group_by_(annotated_regions, .dots=c('annot.type'))
        )
    }

    return(agg)
}

#' Summarize numerical data over groupings of annotated regions
#'
#' Given a \code{GRanges} of annotated regions, summarize numerical data columns based on a grouping.
#'
#' NOTE: We do not take the distinct values of \code{seqnames}, \code{start}, \code{end}, \code{annot.type} as in the other \code{summarize_*()} functions because in the case of a region that intersected two distinct exons, using \code{distinct_()} would destroy the information of the mean of the numerical column over one of the exons, which is not desirable.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}.
#' @param by A character vector of the columns of \code{as.data.frame(annotated_regions)} to group over. Default is \code{c(annot.type, annot.id)}.
#' @param over A character vector of the numerical columns in \code{as.data.frame(annotated_regions)} to \code{count}, take the \code{mean}, and take the \code{sd} over after grouping according to the \code{by} column. NOTE: If more than one value is used, the naming scheme for the resuling \code{dplyr::tbl} summary columns are \code{COLNAME_n}, \code{COLNAME_mean}, \code{COLNAME_sd}. If \code{over} has length one, then the column names are \code{n}, \code{mean}, \code{sd}.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A grouped \code{dplyr::tbl_df}, and the \code{count}, \code{mean}, and \code{sd} of the \code{cols} \code{by} the groupings.
#'
#' @examples
#' ### Test on a very simple bed file to demonstrate different options
#'
#' # Select and build annotations
#' annots = c('hg19_cpg_islands', 'hg19_genes_promoters')
#' annotations = build_annotations(genome = 'hg19', annotations = annots)
#'
#' r_file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
#' extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
#' r = read_regions(con = r_file, genome = 'hg19', extraCols = extraCols, rename_score = 'coverage')
#'
#' a = annotate_regions(
#'        regions = r,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#' # Testing over normal by
#' sn1 = summarize_numerical(
#'        annotated_regions = a,
#'        by = c('annot.type', 'annot.id'),
#'        over = c('coverage', 'mu1', 'mu0'),
#'        quiet = FALSE)
#'
#' # Testing over a different by
#' sn2 = summarize_numerical(
#'        annotated_regions = a,
#'        by = c('diff_exp'),
#'        over = c('coverage', 'mu1', 'mu0'))
#'
#' @export
summarize_numerical = function(annotated_regions, by = c('annot.type', 'annot.id'), over, quiet = FALSE) {
    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_regions = as.data.frame(annotated_regions)

    if(missing(over)) {
        stop("Error: over cannot be missing.")
    }

    if(!quiet) {
        message(sprintf('Grouping regions by %s, and summarizing numerical data over %s',
            paste(by, collapse=' & '), paste(over, collapse=' & ')))
    }
    agg = dplyr::summarize_each_(
        dplyr::group_by_(annotated_regions, .dots = by),
        dplyr::funs(n(), 'mean', 'sd'),
        over)

    return(agg)
}

#' Summarize categorical data over groupings of annotated regions
#'
#' Given a \code{GRanges} of annotated regions, count the number of regions when the annotations are grouped \code{by} categorical columns.
#'
#' If a region is annotated to multiple annotations of the same \code{annot.type}, the region will only be counted once. For example, if a region were annotated to multiple exons, it would only count once toward the exons, but if it were annotated to an exon and an intron, it would count towards both.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}.
#' @param by A character vector to group the data in \code{as.data.frame(annotated_regions)} by and tally over. Default is \code{c('annot.type', 'annot.id')}.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A grouped \code{dplyr::tbl_df} of the counts of groupings according to the \code{by} vector.
#'
#' @examples
#'
#'    # Select and build annotations
#'    annots = c('hg19_cpg_islands', 'hg19_genes_promoters')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
#'
#'    r_file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
#'    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
#'    r = read_regions(con = r_file, genome = 'hg19', extraCols = extraCols, rename_score = 'coverage')
#'
#'    a = annotate_regions(
#'        regions = r,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#'    sc = summarize_categorical(
#'        annotated_regions = a,
#'        by = c('annot.type', 'name'),
#'        quiet = FALSE)
#'
#' @export
summarize_categorical = function(annotated_regions, by = c('annot.type', 'annot.id'), quiet = FALSE) {
    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_regions = as.data.frame(annotated_regions)

    ########################################################################
    # If a region has multiple annotation types that are the same, count only one
    # from each type of annotation
    annotated_regions = dplyr::distinct_(
        dplyr::ungroup(annotated_regions),
        .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)

    if(!quiet) {
        message(sprintf('Grouping regions by %s, and tallying',
            paste(by, collapse=' & ')))
    }

    agg = dplyr::tally(
        dplyr::group_by_(annotated_regions, .dots = by))

    return(agg)
}
