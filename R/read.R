#' Read genomic regions in BEDX+Y format
#'
#' \code{read_regions()} reads genomic regions by calling the \code{rtracklayer::import()} function. This function can automatically deal with BEDX files from BED3 to BED6. For BED6+Y, the \code{extraCols} argument should be used to correctly interpret the extra columns.
#'
#' NOTE: The \code{name} (4th) and \code{score} (5th) columns are so named. If these columns have a particular meaning for your data, they should be renamed with the \code{rename_name} and/or \code{rename_score} parameters.
#'
#' @param con A path, URL, connection or BEDFile object. See \code{rtracklayer::import()} documentation.
#' @param genome From \code{rtracklayer::import()}: The identifier of a genome, or NA if unknown. Typically, this is a UCSC identifier like 'hg19'. An attempt will be made to derive the \code{seqinfo} on the return value using either an installed BSgenome package or UCSC, if network access is available.
#' @param format From \code{rtracklayer::import()}: The format of the output. If not missing, should be one of 'bed', 'bed15', 'bedGraph' or 'bedpe'. If missing and 'con' is a filename, the format is derived from the file extension. This argument is unnecessary when 'con' is a derivative of 'RTLFile'.
#' @param extraCols From \code{rtracklayer::import()}: A character vector in the same form as 'colClasses' from 'read.table'.  It should indicate the name and class of each extra/special column to read from the BED file. As BED does not encode column names, these are assumed to be the last columns in the file. This enables parsing of the various BEDX+Y formats.
#' @param rename_name A string to rename the name column of the BED file. For example, if the name column actually contains a categorical variable.
#' @param rename_score A string to rename the score column of the BED file. For example, if the score column represents a quantity about the data besides the score in the BED specification.
#' @param ... Parameters to pass onto the format-specific method of \code{rtracklayer::import()}.
#'
#' @return A \code{GRanges} object.
#'
#' @examples
#'
#'    # Example of reading a BED3 file
#'    file = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')
#'    gr = read_regions(con = file, genome = 'hg19')
#'
#'    # Example of reading a BED6+3 file where the last 3 columns are non-standard
#'    file = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#'    extraCols = c(diff_meth = 'numeric', mu0 = 'numeric', mu1 = 'numeric')
#'    gr = read_regions(con = file, genome = 'hg19', extraCols = extraCols, format = 'bed',
#'        rename_name = 'DM_status', rename_score = 'pval')
#'
#' @export
read_regions = function(con, genome = NA, format, extraCols = character(), rename_name, rename_score, ...) {

    if(!missing(format)) {
        gr = rtracklayer::import(con = con, genome = genome, format = format, extraCols = extraCols, ...)
    } else {
        gr = rtracklayer::import(con = con, genome = genome, extraCols = extraCols, ...)
    }

    # Rename name and score columns if the user desires
    if(!missing(rename_name)) {
        if(any(colnames(GenomicRanges::mcols(gr)) == 'name')) {
            colnames(GenomicRanges::mcols(gr))[which(colnames(GenomicRanges::mcols(gr)) == 'name')] = rename_name
        } else {
            warning('Ignoring rename_name parameter because con has no name column.')
        }
    }
    if(!missing(rename_score)) {
        if(any(colnames(GenomicRanges::mcols(gr)) == 'score')) {
            colnames(GenomicRanges::mcols(gr))[which(colnames(GenomicRanges::mcols(gr)) == 'score')] = rename_score
        } else {
            warning('Ignoring rename_score parameter because con has no score column.')
        }
    }

    return(gr)
}

#' Read custom annotations
#'
#' \code{read_annotations()} is a wrapper for the \code{rtracklayer::import()} function that creates a \code{GRanges} object matching the structure of annotations built with \code{build_annotations()}. The structure is defined by \code{GRanges}, with the \code{mcols()} with names \code{c('id','gene_id','symbol','type')}.
#'
#' @param con A path, URL, connection or BEDFile object. See \code{rtracklayer::import.bed()} documentation.
#' @param name A string for the name of the annotations to be used in the name of the object, [genome]_custom_[name]
#' @param genome From \code{rtracklayer::import()}: The identifier of a genome, or NA if unknown. Typically, this is a UCSC identifier like 'hg19'. An attempt will be made to derive the \code{seqinfo} on the return value using either an installed BSgenome package or UCSC, if network access is available.
#' @param format From \code{rtracklayer::import()}: The format of the output. If not missing, should be one of 'bed', 'bed15', 'bedGraph' or 'bedpe'. If missing and 'con' is a filename, the format is derived from the file extension. This argument is unnecessary when 'con' is a derivative of 'RTLFile'.
#' @param extraCols From \code{rtracklayer::import.bed()}: A character vector in the same form as 'colClasses' from 'read.table'.  It should indicate the name and class of each extra/special column to read from the BED file. As BED does not encode column names, these are assumed to be the last columns in the file. This enables parsing of the various BEDX+Y formats.
#' @param ... Parameters to pass onto the format-specific method of \code{rtracklayer::import()}.
#'
#' @return A \code{GRanges} object stored in \code{annotatr_cache}. To view a custom annotation, do \code{annotatr_cache$get(name)}. To add a custom annotation to the set of annotations, include \code{'[genome]_custom_[name]'} in the call to \code{build_annotations()}. See example below.
#'
#' @examples
#'
#'  # Read in a BED3 file as a custom annotation
#'  file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
#'  read_annotations(con = file, name = 'test', genome = 'hg19')
#'  annotations = build_annotations(genome = 'hg19', annotations = 'hg19_custom_test')
#'
#' @export
read_annotations = function(con, name, genome = NA, format, extraCols = character(), ...) {

    if(missing(name)) {
        name = 'annotations'
    }
    if(is.na(genome)) {
        genome = 'genome'
    }

    protected_extraCols = c('gene_id','symbol','tx_id')

    if(!missing(format)) {
        gr = rtracklayer::import(con = con, genome = genome, format = format, extraCols = extraCols, ...)
    } else {
        gr = rtracklayer::import(con = con, genome = genome, extraCols = extraCols, ...)
    }

    # Determine whether gene_id or symbol are missing from extraCols
    missing_extraCols = base::setdiff(protected_extraCols, names(extraCols))

    if(any(missing_extraCols == 'gene_id')) {
        GenomicRanges::mcols(gr)$gene_id = NA
    }
    if(any(missing_extraCols == 'symbol')) {
        GenomicRanges::mcols(gr)$symbol = NA
    }
    if(any(missing_extraCols == 'tx_id')) {
        GenomicRanges::mcols(gr)$tx_id = NA
    }

    GenomicRanges::mcols(gr)$id = paste0(name,':',seq_along(gr))
    GenomicRanges::mcols(gr)$type = sprintf('%s_custom_%s', genome, name)

    # Make sure only the desired mcols make it out
    GenomicRanges::mcols(gr) = GenomicRanges::mcols(gr)[,c('id','tx_id','gene_id','symbol','type')]

    ########################################################
    # Write the object named [genome]_custom_[name] to the annotatr_cache
    annotatr_cache$set(sprintf('%s_custom_%s', genome, name), gr)
}
