#' Summarize annotation counts
#'
#' Given a \code{GRanges} of annotated regions, count the number of regions in each annotation type. If \code{annotated_random} is not missing, then the same is computed for the background regions, labeled "Background" in the \code{data_type} column.
#'
#' If a region is annotated to multiple annotations of the same \code{annot.type}, the region will only be counted once. For example, if a region were annotated to multiple exons, it would only count once toward the exons, but if it were annotated to an exon and an intron, it would count towards both.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}.
#' @param annotated_random The \code{GRanges} result of \code{annotate_regions()} on a background of regions, i.e. the regions your data could have come from (e.g. all tested CpGs when \code{annotated_regions} are differentially methylated CpGs). Despite the name, this should not be randomized regions; \code{randomize_regions()} is deprecated.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{tbl_df} of the number of regions per annotation type.
#'
#' @examples
#'    ### An example of differentially methylated (DM) regions compared to
#'    ### all regions tested for differential methylation
#'
#'    # Get premade CpG annotations
#'    data('annotations', package = 'annotatr')
#'
#'    file = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#'    extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
#'    r = read_regions(con = file, genome = 'hg19', extraCols = extraCols,
#'        rename_score = 'pval', rename_name = 'DM_status', format = 'bed')
#'
#'    # Annotate all tested regions, which are the background
#'    tested_annots = annotate_regions(
#'        regions = r,
#'        annotations = annotations,
#'        ignore.strand = TRUE,
#'        quiet = FALSE)
#'
#'    # The data are the DM regions
#'    dm_annots = tested_annots[tested_annots$DM_status != 'none']
#'
#'    # Summarize the annotated DM regions
#'    s = summarize_annotations(annotated_regions = dm_annots)
#'
#'    # Summarize the annotated DM regions and the background
#'    s_bg = summarize_annotations(
#'        annotated_regions = dm_annots,
#'        annotated_random = tested_annots)
#'
#' @export
summarize_annotations = function(annotated_regions, annotated_random, quiet = FALSE) {
    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_regions = as.data.frame(annotated_regions, row.names = NULL)

    ########################################################################
    # If a region has multiple annotation types that are the same, count only one
    # from each type of annotation
    annotated_regions = dplyr::distinct(
        dplyr::ungroup(annotated_regions),
        across(dplyr::all_of(c('seqnames', 'start', 'end', 'annot.type'))), .keep_all=TRUE)

    # Tally over data and background regions if annotated_random isn't missing,
    # otherwise tally over data only
    if(!missing(annotated_random)) {
        # Tidy the GRanges into a tbl_df for use with dplyr functions
        annotated_random = as.data.frame(annotated_random, row.names = NULL)

        # If a region has multiple annotation types that are the same, count only one
        # from each type of annotation
        annotated_random = dplyr::distinct(
            dplyr::ungroup(annotated_random),
            across(dplyr::all_of(c('seqnames', 'start', 'end', 'annot.type'))), .keep_all=TRUE)

        if(!quiet) {
            message('Counting annotation types in data and background regions')
        }

        combined_annots = dplyr::bind_rows('Data' = annotated_regions, 'Background' = annotated_random, .id = 'data_type')

        agg = dplyr::tally(
            dplyr::group_by(combined_annots, across(dplyr::all_of(c('data_type', 'annot.type'))))
        )
    } else {
        if(!quiet) {
            message('Counting annotation types')
        }

        # Tally over the normal data
        agg = dplyr::tally(
            dplyr::group_by(annotated_regions, across(dplyr::all_of(c('annot.type'))))
        )
    }

    return(agg)
}

#' Summarize numerical data over groupings of annotated regions
#'
#' Given a \code{GRanges} of annotated regions, summarize numerical data columns based on a grouping.
#'
#' NOTE: We do not take the distinct values of \code{seqnames}, \code{start}, \code{end}, \code{annot.type} as in the other \code{summarize_*()} functions because in the case of a region that intersected two distinct exons, using \code{distinct()} would destroy the information of the mean of the numerical column over one of the exons, which is not desirable.
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
#' # Get premade CpG annotations
#' data('annotations', package = 'annotatr')
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
    annotated_regions = as.data.frame(annotated_regions, row.names = NULL)

    if(missing(over)) {
        stop("Error: over cannot be missing.")
    }

    if(!quiet) {
        message(sprintf('Grouping regions by %s, and summarizing numerical data over %s',
            paste(by, collapse=' & '), paste(over, collapse=' & ')))
    }
    # The columns are n, mean, and sd for one over column, or [column]_n,
    # [column]_mean, and [column]_sd for several
    if(length(over) == 1) {
        names_pattern = '{.fn}'
    } else {
        names_pattern = '{.col}_{.fn}'
    }
    agg = dplyr::summarize(
        dplyr::group_by(annotated_regions, across(dplyr::all_of(by))),
        across(dplyr::all_of(over), list(n = length, mean = mean, sd = stats::sd), .names = names_pattern))
    # Order the columns by statistic, then by over column
    if(length(over) > 1) {
        agg = dplyr::select(agg, dplyr::all_of(c(by, paste0(over, '_n'), paste0(over, '_mean'), paste0(over, '_sd'))))
    }

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
#'    # Get premade CpG annotations
#'    data('annotations', package = 'annotatr')
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
    annotated_regions = as.data.frame(annotated_regions, row.names = NULL)

    ########################################################################
    # If a region has multiple annotation types that are the same, count only one
    # from each type of annotation
    annotated_regions = dplyr::distinct(
        dplyr::ungroup(annotated_regions),
        across(dplyr::all_of(c('seqnames', 'start', 'end', by))), .keep_all=TRUE)

    if(!quiet) {
        message(sprintf('Grouping regions by %s, and tallying',
            paste(by, collapse=' & ')))
    }

    agg = dplyr::tally(
        dplyr::group_by(annotated_regions, across(dplyr::all_of(by))))

    return(agg)
}

#' Summarize annotated regions by gene
#'
#' Given a \code{GRanges} of annotated regions, summarize the regions annotated to each gene, with one row per gene. Only gene annotations with a gene ID count (e.g. promoters, 1-5kb upstream, UTRs, exons, and introns), so CpG, intergenic, enhancer, and chromatin annotations are left out.
#'
#' A region counts once toward a gene, however many of the gene's annotations it overlaps, and once toward each annotation type of the gene. A region annotated to more than one gene counts toward each of them.
#'
#' MANE annotations (\code{hg38_mane_*}) and annotations from \code{build_txdb_annotations()} (e.g. \code{mm39_gencode_*}) are summarized the same way. With annotations from more than one group, the annotation types of groups other than \code{genes} are prefixed with the group, e.g. \code{n_promoters} and \code{n_mane_promoters}. \code{hg38_genes_*} and \code{hg38_mane_*} both use Entrez gene IDs, so a gene has one row. Groups with different kinds of gene IDs (e.g. Entrez and Ensembl) give separate rows for the same gene.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}, with gene annotations such as \code{[genome]_basicgenes} or \code{hg38_basicmane}.
#' @param over A character vector of numerical data columns to summarize with the mean, median, and standard deviation over each gene's regions. Default \code{NULL}, no numerical summaries.
#' @param by A single categorical data column to count the categories of over each gene's regions, e.g. \code{'DM_status'}. Default \code{NULL}, no category counts.
#' @param format Either \code{'wide'} (the default) for one row per gene with a count column per annotation type, or \code{'long'} for one row per gene and annotation type.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{tbl_df} with columns \code{gene_id} and \code{symbol}, then for \code{format = 'wide'}, \code{n_regions} and \code{n_[type]} for each annotation type (e.g. \code{n_promoters}), or for \code{format = 'long'}, \code{annot.type} and \code{n}. These are followed by \code{n_[category]} for each category in \code{by}, and \code{[column]_mean}, \code{[column]_median}, and \code{[column]_sd} for each column in \code{over}. Genes with the most regions come first.
#'
#' @examples
#'  if(requireNamespace('TxDb.Hsapiens.UCSC.hg19.knownGene', quietly = TRUE) &&
#'      requireNamespace('org.Hs.eg.db', quietly = TRUE)) {
#'    # Build hg19 promoter annotations
#'    annotations = build_annotations(genome = 'hg19', annotations = 'hg19_genes_promoters')
#'
#'    dm_file = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#'    extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
#'    dm_regions = read_regions(con = dm_file, extraCols = extraCols, genome = 'hg19',
#'        rename_score = 'pval', rename_name = 'DM_status', format = 'bed')
#'
#'    dm_annots = annotate_regions(
#'        regions = dm_regions,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#'    # One row per gene, with the mean methylation difference and the
#'    # number of hyper- and hypomethylated regions
#'    genes = summarize_genes(annotated_regions = dm_annots, over = 'diff_meth', by = 'DM_status')
#'
#'    # One row per gene and annotation type
#'    genes_long = summarize_genes(annotated_regions = dm_annots, over = 'diff_meth', format = 'long')
#'  }
#'
#' @export
summarize_genes = function(annotated_regions, over = NULL, by = NULL, format = c('wide', 'long'), quiet = FALSE) {
    format = match.arg(format)

    if(!is.null(by) && length(by) != 1) {
        stop('Error: by must be a single column name.')
    }

    # Tidy the GRanges into a data.frame for use with dplyr functions
    tbl = as.data.frame(annotated_regions, row.names = NULL)

    missing_cols = setdiff(c(over, by), colnames(tbl))
    if(length(missing_cols) > 0) {
        stop(sprintf('Error: %s not column(s) in annotated_regions.', paste(missing_cols, collapse = ', ')))
    }

    # Keep the gene annotations with a gene ID: the genes and mane groups, and
    # the groups from build_txdb_annotations(), [genome]_[group]_[type]
    tokens = strsplit(as.character(tbl$annot.type), '_')
    group = vapply(tokens, function(t) if(length(t) == 3) t[2] else NA_character_, character(1))
    gene_type = vapply(tokens, function(t) if(length(t) == 3) t[3] else NA_character_, character(1))
    is_gene = !is.na(group) & (group %in% c('genes', 'mane') | !(group %in% BUILTIN_GROUPS)) & gene_type %in% GENE_TYPES
    keep = !is.na(tbl$annot.gene_id) & is_gene
    tbl = tbl[keep, , drop = FALSE]
    group = group[keep]
    gene_type = gene_type[keep]
    if(nrow(tbl) == 0) {
        stop('Error: No regions are annotated to genes. Include gene annotations, e.g. [genome]_basicgenes or hg38_basicmane, in build_annotations().')
    }

    if(!quiet) {
        message(sprintf('Summarizing %s regions over %s genes', nrow(dplyr::distinct(tbl, .data$seqnames, .data$start, .data$end)), length(unique(tbl$annot.gene_id))))
    }

    tbl$gene_id = as.character(tbl$annot.gene_id)
    tbl$region = paste(tbl$seqnames, tbl$start, tbl$end, sep = ':')
    # Drop the genome and group prefix, e.g. hg19_genes_promoters to promoters.
    # With more than one group, keep the prefix of groups other than genes,
    # e.g. promoters and mane_promoters.
    groups = unique(group)
    groups = c(intersect('genes', groups), setdiff(groups, 'genes'))
    if(length(groups) == 1) {
        tbl$annot.type = gene_type
    } else {
        tbl$annot.type = ifelse(group == 'genes', gene_type, paste(group, gene_type, sep = '_'))
    }
    type_levels = unlist(lapply(groups, function(g) {
        if(g == 'genes' || length(groups) == 1) GENE_TYPES else paste(g, GENE_TYPES, sep = '_')
    }))

    # A gene's symbol, from any of its annotations
    symbols = dplyr::summarize(
        dplyr::group_by(tbl, .data$gene_id),
        symbol = dplyr::first(stats::na.omit(.data$annot.symbol), default = NA_character_))

    # Count each region once per group (gene, or gene and annotation type), and
    # summarize the data columns
    summarize_group = function(regions, keys) {
        regions = dplyr::distinct(regions, dplyr::across(dplyr::all_of(c(keys, 'region'))), .keep_all = TRUE)
        grouped = dplyr::group_by(regions, dplyr::across(dplyr::all_of(keys)))

        agg = dplyr::summarize(grouped, n = dplyr::n(), .groups = 'drop')
        if(!is.null(by)) {
            counts = dplyr::count(regions, dplyr::across(dplyr::all_of(c(keys, by))))
            counts = reshape2::dcast(counts, stats::as.formula(sprintf('%s ~ `%s`', paste(keys, collapse = ' + '), by)),
                value.var = 'n', fill = 0L)
            category_cols = setdiff(colnames(counts), keys)
            counts[category_cols] = lapply(counts[category_cols], as.integer)
            colnames(counts)[match(category_cols, colnames(counts))] = paste0('n_', category_cols)
            agg = dplyr::left_join(agg, counts, by = keys)
        }
        if(length(over) > 0) {
            stats = dplyr::summarize(grouped,
                dplyr::across(dplyr::all_of(over),
                    list(mean = ~ mean(.x, na.rm = TRUE), median = ~ stats::median(.x, na.rm = TRUE), sd = ~ stats::sd(.x, na.rm = TRUE)),
                    .names = '{.col}_{.fn}'),
                .groups = 'drop')
            agg = dplyr::left_join(agg, stats, by = keys)
        }

        return(agg)
    }

    if(format == 'wide') {
        agg = summarize_group(tbl, 'gene_id')
        colnames(agg)[colnames(agg) == 'n'] = 'n_regions'

        # One count column per annotation type, in genomic order
        type_counts = dplyr::count(dplyr::distinct(tbl, .data$gene_id, .data$annot.type, .data$region), .data$gene_id, .data$annot.type)
        type_counts = reshape2::dcast(type_counts, gene_id ~ annot.type, value.var = 'n', fill = 0L)
        types = c(intersect(type_levels, colnames(type_counts)), setdiff(colnames(type_counts), c('gene_id', type_levels)))
        type_counts = type_counts[, c('gene_id', types), drop = FALSE]
        type_counts[types] = lapply(type_counts[types], as.integer)
        colnames(type_counts) = c('gene_id', paste0('n_', types))

        agg = dplyr::left_join(agg, type_counts, by = 'gene_id')
        agg = agg[, c('gene_id', 'n_regions', paste0('n_', types), setdiff(colnames(agg), c('gene_id', 'n_regions', paste0('n_', types))))]
        count_col = 'n_regions'
    } else {
        agg = summarize_group(tbl, c('gene_id', 'annot.type'))
        agg$annot.type = factor(agg$annot.type, levels = c(intersect(type_levels, agg$annot.type), setdiff(agg$annot.type, type_levels)))
        count_col = 'n'
    }

    agg = dplyr::left_join(symbols, agg, by = 'gene_id')
    if(format == 'wide') {
        agg = dplyr::arrange(agg, dplyr::desc(.data[[count_col]]), .data$symbol, .data$gene_id)
    } else {
        # Genes with the most regions first, then annotation types in genomic order
        gene_n = dplyr::count(dplyr::distinct(tbl, .data$gene_id, .data$region), .data$gene_id, name = 'gene_n')
        agg = dplyr::left_join(agg, gene_n, by = 'gene_id')
        agg = dplyr::arrange(agg, dplyr::desc(.data$gene_n), .data$symbol, .data$gene_id, .data$annot.type)
        agg$gene_n = NULL
        agg$annot.type = as.character(agg$annot.type)
    }

    return(dplyr::as_tibble(agg))
}
