#' Plot the number of regions per annotation
#'
#' Given a \code{GRanges} of annotated regions, plot the number of regions with the corresponding genomic annotations used in \code{annotation_order}. If a region is annotated to multiple annotations of the same \code{annot.type}, the region will only be counted once in the corresponding bar plot. For example, if a region were annotated to multiple exons, it would only count once toward the exon bar in the plot, but if it were annotated to an exon and an intron, it would count towards both.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}.
#' @param annotated_random The \code{GRanges} result of \code{annotate_regions()} on the randomized regions created from \code{randomize_regions()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for the plot as well as the ordering. If \code{NULL}, all annotations are displayed.
#' @param plot_title A string used for the title of the plot. If missing, no title is displayed.
#' @param x_label A string used for the x-axis label. If missing, no x-axis label is displayed.
#' @param y_label A string used for the y-axis label. If missing, no y-axis label is displayed.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{ggplot} object which can be viewed by calling it, saved with \code{ggplot2::ggsave}, or edited.
#'
#' @examples
#'    ########################################################################
#'    # An example of ChIP-seq peaks with signalValue used for score
#'
#'    # Select and build annotations
#'    annots = c('hg19_cpg_islands','hg19_cpg_shores')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
#'
#'    chip_bed = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')
#'    chip_regions = read_regions(con = chip_bed, genome = 'hg19')
#'
#'    chip_rnd = randomize_regions(regions = chip_regions)
#'
#'    chip_annots = annotate_regions(
#'        regions = chip_regions,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#'    chip_rnd_annots = annotate_regions(
#'        regions = chip_rnd,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#'    annots_order = c(
#'        'hg19_cpg_islands',
#'        'hg19_cpg_shores')
#'
#'    p_annots = plot_annotation(annotated_regions = chip_annots,
#'        annotation_order = annots_order)
#'    p_annots_rnd = plot_annotation(annotated_regions = chip_annots,
#'        annotated_random = chip_rnd_annots, annotation_order = annots_order)
#'
#' @export
plot_annotation = function(annotated_regions, annotated_random, annotation_order = NULL,
    plot_title, x_label, y_label, quiet = FALSE) {

    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_regions = as.data.frame(annotated_regions)

    ########################################################################
    # Order and subset the annotations
    annotated_regions = subset_order_tbl(tbl = annotated_regions, col='annot.type', col_order=annotation_order)

    ########################################################################
    # If a region has multiple annotation types that are the same, count only one
    # from each type of annotation
    annotated_regions = dplyr::distinct_(
        dplyr::ungroup(annotated_regions),
        .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)

    # Do particular things if annotated_random isn't NULL
    if(!missing(annotated_random)) {
        # Tidy the GRanges into a tbl_df for use with dplyr functions
        annotated_random = as.data.frame(annotated_random)

        # Order and subset the randomized annotations
        annotated_random = subset_order_tbl(tbl = annotated_random, col='annot.type', col_order=annotation_order)

        # If a region has multiple annotation types that are the same, count only one
        # from each type of annotation
        annotated_random = dplyr::distinct_(
            dplyr::ungroup(annotated_random),
            .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)

        # Combine the tbl_dfs in preparation for visualization
        annotated_regions = dplyr::bind_rows("Data" = annotated_regions, "Random Regions" = annotated_random, .id = 'data_type')
    }

    ########################################################################
    # Construct the plot

    # Make the base ggplot
    # NOTE: binwidth may need to be a parameter
    if(missing(annotated_random)) {
        plot =
        ggplot(annotated_regions, aes_string(x='annot.type')) +
            geom_bar() +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = 'white'))
    } else {
        plot =
            ggplot(annotated_regions, aes_string(x='annot.type')) +
            geom_bar(aes_string(fill = 'data_type'), position='dodge') +
            theme_bw() +
            scale_fill_grey() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = 'white'))
    }

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!missing(plot_title)) {
        plot = plot + ggtitle(plot_title)
    }
    if(!missing(x_label)) {
        plot = plot + xlab(x_label)
    }
    if(!missing(y_label)) {
        plot = plot + ylab(y_label)
    }

    return(plot)
}

#' Plot pair-wise annotations across regions
#'
#' All co-occurring annotations associated with a region are computed and displayed as a heatmap.
#'
#' As with \code{plot_annotation()}, the number in each cell is the number of unique regions annotated to the pair of annotations.
#'
#' For example, if a region is annotated to both a CpG shore and to two different exons simultaneously, the region will only be counted once in the CpG shore / exon cell. NOTE, this same region will count once in both the CpG shore and exon cells on the diagonal.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}.
#' @param annotation_order A character vector which doubles as the subset of annotations desired for plot as well as the ordering. If \code{NULL}, all annotations are displayed.
#' @param plot_title A string used for the title of the plot. If missing, no plot title label is displayed.
#' @param axes_label A string used for the axis labels. If missing, corresponding variable name used.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{ggplot} object which can be viewed by calling it, saved with \code{ggplot2::ggsave}, or edited.
#'
#' @examples
#'    # Select and build annotations
#'    annots = c('hg19_cpgs')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
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
#'    all_order = c(
#'        'hg19_cpg_islands',
#'        'hg19_cpg_shores',
#'        'hg19_cpg_shelves',
#'        'hg19_cpg_inter')
#'
#'    dm_vs_ca = plot_coannotations(
#'        annotated_regions = dm_annots,
#'        annotation_order = all_order,
#'        axes_label = 'Annotations',
#'        plot_title = 'Co-occurrence of Annotations')
#'
#' @export
plot_coannotations = function(annotated_regions, annotation_order = NULL,
    plot_title, axes_label, quiet = FALSE) {

    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_regions = as.data.frame(annotated_regions)

    ########################################################################
    # Order and subset the annotations
    annotated_regions = subset_order_tbl(tbl = annotated_regions, col='annot.type', col_order=annotation_order)

    ########################################################################
    # Find the co-annotations

    annotation_pairs_by_region = dplyr::do(
        dplyr::group_by_(annotated_regions, .dots=c('seqnames', 'start', 'end')),
        expand.grid(.$annot.type, .$annot.type, stringsAsFactors = FALSE))

    annotation_pairs_by_region = dplyr::distinct_(dplyr::ungroup(annotation_pairs_by_region),
        .dots=c('seqnames', 'start', 'end', 'Var1', 'Var2'), .keep_all=TRUE)

    pairwise_annotation_counts = table(annotation_pairs_by_region[['Var1']], annotation_pairs_by_region[['Var2']])

    pac_m = reshape2::melt(pairwise_annotation_counts, value.name = 'Counts')

    ########################################################################
    # Construct the plot

    # Make the base ggplot
    # NOTE: binwidth may need to be a parameter
    plot = ggplot(pac_m, aes_string('Var1', 'Var2')) +
        geom_raster(aes_string(fill = 'Counts')) +
        geom_text(aes_string(fill = 'Counts', label = 'Counts')) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(angle = 30, hjust = 1))

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!missing(plot_title)) {
        plot = plot + ggtitle(plot_title)
    }
    if(!missing(axes_label)) {
        plot = plot + xlab(axes_label)
        plot = plot + ylab(axes_label)
    }

    return(plot)
}

#' Plot numerical data over regions or regions summarized over annotations
#'
#' This function produces either histograms over \code{facet}, or x-y scatterplots over \code{facet}. In the case of histograms over facets, the All distribution (hollow histogram with red outline) is the distribution of \code{x} over all the regions in the data. The facet specific distributions (solid gray) are the distribution of \code{x} over the regions in each facet. For example, a CpG with associated percent methylation annotated to a CpG island and a promoter will count once in the All distribution, but will count once each in the CpG island and promoter facet distributions.
#'
#' @param annotated_regions A \code{GRanges} returned from \code{annotate_regions()}. If the data is not summarized, the data is at the region level. If it is summarized, it represents the average or standard deviation of the regions by the character vector used for \code{by} in \code{summarize_numerical()}.
#' @param x A string indicating the column of the \code{GRanges} to use for the x-axis.
#' @param y A string indicating the column of the \code{GRanges} to use for the y-axis. If missing, a a histogram over \code{x} will be plotted. If not missing, a scatterplot is plotted.
#' @param facet A string indicating which categorical variable in the \code{GRanges} to make \code{ggplot2} facets over. Default is \code{annot.type}.
#' @param facet_order A character vector which give the order of the facets, and can be used to subset the column in the \code{GRanges} used for the \code{facet}. For example, if \code{facet = 'annot.type'}, then the annotations maybe subsetted to just CpG annotations. Default is \code{NULL}, meaning all annotations in their default order are used.
#' @param bin_width An integer indicating the bin width of the histogram used for score. Default 10. Select something appropriate for the data. NOTE: This is only used if \code{y} is \code{NULL}.
#' @param plot_title A string used for the title of the plot. If missing, no title is displayed.
#' @param x_label A string used for the x-axis label. If missing, no x-axis label is displayed.
#' @param y_label A string used for the y-axis label. If missing, no y-axis label is displayed.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#'    # An example with multi-columned data
#'
#'    # Select and build annotations
#'    annots = c('hg19_cpgs')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
#'
#'    dm_file = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
#'    extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
#'    dm_regions = read_regions(con = dm_file, extraCols = extraCols, genome = 'hg19',
#'        rename_score = 'pval', rename_name = 'DM_status', format = 'bed')
#'
#'    # Annotate the regions
#'    dm_annots = annotate_regions(
#'        regions = dm_regions,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#'    # Plot histogram of group 1 methylation rates across the CpG annotations.
#'    # NOTE: Overall distribution (everything in \code{facet_order})
#'    # is plotted in each facet for comparison.
#'    dm_vs_regions_mu1 = plot_numerical(
#'        annotated_regions = dm_annots,
#'        x = 'mu1',
#'        facet = 'annot.type',
#'        facet_order = c('hg19_cpg_islands','hg19_cpg_shores',
#'            'hg19_cpg_shelves','hg19_cpg_inter'),
#'        bin_width = 5,
#'        plot_title = 'Group 1 Methylation over CpG Annotations',
#'        x_label = 'Group 1 Methylation')
#'
#'    # Can also use the result of annotate_regions() to plot two numerical
#'    # data columns against each other for each region, and facet by annotations.
#'    dm_vs_regions_annot = plot_numerical(
#'        annotated_regions = dm_annots,
#'        x = 'mu0',
#'        y = 'mu1',
#'        facet = 'annot.type',
#'        facet_order = c('hg19_cpg_islands','hg19_cpg_shores',
#'            'hg19_cpg_shelves','hg19_cpg_inter'),
#'        plot_title = 'Region Methylation: Group 0 vs Group 1',
#'        x_label = 'Group 0',
#'        y_label = 'Group 1')
#'
#'    # Another example, but using differential methylation status as the facets.
#'    dm_vs_regions_name = plot_numerical(
#'        annotated_regions = dm_annots,
#'        x = 'mu0',
#'        y = 'mu1',
#'        facet = 'DM_status',
#'        facet_order = c('hyper','hypo','none'),
#'        plot_title = 'Region Methylation: Group 0 vs Group 1',
#'        x_label = 'Group 0',
#'        y_label = 'Group 1')
#'
#' @export
plot_numerical = function(annotated_regions, x, y, facet = 'annot.type', facet_order = NULL, bin_width=10,
    plot_title, x_label, y_label, quiet = FALSE) {

    # Tidy the GRanges into a tbl_df for use with dplyr functions
    tbl = as.data.frame(annotated_regions)

    ########################################################################
    # Order and subset the annotations
    sub_tbl = subset_order_tbl(tbl = tbl, col = facet, col_order = facet_order)

    ########################################################################
    # Create data objects for plots
    facet_data = dplyr::distinct_(dplyr::ungroup(sub_tbl), .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)
    all_data = dplyr::distinct_(dplyr::select(dplyr::ungroup(tbl), -matches(facet)), .dots=c('seqnames', 'start', 'end'), .keep_all=TRUE)

    ########################################################################
    # Construct the plot
    # Note, data must be dplyr::ungroup()-ed before hand for the proper
    # display of the overall distribution.

    if(missing(y)) {
        legend_facet = sprintf('%s in %s', x, facet)
        legend_cum = sprintf('All %s', x)
        fill_man = c(NA, 'gray')
        names(fill_man) = c(legend_cum, legend_facet)

        # Make the base histogram ggplot
        plot =
            # Facet hists are plotted with distinct (seqnames, start, end, annot.type) combinations
            ggplot(
                data = facet_data,
                aes_string(x=x, y='..density..')) +
            geom_histogram(binwidth=bin_width, aes(fill = legend_facet)) +
            facet_wrap( stats::as.formula(paste("~", facet)) ) + # Over the facets
            # All hist is plotted with distinct (seqnames, start, end) combinations
            geom_histogram(
                data = all_data,
                binwidth=bin_width, aes(fill = legend_cum, color = 'red')) + # All the data
            theme_bw() +
            scale_fill_manual(values = fill_man) +
            guides(color = FALSE) +
            theme(legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = c('red','white')))
    } else {
        # Make the base scatter ggplot
        plot = ggplot(facet_data, aes_string(x=x, y=y)) +
            geom_point(alpha = 1/8, size = 1) +
            facet_wrap( stats::as.formula(paste("~", facet)) ) +
            theme_bw()
    }

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!missing(plot_title)) {
        plot = plot + ggtitle(plot_title)
    }
    if(!missing(x_label)) {
        plot = plot + xlab(x_label)
    }
    if(!missing(y_label)) {
        plot = plot + ylab(y_label)
    }

    return(plot)
}

#' Plot numerical data occurring in pairs of annotations
#'
#' Plot numerical data associated with regions occurring in \code{annot1}, \code{annot2} and in both. As with \code{plot_numerical()}, the result is a plot of histograms or x-y scatterplots.
#'
#' For example, a CpG with associated percent methylation annotated to a CpG island and a promoter will count once in the All distribution and once in the CpG island / promoter facet distribution. However, a CpG associated only with a promoter will count once in the All distribution and once in the promoter / promoter distribution.
#'
#' @param annotated_regions A \code{GRanges} returned from \code{annotate_regions()}.
#' @param x A string indicating the column of the \code{GRanges} to use for the x-axis.
#' @param y A string indicating the column of the \code{GRanges} to use for the y-axis. If missing, a histogram over \code{x} will be plotted. If not missing, a scatterplot is plotted.
#' @param annot1 A string indicating the first annotation type.
#' @param annot2 A string indicating the second annotation type.
#' @param bin_width An integer indicating the bin width of the histogram used for score. Default 10. Select something appropriate for the data. NOTE: This is only used if \code{y} is \code{NULL}.
#' @param plot_title A string used for the title of the plot. If missing, no title is displayed.
#' @param x_label A string used for the x-axis label. If missing, no x-axis label is displayed.
#' @param y_label A string used for the y-axis label. If missing, no y-axis label is displayed.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#'
#'    # Select and build annotations
#'    annots = c('hg19_cpg_islands','hg19_genes_promoters')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
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
#'    dm_vs_num_co = plot_numerical_coannotations(
#'        annotated_regions = dm_annots,
#'        x = 'mu0',
#'        annot1 = 'hg19_cpg_islands',
#'        annot2 = 'hg19_genes_promoters',
#'        bin_width = 5,
#'        plot_title = 'Group 0 Perc. Meth. in CpG Islands and Promoters',
#'        x_label = 'Percent Methylation')
#'
#' @export
plot_numerical_coannotations = function(annotated_regions, x, y, annot1, annot2, bin_width=10,
    plot_title, x_label, y_label, quiet = FALSE) {

    # Tidy the GRanges into a tbl_df for use with dplyr functions
    tbl = as.data.frame(annotated_regions)

    ########################################################################
    # Order and subset the annotations
    annotation_order = c(annot1,annot2)
    sub_tbl = subset_order_tbl(tbl = tbl, col='annot.type', col_order=annotation_order)

    ########################################################################
    # Find the co-annotations

    # Use combn instead of expand.grid because we do not want regions annotated to
    # a CpG island and a promoter having their data value count in the island / island
    # facet as well as the promoter / promoter facet. We want it *ONLY* in the
    # island / promoter facet. Note, sorting ensures island / promoter and promoter / island
    # are aggregated
    pairs_by_region = dplyr::do(
        dplyr::group_by_(sub_tbl, .dots=c('seqnames', 'start', 'end')),
        if(nrow(.) == 1) {
            as.data.frame(
                t(
                    utils::combn(
                        rep.int(as.character(.$annot.type), 2)
                    , 2))
                , stringsAsFactors = FALSE)
        } else {
            as.data.frame(
                t(
                    utils::combn(
                        sort(as.character(.$annot.type))
                    , 2))
            , stringsAsFactors = FALSE)
        }
    )

    # Join on the data chromosome locations
    pairs_by_region = dplyr::inner_join(x = pairs_by_region, y = sub_tbl, by = c('seqnames','start','end'))

    ########################################################################
    # Create data objects for plots
    facet_data = dplyr::distinct_(dplyr::ungroup(pairs_by_region),
        .dots=c('seqnames', 'start', 'end', 'V1', 'V2'), .keep_all=TRUE)
    all_data = dplyr::distinct_(dplyr::ungroup(tbl), .dots=c('seqnames', 'start', 'end'), .keep_all=TRUE)

    ########################################################################
    # Construct the plot
    # Note, data must be dplyr::ungroup()-ed before hand for the proper
    # display of the overall distribution.

    if(missing(y)) {
        legend_facet = sprintf('%s in %s', x, 'annot pair')
        legend_cum = sprintf('All %s', x)
        fill_man = c(NA, 'gray')
        names(fill_man) = c(legend_cum, legend_facet)

        # Make the base histogram ggplot
        plot =
            # Facet hists are plotted with distinct (seqnames, start, end, annot1, annot2) combinations
            ggplot(
                data = facet_data,
                aes_string(x=x, y='..density..')) +
            geom_histogram(binwidth=bin_width, aes(fill = legend_facet)) +
            facet_wrap( V1 ~ V2 ) + # Over the facets
            # All hist is plotted with distinct (seqnames, start, end) combinations
            geom_histogram(
                data = all_data,
                binwidth=bin_width, aes(fill = legend_cum, color = 'red')) + # All the data
            theme_bw() +
            scale_fill_manual(values = fill_man) +
            guides(color = FALSE) +
            theme(legend.title=element_blank(), legend.position="bottom", legend.key = element_rect(color = c('red','white')))
    } else {
        # Make the base scatter ggplot
        plot = ggplot(pairs_by_region, aes_string(x=x, y=y)) +
            geom_point(alpha = 1/8, size = 1) +
            facet_wrap( V1 ~ V2 ) +
            theme_bw()
    }

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!missing(plot_title)) {
        plot = plot + ggtitle(plot_title)
    }
    if(!missing(x_label)) {
        plot = plot + xlab(x_label)
    }
    if(!missing(y_label)) {
        plot = plot + ylab(y_label)
    }

    return(plot)
}

#' Plot a categorical data variable over another
#'
#' Given a \code{GRanges} of annotated regions from \code{annotate_regions()}, visualize the the distribution of categorical data \code{fill} in categorical data \code{x}. A bar representing the distribution of all \code{fill} in \code{x} will be added according to the contents of \code{fill}. This is the distribution over all values of \code{x}. Additionally, when \code{annotated_random} is not missing, a "Random Regions" bar shows the distribution of random regions over \code{fill}.
#'
#' For example, if a differentially methylated region has the categorical label hyper, and is annotated to a promoter, a 5UTR, two exons, and an intron. Each annotation will appear in the All bar once. Likewise for the hyper bar if the differential methylation status is chosen as \code{x} with \code{annot.type} chosen as \code{fill}.
#'
#' @param annotated_regions The \code{GRanges} result of \code{annotate_regions()}.
#' @param annotated_random The \code{GRanges} result of \code{annotate_regions()} on the randomized regions created from \code{randomize_regions()}. Random regions can only be used with \code{fill == 'annot.type'}.
#' @param x One of 'annot.type' or a categorical data column, indicating whether annotation classes or data classes will appear on the x-axis.
#' @param fill One of 'annot.type', a categorical data column, or \code{NULL}, indicating whether annotation classes or data classes will fill the bars. If \code{NULL} then the bars will be the total counts of the x classes.
#' @param x_order A character vector that subsets and orders the x classes. Default \code{NULL}, uses existing values.
#' @param fill_order A character vector that subsets and orders the fill classes. Default \code{NULL}, uses existing values.
#' @param position A string which has the same possible values as in \code{ggplot2::geom_bar(..., position)}, i.e., 'stack', 'fill', 'dodge', etc.
#' @param plot_title A string used for the title of the plot. If missing, no title is displayed.
#' @param legend_title A string used for the legend title to describe fills (if fill is not \code{NULL}). Default displays corresponding variable name.
#' @param x_label A string used for the x-axis label. If missing, corresponding variable name used.
#' @param y_label A string used for the y-axis label. If missing, corresponding variable name used.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#' @return A \code{ggplot} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#'    # Select and build annotations
#'    annots = c('hg19_cpgs')
#'    annotations = build_annotations(genome = 'hg19', annotations = annots)
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
#'    dm_order = c(
#'        'hyper',
#'        'hypo')
#'    cpg_order = c(
#'        'hg19_cpg_islands',
#'        'hg19_cpg_shores',
#'        'hg19_cpg_shelves',
#'        'hg19_cpg_inter')
#'
#'    dm_vn = plot_categorical(
#'        annotated_regions = dm_annots,
#'        x = 'DM_status',
#'        fill = 'annot.type',
#'        x_order = dm_order,
#'        fill_order = cpg_order,
#'        position = 'fill',
#'        legend_title = 'knownGene Annotations',
#'        x_label = 'DM status',
#'        y_label = 'Proportion')
#'
#'    # Create randomized regions
#'    dm_rnd_regions = randomize_regions(regions = dm_regions)
#'    dm_rnd_annots = annotate_regions(
#'        regions = dm_rnd_regions,
#'        annotations = annotations,
#'        ignore.strand = TRUE)
#'
#'    dm_vn_rnd = plot_categorical(
#'        annotated_regions = dm_annots,
#'        annotated_random = dm_rnd_annots,
#'        x = 'DM_status',
#'        fill = 'annot.type',
#'        x_order = dm_order,
#'        fill_order = cpg_order,
#'        position = 'fill',
#'        legend_title = 'knownGene Annotations',
#'        x_label = 'DM status',
#'        y_label = 'Proportion')
#'
#' @export
plot_categorical = function(annotated_regions, annotated_random, x, fill=NULL, x_order=NULL, fill_order=NULL,
    position = 'stack', plot_title, legend_title, x_label, y_label, quiet = FALSE) {

    ########################################################################
    # Argument parsing and error handling

    # Tidy the GRanges into a tbl_df for use with dplyr functions
    annotated_regions = as.data.frame(annotated_regions)

    # Ensure the value of x is a column name in summarized_cats
    if( !(x %in% colnames(annotated_regions)) ) {
        stop('The column name used for x does not exist in annotated_regions.')
    }

    # Ensure the value of fill is a column name in summarized_cats if it isn't NULL
    # Also ensure fill != x
    if( !is.null(fill) ) {
        if( !(fill %in% colnames(annotated_regions)) ) {
            stop('The column name used for fill does not exist in annotated_regions.')
        }
        if( x == fill ) {
            stop('Error: x cannot equal fill')
        }
    }

    # If !is.null(annotated_random), check that fill = 'annot.type'. This is the
    # only situation where random regions can be used, because the data from the
    # original regions is not transferred to the random ones.
    if(!missing(annotated_random) && fill != 'annot.type') {
        stop('Error: Random regions can only be used in plot_categorical() when fill == "annot.type" since data from the original regions are not transferred to the random regions.')
    }

    # Check valid position argument
    if(position != 'stack' && position != 'fill' && position != 'dodge') {
        stop('Error: position must be one of "stack", "fill", or "dodge"')
    }

    ########################################################################
    # Order and subset based on fill_order
    annotated_regions = subset_order_tbl(tbl = annotated_regions, col = fill, col_order = fill_order)

    # Take the distinct annotation types per unique data region
    annotated_regions = dplyr::distinct_(dplyr::ungroup(annotated_regions), .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)

    ########################################################################
    # Order and subset based on x_order
    if(is.null(x_order)) {
        x_order = unique(annotated_regions[[x]])
    }
    sub_annot_regions = subset_order_tbl(tbl = annotated_regions, col = x, col_order = x_order)

    # Do particular things if annotated_random isn't NULL
    if(!missing(annotated_random)) {
        # Tidy the GRanges into a tbl_df for use with dplyr functions
        annotated_random = as.data.frame(annotated_random)

        # Order and subset the randomized annotations
        annotated_random = subset_order_tbl(tbl = annotated_random, col=fill, col_order=fill_order)

        # Take the distinct annotation types per unique random data region
        annotated_random = dplyr::distinct_(dplyr::ungroup(annotated_random), .dots=c('seqnames', 'start', 'end', 'annot.type'), .keep_all=TRUE)

        # Combine the tbl_dfs in preparation for visualization
        annotated_regions = dplyr::bind_rows("All" = annotated_regions, "Random Regions" = annotated_random, .id = 'data_type')
    }

    ########################################################################
    # Construct the plot

    # Make base ggplot
    if(!missing(annotated_random)) {
        plot =
            ggplot(annotated_regions, aes_string(x='data_type')) +
            geom_bar(aes_string(fill=fill), position=position, width=0.5) + # The All bar
            geom_bar(data = sub_annot_regions, aes_string(x=x, fill=fill), position=position, width=0.5) + # The subsets bars
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
    } else {
        plot =
            ggplot(annotated_regions, aes(x='All')) +
            geom_bar(aes_string(fill=fill), position=position, width=0.5) + # The All bar
            geom_bar(data = sub_annot_regions, aes_string(x=x, fill=fill), position=position, width=0.5) + # The subsets bars
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }

    # Change the fill scale and name if legend_title isn't null
    if(!missing(legend_title)) {
        plot = plot + scale_fill_hue(name=legend_title)
    } else {
        plot = plot + scale_fill_hue()
    }

    # Deal with the x-axis labels to make sure the order is correct
    if(!missing(annotated_random)) {
        plot = plot + scale_x_discrete(limits = c('All', x_order, 'Random Regions'))
    } else {
        if(x == 'annot.type') {
            plot = plot + scale_x_discrete(limits = c('All', names(tidy_annotations(x_order))))
        } else {
            plot = plot + scale_x_discrete(limits = c('All', x_order))
        }
    }

    # Add any user defined labels to the plot if their values are not NULL
    # if they are NULL, ggplot() will use defaults
    if(!missing(plot_title)) {
        plot = plot + ggtitle(plot_title)
    }
    if(!missing(x_label)) {
        plot = plot + xlab(x_label)
    }
    if(!missing(y_label)) {
        plot = plot + ylab(y_label)
    }

    return(plot)
}
