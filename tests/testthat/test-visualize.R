dm_order = c('hyper', 'hypo', 'none')

################################################################################
# plot_annotation()

test_that('plot_annotation() warns about annotations not in the data', {
    expect_warning(
        plot_annotation(annotated_regions = annotate_dm_regions(), annotation_order = c('hypor', 'hype', '')),
        'elements in col_order that are not present')
})

test_that('plot_annotation() plots counts per annotation', {
    a = annotate_dm_regions()

    expect_plot_builds(plot_annotation(annotated_regions = a))

    p = plot_annotation(
        annotated_regions = a,
        annotation_order = cpg_types,
        plot_title = 'Testing plot title',
        x_label = 'Test x-label',
        y_label = 'Test y-label')
    built = expect_plot_builds(p)

    # One bar per annotation type, with the counts from summarize_annotations()
    s = summarize_annotations(annotated_regions = a, quiet = TRUE)
    expect_equal(sort(built$data[[1]]$count), sort(s$n))
    expect_equal(p$labels$title, 'Testing plot title')
})

test_that('plot_annotation() plots data and background side by side', {
    a = annotate_dm_regions()
    a_dm = a[a$DM_status != 'none']

    p = plot_annotation(
        annotated_regions = a_dm,
        annotated_random = a,
        annotation_order = cpg_types,
        plot_title = 'Testing dodged bars',
        x_label = 'Annotation Type',
        y_label = 'Count')
    built = expect_plot_builds(p)

    # Two bars (data and background) per annotation type
    expect_equal(nrow(built$data[[1]]), 2 * length(cpg_types))
})

################################################################################
# plot_coannotations()

test_that('plot_coannotations() plots pairs of annotations', {
    p = plot_coannotations(
        annotated_regions = annotate_dm_regions(),
        annotation_order = cpg_types,
        axes_label = 'Annotations',
        plot_title = 'Co-occurrence of Annotations')

    expect_plot_builds(p)
})

################################################################################
# plot_numerical()

test_that('plot_numerical() plots histograms and scatterplots over facets', {
    a = annotate_dm_regions()

    expect_plot_builds(plot_numerical(
        annotated_regions = a,
        x = 'mu1',
        facet = 'annot.type',
        facet_order = cpg_types,
        bin_width = 5,
        plot_title = 'Group 1 Methylation over CpG Annotations',
        x_label = 'Group 1 Methylation',
        legend_facet_label = 'Group 1 Methylation Rate in Annotation',
        legend_cum_label = 'Overall Group 1 Methylation Rate'))

    expect_plot_builds(plot_numerical(
        annotated_regions = a,
        x = 'mu0',
        y = 'mu1',
        facet = 'annot.type',
        facet_order = cpg_types,
        plot_title = 'Region Methylation: Group 0 vs Group 1',
        x_label = 'Group 0',
        y_label = 'Group 1'))

    expect_plot_builds(plot_numerical(
        annotated_regions = a,
        x = 'mu0',
        y = 'mu1',
        facet = 'DM_status',
        facet_order = dm_order,
        plot_title = 'Region Methylation: Group 0 vs Group 1',
        x_label = 'Group 0',
        y_label = 'Group 1'))
})

test_that('plot_numerical() facets over two variables', {
    a = annotate_dm_regions()

    built = expect_plot_builds(plot_numerical(
        annotated_regions = a,
        x = 'mu1',
        facet = c('annot.type', 'DM_status'),
        facet_order = list(c('hg19_cpg_islands', 'hg19_cpg_shores'), dm_order),
        plot_title = 'Region Methylation: Group 0 vs Group 1',
        x_label = 'Group 0',
        y_label = 'Group 1'))
    # A 2 x 3 grid of facets
    expect_equal(nrow(built$layout$layout), 2 * 3)

    expect_plot_builds(plot_numerical(
        annotated_regions = a,
        x = 'mu0',
        y = 'mu1',
        facet = c('annot.type', 'DM_status'),
        facet_order = list(NULL, dm_order),
        plot_title = 'Region Methylation: Group 0 vs Group 1',
        x_label = 'Group 0',
        y_label = 'Group 1'))
})

################################################################################
# plot_numerical_coannotations()

test_that('plot_numerical_coannotations() plots histograms and scatterplots', {
    a = annotate_dm_regions()

    expect_plot_builds(plot_numerical_coannotations(
        annotated_regions = a,
        x = 'mu0',
        annot1 = 'hg19_cpg_islands',
        annot2 = 'hg19_cpg_shores',
        bin_width = 5,
        plot_title = 'Group 0 Perc. Meth. in CpG Islands and Shores',
        x_label = 'Percent Methylation',
        legend_facet_label = 'Perc. Methylation in annotation pair',
        legend_cum_label = 'Overall Perc. Methylation'))

    expect_plot_builds(plot_numerical_coannotations(
        annotated_regions = a,
        x = 'mu0',
        y = 'mu1',
        annot1 = 'hg19_cpg_islands',
        annot2 = 'hg19_cpg_shores',
        bin_width = 5,
        plot_title = 'Group 0 Perc. Meth. in CpG Islands and Shores',
        x_label = 'Percent Methylation',
        y_label = 'Percent Methylation'))
})

################################################################################
# plot_categorical()

test_that('plot_categorical() errors for bad arguments', {
    a = annotate_dm_regions()

    expect_error(plot_categorical(annotated_regions = a), 'argument "x" is missing')
    expect_error(
        plot_categorical(annotated_regions = a, x = 'testing'),
        'column name used for x does not exist in annotated_regions')
    expect_error(
        plot_categorical(annotated_regions = a, x = 'DM_status', fill = 'testing'),
        'column name used for fill does not exist in annotated_regions')
    expect_error(
        plot_categorical(annotated_regions = a, x = 'DM_status', fill = 'DM_status'),
        'x cannot equal fill')
    expect_error(
        plot_categorical(annotated_regions = a, x = 'DM_status', fill = 'annot.type', position = 'no'),
        'position must be one of "stack", "fill"')
})

test_that('plot_categorical() warns about orders not in the data', {
    a = annotate_dm_regions()

    expect_warning(
        plot_categorical(annotated_regions = a, x = 'DM_status', fill = 'annot.type', x_order = cpg_types),
        'elements in col_order that are not present')
    expect_warning(
        plot_categorical(annotated_regions = a, x = 'DM_status', fill = 'annot.type', fill_order = dm_order),
        'elements in col_order that are not present')
})

test_that('plot_categorical() errors for a background with a data fill', {
    a = annotate_dm_regions()
    a_dm = a[a$DM_status != 'none']

    expect_error(
        plot_categorical(
            annotated_regions = a_dm,
            annotated_random = a,
            x = 'annot.type',
            fill = 'DM_status',
            x_order = cpg_types,
            fill_order = c('hyper', 'hypo')),
        'since the background need not have the data columns')
})

test_that('plot_categorical() plots categories with an All bar', {
    a = annotate_dm_regions()

    expect_plot_builds(plot_categorical(annotated_regions = a, x = 'annot.type'))

    p = plot_categorical(
        annotated_regions = a,
        x = 'DM_status',
        fill = 'annot.type',
        x_order = dm_order,
        fill_order = cpg_types,
        position = 'fill',
        legend_title = 'knownGene Annotations',
        plot_title = 'DM status in knownGene Annots.',
        x_label = 'DM status',
        y_label = 'Proportion')
    expect_plot_builds(p)
    expect_equal(p$labels$title, 'DM status in knownGene Annots.')
})

test_that('plot_categorical() adds a Background bar', {
    a = annotate_dm_regions()
    a_dm = a[a$DM_status != 'none']

    built = expect_plot_builds(plot_categorical(
        annotated_regions = a_dm,
        annotated_random = a,
        x = 'DM_status',
        fill = 'annot.type',
        x_order = c('hyper', 'hypo'),
        fill_order = cpg_types,
        position = 'fill',
        legend_title = 'Annotations',
        plot_title = 'DM status by CpG Annotation Proportions',
        x_label = 'DM status',
        y_label = 'Proportion'))

    expect_equal(built$layout$panel_params[[1]]$x$get_labels(), c('All', 'hyper', 'hypo', 'Background'))
})

test_that('plot_categorical() keeps duplicate regions with different categories', {
    r = dm_regions()
    r$cancer_status = 'Cancer'
    r2 = r
    r2$cancer_status = 'NoCancer'
    a = annotate_dm_regions(c(r, r2))

    expect_plot_builds(plot_categorical(
        annotated_regions = a,
        x = 'cancer_status',
        fill = 'annot.type',
        x_order = c('Cancer', 'NoCancer'),
        fill_order = cpg_types,
        position = 'fill'))
})
