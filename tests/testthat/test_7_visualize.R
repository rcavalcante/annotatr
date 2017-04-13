context('Test plot module')

################################################################################
# Setup annotation objects
data('example_annotations', package = 'annotatr')

################################################################################
# Setup objects for plot_categorical()

dm_file = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
dm_regions = suppressMessages(read_regions(con = dm_file, genome = 'hg19', extraCols = extraCols, rename_score = 'pval', rename_name = 'DM_status', format = 'bed'))
dm_regions = dm_regions[1:1000]

dm_random_regions = suppressMessages(randomize_regions(regions = dm_regions))

dm_annots = suppressMessages(annotate_regions(
    regions = dm_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE))

dm_random_annots = suppressMessages(annotate_regions(
    regions = dm_random_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE))

################################################################################
# Setup order vectors and plots that will work

dm_order = c(
    'hyper',
    'hypo',
    'none')
cpgs_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')

################################################################################
# Test plot_annotation()

test_that('Test plot_annotation() errors', {
    expect_warning(
        plot_annotation(annotated_regions = dm_annots, annotation_order = c('hypor','hype','')),
        'elements in col_order that are not present')
})

test_that('Test plot_annotation() success', {
    dm_va_min = plot_annotation(annotated_regions = dm_annots)

    dm_va = plot_annotation(
        annotated_regions = dm_annots,
        annotation_order = cpgs_order,
        plot_title = 'Testing plot title',
        x_label = 'Test x-label',
        y_label = 'Test y-label')

    dm_va_rnd = plot_annotation(
        annotated_regions = dm_annots,
        annotated_random = dm_random_annots,
        annotation_order = NULL,
        plot_title = 'Testing dodged bars',
        x_label = 'Annotation Type',
        y_label = 'Count')

    expect_equal( dplyr::setequal(class(dm_va_min), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_va), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_va_rnd), c('gg','ggplot')), expected = TRUE)
})

################################################################################
# Test plot_coannotations()

test_that('Test plot_coannotations() success', {

    dm_vs_ca = plot_coannotations(
        annotated_regions = dm_annots,
        annotation_order = cpgs_order,
        axes_label = 'Annotations',
        plot_title = 'Co-occurrence of Annotations')

    expect_equal( dplyr::setequal(class(dm_vs_ca), c('gg','ggplot')), expected = TRUE)
})

################################################################################
# Test plot_numerical()

test_that('Test plot_numerical() success', {

    dm_vs_regions_mu1 = plot_numerical(
        annotated_regions = dm_annots,
        x = 'mu1',
        facet = 'annot.type',
        facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
        bin_width = 5,
        plot_title = 'Group 1 Methylation over CpG Annotations',
        x_label = 'Group 1 Methylation')

    dm_vs_regions_annot = plot_numerical(
        annotated_regions = dm_annots,
        x = 'mu0',
        y = 'mu1',
        facet = 'annot.type',
        facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
        plot_title = 'Region Methylation: Group 0 vs Group 1',
        x_label = 'Group 0',
        y_label = 'Group 1')

    dm_vs_regions_name = plot_numerical(
        annotated_regions = dm_annots,
        x = 'mu0',
        y = 'mu1',
        facet = 'DM_status',
        facet_order = c('hyper','hypo','none'),
        plot_title = 'Region Methylation: Group 0 vs Group 1',
        x_label = 'Group 0',
        y_label = 'Group 1')

    expect_equal( dplyr::setequal(class(dm_vs_regions_mu1), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vs_regions_annot), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vs_regions_name), c('gg','ggplot')), expected = TRUE)
})

################################################################################
# Test plot_numerical_coannotations()

test_that('Test plot_numerical_coannotations()', {
  dm_vs_num_co1 = plot_numerical_coannotations(
    annotated_regions = dm_annots,
    x = 'mu0',
    annot1 = 'hg19_cpg_islands',
    annot2 = 'hg19_cpg_shores',
    bin_width = 5,
    plot_title = 'Group 0 Perc. Meth. in CpG Islands and Promoters',
    x_label = 'Percent Methylation')

  dm_vs_num_co2 = plot_numerical_coannotations(
    annotated_regions = dm_annots,
    x = 'mu0',
    y = 'mu1',
    annot1 = 'hg19_cpg_islands',
    annot2 = 'hg19_cpg_shores',
    bin_width = 5,
    plot_title = 'Group 0 Perc. Meth. in CpG Islands and Promoters',
    x_label = 'Percent Methylation',
    y_label = 'Percent Methylation')

    expect_equal( dplyr::setequal(class(dm_vs_num_co1), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vs_num_co2), c('gg','ggplot')), expected = TRUE)
})

################################################################################
# Test plot_categorical()

  test_that('Test plot_categorical() errors', {
    expect_error(
        plot_categorical(
            annotated_regions = dm_annots),
        'argument "x" is missing')

    expect_error(
        plot_categorical(
            annotated_regions = dm_annots,
            x = 'testing'),
        'column name used for x does not exist in annotated_regions')

    expect_error(
        plot_categorical(
            annotated_regions = dm_annots,
            x = 'DM_status',
            fill = 'testing'),
        'column name used for fill does not exist in annotated_regions')

    expect_error(
        plot_categorical(
            annotated_regions = dm_annots,
            x = 'DM_status',
            fill = 'DM_status'),
        'x cannot equal fill')

    expect_error(
        plot_categorical(
            annotated_regions = dm_annots,
            x = 'DM_status',
            fill = 'annot.type',
            position = 'no'),
        'position must be one of "stack", "fill"')

    expect_warning(
        plot_categorical(
            annotated_regions = dm_annots,
            x = 'DM_status',
            fill = 'annot.type',
            x_order = cpgs_order),
        'elements in col_order that are not present')

    expect_warning(
        plot_categorical(
            annotated_regions = dm_annots,
            x = 'DM_status',
            fill = 'annot.type',
            fill_order = dm_order),
        'elements in col_order that are not present')
  })

test_that('Test plot_categorical() error for random regions and non annot fill', {
    expect_error(
        plot_categorical(
            annotated_regions = dm_annots,
            annotated_random = dm_random_annots,
            x = 'annot.type',
            fill = 'DM_status',
            x_order = cpgs_order,
            fill_order = dm_order,
            position = 'fill',
            legend_title = 'Annotations',
            plot_title = 'DM status by CpG Annotation Proportions',
            x_label = 'DM status',
            y_label = 'Proportion'),
        'since data from the original regions are not transferred to the random regions')
    })

test_that('Test plot_categorical() success', {
    dm_vn_min = plot_categorical(
        annotated_regions = dm_annots,
        x = 'annot.type')

    dm_vn = plot_categorical(
        annotated_regions = dm_annots,
        x = 'DM_status',
        fill = 'annot.type',
        x_order = dm_order,
        fill_order = cpgs_order,
        position = 'fill',
        legend_title = 'knownGene Annotations',
        plot_title = 'DM status in knownGene Annots.',
        x_label = 'DM status',
        y_label = 'Proportion')

    dm_vn_rnd = plot_categorical(
        annotated_regions = dm_annots,
        annotated_random = dm_random_annots,
        x = 'DM_status',
        fill = 'annot.type',
        x_order = dm_order,
        fill_order = cpgs_order,
        position = 'fill',
        legend_title = 'Annotations',
        plot_title = 'DM status by CpG Annotation Proportions',
        x_label = 'DM status',
        y_label = 'Proportion')

    expect_equal( dplyr::setequal(class(dm_vn_min), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vn), c('gg','ggplot')), expected = TRUE)
    expect_equal( dplyr::setequal(class(dm_vn_rnd), c('gg','ggplot')), expected = TRUE)
})
