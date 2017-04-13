context('Test intersect/annotate module')

################################################################################
# Test errors

test_that('Test error thrown for non-GRanges regions object in annotate_regions()',{
    annotations = c('hg19_cpg_islands')

    bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
    r = suppressMessages(read_regions(con = bed, format = 'bed'))

    expect_error(
        annotate_regions(
            regions = bed,
            annotations = annotations,
            ignore.strand = TRUE,
            quiet = TRUE),
        "regions object is not GRanges")

    expect_error(
        annotate_regions(
            regions = r,
            annotations = bed,
            ignore.strand = TRUE,
            quiet = TRUE),
        "annotations object is not GRanges")

    a_file = system.file('extdata', 'test_annotation_nooverlap.bed', package = 'annotatr')
    read_annotations(con = a_file, name = 'test')
    annotations = build_annotations(genome = 'hg19', annotations = 'genome_custom_test')
    expect_error(
        annotate_regions(
            regions = r,
            annotations = annotations,
            ignore.strand = TRUE,
            quiet = TRUE),
        "No annotations intersect the regions")
})

################################################################################
# Test annotate_regions()

test_that('Test a la carte annotations in annotate_regions()',{
    # Get premade CpG annotations
    annots = expand_annotations('hg19_cpgs')
    data('example_annotations', package = 'annotatr')

    bed = system.file('extdata', 'test_intersect.bed', package = 'annotatr')
    r = read_regions(con = bed, format = 'bed')

    i = annotate_regions(
        regions = r,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = TRUE)

    expect_true( dplyr::setequal(unique(i$annot$type), annots) )
})

test_that('Test a la carte and shortcut annotations in annotate_regions()',{
    data('example_annotations', package = 'annotatr')

    file = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')
    r = read_regions(con = file, format = 'bed')

    i = annotate_regions(
        regions = r,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = TRUE)

    expect_true( dplyr::setequal(unique(i$annot$type), c('hg19_cpg_islands', 'hg19_cpg_shores', 'hg19_cpg_shelves', 'hg19_cpg_inter')) )
})

test_that('Custom annotations work in annotate_regions()', {
    r_file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
    r = read_regions(con = r_file, extraCols = extraCols, rename_score = 'coverage')

    a_file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
    read_annotations(con = a_file, name = 'TFBS', genome = 'hg19')

    annots = c('hg19_custom_TFBS', 'hg19_cpgs')
    annotations = build_annotations(genome = 'hg19', annotations = annots)

    a = annotate_regions(
        regions = r,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = TRUE)

    expect_equal(length(a) == 10, expected = TRUE)
})

test_that('annotate_regions() works with only custom annotations', {
    r_file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
    r = read_regions(con = r_file, extraCols = extraCols, rename_score = 'coverage')

    a_file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
    read_annotations(con = a_file, name = 'TFBS')
    annotations = build_annotations(genome = 'hg19', annotations = 'genome_custom_TFBS')

    a = annotate_regions(
        regions = r,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)

    expect_equal(length(a) == 5, expected = TRUE)
})

test_that('annotate_regions() uses minoverlap correctly', {
    file = system.file('extdata', 'test_BED3.bed', package = 'annotatr')
    r = read_regions(con = file, format = 'bed')

    a_file = system.file('extdata', 'test_annotations_minoverlap.bed', package='annotatr')
    read_annotations(con = a_file, name = 'TFBS')
    annotations = build_annotations(genome = 'hg19', annotations = 'genome_custom_TFBS')

    a = annotate_regions(
        regions = r,
        annotations = annotations,
        minoverlap = 5)

    expect_equal(length(a) == 2, expected = TRUE)
    expect_true(all(GenomicRanges::start(a) == c(10791,28801)))
})
