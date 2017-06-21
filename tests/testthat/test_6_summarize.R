context('Test summarize module')

data('annotations', package = 'annotatr')

bed = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
r = suppressMessages(read_regions(con = bed, genome = 'hg19', extraCols = extraCols, rename_score = 'pval', rename_name = 'DM_status', format = 'bed'))
r = r[1:1000]

a = suppressMessages(annotate_regions(
    regions = r,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE))

rnd = suppressMessages(randomize_regions(regions = r))

rnd_annot = suppressMessages(annotate_regions(
    regions = rnd,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE))

################################################################################
# Test errors

test_that('Test for error with over=NULL in summarize_numerical()',{
    expect_error(summarize_numerical(annotated_regions = a),
        'over cannot be missing')
})

################################################################################
# Test summarize functions

test_that('Test summarize_annotations()', {
    s = summarize_annotations(annotated_regions = a, quiet = FALSE)

    srand = summarize_annotations(
        annotated_regions = a,
        annotated_random = rnd_annot,
        quiet = FALSE)

    # NOTE: For small data it is possible that the random regions won't
    # intersect all CpG types so the second test may fail. Moreover,
    # if you are going to compute fold changes, corresponding random
    # rows may be missing if the data is too small...
    expect_equal( sum(s[['n']]), expected = 1064)
    expect_equal( nrow(srand), expected = 8)
})

test_that('Test summarize_numerical()', {
    s = summarize_numerical(
        annotated_regions = a,
        by = c('annot.type', 'annot.id'),
        over = 'diff_meth',
        quiet = TRUE)

    expect_equal( mean(s[['diff_meth_mean']]), expected = 2.424537, tolerance = 0.01)
})

test_that('Test summarize_numerical() and summarize_categorical() over small data', {
    # Testing summarize_numerical()
    sn1 = summarize_numerical(
        annotated_regions = a,
        by = c('annot.type', 'annot.id'),
        over = 'diff_meth',
        quiet = FALSE)
    sn2 = summarize_numerical(
        annotated_regions = a,
        by = c('DM_status'),
        over = c('diff_meth', 'mu1', 'mu0'),
        quiet = TRUE)

    # Testing summarize_categorical()
    sc1 = summarize_categorical(
        annotated_regions = a,
        by = c('annot.type', 'DM_status'),
        quiet = FALSE)

    expect_equal( sn1[['diff_meth_mean']][which(sn1[['annot.id']] == 'inter:8599')], expected = -1.0066888, tolerance = 0.01)
    expect_equal( sn2[['mu0_mean']][which(sn2[['DM_status']] == 'hyper')], expected = 16.34614, tolerance = 0.01)

    expect_equal( sc1[['n']][which(sc1[['annot.type']] == 'hg19_cpg_inter' & sc1[,'DM_status'] == 'hyper')], expected = 19)
})
