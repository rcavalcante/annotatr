context('Test summarize module')

annots = c('hg19_cpgs')
annotations = suppressMessages(build_annotations(genome = 'hg19', annotations = annots))

bed = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
r = suppressMessages(read_regions(con = bed, genome = 'hg19', extraCols = extraCols, rename_score = 'pval', rename_name = 'DM_status', format = 'bed'))

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
    expect_equal( sum(s[['n']]), expected = 18656)
    expect_equal( nrow(srand), expected = 8)
})

test_that('Test summarize_numerical()', {
    s = summarize_numerical(
        annotated_regions = a,
        by = c('annot.type', 'annot.id'),
        over = 'diff_meth',
        quiet = TRUE)

    expect_equal( mean(s[['mean']]), expected = 2.745996, tolerance = 0.01)
})

test_that('Test summarize_numerical() and summarize_categorical() over small data', {
    annots = c('hg19_cpg_islands', 'hg19_genes_promoters')
    annotations = build_annotations(genome = 'hg19', annotations = annots)

    r_file = system.file('extdata', 'test_read_multiple_data_nohead.bed', package='annotatr')
    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
    r = read_regions(con = r_file, extraCols = extraCols, rename_name = 'group', rename_score = 'coverage')

    a = annotate_regions(
        regions = r,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = TRUE)

    # Testing summarize_numerical()
    sn1 = summarize_numerical(
        annotated_regions = a,
        by = c('annot.type', 'annot.id'),
        over = 'coverage',
        quiet = FALSE)
    sn2 = summarize_numerical(
        annotated_regions = a,
        by = c('group'),
        over = c('coverage', 'mu1', 'mu0'),
        quiet = TRUE)

    # Testing summarize_categorical()
    sc1 = summarize_categorical(
        annotated_regions = a,
        by = c('annot.type', 'group'),
        quiet = FALSE)
    sc2 = summarize_categorical(
        annotated_regions = a,
        by = c('annot.type', 'diff_exp'),
        quiet = TRUE)
    sc3 = summarize_categorical(
        annotated_regions = a,
        by = c('group','diff_exp'),
        quiet = TRUE)

    expect_equal( sn1[['mean']][which(sn1[['annot.id']] == 'promoter:1')], expected = 66)
    expect_equal( sn1[['mean']][which(sn1[['annot.id']] == 'island:1')], expected = 48)
    expect_equal( sn2[['mu0_mean']][which(sn2[['group']] == 'A')], expected = 30.14, tolerance = 0.01)
    expect_equal( sn2[['mu1_mean']][which(sn2[['group']] == 'B')], expected = 95)

    expect_equal( sc1[['n']][which(sc1[['annot.type']] == 'hg19_genes_promoters' & sc1[,'group'] == 'A')], expected = 2)
    expect_equal( sc2[['n']][which(sc2[,'annot.type'] == 'hg19_cpg_islands' & sc2[,'diff_exp'] == 'Y')], expected = 2)
    expect_equal( sc3[['n']][which(sc3[,'group'] == 'A' & sc3[,'diff_exp'] == 'N')], expected = 1)
})
