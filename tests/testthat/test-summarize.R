# The expected values below were computed from the first 1000 regions of
# IDH2mut_v_NBM_multi_data_chr9.txt.gz annotated to the premade hg19 CpG
# annotations. They guard against unintended changes in the summaries.

test_that('summarize_annotations() counts regions per annotation type', {
    a = annotate_dm_regions()

    s = summarize_annotations(annotated_regions = a, quiet = TRUE)

    expect_setequal(s$annot.type, cpg_types)
    # Regions spanning two CpG annotation types count once for each
    expect_equal(sum(s$n), 1064)
})

test_that('summarize_annotations() counts a region once per annotation type', {
    a = annotate_dm_regions()

    s = summarize_annotations(annotated_regions = a, quiet = TRUE)

    per_type = tapply(
        paste(as.character(GenomicRanges::seqnames(a)), GenomicRanges::start(a)),
        a$annot$type,
        function(x) length(unique(x)))
    expect_equal(s$n[match(names(per_type), s$annot.type)], as.vector(per_type))
})

test_that('summarize_annotations() compares data to a background', {
    a = annotate_dm_regions()
    a_dm = a[a$DM_status != 'none']

    s = summarize_annotations(annotated_regions = a, quiet = TRUE)
    s_bg = summarize_annotations(annotated_regions = a_dm, annotated_random = a, quiet = TRUE)

    expect_setequal(unique(s_bg$data_type), c('Data', 'Background'))
    # The background counts are the counts of all tested regions
    bg = s_bg[s_bg$data_type == 'Background', ]
    expect_equal(bg$n[match(s$annot.type, bg$annot.type)], s$n)
})

test_that('summarize_numerical() requires over', {
    expect_error(summarize_numerical(annotated_regions = annotate_dm_regions()), 'over cannot be missing')
})

test_that('summarize_numerical() summarizes over annotations and data columns', {
    a = annotate_dm_regions()

    s1 = summarize_numerical(annotated_regions = a, by = c('annot.type', 'annot.id'), over = 'diff_meth', quiet = TRUE)
    expect_equal(mean(s1$mean), 2.424537, tolerance = 0.01)
    expect_equal(s1$mean[s1$annot.id == 'inter:8599'], -1.0066888, tolerance = 0.01)

    s2 = summarize_numerical(annotated_regions = a, by = 'DM_status', over = c('diff_meth', 'mu1', 'mu0'), quiet = TRUE)
    expect_contains(names(s2), c('diff_meth_mean', 'mu1_mean', 'mu0_mean'))
    expect_equal(s2$mu0_mean[s2$DM_status == 'hyper'], 16.34614, tolerance = 0.01)
})

test_that('summarize_categorical() counts categories per annotation type', {
    a = annotate_dm_regions()

    s = summarize_categorical(annotated_regions = a, by = c('annot.type', 'DM_status'), quiet = TRUE)
    expect_equal(s$n[s$annot.type == 'hg19_cpg_inter' & s$DM_status == 'hyper'], 19)
})

test_that('summarize_categorical() keeps duplicate regions with different categories', {
    r = dm_regions()
    r$cancer_status = 'Cancer'
    r2 = r
    r2$cancer_status = 'NoCancer'
    a = annotate_dm_regions(c(r, r2))

    s = summarize_categorical(annotated_regions = a, by = c('annot.type', 'cancer_status'), quiet = TRUE)

    cancer = s[s$cancer_status == 'Cancer', ]
    no_cancer = s[s$cancer_status == 'NoCancer', ]
    expect_equal(no_cancer$n[match(cancer$annot.type, no_cancer$annot.type)], cancer$n)
})
