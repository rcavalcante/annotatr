################################################################################
# Errors

test_that('annotate_regions() errors for non-GRanges input', {
    bed = extdata('test_intersect.bed')
    r = read_regions(con = bed, format = 'bed')

    expect_error(
        annotate_regions(regions = bed, annotations = annotatr::annotations, quiet = TRUE),
        'regions object is not GRanges')
    expect_error(
        annotate_regions(regions = r, annotations = bed, quiet = TRUE),
        'annotations object is not GRanges')
})

test_that('annotate_regions() errors when nothing overlaps', {
    r = read_regions(con = extdata('test_intersect.bed'), format = 'bed')

    read_annotations(con = extdata('test_annotation_nooverlap.bed'), name = 'annotnooverlap')
    annotations = build_annotations(genome = 'hg19', annotations = 'genome_custom_annotnooverlap')

    expect_error(
        annotate_regions(regions = r, annotations = annotations, quiet = TRUE),
        'No annotations intersect the regions')
})

################################################################################
# annotate_regions()

test_that('annotate_regions() annotates to the premade CpG annotations', {
    r = read_regions(con = extdata('Gm12878_Stat3_chr2.bed.gz'), format = 'bed')

    a = annotate_regions(regions = r, annotations = annotatr::annotations, ignore.strand = TRUE, quiet = TRUE)

    expect_s4_class(a, 'GRanges')
    expect_true(all(a$annot$type %in% cpg_types))
    # Every annotated region is one of the input regions
    expect_true(all(IRanges::overlapsAny(a, r)))
})

# The four regions in test_read_multiple_data_nohead.bed overlap the three
# regions in test_annotations_3.bed 5 times: regions 1 and 2 overlap
# annotation 1, region 3 overlaps annotations 2 and 3, and region 4 overlaps
# annotation 3.
test_that('annotate_regions() works with only custom annotations', {
    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
    r = read_regions(con = extdata('test_read_multiple_data_nohead.bed'), extraCols = extraCols, rename_score = 'coverage')

    read_annotations(con = extdata('test_annotations_3.bed'), name = 'annotcustom')
    annotations = build_annotations(genome = 'hg19', annotations = 'genome_custom_annotcustom')

    a = annotate_regions(regions = r, annotations = annotations, ignore.strand = TRUE, quiet = TRUE)

    expect_length(a, 5)
    # The data columns of the regions are kept
    expect_contains(names(GenomicRanges::mcols(a)), c('pval', 'mu1', 'mu0', 'diff_exp', 'coverage', 'annot'))
})

test_that('annotate_regions() works with custom and builtin annotations', {
    extraCols = c(pval = 'numeric', mu1 = 'integer', mu0 = 'integer', diff_exp = 'character')
    r = read_regions(con = extdata('test_read_multiple_data_nohead.bed'), extraCols = extraCols, rename_score = 'coverage')

    read_annotations(con = extdata('test_annotations_3.bed'), name = 'annotboth', genome = 'hg19')
    annotations = c(annotatr_cache$get('hg19_custom_annotboth'), annotatr::annotations)

    a = annotate_regions(regions = r, annotations = annotations, ignore.strand = TRUE, quiet = TRUE)

    # The 5 custom overlaps above, and 5 CpG annotation overlaps
    expect_length(a, 10)
    expect_setequal(unique(a$annot$type), c('hg19_custom_annotboth', 'hg19_cpg_islands', 'hg19_cpg_shores', 'hg19_cpg_inter'))
})

# The regions in test_BED3.bed overlap the annotations in
# test_annotations_minoverlap.bed by 5, 1, and 200 bases.
test_that('annotate_regions() uses minoverlap', {
    r = read_regions(con = extdata('test_BED3.bed'), format = 'bed')

    read_annotations(con = extdata('test_annotations_minoverlap.bed'), name = 'annotminoverlap')
    annotations = build_annotations(genome = 'hg19', annotations = 'genome_custom_annotminoverlap')

    a = annotate_regions(regions = r, annotations = annotations, minoverlap = 5, quiet = TRUE)
    expect_equal(GenomicRanges::start(a), c(10791, 28801))

    a = annotate_regions(regions = r, annotations = annotations, minoverlap = 1, quiet = TRUE)
    expect_length(a, 3)
})
