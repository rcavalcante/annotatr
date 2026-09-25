test_that('randomize_regions() is deprecated', {
    r = read_regions(con = extdata('Gm12878_Stat3_chr2.bed.gz'), genome = 'hg19', format = 'bed')

    expect_warning(
        randomize_regions(regions = r, quiet = TRUE),
        'randomize_regions\\(\\) is deprecated')
})

test_that('randomize_regions() errors for bad input', {
    r_nogenome = read_regions(con = extdata('Gm12878_Stat3_chr2.bed.gz'), format = 'bed')

    expect_error(
        suppressWarnings(randomize_regions(regions = 'hello')),
        'regions must have class GRanges')
    expect_error(
        suppressWarnings(randomize_regions(regions = r_nogenome)),
        'GRanges object must have a valid genome')
})

test_that('randomize_regions() keeps region widths and chromosomes', {
    r = read_regions(con = extdata('Gm12878_Stat3_chr2.bed.gz'), genome = 'hg19', format = 'bed')

    random_r = suppressWarnings(randomize_regions(regions = r, allow.overlaps = TRUE, per.chromosome = TRUE, quiet = TRUE))

    expect_s4_class(random_r, 'GRanges')
    expect_length(random_r, length(r))
    expect_setequal(GenomicRanges::width(random_r), GenomicRanges::width(r))
    expect_equal(unique(as.character(GenomicRanges::seqnames(random_r))), 'chr2')
})
