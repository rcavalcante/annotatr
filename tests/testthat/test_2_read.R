context('Test read module')

################################################################################
# Test warnings in read_regions()

test_that('Test rename_* warnings', {
    file = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')

    expect_warning(
        read_regions(con = file, format = 'bed', rename_name = 'hello'),
        'Ignoring rename_name parameter because')

    expect_warning(
        read_regions(con = file, format = 'bed', rename_score = 'score'),
        'Ignoring rename_score parameter because')
})

################################################################################
# Test BED3-6+ and bedGraph

test_that('Test BED3', {
    file = system.file('extdata', 'test_BED3.bed', package = 'annotatr')
    r = read_regions(con = file, format = 'bed')

    expect_true(is(r, 'GRanges'))
})

test_that('Test BED4', {
    file = system.file('extdata', 'test_BED4.bed', package = 'annotatr')
    r = read_regions(con = file, format = 'bed')

    expect_true(is(r, 'GRanges'))
})

test_that('Test BED5', {
    file = system.file('extdata', 'test_BED5.bed', package = 'annotatr')
    r = read_regions(con = file, format = 'bed')

    expect_true(is(r, 'GRanges'))
})

test_that('Test BED6', {
    file = system.file('extdata', 'test_BED6.bed', package = 'annotatr')
    r = read_regions(con = file, format = 'bed')

    expect_true(is(r, 'GRanges'))
})

test_that('Test BED6+ with renaming', {
    file = system.file('extdata', 'IDH2mut_v_NBM_multi_data_chr9.txt.gz', package = 'annotatr')
    extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
    r = read_regions(con = file, extraCols = extraCols, rename_score = 'pval', rename_name = 'DM_status', format = 'bed')

    expect_true(is(r, 'GRanges'))
})

test_that('Test bedGraph', {
    file = system.file('extdata', 'test_bedGraph.bedGraph', package = 'annotatr')
    r = read_regions(con = file, format = 'bedGraph')

    expect_true(is(r, 'GRanges'))
})

################################################################################
# Test

test_that('Test custom BED3 with no genome and a name', {
    file = system.file('extdata', 'test_annotations_3.bed', package = 'annotatr')
    read_annotations(con = file, name = 'test', format = 'bed')

    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_test'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED3 with no name and a genome', {
    file = system.file('extdata', 'test_annotations_3.bed', package = 'annotatr')
    read_annotations(con = file, genome = 'hg19', format = 'bed')

    expect_true( all(colnames(mcols(annotatr_cache$get('hg19_custom_annotations'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED3 with no name or genome', {
    file = system.file('extdata', 'test_annotations_3.bed', package = 'annotatr')
    read_annotations(con = file, format = 'bed')
    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_annotations'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED4', {
    file = system.file('extdata', 'test_annotations_4.bed', package = 'annotatr')
    read_annotations(con = file, format = 'bed')

    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_test'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED5', {
    file = system.file('extdata', 'test_annotations_5.bed', package = 'annotatr')
    read_annotations(con = file, format = 'bed')

    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_test'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED6', {
    file = system.file('extdata', 'test_annotations_6.bed', package = 'annotatr')
    read_annotations(con = file, name = 'six', format = 'bed')

    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_test'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED6 with gene_id', {
    file = system.file('extdata', 'test_annotations_6_gene.bed', package = 'annotatr')
    extraCols = c(gene_id = 'character')
    read_annotations(con = file, name = 'geneid', format = 'bed', extraCols = extraCols)

    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_test'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED6 with symbol', {
    file = system.file('extdata', 'test_annotations_6_symbol.bed', package = 'annotatr')
    extraCols = c(symbol = 'character')
    read_annotations(con = file, name = 'symbol', format = 'bed', extraCols = extraCols)

    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_test'))) == c('id','gene_id','symbol','type')) )
})

test_that('Test custom BED6 with symbol nad gene_id', {
    file = system.file('extdata', 'test_annotations_6_gene_symbol.bed', package = 'annotatr')
    extraCols = c(gene_id = 'character', symbol = 'character')
    read_annotations(con = file, name = 'genesymbol', format = 'bed', extraCols = extraCols)

    expect_true( all(colnames(mcols(annotatr_cache$get('genome_custom_test'))) == c('id','gene_id','symbol','type')) )
})
