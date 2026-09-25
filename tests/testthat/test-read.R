################################################################################
# read_regions()

test_that('read_regions() warns when rename_* columns are absent', {
    file = extdata('Gm12878_Stat3_chr2.bed.gz')

    expect_warning(
        read_regions(con = file, format = 'bed', rename_name = 'hello'),
        'Ignoring rename_name parameter because')
    expect_warning(
        read_regions(con = file, format = 'bed', rename_score = 'score'),
        'Ignoring rename_score parameter because')
})

test_that('read_regions() reads BED3, BED4, BED5, BED6, and bedGraph', {
    # rtracklayer converts the 0-based BED starts to 1-based
    for(file in c('test_BED3.bed', 'test_BED4.bed', 'test_BED5.bed', 'test_BED6.bed')) {
        r = read_regions(con = extdata(file), format = 'bed')

        expect_s4_class(r, 'GRanges')
        expect_length(r, 3)
        expect_equal(GenomicRanges::end(r), c(if(file == 'test_BED3.bed') 10805 else 11000, 28000, 29000))
    }

    r = read_regions(con = extdata('test_BED6.bed'), format = 'bed')
    expect_equal(as.character(GenomicRanges::strand(r)), c('+', '-', '-'))

    r = read_regions(con = extdata('test_bedGraph.bedGraph'), format = 'bedGraph')
    expect_s4_class(r, 'GRanges')
    expect_equal(r$score, c(31, 36, 83))
})

test_that('read_regions() reads BED6+ and renames name and score', {
    extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
    r = read_regions(con = extdata('IDH2mut_v_NBM_multi_data_chr9.txt.gz'), extraCols = extraCols,
        rename_score = 'pval', rename_name = 'DM_status', format = 'bed')

    expect_s4_class(r, 'GRanges')
    expect_named(GenomicRanges::mcols(r), c('DM_status', 'pval', 'diff_meth', 'mu1', 'mu0'))
    expect_setequal(unique(r$DM_status), c('hyper', 'hypo', 'none'))
})

################################################################################
# read_annotations()

annotation_cols = c('id', 'tx_id', 'gene_id', 'symbol', 'entrez_id', 'ensembl_id', 'type')

test_that('read_annotations() names the cache entry from genome and name', {
    file = extdata('test_annotations_3.bed')

    read_annotations(con = file, name = 'readname', format = 'bed')
    a = annotatr_cache$get('genome_custom_readname')
    expect_named(GenomicRanges::mcols(a), annotation_cols)
    expect_equal(a$type, rep('genome_custom_readname', 3))
    expect_equal(a$id, paste0('readname:', 1:3))

    # Without a name, the name is 'annotations'
    read_annotations(con = file, genome = 'hg19', format = 'bed')
    a = annotatr_cache$get('hg19_custom_annotations')
    expect_named(GenomicRanges::mcols(a), annotation_cols)
    expect_equal(a$type, rep('hg19_custom_annotations', 3))
})

test_that('read_annotations() reads BED3 through BED6', {
    for(n in 3:6) {
        name = sprintf('readbed%s', n)
        read_annotations(con = extdata(sprintf('test_annotations_%s.bed', n)), name = name, format = 'bed')
        a = annotatr_cache$get(sprintf('genome_custom_%s', name))

        expect_named(GenomicRanges::mcols(a), annotation_cols)
        expect_length(a, 3)
        # Without extraCols, the gene columns are empty
        expect_true(all(is.na(a$tx_id) & is.na(a$gene_id) & is.na(a$symbol)))
    }
})

test_that('read_annotations() keeps gene_id, symbol, and tx_id extraCols', {
    read_annotations(con = extdata('test_annotations_6_gene.bed'), name = 'readgene',
        format = 'bed', extraCols = c(gene_id = 'character'))
    a = annotatr_cache$get('genome_custom_readgene')
    expect_named(GenomicRanges::mcols(a), annotation_cols)
    expect_equal(a$gene_id, c('324', '4624', '3447'))

    read_annotations(con = extdata('test_annotations_6_symbol.bed'), name = 'readsymbol',
        format = 'bed', extraCols = c(symbol = 'character'))
    a = annotatr_cache$get('genome_custom_readsymbol')
    expect_equal(a$symbol, c('BRCA', 'TP53', 'HOX1A'))

    read_annotations(con = extdata('test_annotations_6_tx_gene_symbol.bed'), name = 'readtxgenesymbol',
        format = 'bed', extraCols = c(gene_id = 'character', symbol = 'character', tx_id = 'character'))
    a = annotatr_cache$get('genome_custom_readtxgenesymbol')
    expect_named(GenomicRanges::mcols(a), annotation_cols)
    expect_equal(a$gene_id, c('351236', '4624', '3447'))
    expect_equal(a$symbol, c('BRCA', 'TP53', 'HOX1A'))
    expect_equal(a$tx_id, c('ENST00000473358.1', 'ENST00000607096.1', 'ENST00000496488.1'))
})
