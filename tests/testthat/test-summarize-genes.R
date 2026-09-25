# A small example whose summaries can be worked out by hand:
#
#   region  position  diff_meth  DM_status  annotated to
#   r1      100-199       10     hyper      geneA promoter, geneA exon
#   r2      300-399      -20     hypo       geneA promoter, geneB promoter
#   r3      500-599       30     hyper      geneA exon, geneB promoter
#   r4      900-999        5     none       a CpG island only
gene_example = function() {
    regions = GenomicRanges::GRanges(
        seqnames = 'chr1',
        ranges = IRanges::IRanges(start = c(100, 300, 500, 900), end = c(199, 399, 599, 999)),
        diff_meth = c(10, -20, 30, 5),
        DM_status = c('hyper', 'hypo', 'hyper', 'none'))

    annotations = GenomicRanges::GRanges(
        seqnames = 'chr1',
        ranges = IRanges::IRanges(start = c(50, 150, 520, 380, 900), end = c(350, 160, 540, 600, 950)),
        id = c('promoter:1', 'exon:1', 'exon:2', 'promoter:2', 'island:1'),
        tx_id = c('tx1', 'tx1', 'tx1', 'tx2', NA),
        gene_id = c('1', '1', '1', '2', NA),
        symbol = c('geneA', 'geneA', 'geneA', 'geneB', NA),
        type = c('hg19_genes_promoters', 'hg19_genes_exons', 'hg19_genes_exons', 'hg19_genes_promoters', 'hg19_cpg_islands'))

    annotate_regions(regions = regions, annotations = annotations, quiet = TRUE)
}

test_that('summarize_genes() gives one row per gene', {
    g = summarize_genes(gene_example(), over = 'diff_meth', by = 'DM_status', quiet = TRUE)

    expect_s3_class(g, 'tbl_df')
    expect_equal(g$gene_id, c('1', '2'))
    expect_equal(g$symbol, c('geneA', 'geneB'))
    expect_named(g, c('gene_id', 'symbol', 'n_regions', 'n_promoters', 'n_exons',
        'n_hyper', 'n_hypo', 'diff_meth_mean', 'diff_meth_median', 'diff_meth_sd'))

    # geneA has r1, r2, and r3. r1 overlaps two of its annotations but counts once.
    expect_equal(g$n_regions, c(3L, 2L))
    expect_equal(g$n_promoters, c(2L, 2L))
    expect_equal(g$n_exons, c(2L, 0L))
    expect_equal(g$n_hyper, c(2L, 1L))
    expect_equal(g$n_hypo, c(1L, 1L))
    expect_equal(g$diff_meth_mean, c(mean(c(10, -20, 30)), mean(c(-20, 30))))
    expect_equal(g$diff_meth_median, c(10, 5))
    expect_equal(g$diff_meth_sd, c(sd(c(10, -20, 30)), sd(c(-20, 30))))
})

test_that('summarize_genes() leaves out non-gene annotations', {
    g = summarize_genes(gene_example(), by = 'DM_status', quiet = TRUE)

    # r4 is only in a CpG island, so there is no 'none' category
    expect_false('n_none' %in% names(g))
    expect_false(any(grepl('cpg', names(g))))
})

test_that('summarize_genes() gives one row per gene and annotation type in long format', {
    l = summarize_genes(gene_example(), over = 'diff_meth', by = 'DM_status', format = 'long', quiet = TRUE)

    expect_named(l, c('gene_id', 'symbol', 'annot.type', 'n', 'n_hyper', 'n_hypo', 'diff_meth_mean', 'diff_meth_median', 'diff_meth_sd'))
    expect_equal(l$gene_id, c('1', '1', '2'))
    expect_equal(l$annot.type, c('promoters', 'exons', 'promoters'))
    expect_equal(l$n, c(2L, 2L, 2L))
    expect_equal(l$diff_meth_mean, c(mean(c(10, -20)), mean(c(10, 30)), mean(c(-20, 30))))
})

test_that('summarize_genes() works without over or by', {
    g = summarize_genes(gene_example(), quiet = TRUE)

    expect_named(g, c('gene_id', 'symbol', 'n_regions', 'n_promoters', 'n_exons'))
})

test_that('summarize_genes() summarizes MANE annotations, alone or with genes', {
    a = gene_example()
    mane = a
    mane$annot$type = sub('hg19_genes_', 'hg38_mane_', mane$annot$type)

    g = summarize_genes(mane, quiet = TRUE)
    expect_named(g, c('gene_id', 'symbol', 'n_regions', 'n_promoters', 'n_exons'))
    expect_equal(g$n_regions, c(3L, 2L))

    # With both groups, a gene has one row, and the MANE types are prefixed
    both = c(a, mane[mane$annot$type == 'hg38_mane_promoters'])
    g = summarize_genes(both, quiet = TRUE)
    expect_named(g, c('gene_id', 'symbol', 'n_regions', 'n_promoters', 'n_exons', 'n_mane_promoters'))
    expect_equal(g$n_mane_promoters, c(2L, 2L))

    l = summarize_genes(both, format = 'long', quiet = TRUE)
    expect_equal(l$annot.type, c('promoters', 'exons', 'mane_promoters', 'promoters', 'mane_promoters'))
})

test_that('summarize_genes() errors for bad arguments', {
    a = gene_example()

    expect_error(summarize_genes(a, over = 'nonsense', quiet = TRUE), 'nonsense not column')
    expect_error(summarize_genes(a, by = c('DM_status', 'diff_meth'), quiet = TRUE), 'by must be a single column')
    expect_error(summarize_genes(a, format = 'tall', quiet = TRUE), "should be one of")

    # Regions annotated only to CpG annotations have no genes
    expect_error(summarize_genes(annotate_dm_regions(), quiet = TRUE), 'No regions are annotated to genes')
})
