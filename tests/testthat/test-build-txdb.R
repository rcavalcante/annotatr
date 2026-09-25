# A tiny TxDb whose annotations can be worked out by hand:
#
#   gene  tx  strand  exons               CDS
#   g1    t1  +       100-200, 300-400    150-200, 300-350
#   g2    t2  -       1000-1100           1020-1080
tiny_txdb = function(gene_ids = c('g1', 'g2'), genome = NA) {
    # GTF-style features, as from rtracklayer::import() of a GTF
    gr = GenomicRanges::GRanges(
        seqnames = 'chr1',
        ranges = IRanges::IRanges(
            start = c(100, 100, 300, 150, 300, 1000, 1000, 1020),
            end = c(400, 200, 400, 200, 350, 1100, 1100, 1080)),
        strand = rep(c('+', '-'), c(5, 3)),
        type = c('transcript', 'exon', 'exon', 'CDS', 'CDS', 'transcript', 'exon', 'CDS'),
        phase = c(NA, NA, NA, 0L, 0L, NA, NA, 0L),
        gene_id = rep(gene_ids, c(5, 3)),
        transcript_id = rep(c('t1', 't2'), c(5, 3)),
        seqinfo = Seqinfo::Seqinfo('chr1', 5000, FALSE, genome))
    suppressWarnings(txdbmaker::makeTxDbFromGRanges(gr))
}

test_that('build_txdb_annotations() builds gene annotations from a TxDb', {
    skip_if_not_installed('txdbmaker')

    a = build_txdb_annotations(tiny_txdb(), genome = 'hg19', group = 'mytx')

    expect_s4_class(a, 'GRanges')
    expect_named(GenomicRanges::mcols(a), c('id', 'tx_id', 'gene_id', 'symbol', 'entrez_id', 'ensembl_id', 'type'))
    expect_setequal(unique(a$type), sprintf('hg19_mytx_%s', c('1to5kb', 'promoters', '5UTRs', 'exons', 'introns', '3UTRs')))
    expect_equal(unique(Seqinfo::genome(a)), 'hg19')

    exons = a[a$type == 'hg19_mytx_exons']
    expect_equal(IRanges::start(exons), c(100, 300, 1000))
    expect_equal(exons$tx_id, c('t1', 't1', 't2'))
    expect_equal(exons$gene_id, c('g1', 'g1', 'g2'))
    expect_true(all(is.na(exons$symbol)))
    # Without an OrgDb, the kind of gene ID is unknown
    expect_true(all(is.na(exons$entrez_id)))
    expect_true(all(is.na(exons$ensembl_id)))

    introns = a[a$type == 'hg19_mytx_introns']
    expect_equal(as.character(IRanges::ranges(introns)), '201-299')

    # Promoters are 1kb upstream of the TSS, trimmed at the chromosome start
    promoters = a[a$type == 'hg19_mytx_promoters']
    expect_equal(as.character(IRanges::ranges(promoters)), c('1-99', '1101-2100'))
})

test_that('build_txdb_annotations() builds the requested types, including intergenic', {
    skip_if_not_installed('txdbmaker')

    a = build_txdb_annotations(tiny_txdb(), genome = 'hg19', group = 'mytx', annotations = c('cds', 'intergenic'))
    expect_setequal(unique(a$type), c('hg19_mytx_cds', 'hg19_mytx_intergenic'))
    expect_equal(sum(a$type == 'hg19_mytx_cds'), 3)
    expect_true(all(is.na(a$gene_id[a$type == 'hg19_mytx_intergenic'])))
})

test_that('build_txdb_annotations() gets symbols from an OrgDb', {
    skip_if_not_installed('txdbmaker')
    skip_if_not_installed('org.Hs.eg.db')
    orgdb = org.Hs.eg.db::org.Hs.eg.db

    a = build_txdb_annotations(tiny_txdb(c('1', '2')), genome = 'hg19', group = 'mytx', annotations = 'exons', orgdb = orgdb)
    expect_equal(a$symbol, c('A1BG', 'A1BG', 'A2M'))
    expect_equal(a$entrez_id, c('1', '1', '2'))
    expect_equal(a$ensembl_id, c('ENSG00000121410', 'ENSG00000121410', 'ENSG00000175899'))

    # Ensembl gene IDs, with versions as in GENCODE
    txdb = tiny_txdb(c('ENSG00000121410.14', 'ENSG00000175899.17'))
    a = build_txdb_annotations(txdb, genome = 'hg19', group = 'gencode', annotations = 'exons', orgdb = orgdb, keytype = 'ENSEMBL')
    expect_equal(a$gene_id, c('ENSG00000121410.14', 'ENSG00000121410.14', 'ENSG00000175899.17'))
    expect_equal(a$symbol, c('A1BG', 'A1BG', 'A2M'))
    expect_equal(a$entrez_id, c('1', '1', '2'))
    expect_equal(a$ensembl_id, c('ENSG00000121410', 'ENSG00000121410', 'ENSG00000175899'))

    expect_error(build_txdb_annotations(txdb, genome = 'hg19', group = 'gencode', orgdb = orgdb, keytype = 'NONSENSE'), 'not a keytype of orgdb')
    expect_error(build_txdb_annotations(txdb, genome = 'hg19', group = 'gencode', orgdb = 'org.Hs.eg.db'), 'orgdb must be an OrgDb')
})

test_that('build_txdb_annotations() checks its arguments', {
    skip_if_not_installed('txdbmaker')
    txdb = tiny_txdb()

    expect_error(build_txdb_annotations('txdb', genome = 'hg19', group = 'mytx'), 'txdb must be a TxDb or EnsDb')
    expect_error(build_txdb_annotations(txdb, group = 'mytx'), 'genome must be')
    expect_error(build_txdb_annotations(txdb, genome = 'hg_19', group = 'mytx'), 'genome must be')
    expect_error(build_txdb_annotations(txdb, genome = 'hg19'), 'group must be')
    expect_error(build_txdb_annotations(txdb, genome = 'hg19', group = 'my_tx'), 'group must be')
    expect_error(build_txdb_annotations(txdb, genome = 'hg19', group = 'genes'), 'built-in group')
    expect_error(build_txdb_annotations(txdb, genome = 'hg19', group = 'Custom'), 'built-in group')
    expect_error(build_txdb_annotations(txdb, genome = 'hg19', group = 'mytx', annotations = c('exons', 'enhancers')), '"enhancers" is\\(are\\) not gene annotation types')

    # The genome must agree with the TxDb's genome, if it has one
    expect_error(build_txdb_annotations(tiny_txdb(genome = 'mm39'), genome = 'mm10', group = 'mytx'), 'genome is mm10, but txdb is from genome mm39')
    a = build_txdb_annotations(tiny_txdb(genome = 'mm39'), genome = 'mm39', group = 'mytx', annotations = 'exons')
    expect_equal(unique(Seqinfo::genome(a)), 'mm39')
})

test_that('Annotations from build_txdb_annotations() work with annotate_regions() and summarize_genes()', {
    skip_if_not_installed('txdbmaker')

    annots = build_txdb_annotations(tiny_txdb(), genome = 'hg19', group = 'mytx', annotations = c('promoters', 'exons', 'intergenic'))
    regions = GenomicRanges::GRanges('chr1', IRanges::IRanges(c(50, 150, 1050, 3000), width = 10),
        seqinfo = Seqinfo::Seqinfo('chr1', 5000, FALSE, 'hg19'))
    annotated = annotate_regions(regions, annots, quiet = TRUE)

    g = summarize_genes(annotated, quiet = TRUE)
    expect_equal(g$gene_id, c('g1', 'g2'))
    expect_named(g, c('gene_id', 'symbol', 'entrez_id', 'ensembl_id', 'n_regions', 'n_promoters', 'n_exons'))
    expect_equal(g$n_promoters, c(1L, 0L))
    expect_equal(g$n_exons, c(1L, 1L))

    # Alongside another group, the types are prefixed with their group
    other = annots
    other$type = sub('_mytx_', '_other_', other$type)
    g = summarize_genes(annotate_regions(regions, c(annots, other), quiet = TRUE), quiet = TRUE)
    expect_named(g, c('gene_id', 'symbol', 'entrez_id', 'ensembl_id', 'n_regions', 'n_mytx_promoters', 'n_mytx_exons', 'n_other_promoters', 'n_other_exons'))
})
