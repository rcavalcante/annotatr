test_that('get_txdb_name() and get_orgdb_name() map genomes to packages', {
    expect_equal(get_txdb_name('hg19'), 'TxDb.Hsapiens.UCSC.hg19.knownGene')
    expect_equal(get_txdb_name('mm10'), 'TxDb.Mmusculus.UCSC.mm10.knownGene')
    expect_equal(get_orgdb_name('hg38'), 'Hs')
    expect_equal(get_orgdb_name('rn6'), 'Rn')
})

test_that('get_txdb_name() and get_orgdb_name() reject unsupported genomes', {
    expect_error(get_txdb_name(genome = 'hg18'), 'should be one of')
    expect_error(get_orgdb_name(genome = 'hg18'), 'should be one of')
})

test_that('GenArk genomes are builtin', {
    expect_true(all(GENARK$genome %in% builtin_genomes()))
    expect_contains(builtin_annotations(), c('oviariramb2_basicgenes', 'oviariramb2_cpgs'))

    expect_equal(
        get_genark_url('oviariramb2', 'GCF_016772045.1.chromAlias.txt'),
        'https://hgdownload.soe.ucsc.edu/hubs/GCF/016/772/045/GCF_016772045.1/GCF_016772045.1.chromAlias.txt')
})

test_that('builtin_annotations() has CpG annotations for all genomes but fly', {
    annots = builtin_annotations()

    expect_contains(annots, sprintf('%s_cpgs', setdiff(builtin_genomes(), c('dm3', 'dm6'))))
    expect_false(any(c('dm3_cpgs', 'dm6_cpgs') %in% annots))
})

test_that('builtin_annotations() has canonical annotations, with intergenic, for current assemblies', {
    annots = builtin_annotations()
    canonical = grep('_canonical_|basiccanonical', annots, value = TRUE)

    expect_setequal(unique(sub('_.*', '', canonical)), c('hg38', 'mm39', 'rn7', 'danRer11', 'dm6', 'oviariramb2'))
    expect_contains(canonical, c('mm39_basiccanonical', 'mm39_canonical_promoters', 'mm39_canonical_intergenic'))
    expect_false('hg19_canonical_promoters' %in% annots)

    expect_setequal(expand_annotations('mm39_basiccanonical'), sprintf('mm39_canonical_%s', c('1to5kb', 'promoters', '5UTRs', 'exons', 'introns', '3UTRs')))
    expect_named(tidy_annotations(c('mm39_canonical_promoters', 'mm39_canonical_firstexons')), c('canonical promoters', 'canonical first exons'))
})

test_that('builtin_annotations() has MANE annotations for hg38 only, without intergenic', {
    annots = builtin_annotations()
    mane = grep('_mane_|basicmane', annots, value = TRUE)

    expect_contains(mane, c('hg38_basicmane', 'hg38_mane_promoters', 'hg38_mane_exons', 'hg38_mane_cds'))
    expect_true(all(startsWith(mane, 'hg38_')))
    expect_false('hg38_mane_intergenic' %in% annots)
    expect_length(grep('_mane_', mane), length(grep('hg38_genes_', annots)) - 1)
})

test_that('tidy_annotations() gives readable names', {
    hg19_annots = c('hg19_cpg_islands', 'hg19_cpg_inter', 'hg19_genes_firstexons',
        'hg19_genes_intronexonboundaries', 'hg19_genes_exonintronboundaries',
        'hg19_lncrna_gencode', 'hg19_chromatin_Gm12878-ActivePromoter')
    mm9_annots = c('mm9_cpg_islands', 'mm9_genes_exonsCDSs', 'mm9_cpg_inter')
    rn4_custom_annots = c('rn4_custom_cpgislands', 'rn4_custom_TFBS')

    expect_named(tidy_annotations(hg19_annots), c('CpG islands', 'interCGI', 'first exons',
        'intron/exon boundaries', 'exon/intron boundaries', 'GENCODE lncRNA', 'Gm12878-ActivePromoter'))
    expect_named(tidy_annotations(mm9_annots), c('CpG islands', 'exonsCDSs', 'interCGI'))
    expect_named(tidy_annotations(rn4_custom_annots), c('cpgislands', 'TFBS'))
    expect_named(tidy_annotations(c('hg38_genes_promoters', 'hg38_mane_promoters', 'hg38_mane_firstexons')),
        c('promoters', 'MANE promoters', 'MANE first exons'))
    expect_named(tidy_annotations(c('mm39_gencode_promoters', 'mm39_gencode_exonintronboundaries', 'mm39_custom_TFBS')),
        c('gencode promoters', 'gencode exon/intron boundaries', 'TFBS'))

    # The values map back to the original annotation codes
    expect_equal(unname(unlist(tidy_annotations(mm9_annots))), mm9_annots)
})

test_that('check_annotations() accepts valid annotations', {
    expect_null(check_annotations(c('hg19_genes_promoters', 'hg19_cpgs')))
    expect_null(check_annotations(c('hg19_cpgs', 'hg19_custom_TFBS')))
    expect_null(check_annotations(c('hg19_custom_TFBS')))
})

test_that('check_annotations() rejects invalid annotations', {
    expect_error(check_annotations(c('hg17_genes_promoters', 'hg19_cpgs')), 'not supported. See builtin_annotations()')
    expect_error(check_annotations(c('hello', 'hg19_genes_promoters', 'hg19_cpgs')), 'not supported. See builtin_annotations()')
    expect_error(check_annotations(c('hg19_genes_promoters', 'mm9_cpg_islands')), 'genome prefix on all annotations must be the same')
    expect_error(check_annotations('hg19_mane_promoters'), 'not supported. See builtin_annotations()')
    expect_error(check_annotations('hg38_mane_intergenic'), 'not supported. See builtin_annotations()')
})

test_that('expand_annotations() expands shortcuts', {
    # No shortcuts
    annots = c('hg19_genes_promoters', 'hg19_genes_exons')
    expect_setequal(expand_annotations(annots), annots)

    expect_setequal(expand_annotations(c('mm9_basicgenes', 'mm9_cpgs')), c(
        'mm9_cpg_islands', 'mm9_cpg_shores', 'mm9_cpg_shelves', 'mm9_cpg_inter',
        'mm9_genes_1to5kb', 'mm9_genes_promoters', 'mm9_genes_5UTRs', 'mm9_genes_exons',
        'mm9_genes_introns', 'mm9_genes_3UTRs'))

    expect_setequal(expand_annotations(c('hg38_basicmane', 'hg38_genes_promoters')), c(
        'hg38_genes_promoters', 'hg38_mane_1to5kb', 'hg38_mane_promoters', 'hg38_mane_5UTRs',
        'hg38_mane_exons', 'hg38_mane_introns', 'hg38_mane_3UTRs'))

    # A shortcut overlapping an a la carte annotation doesn't duplicate it
    expect_setequal(expand_annotations(c('hg19_cpg_shores', 'hg19_cpgs')), cpg_types)

    expect_setequal(expand_annotations('hg19_Hepg2-chromatin'), sprintf('hg19_chromatin_Hepg2-%s', c(
        'ActivePromoter', 'WeakPromoter', 'PoisedPromoter', 'StrongEnhancer', 'WeakEnhancer',
        'Insulator', 'TxnTransition', 'TxnElongation', 'WeakTxn', 'Repressed',
        'Heterochrom/lo', 'Repetitive/CNV')))
})

test_that('get_chrom_aliases() maps Ensembl and other names to UCSC names', {
    skip_network()

    # UCSC database genome, with a 'ucsc' column in the header
    aliases = get_chrom_aliases('mm39', cache = FALSE)
    expect_equal(unname(aliases[c('1', 'MT', 'chr1', 'NC_000067.7')]), c('chr1', 'chrM', 'chr1', 'chr1'))

    # hg38's header has no 'ucsc' column, so the first column is the UCSC name
    aliases = get_chrom_aliases('hg38', cache = FALSE)
    expect_equal(unname(aliases[c('1', 'X', 'KI270728.1')]), c('chr1', 'chrX', 'chr16_KI270728v1_random'))

    # GenArk genome
    aliases = get_chrom_aliases('oviariramb2', cache = FALSE)
    expect_equal(unname(aliases['MT']), 'chrM')
})

test_that('ucsc_seqlevels() renames to UCSC names, and drops sequences without one', {
    skip_network()

    gr = GenomicRanges::GRanges(c('1', 'MT', 'not_a_chromosome'), IRanges::IRanges(1, 10))
    gr = ucsc_seqlevels(gr, 'mm39', cache = FALSE)

    expect_equal(as.character(GenomicRanges::seqnames(gr)), c('chr1', 'chrM'))
    expect_equal(unique(Seqinfo::genome(gr)), 'mm39')
    expect_equal(unname(Seqinfo::seqlengths(gr)['chr1']), 195154279)
})

test_that('standardize_mcols() gives annotations the same columns', {
    gr = GenomicRanges::GRanges('chr1', IRanges::IRanges(1, 10), type = 'hg19_custom_x', id = 'x:1', score = 5)
    gr = standardize_mcols(gr)

    expect_named(GenomicRanges::mcols(gr), c('id', 'tx_id', 'gene_id', 'symbol', 'entrez_id', 'ensembl_id', 'type'))
    expect_equal(gr$id, 'x:1')
    expect_true(is.na(gr$entrez_id))
})

test_that('get_gene_table() maps gene IDs, using the first of several matches', {
    skip_if_not_installed('org.Hs.eg.db')
    orgdb = org.Hs.eg.db::org.Hs.eg.db

    t = get_gene_table(c('1', '2', '1', NA, 'nonsense'), orgdb = orgdb, keytype = 'ENTREZID')
    expect_named(t, c('gene_id', 'symbol', 'entrez_id', 'ensembl_id'))
    expect_equal(t$gene_id, c('1', '2', 'nonsense'))
    expect_equal(t$symbol, c('A1BG', 'A2M', NA))
    expect_equal(t$ensembl_id, c('ENSG00000121410', 'ENSG00000175899', NA))

    # Ensembl IDs with versions
    t = get_gene_table(c('ENSG00000121410.14', 'ENSG00000175899.17'), orgdb = orgdb, keytype = 'ENSEMBL')
    expect_equal(t$gene_id, c('ENSG00000121410.14', 'ENSG00000175899.17'))
    expect_equal(t$ensembl_id, c('ENSG00000121410', 'ENSG00000175899'))
    expect_equal(t$entrez_id, c('1', '2'))

    # Without an OrgDb, only the IDs themselves
    t = get_gene_table(c('1', '2'), keytype = 'ENTREZID')
    expect_equal(t$entrez_id, c('1', '2'))
    expect_true(all(is.na(t$symbol)))
    expect_equal(nrow(get_gene_table(character(0))), 0)
})

test_that('set_genome_seqinfo() gives ranges the seqinfo of a genome', {
    gr = GenomicRanges::GRanges('chr1', IRanges::IRanges(1, 1e9))

    # Unknown genomes get only the genome
    unknown = set_genome_seqinfo(gr, 'notagenome')
    expect_equal(unname(Seqinfo::genome(unknown)), 'notagenome')

    skip_network()
    hg19 = set_genome_seqinfo(gr, 'hg19')
    expect_equal(unique(unname(Seqinfo::genome(hg19))), 'hg19')
    expect_equal(unname(Seqinfo::seqlengths(hg19)['chr1']), 249250621)
    # Trimmed to the end of chr1
    expect_equal(GenomicRanges::end(hg19), 249250621)
})
