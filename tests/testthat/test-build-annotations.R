################################################################################
# Errors

test_that('build_annotations() errors for custom annotations not in the cache', {
    expect_error(
        build_annotations(genome = 'hg19', annotations = 'hg19_custom_notreadyet'),
        'not in annotatr_cache')
})

test_that('build_annotations() errors for unsupported annotations', {
    expect_error(
        build_annotations(genome = 'hg19', annotations = 'hg19_genes_nonsense'),
        'not supported')
})

test_that('build_annotations() returns custom annotations without downloading', {
    read_annotations(con = extdata('test_annotations_3.bed'), name = 'buildcustom', genome = 'hg19', format = 'bed')
    a = build_annotations(genome = 'hg19', annotations = 'hg19_custom_buildcustom')

    expect_s4_class(a, 'GRanges')
    expect_length(a, 3)
    expect_equal(unique(a$type), 'hg19_custom_buildcustom')
})

################################################################################
# Light network tests: small downloads from each kind of data source, run
# everywhere unless offline. They use cache = FALSE so they test downloading and
# building, and catch changes in the sources.

expect_built = function(annotations, annots) {
    expect_s4_class(annotations, 'GRanges')
    expect_setequal(unique(annotations$type), expand_annotations(annots))
    expect_named(GenomicRanges::mcols(annotations), c('id', 'tx_id', 'gene_id', 'symbol', 'entrez_id', 'ensembl_id', 'type'))
}

test_that('CpG annotations build from UCSC', {
    skip_network()

    a = suppressMessages(build_annotations(genome = 'hg38', annotations = 'hg38_cpgs', cache = FALSE))
    expect_built(a, 'hg38_cpgs')
    expect_equal(unique(Seqinfo::genome(a)), 'hg38')

    # UCSC tables are 0-based, so islands have the widths in UCSC's length
    # column, e.g. the first is chr1:28736-29737 (issue #54)
    islands = a[a$type == 'hg38_cpg_islands']
    ucsc = readr::read_tsv(download_annotation_file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz', genome = 'hg38'),
        col_names = c('chr', 'start', 'end', 'length'), col_types = '-cii-i-----')
    ucsc_gr = GenomicRanges::GRanges(ucsc$chr, IRanges::IRanges(ucsc$start + 1L, ucsc$end), length = ucsc$length)
    m = S4Vectors::match(ucsc_gr, islands, ignore.strand = TRUE)
    expect_false(anyNA(m))
    expect_equal(GenomicRanges::width(ucsc_gr), ucsc_gr$length)

    # interCGI is only on chromosomes with CpG islands
    inter = a[a$type == 'hg38_cpg_inter']
    expect_true(all(as.character(GenomicRanges::seqnames(inter)) %in% as.character(GenomicRanges::seqnames(islands))))
})

test_that('FANTOM5 enhancers build', {
    skip_network()

    a = suppressMessages(build_annotations(genome = 'hg19', annotations = 'hg19_enhancers_fantom', cache = FALSE))
    expect_built(a, 'hg19_enhancers_fantom')
})

test_that('GENCODE lncRNA annotations build', {
    skip_network()
    skip_if_not_installed('org.Hs.eg.db')

    a = suppressMessages(build_annotations(genome = 'hg19', annotations = 'hg19_lncrna_gencode', cache = FALSE))
    expect_built(a, 'hg19_lncrna_gencode')
    # Every lncRNA comes with a gene symbol from GENCODE
    expect_false(anyNA(a$symbol))

    skip_if_not_installed('org.Mm.eg.db')
    a = suppressMessages(build_annotations(genome = 'mm10', annotations = 'mm10_lncrna_gencode', cache = FALSE))
    expect_built(a, 'mm10_lncrna_gencode')
})

test_that('MANE annotations build, with one transcript per gene', {
    skip_network()
    skip_if_not_installed('txdbmaker')

    a = suppressMessages(build_annotations(genome = 'hg38', annotations = c('hg38_mane_promoters', 'hg38_mane_exons'), cache = FALSE))
    expect_built(a, c('hg38_mane_promoters', 'hg38_mane_exons'))

    # MANE Select only: one promoter per gene, and no MANE Plus Clinical, e.g.
    # the second BRAF transcript
    promoters = a[a$type == 'hg38_mane_promoters']
    expect_false(anyDuplicated(promoters$gene_id) > 0)
    expect_gt(length(promoters), 19000)
    expect_equal(unique(a$tx_id[a$symbol %in% 'BRAF']), 'ENST00000646891.2')

    # Entrez gene IDs and symbols, as for hg38_genes_*
    expect_equal(unique(a$gene_id[a$symbol %in% 'BRAF']), '673')
    expect_equal(unique(a$entrez_id[a$symbol %in% 'BRAF']), '673')
    expect_equal(unique(a$ensembl_id[a$symbol %in% 'BRAF']), 'ENSG00000157764')
    expect_false(anyNA(a$symbol))
    expect_equal(unique(Seqinfo::genome(a)), 'hg38')
})

################################################################################
# Full network tests: every builtin genome, and the large downloads. These take
# a while, so they run only with ANNOTATR_FULL_TESTS=true.

# Each genome differs only in where its gene models (TxDb and org packages, or
# an EnsDb) and CpG islands come from, so basicgenes and cpgs test each genome.
# Building every annotation type for every genome uses more memory than an
# 8 GB Docker VM has.
for(genome in builtin_genomes()) {
    test_that(sprintf('Gene and CpG annotations build for %s', genome), {
        skip_if_not_full_tests()
        skip_network()
        if(genome %in% GENARK$genome) {
            skip_if_not_installed('ensembldb')
        } else {
            skip_if_not_installed(get_txdb_name(genome))
            skip_if_not_installed(sprintf('org.%s.eg.db', get_orgdb_name(genome)))
        }

        annots = sprintf('%s_basicgenes', genome)
        if(!(genome %in% c('dm3', 'dm6'))) {
            annots = c(annots, sprintf('%s_cpgs', genome))
        }

        a = suppressMessages(build_annotations(genome = genome, annotations = annots, cache = FALSE))
        expect_built(a, annots)

        # Most gene annotations have a symbol, Entrez ID, and Ensembl ID,
        # whatever kind of gene ID the source has (e.g. FlyBase for dm6)
        genic = a[grepl('_genes_', a$type) & !is.na(a$gene_id)]
        expect_gt(mean(!is.na(genic$symbol)), 0.5)
        expect_gt(mean(!is.na(genic$entrez_id)), 0.5)
        expect_gt(mean(!is.na(genic$ensembl_id)), 0.5)
    })
}

# ENCODE cCREs are large downloads (129 MB for hg38)
for(genome in names(CCRE$files)) {
    test_that(sprintf('ENCODE cCRE annotations build for %s', genome), {
        skip_if_not_full_tests()
        skip_network()

        annots = sprintf('%s_ccres', genome)
        a = suppressMessages(build_annotations(genome = genome, annotations = annots, cache = FALSE))
        expect_built(a, annots)
        expect_equal(unique(Seqinfo::genome(a)), genome)

        # The id is the cCRE accession, e.g. EH38E2776516 or EM10E0932225
        expect_false(anyDuplicated(a$id) > 0)
        expect_true(all(grepl('^E[HM][0-9]{2}E[0-9]+$', a$id)))
        # cCREs are 150 to 350 bp, so 1-based widths are too
        expect_true(all(GenomicRanges::width(a) >= 150 & GenomicRanges::width(a) <= 350))
    })
}

# Ensembl canonical transcripts come from a different EnsDb for each genome
for(genome in CANONICAL$genome) {
    test_that(sprintf('Canonical annotations build for %s, with one transcript per gene', genome), {
        skip_if_not_full_tests()
        skip_network()
        skip_if_not_installed('ensembldb')

        annots = c(sprintf('%s_canonical_promoters', genome), sprintf('%s_canonical_intergenic', genome))
        a = suppressMessages(build_annotations(genome = genome, annotations = annots, cache = FALSE))
        expect_built(a, annots)
        expect_equal(unique(Seqinfo::genome(a)), genome)
        expect_true(all(startsWith(Seqinfo::seqlevelsInUse(a), 'chr')))
        # No second copies of genes on alternate haplotypes or fix patches
        expect_false(any(grepl('_(alt|fix)$', Seqinfo::seqlevelsInUse(a))))

        promoters = a[a$type == sprintf('%s_canonical_promoters', genome)]
        expect_false(anyDuplicated(promoters$gene_id) > 0)
        expect_equal(promoters$ensembl_id, promoters$gene_id)
        expect_gt(mean(!is.na(promoters$symbol)), 0.5)
    })
}

# The remaining gene annotation types use the same code for every genome
test_that('All gene annotation types build', {
    skip_if_not_full_tests()
    skip_network()
    skip_if_not_installed('TxDb.Hsapiens.UCSC.hg19.knownGene')
    skip_if_not_installed('org.Hs.eg.db')

    annots = sprintf('hg19_genes_%s', c('intergenic', 'cds', 'firstexons',
        'intronexonboundaries', 'exonintronboundaries'))
    a = suppressMessages(build_annotations(genome = 'hg19', annotations = annots, cache = FALSE))
    expect_built(a, annots)

    # Intergenic is only on chromosomes with genes, so no contig is entirely intergenic
    intergenic = a[a$type == 'hg19_genes_intergenic']
    genic = a[a$type != 'hg19_genes_intergenic']
    expect_true(all(as.character(GenomicRanges::seqnames(intergenic)) %in% as.character(GenomicRanges::seqnames(genic))))
})

test_that('hg38 GENCODE lncRNA annotations build', {
    skip_if_not_full_tests()
    skip_network()
    skip_if_not_installed('org.Hs.eg.db')

    a = suppressMessages(build_annotations(genome = 'hg38', annotations = 'hg38_lncrna_gencode', cache = FALSE))
    expect_built(a, 'hg38_lncrna_gencode')
})

test_that('FANTOM5 enhancers build for hg38, mm9, and mm10', {
    skip_if_not_full_tests()
    skip_network()

    for(genome in c('hg38', 'mm9', 'mm10')) {
        annot = sprintf('%s_enhancers_fantom', genome)
        a = suppressMessages(build_annotations(genome = genome, annotations = annot, cache = FALSE))
        expect_built(a, annot)
        # Enhancers lifted over from hg19 and mm9 have the seqinfo of their genome
        expect_equal(unique(Seqinfo::genome(a)), genome)
        expect_false(anyNA(Seqinfo::seqlengths(a)))
    }
})

test_that('chromHMM chromatin state annotations build', {
    skip_if_not_full_tests()
    skip_network()

    a = suppressMessages(build_annotations(genome = 'hg19', annotations = 'hg19_Gm12878-chromatin', cache = FALSE))
    expect_built(a, 'hg19_Gm12878-chromatin')

    # The states tile the genome, so with 1-based coordinates they don't
    # overlap. Reading UCSC's 0-based starts as 1-based overlapped them by 1 bp.
    expect_equal(sum(as.numeric(GenomicRanges::width(GenomicRanges::reduce(a)))), sum(as.numeric(GenomicRanges::width(a))))
})
