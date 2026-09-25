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
# everywhere unless offline

expect_built = function(annotations, annots) {
    expect_s4_class(annotations, 'GRanges')
    expect_setequal(unique(annotations$type), expand_annotations(annots))
    expect_named(GenomicRanges::mcols(annotations), c('id', 'tx_id', 'gene_id', 'symbol', 'type'))
}

test_that('CpG annotations build from UCSC', {
    skip_network()

    a = suppressMessages(build_annotations(genome = 'hg38', annotations = 'hg38_cpgs'))
    expect_built(a, 'hg38_cpgs')
    expect_equal(unique(Seqinfo::genome(a)), 'hg38')
})

test_that('FANTOM5 enhancers build', {
    skip_network()

    a = suppressMessages(build_annotations(genome = 'hg19', annotations = 'hg19_enhancers_fantom'))
    expect_built(a, 'hg19_enhancers_fantom')
})

test_that('GENCODE lncRNA annotations build', {
    skip_network()
    skip_if_not_installed('org.Hs.eg.db')

    a = suppressMessages(build_annotations(genome = 'hg19', annotations = 'hg19_lncrna_gencode'))
    expect_built(a, 'hg19_lncrna_gencode')
    # Every lncRNA comes with a gene symbol from GENCODE
    expect_false(anyNA(a$symbol))

    skip_if_not_installed('org.Mm.eg.db')
    a = suppressMessages(build_annotations(genome = 'mm10', annotations = 'mm10_lncrna_gencode'))
    expect_built(a, 'mm10_lncrna_gencode')
})

################################################################################
# Full network tests: every builtin genome, and the large downloads. These take
# a while, so they run only with ANNOTATR_FULL_TESTS=true.

for(genome in builtin_genomes()) {
    test_that(sprintf('All gene and CpG annotations build for %s', genome), {
        skip_if_not_full_tests()
        skip_network()
        if(genome %in% GENARK$genome) {
            skip_if_not_installed('ensembldb')
        } else {
            skip_if_not_installed(get_txdb_name(genome))
            skip_if_not_installed(sprintf('org.%s.eg.db', get_orgdb_name(genome)))
        }

        annots = sprintf('%s_%s', genome, c('basicgenes', 'genes_intergenic', 'genes_cds',
            'genes_firstexons', 'genes_intronexonboundaries', 'genes_exonintronboundaries'))
        if(!(genome %in% c('dm3', 'dm6'))) {
            annots = c(annots, sprintf('%s_cpgs', genome))
        }

        a = suppressMessages(build_annotations(genome = genome, annotations = annots))
        expect_built(a, annots)
    })
}

test_that('hg38 GENCODE lncRNA annotations build', {
    skip_if_not_full_tests()
    skip_network()
    skip_if_not_installed('org.Hs.eg.db')

    a = suppressMessages(build_annotations(genome = 'hg38', annotations = 'hg38_lncrna_gencode'))
    expect_built(a, 'hg38_lncrna_gencode')
})

test_that('FANTOM5 enhancers build for hg38, mm9, and mm10', {
    skip_if_not_full_tests()
    skip_network()

    for(genome in c('hg38', 'mm9', 'mm10')) {
        annot = sprintf('%s_enhancers_fantom', genome)
        a = suppressMessages(build_annotations(genome = genome, annotations = annot))
        expect_built(a, annot)
    }
})

test_that('chromHMM chromatin state annotations build', {
    skip_if_not_full_tests()
    skip_network()

    a = suppressMessages(build_annotations(genome = 'hg19', annotations = 'hg19_Gm12878-chromatin'))
    expect_built(a, 'hg19_Gm12878-chromatin')
})
