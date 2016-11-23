context('Test build annotations module')

################################################################################
# Test errors

test_that('Test error for genome without enhancers', {
    expect_error(
        build_annotations(genome = 'hg38', annotations = 'hg38_enhancer_fantom'),
        'not supported'
    )
})

test_that('Test error for non-existent custom annotations', {
    expect_error(
        build_annotations(genome = 'hg19', annotations = 'hg19_custom_ezh2'),
        'not in annotatr_cache'
    )
})

################################################################################
# Test annotations that aren't otherwise tested
# intergenic, cds, firstexons, and both boundaries

test_that('Test otherwise untested annotations', {
    annots = c('hg19_basicgenes', 'hg19_cpgs', 'hg19_genes_intergenic', 'hg19_genes_cds', 'hg19_genes_firstexons', 'hg19_genes_intronexonboundaries', 'hg19_genes_exonintronboundaries', 'hg19_lncrna_gencode')
    annotations = build_annotations(genome = 'hg19', annotations = annots)
    expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )

    annotations = build_annotations(genome = 'mm9', annotations = 'mm9_enhancers_fantom')
    expect_true( dplyr::setequal(unique(annotations$type), 'mm9_enhancers_fantom'))

    annotations = build_annotations(genome = 'mm10', annotations = 'mm10_lncrna_gencode')
    expect_true( dplyr::setequal(unique(annotations$type), 'mm10_lncrna_gencode'))

    annotations = build_annotations(genome = 'hg38', annotations = 'hg38_lncrna_gencode')
    expect_true( dplyr::setequal(unique(annotations$type), 'hg38_lncrna_gencode'))

    annotations = build_annotations(genome = 'hg38', annotations = 'hg38_cpgs')
    expect_true( dplyr::setequal(unique(annotations$type), c('hg38_cpg_islands','hg38_cpg_shores','hg38_cpg_shelves','hg38_cpg_inter')) )

    annotations = build_annotations(genome = 'mm10', annotations = 'mm10_cpgs')
    expect_true( dplyr::setequal(unique(annotations$type), c('mm10_cpg_islands','mm10_cpg_shores','mm10_cpg_shelves','mm10_cpg_inter')) )

    annotations = build_annotations(genome = 'rn6', annotations = 'rn6_cpgs')
    expect_true( dplyr::setequal(unique(annotations$type), c('rn6_cpg_islands','rn6_cpg_shores','rn6_cpg_shelves','rn6_cpg_inter')) )
})
