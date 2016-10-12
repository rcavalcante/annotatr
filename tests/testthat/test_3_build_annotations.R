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

test_that('Test error for genome with no CpGs in AnnotationHub', {
    expect_error(
        build_annotations(genome = 'hg38', annotations = 'hg38_cpg_islands'),
        'AnnotationHub does not contain CpG Island annotations'
    )
})

################################################################################
# Test annotations that aren't otherwise tested
# intergenic, cds, firstexons, and both boundaries

test_that('Test otherwise untested annotations', {
    annots = c('hg19_basicgenes', 'hg19_cpgs', 'hg19_genes_intergenic', 'hg19_genes_cds', 'hg19_genes_firstexons', 'hg19_genes_intronexonboundaries', 'hg19_genes_exonintronboundaries')
    annotations = build_annotations(genome = 'hg19', annotations = annots)
    expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )

    annotations = build_annotations(genome = 'mm9', annotations = 'mm9_enhancers_fantom')
    expect_true( dplyr::setequal(unique(annotations$type), 'mm9_enhancers_fantom'))
})
