context('Test build annotations module')

################################################################################
# Test errors

test_that('Test error for non-existent custom annotations', {
    expect_error(
        build_annotations(genome = 'hg19', annotations = 'hg19_custom_ezh2'),
        'not in annotatr_cache'
    )
})

################################################################################
# Test annotations that aren't otherwise tested
# intergenic, cds, firstexons, and both boundaries

# test_that('Test all annotations', {
#     annots = c('hg19_basicgenes', 'hg19_cpgs', 'hg19_genes_intergenic', 'hg19_genes_cds', 'hg19_genes_firstexons', 'hg19_genes_intronexonboundaries', 'hg19_genes_exonintronboundaries', 'hg19_lncrna_gencode', 'hg19_Gm12878-chromatin', 'hg19_H1hesc-chromatin', 'hg19_Hepg2-chromatin', 'hg19_Hmec-chromatin', 'hg19_Hsmm-chromatin', 'hg19_Huvec-chromatin', 'hg19_K562-chromatin', 'hg19_Nhek-chromatin', 'hg19_Nhlf-chromatin')
#     annotations = build_annotations(genome = 'hg19', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('hg38_basicgenes', 'hg38_cpgs', 'hg38_genes_intergenic', 'hg38_genes_cds', 'hg38_genes_firstexons', 'hg38_genes_intronexonboundaries', 'hg38_genes_exonintronboundaries', 'hg38_lncrna_gencode', 'hg38_enhancers_fantom')
#     annotations = build_annotations(genome = 'hg38', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('mm10_basicgenes', 'mm10_cpgs', 'mm10_genes_intergenic', 'mm10_genes_cds', 'mm10_genes_firstexons', 'mm10_genes_intronexonboundaries', 'mm10_genes_exonintronboundaries', 'mm10_lncrna_gencode', 'mm10_enhancers_fantom')
#     annotations = build_annotations(genome = 'mm10', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('mm9_basicgenes', 'mm9_cpgs', 'mm9_genes_intergenic', 'mm9_genes_cds', 'mm9_genes_firstexons', 'mm9_genes_intronexonboundaries', 'mm9_genes_exonintronboundaries', 'mm9_enhancers_fantom')
#     annotations = build_annotations(genome = 'mm9', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('rn4_basicgenes', 'rn4_cpgs', 'rn4_genes_intergenic', 'rn4_genes_cds', 'rn4_genes_firstexons', 'rn4_genes_intronexonboundaries', 'rn4_genes_exonintronboundaries')
#     annotations = build_annotations(genome = 'rn4', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('rn5_basicgenes', 'rn5_cpgs', 'rn5_genes_intergenic', 'rn5_genes_cds', 'rn5_genes_firstexons', 'rn5_genes_intronexonboundaries', 'rn5_genes_exonintronboundaries')
#     annotations = build_annotations(genome = 'rn5', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('rn6_basicgenes', 'rn6_cpgs', 'rn6_genes_intergenic', 'rn6_genes_cds', 'rn6_genes_firstexons', 'rn6_genes_intronexonboundaries', 'rn6_genes_exonintronboundaries')
#     annotations = build_annotations(genome = 'rn6', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('dm3_basicgenes', 'dm3_genes_intergenic', 'dm3_genes_cds', 'dm3_genes_firstexons', 'dm3_genes_intronexonboundaries', 'dm3_genes_exonintronboundaries')
#     annotations = build_annotations(genome = 'dm3', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
#
#     annots = c('dm6_basicgenes', 'dm6_genes_intergenic', 'dm6_genes_cds', 'dm6_genes_firstexons', 'dm6_genes_intronexonboundaries', 'dm6_genes_exonintronboundaries')
#     annotations = build_annotations(genome = 'dm6', annotations = annots)
#     expect_true( dplyr::setequal(unique(annotations$type), expand_annotations(annots)) )
# })
