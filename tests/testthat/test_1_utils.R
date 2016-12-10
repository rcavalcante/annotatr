context('Test utility functions')

test_that('Test get_*_name() functions', {
    expect_error(
        get_txdb_name(genome = 'hg18'),
        'should be one of')

    expect_error(
        get_orgdb_name(genome = 'hg18'),
        'should be one of')
})

test_that('Test tidy_annotations()', {
    hg19_annots = c('hg19_cpg_islands', 'hg19_cpg_inter', 'hg19_genes_firstexons', 'hg19_genes_intronexonboundaries', 'hg19_genes_exonintronboundaries', 'hg19_lncrna_gencode', 'hg19_chromatin_Gm12878-ActivePromoter')
    mm9_annots = c('mm9_cpg_islands','mm9_genes_exonsCDSs','mm9_cpg_inter')
    rn4_custom_annots = c('rn4_custom_cpgislands','rn4_custom_TFBS')

    hg19_tidy_annots = tidy_annotations(hg19_annots)
    mm9_tidy_annots = tidy_annotations(mm9_annots)
    rn4_tidy_annots = tidy_annotations(rn4_custom_annots)

    expect_equal( all(names(hg19_tidy_annots) == c('CpG islands', 'interCGI', 'first exons', 'intron/exon boundaries', 'exon/intron boundaries', 'GENCODE lncRNA', 'Gm12878-ActivePromoter')), expected = TRUE)
    expect_equal( all(names(mm9_tidy_annots) == c('CpG islands', 'exonsCDSs', 'interCGI')), expected = TRUE)
    expect_equal( all(names(rn4_tidy_annots) == c('cpgislands', 'TFBS')), expected = TRUE)
    expect_equal
})

test_that('Test check_annotations()', {
    annots1 = c('hg17_genes_promoters','hg19_cpgs')
    annots2 = c('hello','hg19_genes_promoters','hg19_cpgs')
    annots3 = c('hg19_genes_promoters', 'mm9_cpg_islands')

    expect_error( check_annotations(annots1), 'not supported. See supported_annotations()' )
    expect_error( check_annotations(annots2), 'not supported. See supported_annotations()' )
    expect_error( check_annotations(annots3), 'genome prefix on all annotations must be the same' )
})

test_that('Test expand_annotations()', {
    annots1 = c('hg19_genes_promoters', 'hg19_genes_exons')

    annots2 = c('mm9_basicgenes', 'mm9_cpgs')
    expanded_annots2 = c('mm9_cpg_islands', 'mm9_cpg_shores', 'mm9_cpg_shelves', 'mm9_cpg_inter', 'mm9_genes_1to5kb', 'mm9_genes_promoters', 'mm9_genes_5UTRs', 'mm9_genes_exons', 'mm9_genes_introns', 'mm9_genes_3UTRs')

    annots3 = c('hg19_cpg_shores', 'hg19_cpgs')
    expanded_annots3 = c('hg19_cpg_islands', 'hg19_cpg_shores', 'hg19_cpg_shelves','hg19_cpg_inter')

    annots4 = c('hg19_Hepg2-chromatin')
    expanded_annots4 = c('hg19_chromatin_Hepg2-ActivePromoter','hg19_chromatin_Hepg2-WeakPromoter','hg19_chromatin_Hepg2-PoisedPromoter','hg19_chromatin_Hepg2-StrongEnhancer','hg19_chromatin_Hepg2-WeakEnhancer','hg19_chromatin_Hepg2-Insulator','hg19_chromatin_Hepg2-TxnTransition','hg19_chromatin_Hepg2-TxnElongation','hg19_chromatin_Hepg2-WeakTxn','hg19_chromatin_Hepg2-Repressed','hg19_chromatin_Hepg2-Heterochrom/lo','hg19_chromatin_Hepg2-Repetitive/CNV')

    expect_equal( dplyr::setequal(expand_annotations(annots1), annots1), expected = TRUE )
    expect_equal( dplyr::setequal(expand_annotations(annots2), expanded_annots2), expected = TRUE )
    expect_equal( dplyr::setequal(expand_annotations(annots3), expanded_annots3), expected = TRUE )
    expect_equal( dplyr::setequal(expand_annotations(annots4),
    expanded_annots4), expected = TRUE)
})
