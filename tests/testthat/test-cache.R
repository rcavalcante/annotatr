################################################################################
# Offline tests: seed the cache directly with the premade CpG annotations

cached_islands = function() {
    a = annotatr::annotations
    a[a$type == 'hg19_cpg_islands']
}

test_that('build_annotations() loads cached annotations without building', {
    local_empty_cache()
    islands = cached_islands()
    save_cached_annotation(islands, 'hg19_cpg_islands', 'hg19')

    expect_message(
        a <- build_annotations(genome = 'hg19', annotations = 'hg19_cpg_islands'),
        'Loading hg19_cpg_islands from the cache')
    expect_equal(a, islands)
})

test_that('list_cached_annotations() describes cached items', {
    local_empty_cache()
    expect_equal(nrow(list_cached_annotations()), 0)

    save_cached_annotation(cached_islands(), 'hg19_cpg_islands', 'hg19')
    cached = list_cached_annotations()

    expect_equal(nrow(cached), 1)
    expect_equal(cached$type, 'annotation')
    expect_equal(cached$genome, 'hg19')
    expect_equal(cached$name, 'hg19_cpg_islands')
    expect_match(cached$sources, 'annotatr ')
    expect_true(file.exists(cached$path))
})

test_that('Saving an annotation replaces earlier versions of it', {
    local_empty_cache()
    islands = cached_islands()

    save_cached_annotation(islands[1:10], 'hg19_cpg_islands', 'hg19')
    save_cached_annotation(islands, 'hg19_cpg_islands', 'hg19')

    expect_equal(nrow(list_cached_annotations()), 1)
    expect_equal(load_cached_annotation('hg19_cpg_islands', 'hg19'), islands)
})

test_that('A corrupted cached annotation is removed', {
    local_empty_cache()
    save_cached_annotation(cached_islands(), 'hg19_cpg_islands', 'hg19')
    writeLines('not an RDS file', list_cached_annotations()$path)

    expect_message(
        expect_null(load_cached_annotation('hg19_cpg_islands', 'hg19')),
        'Removing hg19_cpg_islands from the cache')
    expect_equal(nrow(list_cached_annotations()), 0)
})

test_that('clear_cached_annotations() removes by genome, or everything', {
    local_empty_cache()
    save_cached_annotation(cached_islands(), 'hg19_cpg_islands', 'hg19')
    save_cached_annotation(cached_islands(), 'mm10_cpg_islands', 'mm10')

    removed = suppressMessages(clear_cached_annotations(genome = 'mm10'))
    expect_equal(removed$name, 'mm10_cpg_islands')
    expect_equal(list_cached_annotations()$genome, 'hg19')

    expect_message(clear_cached_annotations(), 'Removed 1 item')
    expect_equal(nrow(list_cached_annotations()), 0)
})

test_that('Gene annotation cache names include the TxDb and org versions', {
    skip_if_not_installed('TxDb.Hsapiens.UCSC.hg19.knownGene')
    skip_if_not_installed('org.Hs.eg.db')

    rname = get_annotation_rname('hg19_genes_promoters', 'hg19')
    expect_match(rname, '^annotation\\|hg19_genes_promoters\\|annotatr ')
    expect_match(rname, 'TxDb.Hsapiens.UCSC.hg19.knownGene ', fixed = TRUE)
    expect_match(rname, 'org.Hs.eg.db ', fixed = TRUE)

    expect_match(get_annotation_rname('oviariramb2_genes_exons', 'oviariramb2'), 'AH119381', fixed = TRUE)
    expect_match(get_annotation_rname('hg38_mane_exons', 'hg38'), sprintf('MANE %s', MANE$version), fixed = TRUE)
    expect_match(get_annotation_rname('mm39_canonical_exons', 'mm39'), 'AH119358', fixed = TRUE)
    expect_no_match(get_annotation_rname('hg19_cpg_islands', 'hg19'), 'TxDb')
})

################################################################################
# Network tests

test_that('build_annotations() caches what it builds, and cache = FALSE bypasses it', {
    skip_network()
    local_empty_cache()

    a = suppressMessages(build_annotations(genome = 'hg38', annotations = 'hg38_cpgs'))
    cached = list_cached_annotations()
    expect_setequal(cached$name[cached$type == 'annotation'], expand_annotations('hg38_cpgs'))
    expect_equal(sum(cached$type == 'download'), 1)

    messages = capture_messages(a_cached <- build_annotations(genome = 'hg38', annotations = 'hg38_cpgs'))
    expect_length(grep('from the cache', messages), 4)
    expect_equal(a_cached, a)

    local_empty_cache()
    a_uncached = suppressMessages(build_annotations(genome = 'hg38', annotations = 'hg38_cpgs', cache = FALSE))
    expect_equal(a_uncached, a)
    expect_equal(nrow(list_cached_annotations()), 0)
})

test_that('A failed download is not cached', {
    skip_network()
    local_empty_cache()

    expect_error(
        suppressMessages(download_annotation_file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/doesNotExist.txt.gz',
            genome = 'hg38', retries = 2)),
        'Failed to download')
    expect_equal(nrow(list_cached_annotations()), 0)
})
