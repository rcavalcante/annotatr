# Shared fixtures for the tests. testthat sources helper-*.R files before the
# tests, so each test builds the data it needs instead of relying on objects
# created by other test files.

# Keep the tests out of the user's cache
options(annotatr.cache = file.path(tempdir(), 'annotatr-test-cache'))

# Give a test its own empty cache
local_empty_cache = function(env = parent.frame()) {
    withr::local_options(annotatr.cache = withr::local_tempdir(.local_envir = env), .local_envir = env)
}

extdata = function(file) {
    system.file('extdata', file, package = 'annotatr')
}

# Tests that download data are skipped offline. The heavy ones (full builds for
# every genome, AnnotationHub resources) run only when ANNOTATR_FULL_TESTS=true,
# which docker/check.sh sets.
skip_if_not_full_tests = function() {
    testthat::skip_if_not(
        identical(Sys.getenv('ANNOTATR_FULL_TESTS'), 'true'),
        'Set ANNOTATR_FULL_TESTS=true to run the full build tests')
}

# annotatr isn't on CRAN, and Bioconductor's builders don't set NOT_CRAN, so
# don't use testthat::skip_if_offline(), which also skips on CRAN. Otherwise
# the light network tests would never run on the builders.
skip_network = function() {
    testthat::skip_if_not(curl::has_internet(), 'No internet connection')
}

# Regions tested for differential methylation (DM) on chr9, with a DM_status
# column of hyper, hypo, or none
dm_regions = function(n = 1000) {
    extraCols = c(diff_meth = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
    r = suppressMessages(read_regions(
        con = extdata('IDH2mut_v_NBM_multi_data_chr9.txt.gz'),
        genome = 'hg19',
        extraCols = extraCols,
        rename_score = 'pval',
        rename_name = 'DM_status',
        format = 'bed'))

    return(r[seq_len(n)])
}

# dm_regions() annotated to the premade hg19 CpG annotations
annotate_dm_regions = function(regions = dm_regions()) {
    annotate_regions(
        regions = regions,
        annotations = annotatr::annotations,
        ignore.strand = TRUE,
        quiet = TRUE)
}

# ggplot2 does most of its work when a plot is built, so build it to catch errors
expect_plot_builds = function(plot) {
    testthat::expect_s3_class(plot, 'ggplot')
    built = NULL
    testthat::expect_no_error(built <- ggplot2::ggplot_build(plot))
    invisible(built)
}

cpg_types = c('hg19_cpg_islands', 'hg19_cpg_shores', 'hg19_cpg_shelves', 'hg19_cpg_inter')
