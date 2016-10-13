context('Test randomize module')

################################################################################
# Setup objects for plot_categorical()

    file = system.file('extdata', 'Gm12878_Stat3_chr2.bed.gz', package = 'annotatr')
    regions_genome = read_regions(con = file, genome = 'hg19', format = 'bed')
    regions_nogenome = read_regions(con = file, format = 'bed')

################################################################################
# Test errors

    test_that('Test errors', {
        expect_error(
            randomize_regions(regions = 'hello', allow.overlaps = TRUE, per.chromosome = TRUE),
            'regions must have class GRanges')
        expect_error(
            randomize_regions(regions = regions_nogenome),
            'GRanges object must have a valid genome'
            )
    })

################################################################################
# Test randomize_regions()

    test_that('Test randomized regions', {
        random_regions = randomize_regions(
            regions = regions_genome,
            allow.overlaps = TRUE,
            per.chromosome = TRUE)

        expect_equal(class(random_regions)[1], expected = 'GRanges')
        expect_equal(length(random_regions), expected = length(regions_genome))
    })
