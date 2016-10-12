#' A global-variable to hold custom annotations loaded in an R session
#'
#' Code thanks to Martin Morgan. This is a global variable that will store custom
#' annotations that a user reads in during a session in which annotatr is loaded.
#'
#' @return An environment to contain custom annotations from \code{read_annotations}.
#'
#' @examples
#'  # Example usage
#'  annotatr_cache$set("foo", 1:10)
#'  annotatr_cache$get("foo")
#'
#'  # Read in a BED3 file as a custom annotation
#'  file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
#'  # The custom annotation is added to the annotatr_cache environment in this function
#'  read_annotations(con = file, name = 'test', genome = 'hg19')
#'  # The result of read_annotations() is not visible in .GlobalEnv, instead
#'  # need to use the get method
#'  print(annotatr_cache$get('hg19_custom_test'))
#'
#' @export
annotatr_cache <- local({
    env = new.env(parent = emptyenv())
    list(set=function(key, value) {
        env[[key]] = value
        invisible(value)
    }, get = function(key) {
        if (!exists(key, env))
            stop(sQuote(key), " not in annotatr_cache", call.=FALSE)
        env[[key]]
    })
})

#' A function to build annotations from TxDb.* and AnnotationHub resources
#'
#' Create a \code{GRanges} object consisting of all the desired \code{annotations}. Supported annotation codes are listed by \code{supported_annotations()}. The basis for enhancer annotations are FANTOM5 data, the basis for CpG related annotations are CpG island tracks from \code{AnnotationHub}, and the basis for genic annotations are from the \code{TxDb.*} and \code{org.db} group of packages.
#'
#' @param genome The genome assembly.
#' @param annotations A character vector of annotations to build. Valid annotation codes are listed with \code{supported_annotations()}. The "basicgenes" shortcut builds the following regions: 1-5Kb upstream of TSSs, promoters, 5UTRs, exons, introns, and 3UTRs. The "cpgs" shortcut builds the following regions: CpG islands, shores, shelves, and interCGI regions. NOTE: Shortcuts need to be appended by the genome, e.g. \code{hg19_basicgenes}.
#' Custom annotations whose names are of the form \code{[genome]_custom_[name]} should also be included. Custom annotations should be read in and converted to \code{GRanges} with \code{read_annotations()}. They can be for a \code{supported_genome()}, or for an unsupported genome.
#'
#' @return A \code{GRanges} object of all the \code{annotations} combined. The \code{mcols} are \code{id, gene_id, symbol, type}.
#'
#' @examples
#' # Example with hg19
#' annots = c('hg19_cpg_islands','hg19_cpg_shores','hg19_genes_promoters')
#' annots_gr = build_annotations(genome = 'hg19', annotations = annots)
#'
#' # Example with a custom annotation
#' file = system.file('extdata', 'test_annotations_3.bed', package='annotatr')
#' read_annotations(con = file, name = 'test', genome = 'hg19')
#' annots = c('hg19_genes_promoters','hg19_custom_test')
#' annots_gr = build_annotations(genome = 'hg19', annotations = annots)
#'
#' @export
build_annotations = function(genome, annotations) {
    # Check annotations and expand any shortcuts
    check_annotations(annotations)
    annotations = expand_annotations(annotations)

    enh_annotations = grep('_enhancers_', annotations, value=TRUE)
    gene_annotations = grep('_genes_', annotations, value=TRUE)
    cpg_annotations = grep('_cpg_', annotations, value=TRUE)
    custom_annotations = grep('_custom_', annotations, value=TRUE)

    annots_grl = GenomicRanges::GRangesList()

    if(length(enh_annotations) != 0) {
        annots_grl = c(annots_grl, GenomicRanges::GRangesList(enhancers_fantom = suppressWarnings(build_enhancer_annots(genome = genome))))
    }
    if(length(gene_annotations) != 0) {
        annots_grl = c(annots_grl, suppressWarnings(build_gene_annots(genome = genome, annotations = gene_annotations)))
    }
    if(length(cpg_annotations) != 0) {
        annots_grl = c(annots_grl, suppressWarnings(build_cpg_annots(genome = genome, annotations = cpg_annotations)))
    }
    if(length(custom_annotations) > 0) {
        annots_grl = c(annots_grl, GenomicRanges::GRangesList(sapply(custom_annotations, function(ca){annotatr_cache$get(ca)})))
    }

    return(unlist(annots_grl, use.names=FALSE))
}

#' A helper function to build enhancer annotations for hg19 and mm10 from FANTOM5.
#'
#' @param genome The genome assembly.
#'
#' @return A \code{GRanges} object.
build_enhancer_annots = function(genome = supported_genomes()) {
    # Ensure valid arguments
    genome = match.arg(genome)

    message('Building enhancers...')

    # NOTE: Since something like hg38_enhancers_fantom will be caught by check_annotations() there is no need to worry about anything other than hg19 or mm10 will get to this point.

    # Get the enhancer annotations from FANTOM5
    if(genome == 'hg19') {
        enhancers = rtracklayer::import.bed('http://fantom.gsc.riken.jp/5/datafiles/phase2.0/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz', genome = 'hg19')
    } else if (genome == 'mm9') {
        enhancers = rtracklayer::import.bed('http://fantom.gsc.riken.jp/5/datafiles/phase2.0/extra/Enhancers/mouse_permissive_enhancers_phase_1_and_2.bed.gz', genome = 'mm9')
    }

    enhancers = GenomicRanges::granges(enhancers)
    GenomicRanges::mcols(enhancers)$id = paste0('enhancer:', seq_along(enhancers))
    GenomicRanges::mcols(enhancers)$gene_id = NA
    GenomicRanges::mcols(enhancers)$symbol = NA
    GenomicRanges::mcols(enhancers)$type = sprintf('%s_enhancers_fantom', genome)

    return(enhancers)
}

#' A helper function to build CpG related annotations.
#'
#' Using the \code{AnnotationHub} package, extract CpG island track for the appropriate \code{genome} and construct the shores, shelves, and interCGI annotations as desired.
#'
#' @param genome The genome assembly.
#' @param annotations A character vector with entries of the form \code{[genome]_cpg_{islands,shores,shelves,inter}}.
#'
#' @return A list of \code{GRanges} objects.
build_cpg_annots = function(genome = supported_genomes(), annotations = supported_annotations()) {
    # Ensure valid arguments
    genome = match.arg(genome)
    annotations = match.arg(annotations, several.ok = TRUE)

    annot_codes = data.frame(
        code = c(sprintf('%s_cpg_islands', genome),
            sprintf('%s_cpg_shores', genome),
            sprintf('%s_cpg_shelves', genome),
            sprintf('%s_cpg_inter', genome)),
        var = c('islands','shores','shelves','inter_cgi'),
        stringsAsFactors = FALSE)

    # Create AnnotationHub connection
    ah = AnnotationHub::AnnotationHub()

    # And do the query for available CpG Islands
    query = AnnotationHub::query(ah, c('CpG Islands'))

    if(!(genome %in% query$genome)) {
        stop(sprintf('AnnotationHub does not contain CpG Island annotations for %s', genome))
    }

    # Determine the correct ID to extract data from AnnotationHub
    ID = row.names(GenomicRanges::mcols(query)[GenomicRanges::mcols(query)$genome == genome, ])

    if(any(grepl('islands', annotations)) || any(grepl('shores', annotations)) || any(grepl('shelves', annotations)) || any(grepl('inter', annotations))) {
        message('Building CpG islands...')
        ### Islands
            # Extract and sort the islands
            islands = ah[[ID]]
            GenomicRanges::strand(islands) = '*'
            islands = GenomicRanges::sort(islands)

            # Rename the islands
            GenomicRanges::mcols(islands)$id = paste0('island:', seq_along(islands))

            # Add tx_name, gene_id, and symbol columns
            GenomicRanges::mcols(islands)$gene_id = NA
            GenomicRanges::mcols(islands)$symbol = NA

            # Give it the correct type
            GenomicRanges::mcols(islands)$type = sprintf('%s_cpg_islands', genome)

            GenomicRanges::mcols(islands) = GenomicRanges::mcols(islands)[, c('id','gene_id','symbol','type')]

        if(any(grepl('shores', annotations)) || any(grepl('shelves', annotations)) || any(grepl('inter', annotations))) {
            message('Building CpG shores...')
            ### Shores
                # Construct the shores based on:
                # upstream from the island start and downstream from the island end
                up_shores = GenomicRanges::flank(islands, width = 2000, start = TRUE, both = FALSE)
                down_shores = GenomicRanges::flank(islands, width = 2000, start = FALSE, both = FALSE)

                # Combine, sort, trim, and reduce combined up_shores and down_shores
                shores = c(up_shores, down_shores)
                shores = GenomicRanges::sort(shores)
                shores = GenomicRanges::trim(shores)
                shores = GenomicRanges::reduce(shores)

                # Remove islands from the shores
                shores = GenomicRanges::setdiff(shores, islands)

                # Rename the shores
                GenomicRanges::mcols(shores)$id = paste0('shore:', seq_along(shores))

                # Add tx_name, gene_id, and symbol columns
                GenomicRanges::mcols(shores)$gene_id = NA
                GenomicRanges::mcols(shores)$symbol = NA

                # Give it the correct type
                GenomicRanges::mcols(shores)$type = sprintf('%s_cpg_shores', genome)

                GenomicRanges::mcols(shores) = GenomicRanges::mcols(shores)[, c('id','gene_id','symbol','type')]

            if(any(grepl('shelves', annotations)) || any(grepl('inter', annotations))) {
                message('Building CpG shelves...')
                ### Shelves
                    # Construct the shelves based on:
                    # upstream from the up_shores start and downstream from the down_shores end
                    up_shelves = GenomicRanges::flank(shores, width = 2000, start = TRUE, both = FALSE)
                    down_shelves = GenomicRanges::flank(shores, width = 2000, start = FALSE, both = FALSE)

                    # Combine, sort, trim, and reduce combined up_shelves and down_shelves
                    shelves = c(up_shelves, down_shelves)
                    shelves = GenomicRanges::sort(shelves)
                    shelves = GenomicRanges::trim(shelves)
                    shelves = GenomicRanges::reduce(shelves)

                    # Remove islands and/or shores from the shelves
                    shelves = GenomicRanges::setdiff(shelves, islands)
                    shelves = GenomicRanges::setdiff(shelves, shores)

                    # Rename the shelves
                    GenomicRanges::mcols(shelves)$id = paste0('shelf:', seq_along(shelves))

                    # Add gene_id, and symbol columns
                    GenomicRanges::mcols(shelves)$gene_id = NA
                    GenomicRanges::mcols(shelves)$symbol = NA

                    # Give it the correct type
                    GenomicRanges::mcols(shelves)$type = sprintf('%s_cpg_shelves', genome)

                    GenomicRanges::mcols(shelves) = GenomicRanges::mcols(shelves)[, c('id','gene_id','symbol','type')]

                if(any(grepl('inter', annotations))) {
                    message('Building inter-CpG-islands...')
                    ### interCGI
                        # Take the union of all the objects so far
                        extended_cgi = Reduce(
                            function(x,y){
                                GenomicRanges::union(x,y)
                            },
                            list(islands, shores, shelves)
                        )

                        # A quirk of GenomicRanges::gaps() on unstranded ranges gives the whole
                        # chroms on the + and - strands in addition to the * gaps. Remove +/- gaps.
                        inter_cgi = GenomicRanges::gaps(extended_cgi)
                        inter_cgi = inter_cgi[GenomicRanges::strand(inter_cgi) == '*']

                        # Rename the interCGI
                        GenomicRanges::mcols(inter_cgi)$id = paste0('inter:', seq_along(inter_cgi))

                        # Add gene_id, and symbol columns
                        GenomicRanges::mcols(inter_cgi)$gene_id = NA
                        GenomicRanges::mcols(inter_cgi)$symbol = NA

                        # Give it the correct type
                        GenomicRanges::mcols(inter_cgi)$type = sprintf('%s_cpg_inter', genome)

                        GenomicRanges::mcols(inter_cgi) = GenomicRanges::mcols(inter_cgi)[, c('id','gene_id','symbol','type')]
                }
            }
        }
    }

    ### Put it all together
    mgets = annot_codes[annot_codes$code %in% annotations, 'var']
    cpgs = do.call('GRangesList', mget(mgets))
    names(cpgs) = annotations

    return(cpgs)
}

#' A helper function to build genic annotations.
#'
#' Using the \code{TxDb.*} group of packages, construct genic annotations consisting of any combination of 1-5kb upstream of a TSS, promoters (< 1kb from TSS), 5UTRs, CDS, exons, first exons, introns, intron/exon and exon/intron boundaries, 3UTRs, and intergenic.
#'
#' @param genome The genome assembly.
#' @param annotations A character vector with entries of the form \code{[genome]_genes_{1to5kb,promoters,5UTRs,cds,exons,firstexons,introns,intronexonboundaries,exonintronboundaries,3UTRs,intergenic}}.
#'
#' @return A list of \code{GRanges} objects.
build_gene_annots = function(genome = supported_genomes(), annotations = supported_annotations()) {
    # Ensure valid arguments
    genome = match.arg(genome)
    annotations = match.arg(annotations, several.ok = TRUE)

    annot_codes = data.frame(
        code = c(sprintf('%s_genes_promoters', genome),
            sprintf('%s_genes_1to5kb', genome),
            sprintf('%s_genes_cds', genome),
            sprintf('%s_genes_5UTRs', genome),
            sprintf('%s_genes_exons', genome),
            sprintf('%s_genes_firstexons', genome),
            sprintf('%s_genes_introns', genome),
            sprintf('%s_genes_intronexonboundaries', genome),
            sprintf('%s_genes_exonintronboundaries', genome),
            sprintf('%s_genes_3UTRs', genome),
            sprintf('%s_genes_intergenic', genome)),
        var = c('promoters_gr','onetofive_gr','cds_gr','fiveUTRs_gr','exons_gr',
            'firstexons_gr','introns_gr','intronexon_gr','exonintron_gr',
            'threeUTRs_gr','intergenic_gr'),
        stringsAsFactors = FALSE)

    # Load the appropriate TxDb.* library and get the txdb
    txdb_name = get_txdb_name(genome)
    library(txdb_name, character.only = TRUE)
    txdb = get(txdb_name)

    # Get the org.XX.eg.db mapping from Entrez ID to gene symbol
    # First element returned is package name, second is eg2SYMBOL name
    orgdb_name = get_orgdb_name(genome)
    library(orgdb_name[1], character.only = TRUE)
    x = get(orgdb_name[2])
    mapped_genes = mappedkeys(x)
    eg2symbol = as.data.frame(x[mapped_genes])

    # Build the base transcripts
    tx_gr = transcripts(txdb, columns = c('TXID','GENEID','TXNAME'))
    # Create TSS GRanges for later use with intronexon boundaries
    tss_gr = GenomicRanges::GRanges(
        seqnames = seqnames(tx_gr),
        ranges = IRanges(start = start(tx_gr), end = start(tx_gr)),
        strand = '*'
    )
    GenomicRanges::mcols(tss_gr) = GenomicRanges::mcols(tx_gr)
    seqinfo(tss_gr) = seqinfo(tx_gr)

    # Create TES GRanges for later use with exonintron boundaries
    tes_gr = GenomicRanges::GRanges(
        seqnames = seqnames(tx_gr),
        ranges = IRanges(start = end(tx_gr), end = end(tx_gr)),
        strand = '*'
    )
    GenomicRanges::mcols(tes_gr) = GenomicRanges::mcols(tx_gr)
    seqinfo(tes_gr) = seqinfo(tx_gr)

    # Build tables to map TXID to TXNAME and GENEID
    id_maps = AnnotationDbi::select(txdb, keys = as.character(GenomicRanges::mcols(tx_gr)$TXID), columns = c('TXNAME','GENEID'), keytype = 'TXID')

    # Each annotation should be a GRanges object with the following mcols:
    # tx_name, gene_id, symbol, type

    if(any(grepl('promoters', annotations)) || any(grepl('1to5kb', annotations)) || any(grepl('intergenic', annotations))) {
        message('Building promoters...')
        ### promoters
            promoters_gr = GenomicFeatures::promoters(txdb, upstream = 1000, downstream = 0)
            # Add Entrez ID, symbol, and type
            GenomicRanges::mcols(promoters_gr)$gene_id = id_maps[match(GenomicRanges::mcols(promoters_gr)$tx_id, id_maps$TXID), 'GENEID']
            GenomicRanges::mcols(promoters_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(promoters_gr)$gene_id, eg2symbol$gene_id), 'symbol']
            GenomicRanges::mcols(promoters_gr)$type = sprintf('%s_genes_promoters', genome)

            GenomicRanges::mcols(promoters_gr) = GenomicRanges::mcols(promoters_gr)[, c('tx_name','gene_id','symbol','type')]
            colnames(GenomicRanges::mcols(promoters_gr)) = c('id','gene_id','symbol','type')

        if(any(grepl('1to5kb', annotations)) || any(grepl('intergenic', annotations))) {
            message('Building 1to5kb upstream of TSS...')
            ### 1-5kb
                onetofive_gr = GenomicRanges::flank(promoters_gr, width = 4000, start = TRUE, both = FALSE)
                onetofive_gr = GenomicRanges::trim(onetofive_gr)
                # Add Entrez ID, symbol, and type (all but type are inherited from promoters_gr)
                GenomicRanges::mcols(onetofive_gr)$type = sprintf('%s_genes_1to5kb', genome)

            if(any(grepl('intergenic', annotations))) {
                message('Building intergenic...')
                ### intergenic
                    genic_gr = c(GenomicRanges::granges(tx_gr), GenomicRanges::granges(promoters_gr), GenomicRanges::granges(onetofive_gr))
                    GenomicRanges::strand(genic_gr) = '*'
                    intergenic_gr = GenomicRanges::reduce(genic_gr)
                    intergenic_gr = GenomicRanges::gaps(intergenic_gr)

                    # A quirk in gaps gives the entire + and - strand of a chromosome, ignore those
                    intergenic_gr = intergenic_gr[GenomicRanges::strand(intergenic_gr) == '*']

                    GenomicRanges::mcols(intergenic_gr)$id = paste0('intergenic:', seq_along(intergenic_gr))
                    GenomicRanges::mcols(intergenic_gr)$gene_id = NA
                    GenomicRanges::mcols(intergenic_gr)$symbol = NA
                    GenomicRanges::mcols(intergenic_gr)$type = sprintf('%s_genes_intergenic', genome)
            }
        }
    }

    if(any(grepl('cds', annotations))) {
        message('Building cds...')
        ### cds
            cds_grl = GenomicFeatures::cdsBy(txdb, by = 'tx', use.names = TRUE)
            # Create Rle of the tx_names
            cds_txname_rle = S4Vectors::Rle(names(cds_grl), S4Vectors::elementNROWS(cds_grl))
            cds_txname_vec = as.character(cds_txname_rle)
            # Unlist and add the tx_names
            cds_gr = unlist(cds_grl, use.names = FALSE)
            GenomicRanges::mcols(cds_gr)$tx_name = cds_txname_vec
            # Add Entrez ID, symbol, and type
            GenomicRanges::mcols(cds_gr)$gene_id = id_maps[match(GenomicRanges::mcols(cds_gr)$tx_name, id_maps$TXNAME), 'GENEID']
            GenomicRanges::mcols(cds_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(cds_gr)$gene_id, eg2symbol$gene_id), 'symbol']
            GenomicRanges::mcols(cds_gr)$type = sprintf('%s_genes_cds', genome)

            GenomicRanges::mcols(cds_gr) = GenomicRanges::mcols(cds_gr)[, c('tx_name','gene_id','symbol','type')]
            colnames(GenomicRanges::mcols(cds_gr)) = c('id','gene_id','symbol','type')
    }

    if(any(grepl('5UTR', annotations))) {
        message('Building 5UTRs...')
        ### fiveUTRs
            fiveUTRs_grl = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
            # Create Rle of the tx_names
            fiveUTRs_txname_rle = S4Vectors::Rle(names(fiveUTRs_grl), S4Vectors::elementNROWS(fiveUTRs_grl))
            fiveUTRs_txname_vec = as.character(fiveUTRs_txname_rle)
            # Unlist and add the tx_names
            fiveUTRs_gr = unlist(fiveUTRs_grl, use.names = FALSE)
            GenomicRanges::mcols(fiveUTRs_gr)$tx_name = fiveUTRs_txname_vec
            # Add Entrez ID, symbol, and type
            # NOTE: here we match on the tx_name because the tx_id is not given
            GenomicRanges::mcols(fiveUTRs_gr)$gene_id = id_maps[match(GenomicRanges::mcols(fiveUTRs_gr)$tx_name, id_maps$TXNAME), 'GENEID']
            GenomicRanges::mcols(fiveUTRs_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(fiveUTRs_gr)$gene_id, eg2symbol$gene_id), 'symbol']
            GenomicRanges::mcols(fiveUTRs_gr)$type = sprintf('%s_genes_5UTRs', genome)

            GenomicRanges::mcols(fiveUTRs_gr) = GenomicRanges::mcols(fiveUTRs_gr)[, c('tx_name','gene_id','symbol','type')]
            colnames(GenomicRanges::mcols(fiveUTRs_gr)) = c('id','gene_id','symbol','type')

    }

    if(any(grepl('3UTR', annotations))) {
        message('Building 3UTRs...')
        ### threeUTRs
            threeUTRs_grl = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
            # Create Rle of the tx_names
            threeUTRs_txname_rle = S4Vectors::Rle(names(threeUTRs_grl), S4Vectors::elementNROWS(threeUTRs_grl))
            threeUTRs_txname_vec = as.character(threeUTRs_txname_rle)
            # Unlist and add the tx_names
            threeUTRs_gr = unlist(threeUTRs_grl, use.names = FALSE)
            GenomicRanges::mcols(threeUTRs_gr)$tx_name = threeUTRs_txname_vec
            # NOTE: here we match on the tx_name because the tx_id is not given
            GenomicRanges::mcols(threeUTRs_gr)$gene_id = id_maps[match(GenomicRanges::mcols(threeUTRs_gr)$tx_name, id_maps$TXNAME), 'GENEID']
            GenomicRanges::mcols(threeUTRs_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(threeUTRs_gr)$gene_id, eg2symbol$gene_id), 'symbol']
            GenomicRanges::mcols(threeUTRs_gr)$type = sprintf('%s_genes_3UTRs', genome)

            GenomicRanges::mcols(threeUTRs_gr) = GenomicRanges::mcols(threeUTRs_gr)[, c('tx_name','gene_id','symbol','type')]
            colnames(GenomicRanges::mcols(threeUTRs_gr)) = c('id','gene_id','symbol','type')
    }

    if(any(grepl('exon', annotations)) || any(grepl('intron', annotations))) {

        message('Building exons...')
        ### exons
            exons_grl = GenomicFeatures::exonsBy(txdb, by = 'tx', use.names = TRUE)
            # Create Rle of the tx_names
            exons_txname_rle = S4Vectors::Rle(names(exons_grl), S4Vectors::elementNROWS(exons_grl))
            exons_txname_vec = as.character(exons_txname_rle)
            # Unlist and add the tx_names
            exons_gr = unlist(exons_grl, use.names = FALSE)
            GenomicRanges::mcols(exons_gr)$tx_name = exons_txname_vec
            # Add Entrez ID, symbol, and type
            GenomicRanges::mcols(exons_gr)$gene_id = id_maps[match(GenomicRanges::mcols(exons_gr)$tx_name, id_maps$TXNAME), 'GENEID']
            GenomicRanges::mcols(exons_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(exons_gr)$gene_id, eg2symbol$gene_id), 'symbol']
            GenomicRanges::mcols(exons_gr)$type = sprintf('%s_genes_exons', genome)

            # This needs to be here before we remove the exon_rank mcol in exons_gr
            if(any(grepl('firstexons', annotations))) {
                message('Building first exons...')
                ### first exons
                    firstexons_gr = exons_gr[sapply(exons_gr$exon_rank, function(er){1 %in% er})]
                    # NOTE: The mcol() contents are CharacterLists and IntegerLists, which requires a different approach from previous
                    GenomicRanges::mcols(firstexons_gr)$type = sprintf('%s_genes_firstexons', genome)

                    GenomicRanges::mcols(firstexons_gr) = GenomicRanges::mcols(firstexons_gr)[, c('tx_name','gene_id','symbol','type')]
                    colnames(GenomicRanges::mcols(firstexons_gr)) = c('id','gene_id','symbol','type')
            }

            GenomicRanges::mcols(exons_gr) = GenomicRanges::mcols(exons_gr)[, c('tx_name','gene_id','symbol','type')]
            colnames(GenomicRanges::mcols(exons_gr)) = c('id','gene_id','symbol','type')

        message('Building introns...')
        ### introns
            introns_grl = GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
            # Create Rle of the tx_names
            introns_txname_rle = S4Vectors::Rle(names(introns_grl), S4Vectors::elementNROWS(introns_grl))
            introns_txname_vec = as.character(introns_txname_rle)
            # Unlist and add the tx_names
            introns_gr = unlist(introns_grl, use.names = FALSE)
            GenomicRanges::mcols(introns_gr)$tx_name = introns_txname_vec
            # NOTE: here we match on the tx_name because the tx_id is not given
            GenomicRanges::mcols(introns_gr)$gene_id = id_maps[match(GenomicRanges::mcols(introns_gr)$tx_name, id_maps$TXNAME), 'GENEID']
            GenomicRanges::mcols(introns_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(introns_gr)$gene_id, eg2symbol$gene_id), 'symbol']
            GenomicRanges::mcols(introns_gr)$type = sprintf('%s_genes_introns', genome)

            GenomicRanges::mcols(introns_gr) = GenomicRanges::mcols(introns_gr)[, c('tx_name','gene_id','symbol','type')]
            colnames(GenomicRanges::mcols(introns_gr)) = c('id','gene_id','symbol','type')

        if(any(grepl('intronexonboundaries', annotations))) {
            message('Building intron exon boundaries...')
            ### Intron/exon boundary
                # Intron to exon transition will be the starts of exons_gr
                intronexon_gr = GenomicRanges::GRanges(
                    seqnames = seqnames(exons_gr),
                    ranges = IRanges(start = start(exons_gr), end = start(exons_gr)),
                    strand = GenomicRanges::strand(exons_gr))
                GenomicRanges::mcols(intronexon_gr) = GenomicRanges::mcols(exons_gr)
                seqinfo(intronexon_gr) = seqinfo(exons_gr)

                # Need to remove those boundaries that are TSSs
                tss_idx = unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(intronexon_gr, tss_gr)))
                intronexon_gr = intronexon_gr[-tss_idx]

                # Expand 200bp up and down
                intronexon_gr = GenomicRanges::flank(intronexon_gr, width = 200, both = TRUE)
                GenomicRanges::mcols(intronexon_gr)$type = sprintf('%s_genes_intronexonboundaries', genome)

                GenomicRanges::mcols(intronexon_gr) = GenomicRanges::mcols(intronexon_gr)[, c('id','gene_id','symbol','type')]
        }

        if(any(grepl('exonintronboundaries', annotations))) {
            message('Building exon intron boundaries...')
            ### Exon/intron boundary
                # Exon to intron transition will be the ends of exons_gr
                exonintron_gr = GenomicRanges::GRanges(
                    seqnames = seqnames(exons_gr),
                    ranges = IRanges(start = end(exons_gr), end = end(exons_gr)),
                    strand = GenomicRanges::strand(exons_gr))
                GenomicRanges::mcols(exonintron_gr) = GenomicRanges::mcols(exons_gr)
                seqinfo(exonintron_gr) = seqinfo(exons_gr)

                # Need to remove those boundaries that are TESs
                tes_idx = unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(exonintron_gr, tes_gr)))
                exonintron_gr = exonintron_gr[-tes_idx]

                # Expand 200bp up and down
                exonintron_gr = GenomicRanges::flank(exonintron_gr, width = 200, both = TRUE)
                GenomicRanges::mcols(exonintron_gr)$type = sprintf('%s_genes_exonintronboundaries', genome)

                GenomicRanges::mcols(exonintron_gr) = GenomicRanges::mcols(exonintron_gr)[, c('id','gene_id','symbol','type')]
        }
    }

    ### Put it all together
    mgets = annot_codes[annot_codes$code %in% annotations, 'var']
    genes = do.call('GRangesList', mget(mgets))
    names(genes) = annotations

    return(genes)
}
