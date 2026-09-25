#' Build gene annotations from any TxDb or EnsDb
#'
#' Build the same gene annotations as \code{build_annotations()} (promoters, exons, introns, etc.) from gene models you provide, e.g. GENCODE or RefSeq gene models, or a genome without built-in annotations. The result can be combined with built-in annotations with \code{c()} and used with \code{annotate_regions()}, \code{summarize_genes()}, and the plotting functions.
#'
#' To build from a GTF or GFF3 file, first make a \code{TxDb} with \code{txdbmaker::makeTxDbFromGFF()}. An \code{EnsDb} can come from \code{AnnotationHub} or \code{ensembldb::ensDbFromGtf()}.
#'
#' The annotation codes are \code{[genome]_[group]_[type]}, e.g. \code{mm39_gencode_promoters}, so annotations from different gene models can be used together. The chromosome names are those of the \code{txdb}. An \code{EnsDb} has Ensembl names (e.g. \code{1} rather than \code{chr1}), which can be changed with \code{GenomeInfoDb::seqlevelsStyle(edb) = 'UCSC'} before building, to match regions with UCSC names.
#'
#' Annotations built with this function aren't cached on disk. See \code{\link{cached-annotations}}.
#'
#' @param txdb A \code{TxDb} or \code{EnsDb} object.
#' @param genome A string giving the genome assembly of \code{txdb}, e.g. \code{'mm39'}. If \code{txdb} records a genome, they must agree.
#' @param group A string naming the source of the gene models, e.g. \code{'gencode'} or \code{'refseq'}. It must be letters and numbers only, and can't be a built-in group (\code{genes}, \code{mane}, \code{cpg}, \code{enhancers}, \code{chromatin}, \code{lncrna}, \code{custom}).
#' @param annotations A character vector of annotation types to build: \code{1to5kb}, \code{promoters}, \code{5UTRs}, \code{cds}, \code{exons}, \code{firstexons}, \code{introns}, \code{intronexonboundaries}, \code{exonintronboundaries}, \code{3UTRs}, and \code{intergenic}. The \code{basicgenes} shortcut (the default) builds 1-5kb upstream of TSSs, promoters, 5UTRs, exons, introns, and 3UTRs.
#' @param orgdb An optional \code{OrgDb} object, e.g. \code{org.Mm.eg.db::org.Mm.eg.db}, to get gene symbols, Entrez IDs, and Ensembl IDs for the gene IDs of a \code{TxDb}. An \code{EnsDb} has them already.
#' @param keytype A string giving the type of the gene IDs in \code{txdb}, as a \code{keytype} of \code{orgdb}, e.g. \code{'ENTREZID'} (the default) or \code{'ENSEMBL'}. For \code{'ENSEMBL'}, versions of gene IDs (e.g. the \code{.12} in \code{ENSMUSG00000000001.12}) are ignored when looking up symbols.
#'
#' @return A \code{GRanges} object with \code{mcols} \code{id}, \code{tx_id}, \code{gene_id}, \code{symbol}, \code{entrez_id}, \code{ensembl_id}, and \code{type}, as from \code{build_annotations()}. The \code{tx_id} and \code{gene_id} are the transcript names and gene IDs of \code{txdb}. The \code{symbol}, \code{entrez_id}, and \code{ensembl_id} come from an \code{EnsDb}, or from \code{orgdb} (the first, when a gene ID maps to more than one), and are \code{NA} without either.
#'
#' @seealso \code{\link{build_annotations}}, \code{\link{read_annotations}} for annotations from BED files.
#'
#' @examples
#'  if(requireNamespace('TxDb.Hsapiens.UCSC.hg19.knownGene', quietly = TRUE) &&
#'      requireNamespace('org.Hs.eg.db', quietly = TRUE)) {
#'    # Any TxDb works, e.g. one made from a GTF with txdbmaker::makeTxDbFromGFF()
#'    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#'    annots = build_txdb_annotations(txdb, genome = 'hg19', group = 'knowngene',
#'        annotations = 'promoters', orgdb = org.Hs.eg.db::org.Hs.eg.db)
#'
#'    unique(annots$type)
#'  }
#'
#' @export
build_txdb_annotations = function(txdb, genome, group, annotations = 'basicgenes', orgdb = NULL, keytype = 'ENTREZID') {
    if(!(methods::is(txdb, 'TxDb') || methods::is(txdb, 'EnsDb'))) {
        stop('Error: txdb must be a TxDb or EnsDb object.')
    }

    # Check the genome, and that it agrees with the txdb
    if(missing(genome) || !is.character(genome) || length(genome) != 1 || is.na(genome) || !grepl('^[A-Za-z0-9.]+$', genome)) {
        stop("Error: genome must be a single genome assembly name without underscores, e.g. 'mm39'.")
    }
    txdb_genome = unique(stats::na.omit(Seqinfo::genome(Seqinfo::seqinfo(txdb))))
    if(length(txdb_genome) > 0 && !(genome %in% txdb_genome)) {
        stop(sprintf('Error: genome is %s, but txdb is from genome %s.', genome, paste(txdb_genome, collapse = ', ')))
    }

    # Check the group
    if(missing(group) || !is.character(group) || length(group) != 1 || is.na(group) || !grepl('^[A-Za-z0-9]+$', group)) {
        stop("Error: group must be a single name of letters and numbers, e.g. 'gencode'.")
    }
    if(tolower(group) %in% BUILTIN_GROUPS) {
        stop(sprintf('Error: group can\'t be a built-in group: %s.', paste(BUILTIN_GROUPS, collapse = ', ')))
    }

    # Expand the shortcut, and check the annotation types
    if('basicgenes' %in% annotations) {
        annotations = union(setdiff(annotations, 'basicgenes'), BASIC_GENE_TYPES)
    }
    unsupported = setdiff(annotations, c(GENE_TYPES, 'intergenic'))
    if(length(unsupported) > 0) {
        stop(sprintf('Error: "%s" is(are) not gene annotation types. Use one or more of %s, or basicgenes.',
            paste(unsupported, collapse = ', '), paste(c(GENE_TYPES, 'intergenic'), collapse = ', ')))
    }

    # Map gene IDs to symbols, Entrez IDs, and Ensembl IDs. An EnsDb has them.
    gene_table = NULL
    if(!is.null(orgdb)) {
        if(!methods::is(orgdb, 'OrgDb')) {
            stop('Error: orgdb must be an OrgDb object, e.g. org.Mm.eg.db::org.Mm.eg.db.')
        }
        if(!(keytype %in% AnnotationDbi::keytypes(orgdb))) {
            stop(sprintf('Error: keytype %s is not a keytype of orgdb. See AnnotationDbi::keytypes(orgdb).', keytype))
        }
        if(methods::is(txdb, 'EnsDb')) {
            message('Using the gene names and Entrez IDs in the EnsDb, not orgdb.')
        } else {
            gene_table = get_gene_table(AnnotationDbi::keys(txdb, keytype = 'GENEID'), orgdb = orgdb, keytype = keytype)
        }
    }

    prefix = sprintf('%s_%s', genome, group)
    # Promoters and flanks can extend past the ends of chromosomes before
    # they are trimmed, which warns
    genes = withCallingHandlers(
        build_txdb_gene_annots(txdb, prefix = prefix, annotations = sprintf('%s_%s', prefix, annotations), gene_table = gene_table),
        warning = function(w) {
            if(grepl('out-of-bound range', conditionMessage(w))) {
                invokeRestart('muffleWarning')
            }
        })

    gr = unlist(genes, use.names = FALSE)
    names(gr) = NULL

    # Label the genome, e.g. for a TxDb made from a GTF without one
    if(anyNA(Seqinfo::genome(gr))) {
        Seqinfo::genome(gr) = genome
    }

    return(gr)
}
