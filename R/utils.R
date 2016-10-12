### Constants
# TxDb.* family of packages
TXDBS = c(
    'TxDb.Dmelanogaster.UCSC.dm3.ensGene',
    'TxDb.Dmelanogaster.UCSC.dm6.ensGene',
    'TxDb.Hsapiens.UCSC.hg19.knownGene',
    'TxDb.Hsapiens.UCSC.hg38.knownGene',
    'TxDb.Mmusculus.UCSC.mm9.knownGene',
    'TxDb.Mmusculus.UCSC.mm10.knownGene',
    'TxDb.Rnorvegicus.UCSC.rn4.ensGene',
    'TxDb.Rnorvegicus.UCSC.rn5.refGene',
    'TxDb.Rnorvegicus.UCSC.rn6.refGene')

# org.* family of packages
ORGDBS = data.frame(
    genome = c('dm3','dm6','hg19','hg38','mm9','mm10','rn4','rn5','rn6'),
    orgdb = c('org.Dm.eg.db', 'org.Dm.eg.db', 'org.Hs.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Rn.eg.db', 'org.Rn.eg.db'),
    tosymbol = c('org.Dm.egSYMBOL', 'org.Dm.egSYMBOL', 'org.Hs.egSYMBOL', 'org.Hs.egSYMBOL', 'org.Mm.egSYMBOL', 'org.Mm.egSYMBOL', 'org.Rn.egSYMBOL', 'org.Rn.egSYMBOL', 'org.Rn.egSYMBOL'),
    chrlengths = c('org.Dm.egCHRLENGTHS', 'org.Dm.egCHRLENGTHS', 'org.Hs.egCHRLENGTHS', 'org.Hs.egCHRLENGTHS', 'org.Mm.egCHRLENGTHS', 'org.Mm.egCHRLENGTHS', 'org.Rn.egCHRLENGTHS', 'org.Rn.egCHRLENGTHS', 'org.Rn.egCHRLENGTHS'),
    stringsAsFactors = FALSE)

#' Function listing which annotations are available.
#'
#' This includes the shortcuts. The \code{expand_annotations()} function helps
#' handle the shortcuts.
#'
#' @return A character vector of available annotations.
#'
#' @examples
#' supported_annotations()
#'
#' @export
supported_annotations = function() {
    shortcuts = c('basicgenes','cpgs')

    genes = c('1to5kb', 'promoters', 'cds', '5UTRs', 'exons', 'firstexons', 'introns', 'intronexonboundaries', 'exonintronboundaries', '3UTRs', 'intergenic')
    cpgs = c('islands', 'shores', 'shelves', 'inter')

    combos = rbind(
        expand.grid(supported_genomes(), 'genes', genes, stringsAsFactors = FALSE),
        expand.grid(supported_genomes(), 'cpg', cpgs, stringsAsFactors= FALSE)
    )

    annots = as.character(apply(combos, 1, paste, collapse='_'))
    annots = c(annots, apply(expand.grid(supported_genomes(), shortcuts, stringsAsFactors = FALSE), 1, paste, collapse='_'))
    annots = c(annots, c('hg19_enhancers_fantom','mm9_enhancers_fantom'))
    return(annots)
}

#' Function returning supported TxDb.* genomes
#'
#' @return A character vector of genomes for supported TxDb.* packages
#'
#' @examples
#' supported_genomes()
#'
#' @export
supported_genomes = function() {
    return(c('dm3','dm6','hg19','hg38','mm9','mm10','rn4','rn5','rn6'))
}

#' Function to get correct TxDb.* package name based on genome
#'
#' @param genome A string giving the genome assembly.
#'
#' @return A string giving the name of the correct TxDb.* package name based on \code{genome}.
get_txdb_name = function(genome = supported_genomes()) {
    # Ensure valid arguments
    genome = match.arg(genome)

    db = grep(genome, TXDBS, value = TRUE)

    return(db)
}

#' Function to get correct org.* package name based on genome
#'
#' @param genome A string giving the genome assembly.
#'
#' @return A string giving the name of the correct TxDb.* package name based on \code{genome}.
get_orgdb_name = function(genome = supported_genomes()) {
    # Ensure valid arguments
    genome = match.arg(genome)

    db = ORGDBS[ORGDBS$genome == genome, 'orgdb']
    tosymbol = ORGDBS[ORGDBS$genome == genome, 'tosymbol']
    chrlengths = ORGDBS[ORGDBS$genome == genome, 'chrlengths']

    return(c(db, tosymbol, chrlengths))
}

#' Function to tidy up annotation accessors for visualization
#'
#' @param annotations A character vector of annotations, in the order they are to appear in the visualization.
#'
#' @return A list of mappings from original annotation names to names ready for visualization.
tidy_annotations = function(annotations) {
    tidy = sapply(annotations, function(a){
        tokens = unlist(strsplit(a,'_'))
        if(tokens[2] == 'cpg') {
            if(tokens[3] == 'inter') {
                return('interCGI')
            } else {
                return(paste('CpG', tokens[3]))
            }
        } else if (tokens[2] == 'genes') {
            if(tokens[3] == 'firstexons') {
                return('first exons')
            } else if (tokens[3] == 'intronexonboundaries') {
                return('intron/exon boundaries')
            } else if (tokens[3] == 'exonintronboundaries') {
                return('exon/intron boundaries')
            } else {
                return(tokens[3])
            }
        } else if (tokens[2] == 'enhancers') {
            return('enhancers')
        } else if (tokens[2] == 'custom') {
            return(tokens[3])
        }
    })

    flip_tidy = names(tidy)
    names(flip_tidy) = tidy

    return(as.list(flip_tidy))
}

#' Function to check for valid annotations
#'
#' Gives errors if any annotations are not in supported_annotations() (and they are not in the required custom format), basicgenes are used, or the genome prefixes are not the same for all annotations.
#'
#' @param annotations A character vector of annotations possibly using the shortcuts
#' @return If all the checks on the annotations pass, returns NULL to allow code to move forward.
check_annotations = function(annotations) {
    # Pull out any custom annotations before checking
    custom_annotations = grep('custom', annotations, value = TRUE)
    annotations = base::setdiff(annotations, custom_annotations)

    # Check that the annotations are supported, tell the user which are unsupported
    if( !all(annotations %in% supported_annotations()) ) {
        unsupported = base::setdiff(annotations, supported_annotations())

        stop(sprintf('Error: "%s" is(are) not supported. See supported_annotations().',
            paste(unsupported, collapse=', ')))
    }

    # Recombine annotations and custom_annotations or you get failure when
    # there are only custom annotations
    annotations = c(custom_annotations, annotations)

    genomes = sapply(annotations, function(a){
        unlist(strsplit(a, '_'))[1]
    }, USE.NAMES = FALSE)

    # Check for same genome on all annotations
    if( length(unique(genomes)) != 1 ){
        stop('Error: genome prefix on all annotations must be the same.')
    }

    return(NULL)
}

#' Function to expand annotation shortcuts
#'
#' @param annotations A character vector of annotations, possibly using the shortcut accessors
#'
#' @return A vector of data accession-ized names that are ordered from upstream to downstream in the case of knownGenes and islands to interCGI in the case of cpgs.
expand_annotations = function(annotations) {
    are_basicgenes = any(grepl('basicgenes', annotations))
    are_cpgs = any(grepl('cpgs', annotations))
    which_are_shortcuts = c(which(grepl('basicgenes', annotations)), which(grepl('cpgs', annotations)))

    # expand_shortcuts() will always be run after check_annotations() so we can be
    # sure that the genome prefixes are the same for all annotaitons.
    genome = unique( sapply(annotations, function(a){ unlist(strsplit(a, '_'))[1] }, USE.NAMES = FALSE) )

    if(are_basicgenes || are_cpgs) {

        # Check for shortcut annotation accessors 'cpgs', 'basicgenes'
        # and create the right annotations based on the genome
        new_annotations = c()
        remove_shortcuts = c()
        if(are_cpgs) {
            new_annotations = paste(genome, 'cpg', c('islands','shores','shelves','inter'), sep='_')
        }
        if(are_basicgenes) {
            new_annotations = c(new_annotations, paste(genome, 'genes', c('1to5kb','promoters','5UTRs','exons','introns','3UTRs'), sep='_'))
        }
        annotations = base::setdiff(c(annotations, new_annotations), annotations[which_are_shortcuts])
    }

    return(annotations)
}

#' Function to subset a tbl_df or grouped_df by a column
#'
#' @param tbl A \code{tbl_df} or \code{grouped_df}.
#' @param col A string indicating which column of of \code{tbl} to subset and order
#' @param col_order A character vector indicating the order of \code{col}.
#'
#' @return A modified version of \code{summary} with \code{col} subsetted by \code{col_order}.
subset_order_tbl = function(tbl, col, col_order) {
    if(!is.null(col)) {
        # Collect all types in the column
        all_col_names = unique(tbl[[col]])

        # Inherit col_order from the order in tbl
        if(is.null(col_order)) {
            col_order = all_col_names
        }

        # Check set equality of col in the summary and the col_order
        if( !dplyr::setequal(all_col_names, col_order) ) {
            if( all(col_order %in% all_col_names) ) {
                tbl = subset(tbl, tbl[[col]] %in% col_order)
            } else {
        # Intersect col_order with unique(tbl[[col]]) to deal with possible 0 tallies
        col_order = intersect(col_order, unique(tbl[[col]]))
                warning('There are elements in col_order that are not present in the corresponding column. Check for typos, or this could be a result of 0 tallies.')
            }
        }

        # Convert fill to factor with levels in the correct order
        tbl[[col]] = factor(tbl[[col]], levels = col_order)
        # Also convert the levels to tidy names if fill is annotations
        if(col == 'annot.type') {
            levels(tbl[[col]]) = tidy_annotations(col_order)
        }
    }
    return(tbl)
}
