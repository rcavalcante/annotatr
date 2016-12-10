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
    org = c('Dm','Dm','Hs','Hs','Mm','Mm','Rn','Rn','Rn'),
    stringsAsFactors = FALSE)

HMMCELLLINES = c('Gm12878','H1hesc','Hepg2','Hmec','Hsmm','Huvec','K562','Nhek','Nhlf')

HMMCODES = c('1_Active_Promoter', '2_Weak_Promoter' ,'3_Poised_Promoter' ,'4_Strong_Enhancer', '5_Strong_Enhancer', '6_Weak_Enhancer', '7_Weak_Enhancer', '8_Insulator', '9_Txn_Transition', '10_Txn_Elongation', '11_Weak_Txn', '12_Repressed', '13_Heterochrom/lo', '14_Repetitive/CNV')

#' Function to recode classes from chromHMM type column
#'
#' @return A character vector of chromHMM classes with numbers and underscores removed.
reformat_hmm_codes = function(hmm_codes) {
    new_codes = sapply(hmm_codes,
            function(hmm){paste(unlist(strsplit(hmm,'_'))[-1],collapse='')},
            USE.NAMES=FALSE)
    return(new_codes)
}

#' Function to return cell line from chromatin annotation shortcut
#'
#' @return A string of the cell line used in a chromatin annotation shortcut
get_cellline_from_shortcut = function(shortcut) {
    return(unlist(strsplit(unlist(strsplit(shortcut,'_'))[2], '-'))[1])
}

#' Function to return cell line from chromatin annotation code
#'
#' @return A string of the cell line used in a chromatin annotation code
get_cellline_from_code = function(code) {
    return(unlist(strsplit(unlist(strsplit(code,'_'))[3], '-'))[1])
}

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
    # Create annotation code endings
        shortcut_ends = c('basicgenes','cpgs')

        # Gene codes
        gene_ends = c('1to5kb', 'promoters', 'cds', '5UTRs', 'exons', 'firstexons', 'introns', 'intronexonboundaries', 'exonintronboundaries', '3UTRs', 'intergenic')

        # CpG codes
        cpg_ends = c('islands', 'shores', 'shelves', 'inter')

        # Chromatin state codes
        # Remove numbers, and underscores, and take unique
        chromatin_recode = unique(reformat_hmm_codes(HMMCODES))

        chromatin_ends = apply(
            expand.grid(HMMCELLLINES, chromatin_recode, stringsAsFactors = FALSE),
            1, paste, collapse='-')

        chromatin_shortcut_ends = apply(
            expand.grid(HMMCELLLINES, 'chromatin', stringsAsFactors = FALSE),
            1, paste, collapse='-')

    # Create full annotation codes
        gene_codes = apply(
            expand.grid(supported_genomes(), 'genes', gene_ends, stringsAsFactors = FALSE),
            1, paste, collapse='_')
        cpg_codes = apply(
            expand.grid(supported_genomes(), 'cpg', cpg_ends, stringsAsFactors= FALSE),
            1, paste, collapse='_')
        chromatin_codes = apply(
            expand.grid('hg19', 'chromatin', chromatin_ends, stringsAsFactors=FALSE),
            1, paste, collapse='_')

        enhancer_codes = c('hg19_enhancers_fantom','mm9_enhancers_fantom')
        lncrna_codes = c('hg19_lncrna_gencode','hg38_lncrna_gencode','mm10_lncrna_gencode')

        shortcut_codes = apply(
            expand.grid(supported_genomes(), shortcut_ends, stringsAsFactors = FALSE),
            1, paste, collapse='_')
        chromatin_shortcut_codes = paste('hg19', chromatin_shortcut_ends, sep='_')

    # Create the big vector of supported annotations
    annots = c(gene_codes, cpg_codes, chromatin_codes, enhancer_codes, lncrna_codes,
        shortcut_codes, chromatin_shortcut_codes)

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
#' @return A string giving the correct org for org.db packages. e.g. hg19 -> Hs.
get_orgdb_name = function(genome = supported_genomes()) {
    # Ensure valid arguments
    genome = match.arg(genome)

    org = ORGDBS[ORGDBS$genome == genome, 'org']

    return(org)
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
        } else if (tokens[2] == 'chromatin') {
            return(tokens[3])
        } else if (tokens[2] == 'custom') {
            return(tokens[3])
        } else if (tokens[2] == 'lncrna') {
            return('GENCODE lncRNA')
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
    are_hmms = any(grepl('chromatin', annotations))

    which_are_shortcuts = c(which(grepl('basicgenes', annotations)), which(grepl('cpgs', annotations)), which(grepl('-chromatin', annotations)))

    # expand_shortcuts() will always be run after check_annotations() so we can be
    # sure that the genome prefixes are the same for all annotaitons.
    genome = unique( sapply(annotations, function(a){ unlist(strsplit(a, '_'))[1] }, USE.NAMES = FALSE) )

    if(are_basicgenes || are_cpgs || are_hmms) {

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
        if(are_hmms) {
            # Could conceivably use shortcuts for multiple cell lines
            hmms = grep('-chromatin', annotations, value = TRUE)
            cell_lines = sapply(hmms, get_cellline_from_shortcut, USE.NAMES = FALSE)

            new_hmm_codes = apply(
                expand.grid(cell_lines, unique(reformat_hmm_codes(HMMCODES)), stringsAsFactors = FALSE),
                1, paste, collapse='-')

            new_annotations = c(new_annotations,
                paste(genome, 'chromatin', new_hmm_codes, sep='_'))
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
