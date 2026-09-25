### Constants
# TxDb.* family of packages
TXDBS = c(
    'TxDb.Dmelanogaster.UCSC.dm3.ensGene',
    'TxDb.Dmelanogaster.UCSC.dm6.ensGene',
    'TxDb.Drerio.UCSC.danRer10.refGene',
    'TxDb.Drerio.UCSC.danRer11.refGene',
    'TxDb.Ggallus.UCSC.galGal5.refGene',
    'TxDb.Hsapiens.UCSC.hg19.knownGene',
    'TxDb.Hsapiens.UCSC.hg38.knownGene',
    'TxDb.Mmusculus.UCSC.mm9.knownGene',
    'TxDb.Mmusculus.UCSC.mm10.knownGene',
    'TxDb.Mmusculus.UCSC.mm39.knownGene',
    'TxDb.Rnorvegicus.UCSC.rn4.ensGene',
    'TxDb.Rnorvegicus.UCSC.rn5.refGene',
    'TxDb.Rnorvegicus.UCSC.rn6.refGene',
    'TxDb.Rnorvegicus.UCSC.rn7.refGene')

# org.* family of packages
ORGDBS = data.frame(
    genome = c('dm3','dm6','danRer10','danRer11','galGal5','hg19','hg38','mm9','mm10','mm39','oviariramb2','rn4','rn5','rn6','rn7'),
    org = c('Dm','Dm','Dr','Dr','Gg','Hs','Hs','Mm','Mm','Mm',NA,'Rn','Rn','Rn','Rn'),
    stringsAsFactors = FALSE)

# Genomes without TxDb.* / org.* packages. Gene models come from an AnnotationHub
# EnsDb, and CpG islands, chromosome sizes, and chromosome name aliases come from
# the UCSC GenArk assembly hub. Sequences are renamed to the UCSC-style names.
GENARK = data.frame(
    genome = c('oviariramb2'),
    accession = c('GCF_016772045.1'),
    assembly = c('ARS-UI_Ramb_v2.0'),
    ensdb = c('AH119381'),
    stringsAsFactors = FALSE)

# MANE transcripts (hg38 only), for the hg38_mane_* annotations. The release
# directory is versioned, so a new release needs a change here.
MANE = list(
    version = '1.5',
    url = 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5')

# Gene annotation types in genomic order, for summarize_genes()
GENE_TYPES = c('1to5kb', 'promoters', '5UTRs', 'cds', 'firstexons', 'exons', 'intronexonboundaries', 'introns', 'exonintronboundaries', '3UTRs')

# Gene annotation types of the basicgenes and basicmane shortcuts
BASIC_GENE_TYPES = c('1to5kb', 'promoters', '5UTRs', 'exons', 'introns', '3UTRs')

# The groups of builtin annotation codes, [genome]_[group]_[type]. Gene
# annotations from build_txdb_annotations() use other group names.
BUILTIN_GROUPS = c('genes', 'mane', 'canonical', 'cpg', 'enhancers', 'chromatin', 'lncrna', 'ccre', 'custom')

# ENCODE candidate cis-regulatory elements (cCREs), from the SCREEN registry,
# for the [genome]_ccre_* annotations. The registry directory is versioned, so
# a new version needs a change here.
CCRE = list(
    version = 'V4',
    url = 'https://downloads.wenglab.org/Registry-V4',
    files = c(hg38 = 'GRCh38-cCREs.bed', mm10 = 'mm10-cCREs.bed'))

# cCRE classes, and their readable names for tidy_annotations()
CCRE_CLASSES = c(
    'PLS' = 'promoter-like',
    'pELS' = 'proximal enhancer-like',
    'dELS' = 'distal enhancer-like',
    'CA-H3K4me3' = 'accessible + H3K4me3',
    'CA-CTCF' = 'accessible + CTCF',
    'CA-TF' = 'accessible + TF',
    'CA' = 'accessible only',
    'TF' = 'TF only')

# The builtin groups of gene annotations
GENE_GROUPS = c('genes', 'mane', 'canonical')

# Ensembl canonical transcripts, one per gene, for the hg38_canonical_*
# annotations and the like, from the Ensembl 113 EnsDbs in AnnotationHub
CANONICAL = data.frame(
    genome = c('hg38', 'mm39', 'rn7', 'danRer11', 'dm6', 'oviariramb2'),
    ensdb = c('AH119325', 'AH119358', 'AH119437', 'AH119289', 'AH119285', 'AH119381'),
    stringsAsFactors = FALSE)

# The mcols of every annotation, in order
ANNOTATION_MCOLS = c('id', 'tx_id', 'gene_id', 'symbol', 'entrez_id', 'ensembl_id', 'type')

HMMCELLLINES = c('Gm12878','H1hesc','Hepg2','Hmec','Hsmm','Huvec','K562','Nhek','Nhlf')

HMMCODES = c('1_Active_Promoter', '2_Weak_Promoter' ,'3_Poised_Promoter' ,'4_Strong_Enhancer', '5_Strong_Enhancer', '6_Weak_Enhancer', '7_Weak_Enhancer', '8_Insulator', '9_Txn_Transition', '10_Txn_Elongation', '11_Weak_Txn', '12_Repressed', '13_Heterochrom/lo', '14_Repetitive/CNV')

#' Function to recode classes from chromHMM type column
#'
#' @param hmm_codes in the original form from UCSC Genome Browser track.
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
#' @param shortcut The annotation shortcut, used in \code{build_annotations()}.
#'
#' @return A string of the cell line used in a chromatin annotation shortcut
get_cellline_from_shortcut = function(shortcut) {
    return(unlist(strsplit(unlist(strsplit(shortcut,'_'))[2], '-'))[1])
}

#' Function to return cell line from chromatin annotation code
#'
#' @param code The annotation code, used in \code{build_annotations()}.
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
#' builtin_annotations()
#'
#' @export
builtin_annotations = function() {
    # Create annotation code endings
        shortcut_ends = c('basicgenes','basicmane','basiccanonical','cpgs','ccres')

        # Gene codes
        gene_genomes = annotatr::builtin_genomes()
        gene_ends = c('1to5kb', 'promoters', 'cds', '5UTRs', 'exons', 'firstexons', 'introns', 'intronexonboundaries', 'exonintronboundaries', '3UTRs', 'intergenic')

        # MANE codes (hg38 only), the gene codes without intergenic
        mane_ends = setdiff(gene_ends, 'intergenic')

        # CpG codes
        cpg_genomes = base::setdiff(annotatr::builtin_genomes(),c('dm3','dm6'))
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
            expand.grid(gene_genomes, 'genes', gene_ends, stringsAsFactors = FALSE),
            1, paste, collapse='_')
        mane_codes = paste('hg38', 'mane', mane_ends, sep='_')
        canonical_codes = apply(
            expand.grid(CANONICAL$genome, 'canonical', gene_ends, stringsAsFactors = FALSE),
            1, paste, collapse='_')
        ccre_codes = apply(
            expand.grid(names(CCRE$files), 'ccre', names(CCRE_CLASSES), stringsAsFactors = FALSE),
            1, paste, collapse='_')
        cpg_codes = apply(
            expand.grid(cpg_genomes, 'cpg', cpg_ends, stringsAsFactors= FALSE),
            1, paste, collapse='_')
        chromatin_codes = apply(
            expand.grid('hg19', 'chromatin', chromatin_ends, stringsAsFactors=FALSE),
            1, paste, collapse='_')

        enhancer_codes = c('hg19_enhancers_fantom','hg38_enhancers_fantom','mm9_enhancers_fantom','mm10_enhancers_fantom')
        lncrna_codes = c('hg19_lncrna_gencode','hg38_lncrna_gencode','mm10_lncrna_gencode')

        gene_shortcut_codes = apply(
            expand.grid(gene_genomes, 'basicgenes', stringsAsFactors = FALSE),
            1, paste, collapse='_')
        mane_shortcut_codes = 'hg38_basicmane'
        canonical_shortcut_codes = paste(CANONICAL$genome, 'basiccanonical', sep='_')
        ccre_shortcut_codes = paste(names(CCRE$files), 'ccres', sep='_')
        cpg_shortcut_codes = apply(
            expand.grid(cpg_genomes, 'cpgs', stringsAsFactors = FALSE),
            1, paste, collapse='_')
        chromatin_shortcut_codes = paste('hg19', chromatin_shortcut_ends, sep='_')

    # Create the big vector of supported annotations
    annots = c(gene_codes, mane_codes, canonical_codes, cpg_codes, chromatin_codes, enhancer_codes, lncrna_codes, ccre_codes,
        gene_shortcut_codes, mane_shortcut_codes, canonical_shortcut_codes, cpg_shortcut_codes, chromatin_shortcut_codes, ccre_shortcut_codes)

    return(annots)
}

#' Function returning supported TxDb.* genomes
#'
#' @return A character vector of genomes for supported TxDb.* packages
#'
#' @examples
#' builtin_genomes()
#'
#' @export
builtin_genomes = function() {
    return(ORGDBS$genome)
}

#' Function to get correct TxDb.* package name based on genome
#'
#' @param genome A string giving the genome assembly.
#'
#' @return A string giving the name of the correct TxDb.* package name based on \code{genome}.
get_txdb_name = function(genome = annotatr::builtin_genomes()) {
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
get_orgdb_name = function(genome = annotatr::builtin_genomes()) {
    # Ensure valid arguments
    genome = match.arg(genome)

    org = ORGDBS[ORGDBS$genome == genome, 'org']

    return(org)
}

#' Function to give a GRanges the seqinfo of a genome
#'
#' Some operations, e.g. \code{rtracklayer::liftOver()}, drop the genome and sequence lengths. This restores them from \code{Seqinfo::Seqinfo(genome = genome)}. If that fails, e.g. offline, or if the ranges are on sequences the genome doesn't have, only the genome is set.
#'
#' @param gr A \code{GRanges} object.
#' @param genome A string giving the genome assembly, e.g. \code{'hg38'}.
#'
#' @return \code{gr} with the seqinfo of \code{genome}, trimmed to the ends of its sequences.
set_genome_seqinfo = function(gr, genome) {
    seqinfo = tryCatch(Seqinfo::Seqinfo(genome = genome), error = function(e) NULL)

    if(!is.null(seqinfo) && all(Seqinfo::seqlevelsInUse(gr) %in% Seqinfo::seqlevels(seqinfo))) {
        Seqinfo::seqlevels(gr, pruning.mode = 'coarse') = Seqinfo::seqlevels(seqinfo)
        # Out-of-bound ranges warn here, and are trimmed next
        suppressWarnings(Seqinfo::seqinfo(gr) <- seqinfo)
        gr = GenomicRanges::trim(gr)
    } else {
        Seqinfo::genome(gr) = genome
    }

    return(gr)
}

#' Function to get the URL of a file in the UCSC GenArk hub for a genome
#'
#' @param genome A string giving the genome assembly, one of \code{GENARK$genome}.
#' @param file A string giving the path of the file within the hub directory, e.g. \code{'GCF_016772045.1.chromAlias.txt'}.
#'
#' @return A string giving the URL.
get_genark_url = function(genome, file) {
    accession = GENARK[GENARK$genome == genome, 'accession']

    # GCF_016772045.1 lives at GCF/016/772/045/GCF_016772045.1
    digits = substr(accession, 5, 13)
    hub_dir = paste(substr(accession, 1, 3), substr(digits, 1, 3), substr(digits, 4, 6), substr(digits, 7, 9), accession, sep = '/')

    return(sprintf('https://hgdownload.soe.ucsc.edu/hubs/%s/%s', hub_dir, file))
}

#' Function to map any chromosome alias of a genome to its UCSC-style name
#'
#' Uses the chromAlias file of the UCSC GenArk hub (for \code{GENARK} genomes) or of the UCSC database (\code{goldenPath/[genome]/bigZips}), e.g. to rename the Ensembl sequence names of an \code{EnsDb} (\code{1}, \code{MT}, \code{KI270728.1}) to UCSC names (\code{chr1}, \code{chrM}, \code{chr1_KI270728v1_random}).
#'
#' @param genome A string giving the genome assembly, e.g. \code{'mm39'} or \code{'oviariramb2'}.
#'
#' @param cache A logical stating whether to use the cache on disk for downloads.
#'
#' @return A named character vector whose names are aliases (e.g. Ensembl, GenBank, RefSeq, and UCSC names) and whose values are UCSC-style names.
get_chrom_aliases = function(genome, cache = TRUE) {
    if(genome %in% GENARK$genome) {
        url = get_genark_url(genome, sprintf('%s.chromAlias.txt', GENARK[GENARK$genome == genome, 'accession']))
    } else {
        url = sprintf('https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.chromAlias.txt', genome, genome)
    }
    path = download_annotation_file(url, genome = genome, cache = cache)

    alias_tbl = utils::read.delim(path, header = FALSE, comment.char = '#', colClasses = 'character')

    # The header names the columns, e.g. '# ucsc ensembl genbank refseq'. The
    # UCSC name is the 'ucsc' column, or the first column if there is none
    # (e.g. hg38, '# sequenceName alias names').
    header = strsplit(sub('^#\\s*', '', readLines(path, n = 1)), '\t')[[1]]
    ucsc_col = if('ucsc' %in% header) match('ucsc', header) else 1

    ucsc = alias_tbl[[ucsc_col]]
    aliases = unlist(alias_tbl, use.names = FALSE)
    names(aliases) = aliases
    aliases[] = rep(ucsc, times = ncol(alias_tbl))
    aliases = aliases[names(aliases) != '' & !duplicated(names(aliases))]

    return(aliases)
}

#' Function to get the Seqinfo of a GenArk genome with UCSC-style names
#'
#' @param genome A string giving the genome assembly, one of \code{GENARK$genome}.
#'
#' @param cache A logical stating whether to use the cache on disk for downloads.
#'
#' @return A \code{Seqinfo} object.
get_genark_seqinfo = function(genome, cache = TRUE) {
    sizes = utils::read.delim(download_annotation_file(get_genark_url(genome, sprintf('%s.chrom.sizes.txt', GENARK[GENARK$genome == genome, 'accession'])), genome = genome, cache = cache),
        header = FALSE, col.names = c('chr', 'length'), colClasses = c('character', 'numeric'))
    aliases = get_chrom_aliases(genome, cache = cache)

    seqinfo = Seqinfo::Seqinfo(
        seqnames = unname(aliases[sizes$chr]),
        seqlengths = sizes$length,
        isCircular = unname(aliases[sizes$chr]) == 'chrM',
        genome = genome)

    return(seqinfo)
}

#' Function to rename sequences to UCSC-style names
#'
#' Renames sequences from any alias in the genome's chromAlias file (see \code{get_chrom_aliases()}), e.g. the Ensembl names of an \code{EnsDb}, to UCSC-style names, and sets the genome's full \code{seqinfo}. Sequences without a UCSC name are dropped.
#'
#' @param gr A \code{GRanges} or \code{GRangesList}.
#' @param genome A string giving the genome assembly, e.g. \code{'mm39'} or \code{'oviariramb2'}.
#'
#' @param cache A logical stating whether to use the cache on disk for downloads.
#'
#' @return \code{gr} with UCSC-style sequence names and the genome's \code{seqinfo}.
ucsc_seqlevels = function(gr, genome, cache = TRUE) {
    aliases = get_chrom_aliases(genome, cache = cache)

    # Drop sequences without a UCSC name, keeping the rest of each element of
    # a GRangesList
    mapped = Seqinfo::seqlevels(gr)[Seqinfo::seqlevels(gr) %in% names(aliases)]
    Seqinfo::seqlevels(gr, pruning.mode = 'tidy') = mapped
    Seqinfo::seqlevels(gr) = unname(aliases[Seqinfo::seqlevels(gr)])

    if(genome %in% GENARK$genome) {
        seqinfo = get_genark_seqinfo(genome, cache = cache)
        Seqinfo::seqlevels(gr) = Seqinfo::seqlevels(seqinfo)
        Seqinfo::seqinfo(gr) = seqinfo
    } else {
        gr = set_genome_seqinfo(gr, genome)
    }

    return(gr)
}

#' Function to tidy up annotation accessors for visualization
#'
#' @param annotations A character vector of annotations, in the order they are to appear in the visualization.
#'
#' @return A list of mappings from original annotation names to names ready for visualization.
#' @export
tidy_annotations = function(annotations) {
    tidy = sapply(annotations, function(a){
        tokens = unlist(strsplit(a,'_'))
        if(tokens[2] == 'cpg') {
            if(tokens[3] == 'inter') {
                return('interCGI')
            } else {
                return(paste('CpG', tokens[3]))
            }
        } else if (tokens[2] %in% GENE_GROUPS || (!(tokens[2] %in% BUILTIN_GROUPS) && length(tokens) == 3 && tokens[3] %in% c(GENE_TYPES, 'intergenic'))) {
            if(tokens[3] == 'firstexons') {
                type = 'first exons'
            } else if (tokens[3] == 'intronexonboundaries') {
                type = 'intron/exon boundaries'
            } else if (tokens[3] == 'exonintronboundaries') {
                type = 'exon/intron boundaries'
            } else {
                type = tokens[3]
            }
            # Tell MANE and build_txdb_annotations() groups apart from the
            # builtin genes, e.g. in the same plot
            if(tokens[2] == 'mane') {
                type = paste('MANE', type)
            } else if(tokens[2] != 'genes') {
                type = paste(tokens[2], type)
            }
            return(type)
        } else if (tokens[2] == 'ccre') {
            return(paste('cCRE', CCRE_CLASSES[[tokens[3]]]))
        } else if (tokens[2] == 'enhancers') {
            return('enhancers')
        } else if (tokens[2] == 'chromatin') {
            return(tokens[3])
        } else if (tokens[2] == 'custom') {
            return(tokens[3])
        } else if (tokens[2] == 'lncrna') {
            return('GENCODE lncRNA')
        } else {
            return(sprintf('%s %s', tokens[2], tokens[3]))
        }
    })

    flip_tidy = names(tidy)
    names(flip_tidy) = tidy

    return(as.list(flip_tidy))
}

#' Function to check for valid annotations
#'
#' Gives errors if any annotations are not in builtin_annotations() (and they are not in the required custom format), basicgenes are used, or the genome prefixes are not the same for all annotations.
#'
#' @param annotations A character vector of annotations possibly using the shortcuts
#' @return If all the checks on the annotations pass, returns NULL to allow code to move forward.
check_annotations = function(annotations) {
    # Pull out any custom annotations before checking
    custom_annotations = grep('custom', annotations, value = TRUE)
    annotations = base::setdiff(annotations, custom_annotations)

    # Check that the annotations are supported, tell the user which are unsupported
    if( !all(annotations %in% annotatr::builtin_annotations()) ) {
        unsupported = base::setdiff(annotations, annotatr::builtin_annotations())

        stop(sprintf('Error: "%s" is(are) not supported. See builtin_annotations().',
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
#' @export
expand_annotations = function(annotations) {
    are_basicgenes = any(grepl('basicgenes', annotations))
    are_basicmane = any(grepl('basicmane', annotations))
    are_basiccanonical = any(grepl('basiccanonical', annotations))
    are_ccres = any(grepl('_ccres$', annotations))
    are_cpgs = any(grepl('cpgs', annotations))
    are_hmms = any(grepl('-chromatin', annotations))

    which_are_shortcuts = c(which(grepl('basicgenes', annotations)), which(grepl('basicmane', annotations)), which(grepl('basiccanonical', annotations)), which(grepl('_ccres$', annotations)), which(grepl('cpgs', annotations)), which(grepl('-chromatin', annotations)))

    # expand_shortcuts() will always be run after check_annotations() so we can be
    # sure that the genome prefixes are the same for all annotaitons.
    genome = unique( sapply(annotations, function(a){ unlist(strsplit(a, '_'))[1] }, USE.NAMES = FALSE) )

    if(are_basicgenes || are_basicmane || are_basiccanonical || are_ccres || are_cpgs || are_hmms) {

        # Check for shortcut annotation accessors 'cpgs', 'basicgenes', 'basicmane', 'basiccanonical'
        # and create the right annotations based on the genome
        new_annotations = c()
        remove_shortcuts = c()
        if(are_cpgs) {
            new_annotations = paste(genome, 'cpg', c('islands','shores','shelves','inter'), sep='_')
        }
        if(are_basicgenes) {
            new_annotations = c(new_annotations, paste(genome, 'genes', BASIC_GENE_TYPES, sep='_'))
        }
        if(are_basicmane) {
            new_annotations = c(new_annotations, paste(genome, 'mane', BASIC_GENE_TYPES, sep='_'))
        }
        if(are_basiccanonical) {
            new_annotations = c(new_annotations, paste(genome, 'canonical', BASIC_GENE_TYPES, sep='_'))
        }
        if(are_ccres) {
            new_annotations = c(new_annotations, paste(genome, 'ccre', names(CCRE_CLASSES), sep='_'))
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
#' @export
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

#' Function to give an annotation the standard mcols
#'
#' Adds any missing columns of \code{id}, \code{tx_id}, \code{gene_id}, \code{symbol}, \code{entrez_id}, \code{ensembl_id}, and \code{type} as \code{NA}, and drops any others, so annotations can be combined with \code{c()}.
#'
#' @param gr A \code{GRanges} object of an annotation.
#'
#' @return The \code{GRanges} object with the standard mcols, in order.
standardize_mcols = function(gr) {
    for(col in setdiff(ANNOTATION_MCOLS, colnames(GenomicRanges::mcols(gr)))) {
        GenomicRanges::mcols(gr)[[col]] = rep(NA_character_, length(gr))
    }
    GenomicRanges::mcols(gr) = GenomicRanges::mcols(gr)[, ANNOTATION_MCOLS]

    return(gr)
}

#' Function to remove the version from Ensembl IDs
#'
#' @param ids A character vector of IDs, e.g. \code{'ENSG00000121410.14'}.
#'
#' @return The IDs without versions, e.g. \code{'ENSG00000121410'}. IDs without a version are unchanged.
strip_id_version = function(ids) {
    return(sub('\\.[0-9]+$', '', ids))
}

#' Function to map gene IDs to gene symbols, Entrez IDs, and Ensembl IDs
#'
#' When an ID maps to more than one symbol, Entrez ID, or Ensembl ID, the first is used.
#'
#' @param gene_ids A character vector of gene IDs.
#' @param orgdb An \code{OrgDb} object to map the IDs with, or \code{NULL} to only use the IDs themselves (per \code{keytype}).
#' @param keytype A string giving the type of \code{gene_ids}, as a \code{keytype} of \code{orgdb}, e.g. \code{'ENTREZID'}, \code{'ENSEMBL'}, or \code{'FLYBASE'}. Versions of Ensembl IDs are ignored. If \code{NULL}, the type is unknown.
#'
#' @return A \code{data.frame} with one row per unique gene ID, and columns \code{gene_id}, \code{symbol}, \code{entrez_id}, and \code{ensembl_id}.
get_gene_table = function(gene_ids, orgdb = NULL, keytype = NULL) {
    gene_ids = unique(as.character(gene_ids[!is.na(gene_ids)]))
    lookup = if(identical(keytype, 'ENSEMBL')) strip_id_version(gene_ids) else gene_ids

    table = data.frame(
        gene_id = gene_ids,
        symbol = rep(NA_character_, length(gene_ids)),
        entrez_id = if(identical(keytype, 'ENTREZID')) gene_ids else rep(NA_character_, length(gene_ids)),
        ensembl_id = if(identical(keytype, 'ENSEMBL')) lookup else rep(NA_character_, length(gene_ids)),
        stringsAsFactors = FALSE)

    if(!is.null(orgdb) && !is.null(keytype) && length(gene_ids) > 0) {
        targets = c(symbol = 'SYMBOL', entrez_id = 'ENTREZID', ensembl_id = 'ENSEMBL')
        for(col in names(targets)) {
            if(targets[[col]] == keytype || !(targets[[col]] %in% AnnotationDbi::columns(orgdb))) {
                next
            }
            mapped = tryCatch(
                suppressMessages(AnnotationDbi::mapIds(orgdb, keys = unique(lookup), column = targets[[col]], keytype = keytype, multiVals = 'first')),
                error = function(e) NULL)
            if(!is.null(mapped)) {
                table[[col]] = unname(mapped[lookup])
            }
        }
    }

    return(table)
}
