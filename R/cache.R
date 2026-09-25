#' Cached annotations
#'
#' \code{build_annotations()} saves each annotation it builds, and each file it downloads, in a cache on disk. Later calls load them from the cache instead of building them again, which is much faster (e.g. seconds instead of about half a minute for the hg19 gene annotations). This cache is separate from \code{annotatr_cache}, which holds custom annotations for the current R session only.
#'
#' The cache is in \code{tools::R_user_dir('annotatr', which = 'cache')} and is managed with the \code{BiocFileCache} package. To use another directory, set \code{options(annotatr.cache = '/path/to/cache')}, e.g. in your \code{.Rprofile}.
#'
#' A cached annotation is rebuilt automatically after annotatr is updated, and for gene annotations, after the \code{TxDb.*} or \code{org.*.eg.db} package is updated. Downloaded files are not updated automatically, because their sources rarely change.
#'
#' @section Troubleshooting:
#' \describe{
#'   \item{See what is cached}{\code{list_cached_annotations()} lists each cached annotation and download, with its size and when it was added.}
#'   \item{Annotations seem out of date, or a source has been updated}{Remove the cached annotations and downloads for the genome, e.g. \code{clear_cached_annotations(genome = 'hg19')}, and build again.}
#'   \item{Build once without the cache}{\code{build_annotations(..., cache = FALSE)} builds from scratch, and neither reads nor writes the cache. This is useful to check whether a problem comes from the cache.}
#'   \item{Errors reading a cached file}{A corrupted cached annotation (e.g. from an interrupted build) is removed and rebuilt automatically. If errors continue, remove the genome from the cache with \code{clear_cached_annotations(genome = ...)}.}
#'   \item{A download fails}{Downloads are retried 3 times. A failed download isn't cached, so the next call tries again. Check the network connection, and that the source (e.g. UCSC, GENCODE, FANTOM5) is up.}
#'   \item{The cache uses too much disk space}{Check sizes with \code{list_cached_annotations()}, and remove genomes you no longer need with \code{clear_cached_annotations(genome = ...)}. \code{clear_cached_annotations()} removes everything.}
#'   \item{Home directory is small or read-only (e.g. on a cluster)}{Set \code{options(annotatr.cache = '/path/with/space')}.}
#'   \item{"database is locked" errors}{The cache index is a SQLite database, which can lock when parallel jobs write to it at once. Build the annotations once before starting parallel jobs, so the jobs only read from the cache, or give each job its own cache directory with \code{options(annotatr.cache = ...)}.}
#'   \item{Problems with AnnotationHub}{Some annotations come from AnnotationHub (e.g. hg19 CpG islands, sheep gene models), which has its own cache. See \code{AnnotationHub::hubCache()} for its location, and \code{AnnotationHub::removeResources()} or \code{AnnotationHub::removeCache()} to clear it.}
#' }
#'
#' @return None. This help page describes the cache.
#'
#' @seealso \code{\link{list_cached_annotations}}, \code{\link{clear_cached_annotations}}, \code{\link{build_annotations}}
#'
#' @name cached-annotations
NULL

#' Function to get the directory of the cache
#'
#' @return A string giving the cache directory, from the \code{annotatr.cache} option or \code{tools::R_user_dir()}.
get_cache_dir = function() {
    return(getOption('annotatr.cache', tools::R_user_dir('annotatr', which = 'cache')))
}

#' Function to get the BiocFileCache for the cache
#'
#' @return A \code{BiocFileCache} object.
get_bfc = function() {
    return(BiocFileCache::BiocFileCache(cache = get_cache_dir(), ask = FALSE))
}

#' Function to get the resource IDs of cache entries
#'
#' @param bfc A \code{BiocFileCache} object.
#' @param rname A string giving the resource name.
#' @param exact A logical stating whether to match \code{rname} exactly (TRUE) or as a prefix (FALSE).
#'
#' @return A character vector of resource IDs.
get_cache_rids = function(bfc, rname, exact = TRUE) {
    info = BiocFileCache::bfcinfo(bfc)
    if(exact) {
        matches = info$rname == rname
    } else {
        matches = startsWith(info$rname, rname)
    }

    return(info$rid[matches])
}

#' Function to download a file, with retries, into the cache
#'
#' @param url A string giving the URL.
#' @param genome A string giving the genome the file is for, used by \code{clear_cached_annotations()}.
#' @param cache A logical stating whether to use the cache (TRUE) or a temporary file (FALSE).
#' @param retries The number of times to try the download.
#'
#' @return A string giving the path of the downloaded file. The extension of the URL is kept, so functions that read the file can determine its format.
download_annotation_file = function(url, genome, cache = TRUE, retries = 3) {
    ext = regmatches(basename(url), regexpr('\\.[A-Za-z0-9]+(\\.gz)?$', basename(url)))

    if(cache) {
        bfc = get_bfc()
        rname = sprintf('download|%s|%s', genome, url)
        rid = get_cache_rids(bfc, rname)
        if(length(rid) == 1) {
            return(unname(BiocFileCache::bfcrpath(bfc, rids = rid)))
        }
        path = BiocFileCache::bfcnew(bfc, rname = rname, ext = ext)
    } else {
        path = tempfile(fileext = ext)
    }

    # Large files take longer than R's default timeout of 60 seconds
    old_options = options(timeout = max(600, getOption('timeout')))
    on.exit(options(old_options), add = TRUE)

    for(i in seq_len(retries)) {
        status = tryCatch(
            utils::download.file(url, path, mode = 'wb', quiet = TRUE),
            error = function(e) {
                conditionMessage(e)
            },
            warning = function(w) {
                conditionMessage(w)
            })
        if(identical(status, 0L)) {
            return(unname(path))
        }
        if(i < retries) {
            message(sprintf('Download of %s failed, retrying: %s', url, status))
            Sys.sleep(2)
        }
    }

    # Don't keep a partial download in the cache
    if(cache) {
        BiocFileCache::bfcremove(bfc, names(path))
    }
    stop(sprintf('Failed to download %s after %s tries: %s', url, retries, status))
}

#' Function to get the cache resource name of a built annotation
#'
#' The name includes the annotatr version, and for gene annotations, the versions of the \code{TxDb.*} and \code{org.*.eg.db} packages or the AnnotationHub EnsDb, so updating them causes a rebuild.
#'
#' @param code A string giving the annotation code, e.g. \code{'hg19_genes_promoters'}.
#' @param genome A string giving the genome assembly.
#'
#' @return A string giving the resource name.
get_annotation_rname = function(code, genome) {
    sources = sprintf('annotatr %s', utils::packageVersion('annotatr'))

    if(grepl('_genes_', code)) {
        if(genome %in% GENARK$genome) {
            sources = c(sources, GENARK[GENARK$genome == genome, 'ensdb'])
        } else {
            pkgs = c(get_txdb_name(genome), sprintf('org.%s.eg.db', get_orgdb_name(genome)))
            for(pkg in pkgs) {
                if(requireNamespace(pkg, quietly = TRUE)) {
                    sources = c(sources, sprintf('%s %s', pkg, utils::packageVersion(pkg)))
                }
            }
        }
    }

    return(sprintf('annotation|%s|%s', code, paste(sources, collapse = ', ')))
}

#' Function to load a built annotation from the cache
#'
#' A cached annotation that can't be read is removed so it will be rebuilt.
#'
#' @param code A string giving the annotation code.
#' @param genome A string giving the genome assembly.
#'
#' @return A \code{GRanges} object, or \code{NULL} if the annotation isn't cached.
load_cached_annotation = function(code, genome) {
    bfc = get_bfc()
    rid = get_cache_rids(bfc, get_annotation_rname(code, genome))
    if(length(rid) != 1) {
        return(NULL)
    }

    gr = tryCatch(readRDS(BiocFileCache::bfcrpath(bfc, rids = rid)), error = function(e) {
        message(sprintf('Removing %s from the cache because it could not be read: %s', code, conditionMessage(e)))
        BiocFileCache::bfcremove(bfc, rid)
        NULL
    })

    return(gr)
}

#' Function to save a built annotation in the cache
#'
#' Any earlier versions of the annotation in the cache are removed.
#'
#' @param gr A \code{GRanges} object of the annotation.
#' @param code A string giving the annotation code.
#' @param genome A string giving the genome assembly.
#'
#' @return The path of the cached file, invisibly.
save_cached_annotation = function(gr, code, genome) {
    bfc = get_bfc()

    old_rids = get_cache_rids(bfc, sprintf('annotation|%s|', code), exact = FALSE)
    if(length(old_rids) > 0) {
        BiocFileCache::bfcremove(bfc, old_rids)
    }

    path = BiocFileCache::bfcnew(bfc, rname = get_annotation_rname(code, genome), ext = '.rds')
    saveRDS(gr, path)

    return(invisible(unname(path)))
}

#' List cached annotations and downloads
#'
#' List the annotations and downloaded files in the cache used by \code{build_annotations()}. See \code{\link{cached-annotations}} for more about the cache.
#'
#' @return A \code{data.frame} with one row per cached item, and columns \code{type} (\code{'annotation'} or \code{'download'}), \code{genome}, \code{name} (the annotation code or URL), \code{sources} (the versions an annotation was built with), \code{size_mb}, \code{added}, and \code{path}.
#'
#' @examples
#' # Use a temporary cache so the example doesn't change your cache
#' old_options = options(annotatr.cache = tempfile())
#'
#' list_cached_annotations()
#'
#' options(old_options)
#'
#' @export
list_cached_annotations = function() {
    info = as.data.frame(BiocFileCache::bfcinfo(get_bfc()))
    parts = strsplit(info$rname, '|', fixed = TRUE)

    type = vapply(parts, `[`, character(1), 1)
    name = vapply(parts, `[`, character(1), 3)
    genome = vapply(parts, `[`, character(1), 2)
    sources = rep(NA_character_, length(parts))

    # Annotation names are code|sources, with the genome as the code prefix
    is_annotation = type == 'annotation'
    sources[is_annotation] = name[is_annotation]
    name[is_annotation] = genome[is_annotation]
    genome[is_annotation] = sub('_.*', '', name[is_annotation])

    size_mb = round(file.size(info$rpath) / 1e6, 1)

    return(data.frame(
        type = type,
        genome = genome,
        name = name,
        sources = sources,
        size_mb = size_mb,
        added = info$create_time,
        path = info$rpath,
        stringsAsFactors = FALSE))
}

#' Clear cached annotations and downloads
#'
#' Remove annotations and downloaded files from the cache used by \code{build_annotations()}, so they will be built or downloaded again. See \code{\link{cached-annotations}} for more about the cache.
#'
#' @param genome A character vector of genome assemblies to remove, e.g. \code{'hg19'}. If \code{NULL} (the default), everything is removed.
#'
#' @return A \code{data.frame} of the removed items, as from \code{list_cached_annotations()}, invisibly.
#'
#' @examples
#' # Use a temporary cache so the example doesn't change your cache
#' old_options = options(annotatr.cache = tempfile())
#'
#' # Remove the hg19 annotations and downloads
#' clear_cached_annotations(genome = 'hg19')
#'
#' # Remove everything
#' clear_cached_annotations()
#'
#' options(old_options)
#'
#' @export
clear_cached_annotations = function(genome = NULL) {
    bfc = get_bfc()
    cached = list_cached_annotations()
    info = BiocFileCache::bfcinfo(bfc)

    if(is.null(genome)) {
        remove = rep(TRUE, nrow(cached))
    } else {
        remove = cached$genome %in% genome
    }

    if(any(remove)) {
        BiocFileCache::bfcremove(bfc, info$rid[remove])
    }
    message(sprintf('Removed %s item(s) from the annotatr cache.', sum(remove)))

    return(invisible(cached[remove, , drop = FALSE]))
}
