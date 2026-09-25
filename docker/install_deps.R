# Install the Imports and Suggests of a DESCRIPTION file, plus BiocCheck, that are
# not already installed. Retries because downloads from the package mirrors
# occasionally fail.
#
# Usage: Rscript install_deps.R path/to/DESCRIPTION [--update]
#
#   --update  Also update installed packages that have newer versions, e.g. so
#             a cached library in CI picks up bug fixes.
#
# When run in GitHub Actions, sets the step output 'changed' to whether any
# package was installed or updated, so the workflow saves the cache only then.

args = commandArgs(trailingOnly = TRUE)
update = '--update' %in% args
description = setdiff(args, '--update')[1]

fields = read.dcf(description, fields = c('Depends', 'Imports', 'Suggests'))
pkgs = trimws(sub('\\(.*', '', unlist(strsplit(paste(fields[!is.na(fields)], collapse = ','), ','))))
pkgs = setdiff(c(pkgs[pkgs != ''], 'BiocCheck'), 'R')

# Packages whose loaded version (the first on .libPaths()) is older than the
# repositories'. An update goes to the first library, and leaves the old
# version in a later one, e.g. the container's site library.
outdated_packages = function() {
    available = tryCatch(available.packages(repos = BiocManager::repositories()), error = function(e) NULL)
    if(is.null(available)) {
        return(character(0))
    }
    ip = installed.packages()
    ip = ip[!duplicated(ip[, 'Package']), , drop = FALSE]
    common = intersect(ip[, 'Package'], rownames(available))
    newer = package_version(available[common, 'Version']) > package_version(ip[match(common, ip[, 'Package']), 'Version'])
    return(common[newer])
}

installed_versions = function() {
    ip = installed.packages()
    paste(ip[, 'Package'], ip[, 'Version'], ip[, 'LibPath'])
}
before = installed_versions()

tries = 5
for(i in seq_len(tries)) {
    missing = setdiff(pkgs, rownames(installed.packages()))
    outdated = if(update) outdated_packages() else character(0)
    if(length(missing) == 0 && length(outdated) == 0) {
        break
    }
    if(i > 1) {
        message(sprintf('Retrying in 30 seconds (try %s of %s)', i, tries))
        Sys.sleep(30)
    }
    BiocManager::install(c(missing, outdated), ask = FALSE, update = FALSE)
}

missing = setdiff(pkgs, rownames(installed.packages()))
if(length(missing) != 0) {
    stop('Could not install: ', paste(missing, collapse = ', '))
}

changed = !setequal(before, installed_versions())
message(sprintf('Packages installed or updated: %s', changed))
if(nzchar(Sys.getenv('GITHUB_OUTPUT'))) {
    cat(sprintf('changed=%s\n', tolower(changed)), file = Sys.getenv('GITHUB_OUTPUT'), append = TRUE)
}
