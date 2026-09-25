# Install the Imports and Suggests of a DESCRIPTION file, plus BiocCheck, that are
# not already installed. Retries because downloads from the package mirrors
# occasionally fail.
#
# Usage: Rscript install_deps.R path/to/DESCRIPTION

args = commandArgs(trailingOnly = TRUE)

fields = read.dcf(args[1], fields = c('Depends', 'Imports', 'Suggests'))
pkgs = trimws(sub('\\(.*', '', unlist(strsplit(paste(fields[!is.na(fields)], collapse = ','), ','))))
pkgs = setdiff(c(pkgs[pkgs != ''], 'BiocCheck'), 'R')

for(i in 1:3) {
    missing = setdiff(pkgs, rownames(installed.packages()))
    if(length(missing) == 0) {
        break
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
}

missing = setdiff(pkgs, rownames(installed.packages()))
if(length(missing) != 0) {
    stop('Could not install: ', paste(missing, collapse = ', '))
}
