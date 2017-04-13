annotations = build_annotations(genome = 'hg19', annotations = 'hg19_cpgs')
devtools::use_data(annotations, internal = FALSE, compress = 'xz')
