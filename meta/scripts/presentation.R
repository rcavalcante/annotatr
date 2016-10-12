devtools::load_all()
library(ggplot2)

bed = '../data/2607_mc_hmc_perc_meth_test2.txt'
annotations = c('basic_genes','cpgs')

d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = TRUE)

i = intersect_annotations(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = T)

t = annotate_intersections(
  regions = d,
  intersections = i,
  use.score = TRUE)

s = summarize_score(t)

cpgs_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
v_cpgs = visualize_score(s, cpgs_order)

genes_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_5UTRs',
  'hg19_knownGenes_exons',
  'hg19_knownGenes_introns',
  'hg19_knownGenes_3UTRs')
v_genes = visualize_score(s, genes_order)

ggsave(filename='../presentation/2607_mc_hmc_perc_meth_cpgs.png', plot=v_cpgs, width=8, height=8, dpi=300)
ggsave(filename='../presentation/2607_mc_hmc_perc_meth_basic_genes.png', plot=v_genes, width=8, height=8, dpi=300)

################################################################################

bed = '../data/2607_mc_hmc_perc_meth_test2.txt'
annotations = c('detailed_genes')

d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = TRUE)

i = intersect_annotations(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = T)

t = annotate_intersections(
  regions = d,
  intersections = i,
  use.score = TRUE)

s = summarize_score(t)

genes_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_exons5UTRs',
  'hg19_knownGenes_introns5UTRs',
  'hg19_knownGenes_exonsCDSs',
  'hg19_knownGenes_intronsCDSs',
  'hg19_knownGenes_exons3UTRs',
  'hg19_knownGenes_introns3UTRs')
v_genes = visualize_score(s, genes_order)

ggsave(filename='../presentation/2607_mc_hmc_perc_meth_detailed_genes.png', plot=v_genes, width=12, height=12, dpi=300)

################################################################################

bed = '../data/IDH2mut_v_NBM_DM.txt'
annotations = c('detailed_genes','cpgs')

d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = F)

i = intersect_annotations(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = T)

t = annotate_intersections(
  regions = d,
  intersections = i,
  use.score = F)

s = summarize_name(t)

cpgs_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
data_order = c(
  'DMup',
  'DMdown',
  'noDM')
v_cpgs_props = visualize_name(s, cpgs_order, data_order, fill=T)
v_cpgs_counts = visualize_name(s, cpgs_order, data_order, fill=F)

genes_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_exons5UTRs',
  'hg19_knownGenes_introns5UTRs',
  'hg19_knownGenes_exonsCDSs',
  'hg19_knownGenes_intronsCDSs',
  'hg19_knownGenes_exons3UTRs',
  'hg19_knownGenes_introns3UTRs')
data_order = c(
  'DMup',
  'DMdown',
  'noDM')
v_genes_props = visualize_name(s, genes_order, data_order, fill=T)
v_genes_counts = visualize_name(s, genes_order, data_order, fill=F)

ggsave(filename='../presentation/IDH2mut_v_NBM_cpgs_props.png', plot=v_cpgs_props, width=8, height=8, dpi=300)
ggsave(filename='../presentation/IDH2mut_v_NBM_cpgs_counts.png', plot=v_cpgs_counts, width=8, height=8, dpi=300)
ggsave(filename='../presentation/IDH2mut_v_NBM_details_genes_props.png', plot=v_genes_props, width=8, height=8, dpi=300)
ggsave(filename='../presentation/IDH2mut_v_NBM_details_genes_counts.png', plot=v_genes_counts, width=8, height=8, dpi=300)

################################################################################

bed = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
annotations = c('detailed_genes','cpgs')

d = read_bed(file = bed, genome = 'hg19', stranded = F, use.score = F)

i = intersect_annotations(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = T)

t = annotate_intersections(
  regions = d,
  intersections = i,
  use.score = F)

s = summarize_name(t)

cpgs_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
data_order = c(
  'hyper5mC_5hmC',
  'hypo5mC_5hmC',
  'hyper5mC',
  'hypo5mC',
  'hyper5hmC',
  'hypo5hmC',
  'hyper5mC_hypo5hmC',
  'hypo5mC_hyper5hmC')
v_cpgs_props = visualize_name(s, cpgs_order, data_order, fill=T)
v_cpgs_counts = visualize_name(s, cpgs_order, data_order, fill=F)

genes_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_exons5UTRs',
  'hg19_knownGenes_introns5UTRs',
  'hg19_knownGenes_exonsCDSs',
  'hg19_knownGenes_intronsCDSs',
  'hg19_knownGenes_exons3UTRs',
  'hg19_knownGenes_introns3UTRs')
data_order = c(
  'hyper5mC_5hmC',
  'hypo5mC_5hmC',
  'hyper5mC',
  'hypo5mC',
  'hyper5hmC',
  'hypo5hmC',
  'hyper5mC_hypo5hmC',
  'hypo5mC_hyper5hmC')
v_genes_props = visualize_name(s, genes_order, data_order, fill=T)
v_genes_counts = visualize_name(s, genes_order, data_order, fill=F)

ggsave(filename='../presentation/IDH2mut_v_NBM_class_cpgs_props.png', plot=v_cpgs_props, width=12, height=8, dpi=300)
ggsave(filename='../presentation/IDH2mut_v_NBM_class_cpgs_counts.png', plot=v_cpgs_counts, width=12, height=8, dpi=300)
ggsave(filename='../presentation/IDH2mut_v_NBM_class_details_genes_props.png', plot=v_genes_props, width=12, height=8, dpi=300)
ggsave(filename='../presentation/IDH2mut_v_NBM_class_details_genes_counts.png', plot=v_genes_counts, width=12, height=8, dpi=300)
