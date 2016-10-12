devtools::load_all()
library(ggplot2)

################################################################################
# Section 3.X
# Removed for actual manuscript

# hpv_p = '../data/hpv+_mc_hmc_perc_meth.txt'
# hpv_n = '../data/hpv-_mc_hmc_perc_meth.txt'
#
# rp = read_bed(hpv_p, genome = 'hg19', stranded = T, use.score = T)
# rn = read_bed(hpv_n, genome = 'hg19', stranded = T, use.score = T)
#
# ap = annotate_regions(rp, annotations = c('hg19_cpgs','hg19_detailedgenes'), genome = 'hg19', ignore.strand = TRUE, use.score = TRUE)
# an = annotate_regions(rn, annotations = c('hg19_cpgs','hg19_detailedgenes'), genome = 'hg19', ignore.strand = TRUE, use.score = TRUE)
#
# sp = summarize_numerical(ap)
# sn = summarize_numerical(an)
#
# cpgs_order = c(
#   'hg19_cpg_islands',
#   'hg19_cpg_shores',
#   'hg19_cpg_shelves',
#   'hg19_cpg_inter')
# vp_cpg = visualize_numerical(sp, cpgs_order, bin_width = 5,
#   plot_title = 'Avg. Meth. Over CpG Annotations (HPV+)', x_label = 'Avg. Meth. per CpG Annotations')
# vn_cpg = visualize_numerical(sn, cpgs_order, bin_width = 5,
#   plot_title = 'Avg. Meth. Over CpG Annotations (HPV-)', x_label = 'Avg. Meth. per CpG Annotation')
#
# ggsave(filename='../paper/hpv+_score_over_cpgs.png',width=6, height=6, plot=vp_cpg, dpi=300)
# ggsave(filename='../paper/hpv-_score_over_cpgs.png',width=6, height=6, plot=vn_cpg, dpi=300)
#
# genes_order = c(
#   'hg19_knownGenes_1to5kb',
#   'hg19_knownGenes_promoters',
#   'hg19_knownGenes_exons5UTRs',
#   'hg19_knownGenes_introns5UTRs',
#   'hg19_knownGenes_exonsCDSs',
#   'hg19_knownGenes_intronsCDSs',
#   'hg19_knownGenes_exons3UTRs',
#   'hg19_knownGenes_introns3UTRs')
# vp_genes = visualize_numerical(sp, genes_order, bin_width = 5,
#   plot_title = 'Avg. Meth. Over knownGene Annotations (HPV+)', x_label = 'Avg. Meth. per knownGenes Annotation')
# vn_genes = visualize_numerical(sn, genes_order, bin_width = 5,
#   plot_title = 'Avg. Meth. Over knownGene Annotations (HPV-)', x_label = 'Avg. Meth. per knownGenes Annotation')
#
# ggsave(filename='../paper/hpv+_score_over_detailed_genes.png',width=12, height=12, plot=vp_genes, dpi=300)
# ggsave(filename='../paper/hpv-_score_over_detailed_genes.png',width=12, height=12, plot=vn_genes, dpi=300)

################################################################################
# Section 3.2

dm = '../data/IDH2mut_v_NBM_DM_new.txt'

rdm = read_bed(dm,
  genome = 'hg19',
  col.names= c(
    'chr','start','end','DM_status','pval',
    'strand','diff_meth','mu1','mu0'),
  stranded = F,
  use.score = T)
adm = annotate_regions(
  regions = rdm,
  annotations = c('hg19_cpgs','hg19_detailedgenes'),
  ignore.strand = TRUE,
  use.score = TRUE)

sdm = summarize_categorical(
  annotated_regions = adm,
  by = c('annot_type','DM_status'))

adm[['logpval']] = -log10(adm[['pval']])

vdm_volcano = visualize_numerical(
  tbl = adm,
  x = 'diff_meth',
  y = 'logpval',
  facet = 'annot_type',
  facet_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter'),
  plot_title = 'Volcano Plots over CpG Annotations',
  x_label = 'Methylation Difference',
  y_label = '-log10(pval)'
  )
vdm_volcano = vdm_volcano + geom_hline(yintercept = -log10(0.001), color = 'red')

ggsave(filename='../paper/DM_volcano.png',width=8, height=8, plot=vdm_volcano, dpi=300)

x_order = c('hyper','hypo')
fill_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
vdm_cpgs = visualize_categorical(
  summarized_cats = sdm,
  x = 'DM_status',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'stack',
  plot_title = 'CpG Annotations by DM Status',
  legend_title = 'Annotation',
  x_label = 'Differential Methylation Status',
  y_label = 'Count'
  )

ggsave(filename='../paper/DM_in_cpgs.png',width=6, height=4, plot=vdm_cpgs, dpi=300)

vdm_cpgs_fill = visualize_categorical(
  summarized_cats = sdm,
  x = 'DM_status',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'fill',
  plot_title = 'CpG Annotations by DM Status',
  legend_title = 'Annotation',
  x_label = 'Differential Methylation Status',
  y_label = 'Proportion'
  )

ggsave(filename='../paper/DM_in_cpgs_fill.png',width=6, height=4, plot=vdm_cpgs_fill, dpi=300)

fill_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_exons5UTRs',
  'hg19_knownGenes_introns5UTRs',
  'hg19_knownGenes_exonsCDSs',
  'hg19_knownGenes_intronsCDSs',
  'hg19_knownGenes_exons3UTRs',
  'hg19_knownGenes_introns3UTRs')
x_order = c('hyper','hypo')
vdm_genes = visualize_categorical(
  summarized_cats = sdm,
  x = 'DM_status',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'stack',
  plot_title = 'KnownGenes Annotations by DM Status',
  legend_title = 'Annotation',
  x_label = 'Differential Methylation Status',
  y_label = 'Count'
  )

ggsave(filename='../paper/DM_in_genes.png',width=6, height=4, plot=vdm_genes, dpi=300)

vdm_genes_fill = visualize_categorical(
  summarized_cats = sdm,
  x = 'DM_status',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'fill',
  plot_title = 'KnownGenes Annotations by DM Status',
  legend_title = 'Annotation',
  x_label = 'Differential Methylation Status',
  y_label = 'Proportion'
  )

ggsave(filename='../paper/DM_in_genes_fill.png',width=6, height=4, plot=vdm_genes_fill, dpi=300)

################################################################################
# Section 3.4

# cl = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
#
# rcl = read_bed(cl, genome = 'hg19', stranded = F, use.score = F)
# acl = annotate_regions(rcl, annotations = c('hg19_cpgs','hg19_detailedgenes'), genome = 'hg19', ignore.strand = TRUE, use.score = FALSE)
# scl = summarize_categorical(acl)
#
# x_order = c(
#   'hyper5mC_5hmC',
#   'hypo5mC_5hmC',
#   'hyper5hmC',
#   'hypo5hmC',
#   'hyper5mC',
#   'hypo5mC',
#   'hyper5mC_hypo5hmC',
#   'hypo5mC_hyper5hmC')
# fill_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
# vcl_cpgs = visualize_categorical(
#   summarized_cats = scl,
#   x = 'data_name',
#   fill = 'annot_type',
#   x_order = x_order,
#   fill_order = fill_order,
#   position = 'stack',
#   plot_title = 'CpG Annotations by Methylation Classification',
#   legend_title = 'Annotation',
#   x_label = 'mC/hmC Methylation Classification',
#   y_label = 'Count'
#   )
#
# ggsave(filename='../paper/classes_in_cpgs.png',width=6, height=6, plot=vcl_cpgs, dpi=300)
#
# vcl_cpgs_fill = visualize_categorical(
#   summarized_cats = scl,
#   x = 'data_name',
#   fill = 'annot_type',
#   x_order = x_order,
#   fill_order = fill_order,
#   position = 'fill',
#   plot_title = 'CpG Annotations by Methylation Classification',
#   legend_title = 'Annotation',
#   x_label = 'mC/hmC Methylation Classification',
#   y_label = 'Proportion'
#   )
#
# ggsave(filename='../paper/classes_in_cpgs_fill.png',width=6, height=6, plot=vcl_cpgs_fill, dpi=300)
#
#
# fill_order = c(
#   'hg19_knownGenes_1to5kb',
#   'hg19_knownGenes_promoters',
#   'hg19_knownGenes_exons5UTRs',
#   'hg19_knownGenes_introns5UTRs',
#   'hg19_knownGenes_exonsCDSs',
#   'hg19_knownGenes_intronsCDSs',
#   'hg19_knownGenes_exons3UTRs',
#   'hg19_knownGenes_introns3UTRs')
# x_order = c(
#   'hyper5mC_5hmC',
#   'hypo5mC_5hmC',
#   'hyper5hmC',
#   'hypo5hmC',
#   'hyper5mC',
#   'hypo5mC',
#   'hyper5mC_hypo5hmC',
#   'hypo5mC_hyper5hmC')
# scl_genes = visualize_categorical(
#   summarized_cats = scl,
#   x = 'data_name',
#   fill = 'annot_type',
#   x_order = x_order,
#   fill_order = fill_order,
#   position = 'stack',
#   plot_title = 'knownGenes Annotations by Methylation Classification',
#   legend_title = 'Annotation',
#   x_label = 'mC/hmC Methylation Classification',
#   y_label = 'Count'
#   )
#
# ggsave(filename='../paper/classes_in_genes.png',width=6, height=6, plot=scl_genes, dpi=300)
#
# scl_genes_fill = visualize_categorical(
#   summarized_cats = scl,
#   x = 'data_name',
#   fill = 'annot_type',
#   x_order = x_order,
#   fill_order = fill_order,
#   position = 'fill',
#   plot_title = 'knownGenes Annotations by Methylation Classification',
#   legend_title = 'Annotation',
#   x_label = 'mC/hmC Methylation Classification',
#   y_label = 'Proportion'
#   )
#
# ggsave(filename='../paper/classes_in_genes_fill.png',width=6, height=6, plot=scl_genes_fill, dpi=300)
