file = '~/Desktop/IDH2mut_v_NBM_annotatr.txt'
annotations = c('basic_genes','cpgs')

d = read_bed(
  file,
  col.names=c('chr','start','end','name','pval','strand','diff_meth','mu1','mu0'),
  genome = 'hg19',
  stranded = FALSE,
  use.score = TRUE)

i = annotate_regions(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = T,
  use.score = TRUE)

i_sub = subset(i, name %in% c('hyper','hypo'))

plot = ggplot(i_sub, aes(diff_meth, -log10(pval))) +
  geom_point(alpha = 1/8, size = 1) +
  geom_abline(intercept=-log10(0.05), slope=0, color='red') +
  facet_wrap(~ annot_type) +
  theme_bw()
