library(microbenchmark)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
devtools::load_all()

################################################################################
# ChIPpeakAnno

  ################################
  # ~27K rows
  cpa_27k <- function(){
      annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
      file = '../data/Gm12878_Pol2.narrowPeak.gz'
      myPeakList = toGRanges(file, format="BED", header=FALSE)

      annotatedPeak <- annotatePeakInBatch(myPeakList,
                                           AnnotationData=annoData,
                                           output="overlapping")
  }

  ################################
  # ~365K rows
  cpa_365k <- function(){
      annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
      file = '../data/IDH2mut_v_NBM_DM_test_new.txt'
      myPeakList = toGRanges(file, format="BED", header=FALSE)

      annotatedPeak <- annotatePeakInBatch(myPeakList,
                                           AnnotationData=annoData,
                                           output="overlapping")
  }

  ################################
  # ~4M rows
  cpa_4m <- function(){
      annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
      file = '../data/2607_mc_hmc_perc_meth_test2.txt'
      myPeakList = toGRanges(file, format="BED", header=FALSE)
      annotatedPeak <- annotatePeakInBatch(myPeakList,
                                           AnnotationData=annoData,
                                           output="overlapping")
  }
  ################################
  # ~25.4M rows
  cpa_25m <- function(){
      annoData <- annoGR(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
      file = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
      myPeakList = toGRanges(file, format="BED", header=FALSE)
      annotatedPeak <- annotatePeakInBatch(myPeakList,
                                           AnnotationData=annoData,
                                           output="overlapping")
  }

################################################################################
# annotatr

  ################################
  # ~27K rows
  annotatr_27k <- function(){
      file = '../data/Gm12878_Pol2.narrowPeak.gz'
      d = read_bed(
        file = file,
        col.names = FALSE,
        genome = 'hg19',
        stranded = FALSE,
        use.score = FALSE)
      annotations = c('hg19_cpgs', 'hg19_basicgenes')
      t = annotate_regions(
          regions = d,
          annotations = annotations,
          ignore.strand = T,
          use.score = F)
  }

  ################################
  # ~365K rows
  annotatr_365k <- function(){
      file = '../data/IDH2mut_v_NBM_DM_test_new.txt'
      d = read_bed(
        file = file,
        col.names = FALSE,
        genome = 'hg19',
        stranded = FALSE,
        use.score = FALSE)
      annotations = c('hg19_cpgs', 'hg19_basicgenes')
      t = annotate_regions(
          regions = d,
          annotations = annotations,
          ignore.strand = T,
          use.score = F)
  }

  ################################
  # ~4M rows
  annotatr_4m <- function(){
      file = '../data/2607_mc_hmc_perc_meth_test2.txt'
      d = read_bed(
        file = file,
        col.names = FALSE,
        genome = 'hg19',
        stranded = TRUE,
        use.score = FALSE)
      annotations = c('hg19_cpgs', 'hg19_basicgenes')
      t = annotate_regions(
          regions = d,
          annotations = annotations,
          ignore.strand = T,
          use.score = F)
  }

  ################################
  # ~25.4M rows
  annotatr_25m <- function(){
      file = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
      d = read_bed(
        file = file,
        col.names = FALSE,
        genome = 'hg19',
        stranded = FALSE,
        use.score = FALSE)
      annotations = c('hg19_cpgs', 'hg19_basicgenes')
      t = annotate_regions(
          regions = d,
          annotations = annotations,
          ignore.strand = T,
          use.score = F)
  }

message('Benchmarking 27k file')
test_27k = microbenchmark(cpa_27k(), annotatr_27k(), times = 10)
message('Benchmarking 365k file')
test_365k = microbenchmark(cpa_365k(), annotatr_365k(), times = 10)
message('Benchmarking 4m file')
test_4m = microbenchmark(cpa_4m(), annotatr_4m(), times = 10)
message('Benchmarking 25m file')
test_25m = microbenchmark(cpa_25m(), annotatr_25m(), times = 10)

write.table(print(test_27k), file='../paper/microbench_27k_results.txt', sep='\t', quote=F)
write.table(print(test_365k), file='../paper/microbench_365k_results.txt', sep='\t', quote=F)
write.table(print(test_4m), file='../paper/microbench_4m_results.txt', sep='\t', quote=F)
write.table(print(test_25m), file='../paper/microbench_25m_results.txt', sep='\t', quote=F)
