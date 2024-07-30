
######################## DNA Optimization 500kb
~/software/Ginkgo/ginkgo/uploads/OneDrop_20231020

ginkgo = read.delim("~/software/Ginkgo/ginkgo/uploads/OneDrop_20231020/SegNorm")
raw = read.table("~/software/Ginkgo/ginkgo/uploads/OneDrop_20231020/data", header=TRUE, sep="\t")
ginkgo$feature <- paste0("seg_",rownames(ginkgo))
rownames(ginkgo) <- ginkgo$feature
rownames(raw) <- ginkgo$feature
ginkgo <- dplyr::filter(ginkgo, CHR != "chrY")
raw<- raw[rownames(ginkgo),]

ginkgo_cnv_mat <- apply(ginkgo[,-c(1:3,19)], 2, as.numeric)
scWGS_info <- data.frame(sampleID = colnames(ginkgo_cnv_mat),
                         Counts = colSums(raw),
                         mad = c(apply(ginkgo_cnv_mat, 2, mad)) ,
                         fano = c(apply(ginkgo_cnv_mat, 2, function(x){(sd(x))^2/mean(x)}) ),
                         cv = c(apply(ginkgo_cnv_mat, 2, function(x){sd(x)/mean(x)}) )
                         #gini = c(SCOPE::get_gini(Downsample_human_cnv_mat), SCOPE::get_gini(Downsample_mouse_cnv_mat))
                         )
ggplot(scWGS_info, aes(x=sampleID, y=Counts, fill = sampleID)) + geom_col(lwd=0.3) + scale_fill_igv() + labs(x = NULL, y = "Counts") + theme_classic()
ggplot(scWGS_info, aes(x=sampleID, y=mad, fill = sampleID)) + geom_col(lwd=0.3) + scale_fill_igv() + labs(x = NULL, y = "Counts") + theme_classic()
ggplot(scWGS_info, aes(x=sampleID, y=fano, fill = sampleID)) + geom_col(lwd=0.3) + scale_fill_igv() + labs(x = NULL, y = "Counts") + theme_classic()
ggplot(scWGS_info, aes(x=sampleID, y=cv, fill = sampleID)) + geom_col(lwd=0.3) + scale_fill_igv() + labs(x = NULL, y = "Counts") + theme_classic()

plot_sample = unique(scWGS_info$sampleID)[6:15]
ggplot(scWGS_info %>% filter(sampleID %in% plot_sample), aes(x=sampleID, y=mad, fill = sampleID)) + geom_col(lwd=0.3) + scale_fill_npg() + labs(x = NULL, y = "Counts") + theme_classic()

lorenz_prep <- function(x) {
  nReads <- sum(x)
  uniq <- sort(x)
  lorenz = matrix(0, nrow=length(uniq), ncol=2)
  a = c(length(which(x==0)), tabulate(x, nbins=max(x)))
  b = a*(0:(length(a)-1))
  for (i in 2:length(uniq)) {
    lorenz[i,1] = sum(a[1:uniq[i]]) / length(x)
    lorenz[i,2] = sum(b[2:uniq[i]]) / nReads
  }
  return(lorenz)
}
lorenz_plotdataframe = data.frame()
for(i in colnames(raw)){
  temp = data.frame(lorenz_prep(as.matrix(raw[,i])))
  temp$sample = i
  lorenz_plotdataframe = rbind(lorenz_plotdataframe, temp)
}
head(lorenz_plotdataframe)
perfect_lorenz <- data.frame(`X1` = c(0,0.1,0.2,0.3,0.4,1),`X2` = c(0,0.1,0.2,0.3,0.4,1), sample = "Perfect uniformity")
lorenz_plotdataframe = rbind(lorenz_plotdataframe, perfect_lorenz)
plot_sample = c("Perfect uniformity","Control.WGS.1_human_5000000_sampled","Control.WGS.2_human_5000000_sampled","Control.WGS.3_human_5000000_sampled")
plot_sample = c("Perfect uniformity","CHIP.Cells.DNA.Beads.1_human_dedup_sorted","CHIP.Cells.DNA.Beads.2_human_dedup_sorted","CHIP.Cells.DNA.Beads.3_human_dedup_sorted")
plot_sample = c("Perfect uniformity","CHIP.Cells.DR_DNA.Beads.1_human_dedup_sorted","CHIP.Cells.DR_DNA.Beads.2_human_dedup_sorted")
plot_sample = c("Perfect uniformity","CHIP.DNA.Primer_human_5000000_sampled","CHIP.Mock_DNA.Beads_human_5000000_sampled")
plot_sample = c("Perfect uniformity","Tube.DNA.Primer_TCEP_human_5000000_sampled","Tube.DNA.Primer_TCGI_human_5000000_sampled","Tube.DNA.Primer_RT_GITC_human_5000000_sampled","Tube.DNA.Primer.2S_human_5000000_sampled")
plot_sample = c("Perfect uniformity","Tube.DNA.Primer_TCGI_human_5000000_sampled","Tube.DNA.Beads_TCGI_human_5000000_sampled")

plot_sample = c("Perfect uniformity","Tube.DNA.Primer_RT_GITC_human_5000000_sampled","Tube.DNA.Primer.2S_human_5000000_sampled","CHIP.DNA.Primer_human_5000000_sampled","Control.WGS.1_human_5000000_sampled","CHIP.Cells.DNA.Beads.2_human_dedup_sorted","CHIP.Cells.DR_DNA.Beads.1_human_dedup_sorted")

ggplot(lorenz_plotdataframe %>% filter(sample %in% plot_sample), aes(x=X1,y=X2)) +
  geom_line(aes(color = sample)) +
  scale_color_igv() +
  labs(x = "Cumulative fraction of genome", y = "Cumulative fraction of total reads") + theme_classic()

mergeLevels_multi_segmentation <- function(smoothed_cnv = smooth_cnv,
                                             postionanno = bin_anno,
                                             gamma.param = 5, hg38 = hg38) {
    ## all chr
    allchr_2 <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
    ##
    log.ratio.mats <- smoothed_cnv
    log.ratio.mats[log.ratio.mats <= 0] <- 0.0000001
    log.ratio.mats <- apply(log.ratio.mats, 2, function(x){x / median(x)})
    #log.ratio.mats <- log.ratio.mats / normal
    log.ratio.mats <- log2(log.ratio.mats)
    postionanno <- postionanno %>% dplyr::select(CHR, START)
    postionanno$CHR <- gsub('[chr]', '', postionanno$CHR)
    postionanno$CHR  <- factor(postionanno$CHR, levels = allchr_2)
    sample.matrices <- data.frame(cbind(postionanno, log.ratio.mats))
    sample.matrices <- copynumber::winsorize(sample.matrices,assembly = "hg38")
    segmentations <- multipcf(sample.matrices,gamma = gamma.param,assembly = "hg38",return.est=TRUE,fast = FALSE)
    log.ratio.longs <- segmentations$estimates[,3:ncol(segmentations$estimates)]
    ratio.longs <- 2^log.ratio.longs
    rownames(ratio.longs) <- postionanno$feature
    #colnames(ratio.longs) <- gsub('[.]', '-', colnames(ratio.longs))
    seg.ml.data <- data.frame(row.names = rownames(smoothed_cnv),ratio.longs)
    for(i in colnames(smoothed_cnv)){
      logratio <- log.ratio.mats[,i]
      logratio <- 2^logratio
      seg.mean <- ratio.longs[,i]
      ml <- aCGH::mergeLevels(logratio,seg.mean)
      seg.ml.data[,i] <- ml$vecMerged
    }
    seg.data <- ratio.longs
    output <- list(seg.data = seg.data, seg.ml.data = seg.ml.data)
    return(output)
}
process_benchmark <- function(dir, samples = c()){
  sc_cnv <- read.delim(dir)
  sc_cnv_anno <- dplyr::select(sc_cnv,CHR,START,END)
  sc_cnv_anno$feature <- paste0("seg_",seq(1:length(rownames(sc_cnv))))
  sc_cnv_anno <- sc_cnv_anno[sc_cnv_anno$CHR !="chrY",]
  allchr <-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  sc_cnv <- sc_cnv %>% dplyr::filter(CHR !="chrY") %>% dplyr::select(-CHR,-START,-END)
  rownames(sc_cnv) <- sc_cnv_anno$feature
  #sc_cnv <- sc_cnv[,samples]
  col_demode <- 1 / apply(sc_cnv,2,densMode)
  sc_cnv <- sweep(sc_cnv,2,col_demode,"*")
  sc_cnv_seg <- mergeLevels_multi_segmentation(smoothed_cnv = sc_cnv, postionanno = sc_cnv_anno, gamma.param = 40, hg38 = hg38)
  sc_cnv <- data.frame(CHR = sc_cnv_anno$CHR, START = sc_cnv_anno$START, END = sc_cnv_anno$END,
                       feature = sc_cnv_anno$feature, round(sc_cnv_seg$seg.ml.data * 2)) %>%
    makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T)
  out <- list(sc_cnv_seg, sc_cnv)
  return(out)
}
CNV_calling_accuracy <- function(true_cnv, scCNV) {
  features <- names(true_cnv)
  mymat <- mcols(scCNV)[-1]
  samples <- colnames(mymat)
  ginkgo_mabf <- matrix(0, nrow = length(features), ncol = length(samples))
  rownames(ginkgo_mabf) <- features
  colnames(ginkgo_mabf) <- samples

  olaps1 <- findOverlaps(true_cnv, scCNV)

  for (i in seq(1: length(samples))) {
    mysample <- samples[i]

    mk_df1 <- tibble(feature = names(true_cnv)[queryHits(olaps1)],
                     mbaf = mymat[,mysample][subjectHits(olaps1)]
                     ) %>%
      dplyr::distinct(feature, .keep_all = TRUE)

    ginkgo_mabf <- ginkgo_mabf[mk_df1$feature,]
    ginkgo_mabf[,mysample] <- mk_df1$mbaf
  }
  accuracy_table <- data.frame(sample_id = samples, recall = 0, precision = 0, f1_score=0)
  for (i in seq(1: length(samples))) {
    mysample <- samples[i]
    Actual = factor(mcols(true_cnv)[[1]], levels = c(1,2,3))
    Predicted = factor(ginkgo_mabf[,mysample], levels = c(1,2,3))
    confu_mat <- as.matrix(table(Actual = Actual, Predicted =Predicted))
    accuracy = sum(diag(confu_mat)) / sum(confu_mat)
    precision = diag(confu_mat) / apply(confu_mat, 2, sum)
    recall = diag(confu_mat) / apply(confu_mat, 1, sum)
    f1 = (2 * precision * recall) / (precision + recall)
    accuracy_table$recall[i] <- (sum(recall) - recall[2])/2
    accuracy_table$precision[i] <- (sum(precision) - precision[2])/2
    accuracy_table$f1_score[i] <- (sum(f1) - f1[2])/2
  }
  return(accuracy_table)
}
CNV_calling_rate <- function(true_cnv, scCNV) {
  features <- names(true_cnv)
  mymat <- mcols(scCNV)[-1]
  samples <- colnames(mymat)
  ginkgo_mabf <- matrix(0, nrow = length(features), ncol = length(samples))
  rownames(ginkgo_mabf) <- features
  colnames(ginkgo_mabf) <- samples

  olaps1 <- findOverlaps(true_cnv, scCNV)

  for (i in seq(1: length(samples))) {
    mysample <- samples[i]

    mk_df1 <- tibble(feature = names(true_cnv)[queryHits(olaps1)],
                     mbaf = mymat[,mysample][subjectHits(olaps1)]
    ) %>%
      dplyr::distinct(feature, .keep_all = TRUE)

    ginkgo_mabf <- ginkgo_mabf[mk_df1$feature,]
    ginkgo_mabf[,mysample] <- mk_df1$mbaf
  }
  accuracy_table <- data.frame(sample_id = samples, rate = 0)
  for (i in seq(1: length(samples))) {
    mysample <- samples[i]
    Actual = mcols(true_cnv)[[1]]
    Predicted = ginkgo_mabf[,mysample]

    accuracy_table$rate[i] <- sum(Predicted == Actual) / length(Actual)
  }
  return(accuracy_table)
}
detach("package:copynumber", unload=TRUE)
dir = "~/software/Ginkgo/ginkgo/uploads/HCT116_Bulk/HCT_CNV/SegNorm"
dir = "~/software/Ginkgo/ginkgo/uploads/OneDrop_20231020/SegNorm"
hct_bulk_cnv <- process_benchmark("~/software/Ginkgo/ginkgo/uploads/HCT116_Bulk/HCT_CNV/SegNorm")
mcols(hct_bulk_cnv[[2]]) <- mcols(hct_bulk_cnv[[2]])["ccle_hct"]

Benchmark_500kb <- process_benchmark("~/software/Ginkgo/ginkgo/uploads/OneDrop_20231020/SegNorm")
cnv_accuracy <- CNV_calling_accuracy(true_cnv = hct_bulk_cnv[[2]], scCNV = Benchmark_500kb[[2]])
cnv_accuracy$sampleID = rownames(cnv_accuracy)
CNV_calling_rate(true_cnv = hct_bulk_cnv[[2]], scCNV = Benchmark_500kb[[2]])

plot_sample = unique(cnv_accuracy$sample_id)[c()]
ggplot(cnv_accuracy %>% filter(sample_id %in% plot_sample), aes(x=sample_id, y=precision, fill = sample_id)) + geom_col(lwd=0.3) + scale_fill_npg() + labs(x = NULL, y = "Counts") + theme_classic()
