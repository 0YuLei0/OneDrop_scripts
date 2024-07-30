## Load library
library(Seurat)
library(scCustomize)
library(ggsci)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(singleseqgset)
library(UCell)
## Import data
library(tximport)
dir <- "./k_count"
list.files(dir)
mysample <- list.files(dir)
samples<-as.data.frame(mysample)
tx2gene<-read.csv("/import/home/lyuah/reference/tx2gene/gc_ercc_ebv.tx2gene.csv",header = F)
tx2gene <- tx2gene[,2:3]
files <- file.path(dir,samples[,1],"abundance.h5")
names(files) <- mysample
txi.kallisto <- tximport(files, type = 'kallisto', ignoreTxVersion = F, ignoreAfterBar= T, tx2gene = tx2gene)
head(txi.kallisto$counts)
counts<-txi.kallisto$counts
tpm<-txi.kallisto$abundance

sample_anno = as.data.frame(readxl::read_excel("sample_anno.xlsx"))
rownames(sample_anno) = sample_anno$SampleID
head(sample_anno)
all(colnames(counts) %in% sample_anno$SampleID)

r_counts = counts[,sample_anno$SampleID[sample_anno$Assay %in% c("DR_RNA","RNA","SS3")]]
r_tpm = tpm[,sample_anno$SampleID[sample_anno$Assay %in% c("DR_RNA","RNA","SS3")]]
qc_info_all = data.frame(SampleID = colnames(counts),
                      Log10_counts = log10(colSums(counts)),
                      Gene_num =  colSums(counts > 0),
                      MTRNR_ratio = colSums(counts[rownames(counts) %in% c("MT-RNR1","MT-RNR2","MT-RNR3"),]) / colSums(counts),
                      outrRNA_ratio = colSums(counts[rownames(counts) %in% c("FP236383.5","FP236383.4","FP671120.7"),]) / colSums(counts),
                      rRNA_ratio = colSums(counts[rownames(counts) %in% other_rRNA_ratio,]) / colSums(counts),
                      mtRNA_ratio = colSums(counts[rownames(counts) %in% other_mtRNA_ratio,]) / colSums(counts),
                      RPLandRPS_ratio = colSums(counts[rownames(counts) %in% rownames(counts)[c(grep('^RPL', rownames(counts)),grep('^RPS', rownames(counts)))],]) / colSums(counts) ,
                      coding_ratio = colSums(counts[rownames(counts) %in% coding_chr,]) / colSums(counts),
                      log10GenesPerUMI = log10(colSums(counts > 0)) / log10(colSums(counts))
                      )
# Try downsampling
prop = 3203588/colSums(r_counts)
down_counts = DropletUtils::downsampleMatrix(r_counts, prop = prop)
down_counts = as.matrix(down_counts)
havana = read.delim("HGNC_Mart.txt")
coding_chr = havana$Approved.symbol[((havana$Locus.type == "gene with protein product") & (havana$Chromosome != "mitochondria"))]
other_rRNA_ratio = havana$Approved.symbol[havana$Locus.type == "RNA, ribosomal"][-c(1:2)]
other_mtRNA_ratio = na.omit(havana$Approved.symbol[havana$Chromosome == "mitochondria"])[-c(28,29,30)]
coding_chr = havana$Approved.symbol[((havana$Locus.type == "gene with protein product") & (havana$Chromosome != "mitochondria"))]

qc_info = data.frame(SampleID = colnames(down_counts),
                      Log10_counts = log10(colSums(down_counts)),
                      Gene_num =  colSums(down_counts > 0),
                      MTRNR_ratio = colSums(down_counts[rownames(down_counts) %in% c("MT-RNR1","MT-RNR2","MT-RNR3"),]) / colSums(down_counts),
                      outrRNA_ratio = colSums(down_counts[rownames(down_counts) %in% c("FP236383.5","FP236383.4","FP671120.7"),]) / colSums(down_counts),
                      rRNA_ratio = colSums(down_counts[rownames(down_counts) %in% other_rRNA_ratio,]) / colSums(down_counts),
                      mtRNA_ratio = colSums(down_counts[rownames(down_counts) %in% other_mtRNA_ratio,]) / colSums(down_counts),
                      RPLandRPS_ratio = colSums(down_counts[rownames(down_counts) %in% rownames(down_counts)[c(grep('^RPL', rownames(down_counts)),grep('^RPS', rownames(down_counts)))],]) / colSums(down_counts) ,
                      coding_ratio = colSums(down_counts[rownames(down_counts) %in% coding_chr,]) / colSums(down_counts),
                      log10GenesPerUMI = log10(colSums(down_counts > 0)) / log10(colSums(down_counts))
                      )


qc_info = dplyr::inner_join(qc_info, sample_anno)
colSums(down_counts[rownames(down_counts) %in% havana$Approved.symbol[havana$Locus.type == "RNA, long non-coding"],]) / colSums(down_counts)
colSums(down_counts[rownames(down_counts) %in% havana$Approved.symbol[havana$Locus.type == "gene with protein product"],]) / colSums(down_counts)
colSums(down_counts[rownames(down_counts) %in% havana$Approved.symbol[havana$Locus.type == "RNA, ribosomal"],]) / colSums(down_counts)
colSums(down_counts[rownames(down_counts) %in% havana$Approved.symbol[havana$Locus.type == "readthrough"],]) / colSums(down_counts)
#ggplot(qc_info, aes(x=SampleID, y=Log10_counts, fill = Primer)) + geom_col(lwd=0.3) + scale_fill_npg() + labs(x = NULL, y = "Counts")
#ggplot(qc_info, aes(x=SampleID, y=log10GenesPerUMI, fill = Optimization)) + geom_col(lwd=0.3) + scale_fill_npg() + labs(x = NULL, y = "Counts")

ggplot(qc_info, aes(x=SampleID, y=Gene_num, fill = Primer)) + geom_col(lwd=0.3) + scale_fill_npg() + labs(x = NULL, y = "Counts")
ggplot(qc_info, aes(x=SampleID, y=mtRNA_ratio, fill = platform)) + geom_col(lwd=0.3) + scale_fill_npg() + labs(x = NULL, y = "mtRNA_ratio")
ggplot(qc_info, aes(x=SampleID, y=rRNA_ratio, fill = platform)) + geom_col(lwd=0.3) + scale_fill_npg() + labs(x = NULL, y = "Counts")

ggplot(as.data.frame(r_tpm), aes(x="Control-SS3-2", y="Control-SS3-3")) + geom_point()

corrplot::corrplot(cor(log2(r_tpm[rownames(r_tpm) %in% coding_chr, ] + 1),use = "complete.obs",method = "pearson"), is.corr = FALSE,
                                            type="lower", order="hclust", hclust.method = "average",tl.cex = 0.7, method = 'color',
                                           addCoef.col = "black", # Add coefficient of correlation
                                           tl.col="black", tl.srt=90,number.cex = 0.5, col = rev(COL2('RdBu', 100)))

cor(r_tpm[rownames(r_tpm) %in% coding_chr, c("Control-SS3-2","Control-SS3-3")]
qplot( log2(r_tpm[rownames(r_tpm) %in% coding_chr,"Control-SS3-2"] + 1), log2(r_tpm[rownames(r_tpm) %in% coding_chr,"Control-SS3-3"]+1) )

cor(r_tpm[rownames(r_tpm) %in% coding_chr, c("Tube-RNA-Beads","Tube-RNA-Beads-PK")])
qplot( log2(r_tpm[rownames(r_tpm) %in% coding_chr,"CHIP-RNA-Primer"] + 1), log2(r_tpm[rownames(r_tpm) %in% coding_chr,"Tube-RNA-Primer_old"]+1) )

a = "Tube-RNA-Beads"; b = "Tube-RNA-Beads_old"
a = "Tube-RNA-Beads"; b = "Tube-RNA-Beads-PK"
a = "Tube-RNA-Beads"; b = "Tube-RNA-Beads-TFPA"
a = "Tube-RNA-Beads-TFPA"; b = "Tube-RNA-Beads-NEBPA"
a = "Tube-RNA-Beads-TFPA"; b = "Tube-RNA-Beads-NEBPA"
cor(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)])
cor(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)], method = "spearman")
qplot( log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1) )
a = "Tube-DR_RNA-Primer"; b = "Tube-RNA-Primer_new"
a = "CHIP-RNA-Primer"; b = "Tube-RNA-Primer_new"
cor(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)])
cor(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)], method = "spearman")
qplot( log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1) )
a = "CHIP-RNA-Primer"; b = "Tube-RNA-Primer_new"
a = "CHIP-Mock_RNA-Beads"; b = "Tube-RNA-Beads_new"
a = "CHIP-Mock_RNA-Beads"; b = "Tube-RNA-Beads"
a = "CHIP-DR_RNA-Primer"; b = "CHIP-RNA-Primer"
a = "CHIP-Mock_DR_RNA-Beads"; b = "CHIP-Mock_RNA-Beads"

cor(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)])
cor(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)], method = "spearman")
qplot( log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1) )

a = "Control-SS3-2"; b = "Tube-RNA-Primer_new"
a = "Control-SS3-2"; b = "Tube-RNA-Primer-2S"
cor(log2(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)] + 1))
cor(log2(r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)] + 1), method = "spearman")
cor((r_tpm[rownames(r_tpm) %in% coding_chr, c(a, b)] ), method = "spearman")
qplot( log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1) )

#########################################################################################################################################################################
#########################################################################################################################################################################
########################################## DR VS RNA
a = "Tube-DR_RNA-Primer"; b = "Tube-RNA-Primer_new"
DR_vs_RNA_data = data.frame ("Tube-DR_RNA-Primer" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), "Tube-RNA-Primer_new" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1))
DR_vs_RNA_data$color = (DR_vs_RNA_data[,2] > DR_vs_RNA_data[,1] + 1.5) & (DR_vs_RNA_data[,2] > 4) & (DR_vs_RNA_data[,1] > 4)

primer_plot = ggplot(DR_vs_RNA_data) + geom_point(aes(x= `Tube.DR_RNA.Primer`, y = `Tube.RNA.Primer_new`, color = color)) + theme_classic() + scale_color_npg() + geom_line(data=data.frame(a =c(0,15),b = c(0,15)),aes(x=b,y=a))

a = "Tube-DR_RNA-Beads"; b = "Tube-RNA-Beads"
DR_vs_RNA_beads_data = data.frame ("DR_RNA_Beads" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), "RNA_Beads" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1))
DR_vs_RNA_beads_data$color = (DR_vs_RNA_beads_data[,2] > DR_vs_RNA_beads_data[,1] + 1.5) & (DR_vs_RNA_beads_data[,2] > 4) & (DR_vs_RNA_beads_data[,1] > 4)

beads_plot = ggplot(DR_vs_RNA_beads_data) + geom_point(aes(x= `DR_RNA_Beads`, y = `RNA_Beads`, color = color)) + theme_classic() + scale_color_npg() + geom_line(data=data.frame(a =c(0,15),b = c(0,15)),aes(x=b,y=a))
table(DR_vs_RNA_data$color)
table(DR_vs_RNA_beads_data$color)
intersect(rownames(DR_vs_RNA_data)[DR_vs_RNA_data$color], rownames(DR_vs_RNA_beads_data)[DR_vs_RNA_beads_data$color])

library(EnrichmentBrowser)
library("EDASeq")
intersected_genes = filter(havana, Approved.symbol %in% intersect(rownames(DR_vs_RNA_data)[DR_vs_RNA_data$color], rownames(DR_vs_RNA_beads_data)[DR_vs_RNA_beads_data$color])) %>% select("Approved.symbol","Ensembl.gene.ID")
gc_lenght = as.data.frame(getGeneLengthAndGCContent(id=intersected_genes$Ensembl.gene.ID, org="hsa", mode="org.db"))


randome_genes = filter(havana, Approved.symbol %in% sample(coding_chr, 250)) %>% select("Approved.symbol","Ensembl.gene.ID")
randome_genes = randome_genes[!apply(randome_genes == "", 1, any),]
gc_lenght_ran1 = as.data.frame(getGeneLengthAndGCContent(id=randome_genes$Ensembl.gene.ID, org="hsa", mode="org.db"))
gc_lenght_ran2 = as.data.frame(getGeneLengthAndGCContent(id=randome_genes$Ensembl.gene.ID, org="hsa", mode="org.db"))
gc_lenght_ran3 = as.data.frame(getGeneLengthAndGCContent(id=randome_genes$Ensembl.gene.ID, org="hsa", mode="org.db"))

gc_lenght$intersect = "bias"
gc_lenght_ran1$intersect = "random1"
gc_lenght_ran2$intersect = "random2"
gc_lenght_ran3$intersect = "random3"

gc_lenght_all = do.call("rbind", list(gc_lenght, gc_lenght_ran1, gc_lenght_ran2,gc_lenght_ran3))
ggplot(gc_lenght_all,aes(x= intersect, y = log10(length), fill = intersect)) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.2)+ theme_classic() + scale_fill_npg() +
ggplot(gc_lenght_all,aes(x= intersect, y = (gc), fill = intersect)) + geom_violin(width=1) + geom_boxplot(width=0.5, color="grey", alpha=0.2)+ theme_classic() + scale_fill_npg()

gc_lenght = getGeneLengthAndGCContent(id=intersected_genes$Ensembl.gene.ID, org="hsa")
########################################## Tube VS CHIP
a = "CHIP-RNA-Primer"; b = "Tube-RNA-Primer_new"
CP_TP_data = data.frame ("CHIP_RP" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), "Tube_RPnew" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1))
CP_TP_data$color = (CP_TP_data[,2] > CP_TP_data[,1] + 1.5) & (CP_TP_data[,2] > 4) & (CP_TP_data[,1] > 4)
ggplot(CP_TP_data,aes(x= `CHIP_RP`, y = `Tube_RPnew`)) + geom_point(aes(color = color)) + theme_classic() + scale_color_npg() +geom_smooth(method='lm')


a = "CHIP-Mock_RNA-Beads"; b = "Tube-RNA-Beads_new"
CB_TB_data = data.frame ("CHIP_RB" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), "Tube_RBnew" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1))
CB_TB_data$color = (CB_TB_data[,2] > CB_TB_data[,1] + 1.5) & (CB_TB_data[,2] > 4) & (CB_TB_data[,1] > 4)
ggplot(CB_TB_data,aes(x= `CHIP_RB`, y = `Tube_RBnew`)) + geom_point(aes(color = color)) + theme_classic() + scale_color_npg() +geom_smooth(method='lm')

a = "CHIP-Mock_RNA-Beads"; b = "Tube-RNA-Beads"
CB_TB_data2 = data.frame ("CHIP_RB" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,a] + 1), "Tube_RB" = log2(r_tpm[rownames(r_tpm) %in% coding_chr,b]+1))
CB_TB_data2$color = (CB_TB_data2[,2] > CB_TB_data2[,1] + 1.5) & (CB_TB_data2[,2] > 4) & (CB_TB_data2[,1] > 4)
ggplot(CB_TB_data2,aes(x= `CHIP_RB`, y = `Tube_RB`)) + geom_point(aes(color = color)) + theme_classic() + scale_color_npg() +geom_smooth(method='lm')

table(CB_TB_data$color)
table(CB_TB_data2$color)
intersect(rownames(CB_TB_data)[CB_TB_data$color], rownames(CB_TB_data2)[CB_TB_data2$color])
intersect(rownames(CP_TP_data)[CP_TP_data$color], intersect(rownames(CB_TB_data)[CB_TB_data$color], rownames(CB_TB_data2)[CB_TB_data2$color]) )
data.frame(a = intersect(rownames(CP_TP_data)[CP_TP_data$color], intersect(rownames(CB_TB_data)[CB_TB_data$color], rownames(CB_TB_data2)[CB_TB_data2$color]) ))

high_genes_1 = filter(havana, Approved.symbol %in% rownames(CP_TP_data)[CP_TP_data$color]) %>% select("Approved.symbol","Ensembl.gene.ID")
gc_lenght_1 = as.data.frame(getGeneLengthAndGCContent(id=high_genes_1$Ensembl.gene.ID, org="hsa", mode="org.db"))
high_genes_2 = filter(havana, Approved.symbol %in% rownames(CB_TB_data)[CB_TB_data$color]) %>% select("Approved.symbol","Ensembl.gene.ID")
gc_lenght_2 = as.data.frame(getGeneLengthAndGCContent(id=high_genes_2$Ensembl.gene.ID, org="hsa", mode="org.db"))
high_genes_3 = filter(havana, Approved.symbol %in% rownames(CB_TB_data2)[CB_TB_data2$color]) %>% select("Approved.symbol","Ensembl.gene.ID")
gc_lenght_3 = as.data.frame(getGeneLengthAndGCContent(id=high_genes_3$Ensembl.gene.ID, org="hsa", mode="org.db"))

gc_lenght_1$intersect = "primer"
gc_lenght_2$intersect = "beads1"
gc_lenght_3$intersect = "beads2"

gc_lenght_chip_tube = do.call("rbind", list(gc_lenght_1,gc_lenght_2,gc_lenght_3, gc_lenght_ran1, gc_lenght_ran2,gc_lenght_ran3))
ggplot(gc_lenght_chip_tube,aes(x= intersect, y = log10(length), fill = intersect)) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.2)+ theme_classic() + scale_fill_npg() + ggplot(gc_lenght_chip_tube,aes(x= intersect, y = (gc), fill = intersect)) + geom_violin(width=1) + geom_boxplot(width=0.5, color="grey", alpha=0.2)+ theme_classic() + scale_fill_npg()

intersect(
intersect(rownames(CP_TP_data)[CP_TP_data$color], intersect(rownames(CB_TB_data)[CB_TB_data$color], rownames(CB_TB_data2)[CB_TB_data2$color]) ),
intersect(rownames(DR_vs_RNA_data)[DR_vs_RNA_data$color], rownames(DR_vs_RNA_beads_data)[DR_vs_RNA_beads_data$color])
)

sort(ss3_average[intersect(
intersect(rownames(CP_TP_data)[CP_TP_data$color], intersect(rownames(CB_TB_data)[CB_TB_data$color], rownames(CB_TB_data2)[CB_TB_data2$color]) ),
intersect(rownames(DR_vs_RNA_data)[DR_vs_RNA_data$color], rownames(DR_vs_RNA_beads_data)[DR_vs_RNA_beads_data$color])
)])

bias_candidate_genes = intersect(
  intersect(rownames(CP_TP_data)[CP_TP_data$color], intersect(rownames(CB_TB_data)[CB_TB_data$color], rownames(CB_TB_data2)[CB_TB_data2$color]) ),
  intersect(rownames(DR_vs_RNA_data)[DR_vs_RNA_data$color], rownames(DR_vs_RNA_beads_data)[DR_vs_RNA_beads_data$color])
)

DR_vs_RNA_data$color2 = rownames(DR_vs_RNA_data) %in% bias_candidate_genes
DR_vs_RNA_beads_data$color2 = rownames(DR_vs_RNA_beads_data) %in% bias_candidate_genes
CP_TP_data$color2 = rownames(CP_TP_data) %in% bias_candidate_genes
CB_TB_data$color2 = rownames(CB_TB_data) %in% bias_candidate_genes
CB_TB_data2$color2 = rownames(CB_TB_data2) %in% bias_candidate_genes

ggplot(DR_vs_RNA_data) + geom_point(aes(x= `Tube.DR_RNA.Primer`, y = `Tube.RNA.Primer_new`, color = color2)) + theme_classic() + scale_color_npg() + geom_line(data=data.frame(a =c(0,15),b = c(0,15)),aes(x=b,y=a))
ggplot(DR_vs_RNA_beads_data) + geom_point(aes(x= `DR_RNA_Beads`, y = `RNA_Beads`, color = color2)) + theme_classic() + scale_color_npg() + geom_line(data=data.frame(a =c(0,15),b = c(0,15)),aes(x=b,y=a))
ggplot(CP_TP_data,aes(x= `CHIP_RP`, y = `Tube_RPnew`)) + geom_point(aes(color = color2)) + theme_classic() + scale_color_npg() +geom_smooth(method='lm')
ggplot(CB_TB_data,aes(x= `CHIP_RB`, y = `Tube_RBnew`)) + geom_point(aes(color = color2)) + theme_classic() + scale_color_npg() +geom_smooth(method='lm')
ggplot(CB_TB_data2,aes(x= `CHIP_RB`, y = `Tube_RB`)) + geom_point(aes(color = color2)) + theme_classic() + scale_color_npg() +geom_smooth(method='lm')

sort(ss3_average[intersect(
intersect(rownames(CP_TP_data)[CP_TP_data$color], intersect(rownames(CB_TB_data)[CB_TB_data$color], rownames(CB_TB_data2)[CB_TB_data2$color]) ),
intersect(rownames(DR_vs_RNA_data)[DR_vs_RNA_data$color], rownames(DR_vs_RNA_beads_data)[DR_vs_RNA_beads_data$color])
)])

sort(r_tpm[,"Tube-RNA-Primer_new"][bias_candidate_genes], decreasing = T)
sort(ss3_average[bias_candidate_genes], decreasing = T)
ss3_average["GAPDH"]
ss3_average["MYC"]
ss3_average["CASP3"]

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
## analysis
bulk_seu <- CreateSeuratObject(as.matrix(counts), meta.data = sample_anno)
bulk_seu$log10GenesPerUMI <- log10(bulk_seu$nFeature_RNA) / log10(bulk_seu$nCount_RNA)
bulk_seu$percent.mt <- PercentageFeatureSet(object = bulk_seu, pattern = "^MT-")
bulk_seu$percent.RPL <- PercentageFeatureSet(object = bulk_seu, pattern = "^RPL")
bulk_seu$percent.RPS <- PercentageFeatureSet(object = bulk_seu, pattern = "^RPS")
bulk_seu = subset(bulk_seu, subset=(Assay %in%  c("Cell_DR_RNA","Cell_RNA","DR_RNA","RNA","SS3")))
VlnPlot(bulk_seu, features = c("nCount_RNA","nFeature_RNA","percent.mt", "log10GenesPerUMI","percent.RPL","percent.RPS" ), group.by = "SampleID", pt.size = 0.1, ncol = 3)

bulk_seu <- NormalizeData(bulk_seu)
bulk_seu <- ScaleData(bulk_seu)
bulk_seu = FindVariableFeatures(bulk_seu)
bulk_seu <- RunPCA(bulk_seu, npcs = 50, verbose = FALSE)

DimPlot(bulk_seu,group.by = "SampleID")

corrplot::corrplot(cor(tpm[,c(3,4,5,6,9,10,11,12,13,14,17:20)],use = "complete.obs",method = "spearman"),
                     col=col(200), type="lower", order="hclust", hclust.method = "average",tl.cex = 0.7, method = 'color',
                     addCoef.col = "black", # Add coefficient of correlation
                     tl.col="black", tl.srt=90,number.cex = 0.5)
