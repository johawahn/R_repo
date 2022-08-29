setwd("~/Desktop/Werner_Lab/20220617_RNAseq")

library(DESeq2) 
library(edgeR)
library(Biobase)
library(GenomicRanges)
library(ggpubr)

WT1_UN = read.table("GSM5243065_WT1_UN.counts.txt")
WT2_UN = read.table("GSM5243066_WT2_UN.counts.txt")
WT1_W = read.table("GSM5243067_WT1_W.counts.txt")
WT2_W = read.table("GSM5243068_WT2_W.counts.txt")

data = data.frame("WT1_UN" = WT1_UN$V2,
                  "WT2_UN" = WT2_UN$V2,
                  "WT1_W"  = WT1_W$V2,
                  "WT2_W"  = WT2_W$V2,
                  row.names = WT1_UN$V1)

nfat5_i = which(rownames(data) == "Nfat5")
se <- SummarizedExperiment(data.matrix(data))
colSums(assay(se))

#######   Normalization of the counts
dds <- DESeqDataSet( se, design = ~ 1 )

#Estimate size factors
dds <- estimateSizeFactors( dds ) #size factors to account for differences in sequencing depth.
sizeFactors(dds)
colSums(counts(dds))

#Plot column sums according to size factor
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#The argument normalized equals true, divides each column by its size factor.
logcounts <- log2( counts(dds, normalized=TRUE) + 1 ) ### NORMALIZED



##### Plotting
#boxplot(logcounts[nfat5_i,1:2], logcounts[nfat5_i,3:4],
#        main = "Nfat5 RNAseq Counts",
#        ylab = "Count Log2 Size Factor Normalized",
#        names = c("WT_unwounded", "WT_wounded"))


#### With ggplot + pvalue
plot_comparison_gene = function(gene){
  log2nfat5 = data.frame("log2counts" = as.numeric(logcounts[gene,]),
                         "Keratinocytes" = c("UN", "UN", "W", "W"))
  
  p <- ggboxplot(log2nfat5, x = "Keratinocytes", y = "log2counts",
                 color = "Keratinocytes", palette = "jco",
                 add = "jitter")
  
  p + ylab("Log2 Counts") +
    scale_x_discrete(labels=c("Unwounded", "Wounded"), name = "Keratinocytes") +
    xlab("") 
  
  #  Add p-value
  p + stat_compare_means()
  # Change method
  p + stat_compare_means(method = "t.test")
}

plot_comparison_gene(nfat5_i)

#### LOOKING FOR THE OTHER GENES
rps29_i = which(rownames(data) == "Rps29")
akr1b3_i = which(rownames(data) == "Akr1b3")
slc6a6_i = which(rownames(data) == "Slc6a6")
cldn1_i = which(rownames(data) == "Cldn1")
cldn3_i = which(rownames(data) == "Cldn3")
cldn4_i = which(rownames(data) == "Cldn4")

plot_comparison_gene(rps29_i)
plot_comparison_gene(akr1b3_i)
plot_comparison_gene(slc6a6_i)
plot_comparison_gene(cldn1_i)
plot_comparison_gene(cldn3_i)
plot_comparison_gene(cldn4_i)

