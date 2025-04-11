# this script abbitates H2AZ peaks from each cluster to genomic features and
# visualize them as barplot

# dependencies
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(GenomicRanges)

# set workdir
setwd("~/work/postdoc/projects/h2az/annotation/cluster_peaks/")

# function to process the annotated peaks and generate a df with detailed feature counts
get_detailed_feature_counts <- function(stage) {
  # get df of interest to obtain peak numbers for genomic elements
  stage_df <- data.frame(stage)
  feature_counts <- table(stage_df$annotation)
  feature_counts <- as.data.frame(feature_counts)
  colnames(feature_counts) <- c("genomic_feature", "peak_count")
  # change all specific exon and intron names into just exon and intron
  feature_counts$genomic_feature <- gsub("Exon .*", "Exon", feature_counts$genomic_feature)
  feature_counts$genomic_feature <- gsub("Intron .*", "Intron", feature_counts$genomic_feature)
#  feature_counts$genomic_feature <- gsub("Promoter .*", "Promoter", feature_counts$genomic_feature)
#  feature_counts$genomic_feature <- gsub("Distal Intergenic", "Intergenic", feature_counts$genomic_feature)
#  feature_counts$genomic_feature <- gsub("Downstream .*", "Intergenic", feature_counts$genomic_feature)
  # merge all exon and intron rows into one and sum the peak counts
  feature_counts <- feature_counts %>%
    group_by(genomic_feature) %>%
      summarize(peak_count = sum(peak_count))
}

# function for filtering for a cluster and annotation
get_clust_anno <- function(df, clust) {
  # filter for clusters
  df <- filter(df, Group.order == clust)
  # convert dfs into peakfiles
  df_gr <- GRanges(
                         seqnames = df$Chromosome,
                         ranges = IRanges(start = df$Start, end = df$End),
                         strand = df$strand,
                         score = df$Mean.of...4.p10B.AZ.250.14_S10_R1_001.quantile.normalised.4.p1,
                         peak_name = df$Random.values
  )
  df_gr_anno <- annotatePeak(df_gr, TxDb = txdb)
  return(df_gr_anno)
}

# function to plot each cluster
plot_clust_anno <- function(clust_dat) {
  # get object name for plotting
  x <- deparse(substitute(clust_dat))
  # get feature counts
  clust_dat <- get_detailed_feature_counts(clust_dat)
  # calculate percentage distribution of counts
  clust_dat <- mutate(clust_dat, percentage = round(peak_count / sum(peak_count) *100))
  # change genomic_feature to a foctor and set by levels to desired order
  clust_dat$genomic_feature <- factor(clust_dat$genomic_feature, 
                                      levels = rev(c("Promoter (2-3kb)", "Promoter (1-2kb)",
                                                     "Promoter (<=1kb)", "5' UTR",
                                                     "Exon", "Intron", "3' UTR", 
                                                     "Downstream (<=300bp)", "Distal Intergenic")))
  clust_plot <- ggplot(clust_dat, aes(x=percentage, y = "", fill = genomic_feature)) +
                geom_bar(stat = "identity", width = 1) +
                xlab("Percentage (%)") +
                ylab(x) +
#                scale_x_discrete(labels = c("0", "25", "50", "75", "100"))
#                coord_polar(theta = "y") +
#                xlim(c(0.2, hole_size + 0.5)) +
                theme_classic()
  pdf(paste0("./", x, "_plot.pdf"), height = 1)
  plot(clust_plot)
  dev.off()
}

# import data
all_peaks <- read.csv("./H2A.Z Peaks mouse Merged all samples.txt", sep = "\t")

# load mm10 genomic features
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# get peak cluster annotations
c1 <- get_clust_anno(all_peaks, 1)
c2 <- get_clust_anno(all_peaks, 2)
c3 <- get_clust_anno(all_peaks, 3)
c4 <- get_clust_anno(all_peaks, 4)
c5 <- get_clust_anno(all_peaks, 5)
c6 <- get_clust_anno(all_peaks, 6)
c7 <- get_clust_anno(all_peaks, 7)
c8 <- get_clust_anno(all_peaks, 8)
c9 <- get_clust_anno(all_peaks, 9)
c10 <- get_clust_anno(all_peaks, 10)

# list of clusters for barplot
#all_clusters <- c(Cluster1=c1, Cluster2=c2, Cluster3=c3, Cluster4=c4, Cluster5=c5,
#                  Cluster6=c6, Cluster7=c7, Cluster8=c8, Cluster9=c9, Cluster10=c10)

# barplot from ChipSeeker
#pdf("./cluster_annotation.pdf")
#plotAnnoBar(all_clusters)
#dev.off()

#plot cluster annotations
plot_clust_anno(c1)
plot_clust_anno(c2)
plot_clust_anno(c3)
plot_clust_anno(c4)
plot_clust_anno(c5)
plot_clust_anno(c6)
plot_clust_anno(c7)
plot_clust_anno(c8)
plot_clust_anno(c9)
plot_clust_anno(c10)


