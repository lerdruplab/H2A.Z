# this script takes the peaks identified in each oocyte/early embryo stage
# and annotates them to genomic features. pie plot and donut plot from all peaks
# were used in the manuscript. barplot annotation for each stage
# was not included in the final version of the manuscript

# dependencies
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)

# set workdir
setwd("~/work/postdoc/projects/h2az/annotation/genome-wide_peaks/")

# functions
# this function process the annotated peaks and generates a df with generic feature counts
get_generic_feature_counts <- function(stage) {
  # get df of interest to obtain peak numbers for genomic elements
  stage_df <- data.frame(peaks_annotated[[stage]]@anno)

  feature_counts <- table(stage_df$annotation)
  feature_counts <- as.data.frame(feature_counts)
  colnames(feature_counts) <- c("genomic_feature", "peak_count")

  # change all specific exon and intron names into just exon and intron
  feature_counts$genomic_feature <- gsub("Exon .*", "Genic", feature_counts$genomic_feature)
  feature_counts$genomic_feature <- gsub("Intron .*", "Genic", feature_counts$genomic_feature)
  feature_counts$genomic_feature <- gsub("Promoter .*", "Promoter", feature_counts$genomic_feature)
  feature_counts$genomic_feature <- gsub("3' UTR", "Genic", feature_counts$genomic_feature)
  feature_counts$genomic_feature <- gsub("5' UTR", "Genic", feature_counts$genomic_feature)
  feature_counts$genomic_feature <- gsub("Distal Intergenic", "Intergenic", feature_counts$genomic_feature)
  feature_counts$genomic_feature <- gsub("Downstream .*", "Intergenic", feature_counts$genomic_feature)
  
  # merge all exon and intron rows into one and sum the peak counts
  feature_counts <- feature_counts %>%
    group_by(genomic_feature) %>%
      summarize(peak_count = sum(peak_count))

}

# this function process the annotated peaks and generates a df with detailed feature counts
get_detailed_feature_counts <- function(stage) {
  # get df of interest to obtain peak numbers for genomic elements
  stage_df <- data.frame(peaks_annotated[[stage]]@anno)

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

# load mm10 genomic features
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# import peaks
p7 <-readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks p7.txt")
p10 <- readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks p10B.txt")
p12 <- readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks from Pooled p12.txt")
NSN <-readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks from Pooled NSN.txt")
SN <-readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks from Pooled SN.txt")
MII <- readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks from Pooled MII.txt")
c8 <- readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks 8c.txt")
bla <-readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks bla.txt")
mor <-readPeakFile("./peakdata/H2A.Z Peaks mouse Peaks mor.txt")
all <- readPeakFile("./peakdata/H2A.Z Peaks mouse Merged all samples.txt")

# prepare samples list
samples <- list(p7 = p7, p10 = p10, p12 = p12,
                NSN = NSN, SN = SN, MII = MII,
                "8-cell" = c8, Mor = mor, Bla = bla,
                All = all)

# annotate peaks to genomic features
peaks_annotated <- lapply(samples, annotatePeak, TxDb = txdb)

# plot genomic feature distribution as a pie chart
plotAnnoPie(peaks_annotated[["All"]])

# barplot
barplot_stages <- plotAnnoBar(peaks_annotated)

# get feature count tables
p7_df <- data.frame(peaks_annotated[["p7"]])
write.table(p7_df, file = "./p7_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
p7_count <- get_generic_feature_counts("p7")
write.table(p7_count, file = "./p7_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

p10_df <- data.frame(peaks_annotated[["p10"]])
write.table(p10_df, file = "./p10_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
p10_count <- get_generic_feature_counts("p10")
write.table(p10_count, file = "./p10_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

p12_df <- data.frame(peaks_annotated[["p12"]])
write.table(p12_df, file = "./p12_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
p12_count <- get_generic_feature_counts("p12")
write.table(p12_count, file = "./p12_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

NSN_df <- data.frame(peaks_annotated[["NSN"]])
write.table(NSN_df, file = "./NSN_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
NSN_count <- get_generic_feature_counts("NSN")
write.table(NSN_count, file = "./NSN_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

SN_df <- data.frame(peaks_annotated[["SN"]])
write.table(SN_df, file = "./SN_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
SN_count <- get_generic_feature_counts("SN")
write.table(SN_count, file = "./SN_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

MII_df <- data.frame(peaks_annotated[["MII"]])
write.table(MII_df, file = "./MII_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
MII_count <- get_generic_feature_counts("MII")
write.table(MII_count, file = "./MII_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

c8_df <- data.frame(peaks_annotated[["8-cell"]])
write.table(c8_df, file = "./8-cell_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
c8_count <- get_generic_feature_counts("8-cell")
write.table(c8_count, file = "./8-cell_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

mor_df <- data.frame(peaks_annotated[["Mor"]])
write.table(mor_df, file = "./Mor_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
mor_count <- get_generic_feature_counts("Mor")
write.table(mor_count, file = "./Mor_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

bla_df <- data.frame(peaks_annotated[["Bla"]])
write.table(bla_df, file = "./Bla_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
bla_count <- get_generic_feature_counts("Bla")
write.table(bla_count, file = "./Bla_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

all_df <- data.frame(peaks_annotated[["All"]])
write.table(all_df, file = "./All_peak_annotation.csv", quote = FALSE, sep = ",",
            row.names = TRUE)
all_generic_count <- get_generic_feature_counts("All")
all_detailed_count <- get_detailed_feature_counts("All")
write.table(all_generic_count, file = "./All_peak_annotation_count.csv", quote = FALSE, sep = ",",
            row.names = FALSE)
write.table(all_detailed_count, file = "./All_peak_annotation_count_detailed.csv", quote = FALSE, sep = ",",
            row.names = FALSE)

# rearrange the order of rows for plotting
#all_detailed_count <- all_detailed_count[c(2,1,3,5,4,6), ]

# change genomic_feature to a foctor and set by levels to desired order
all_generic_count$genomic_feature <- factor(all_generic_count$genomic_feature, 
                                            levels = rev(c("Promoter", "Genic", "Intergenic")))
all_detailed_count$genomic_feature <- factor(all_detailed_count$genomic_feature, 
                                             levels = rev(c("Promoter (2-3kb)", "Promoter (1-2kb)",
                                                            "Promoter (<=1kb)", "5' UTR",
                                                            "Exon", "Intron", "3' UTR", 
                                                            "Downstream (<=300bp)", "Distal Intergenic")))

# pie plot on all peaks for the paper
pie_plot <- ggplot(all_generic_count, aes(x = "", y = peak_count, fill = genomic_feature)) +
         geom_bar(stat = "identity", width = 1) +
         coord_polar(theta = "y") +
         geom_text(aes(label = paste(genomic_feature, round(peak_count/sum(peak_count) * 100), "%")),
                   position = position_stack(vjust = 0.5)) +
         scale_fill_brewer(palette = "Set2") +
         theme_void()


hole_size <- 3
donut_plot <- ggplot(all_detailed_count, aes(x=hole_size, y = peak_count, fill = genomic_feature)) +
              geom_bar(stat = "identity", width = 1) +
              coord_polar(theta = "y") +
              xlim(c(0.2, hole_size + 0.5)) +
              theme_void()

# save plots
pdf("./pie_plot.pdf")
plot(pie_plot)
dev.off()
pdf("./donut_plot.pdf")
plot(donut_plot)
dev.off()
pdf("./barplot_stages.pdf")
plot(barplot_stages)
dev.off()



