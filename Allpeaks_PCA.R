library('tidyverse')
library('plot3D')

# Load in the data
allPeaks <- read.csv("Merged all samples with New Mike quantified and quantile normalized.txt", sep = "\t")
names(allPeaks)
df <- allPeaks %>% select("X6.p7.Az.250.16_S14_R1_001.quantile.normalised", 
                          "X4.p10B.AZ.250.14_S10_R1_001.quantile.normalised", 
                          "Pooled.1.p12.AZ.200_S1_L001_R1_001.quantile.normalised", 
                          "Pooled.3.NSN.AZ.200_S3_L001_R1_001.quantile.normalised", 
                          "Pooled.4.SN.AZ.200_S4_L001_R1_001.quantile.normalised", 
                          "Pooled.5.MII.AZ.300_S5_L001_R1_001.quantile.normalised", 
                          "X8.8c.AZ.460_S8_L001_R1_001.quantile.normalised", 
                          "X11.mor.AZ460_S11_L001_R1_001.quantile.normalised", 
                          "X13.bla.AZ460_S13_L001_R1_001.quantile.normalised",
                          "X3.MC.ChIP1.2AZ.Zyg.200.20250804_S3_L001_R1_001.quantile.normalised",
                          "X4.MC.ChIP1.2AZ.2c.200.20250804_S4_L001_R1_001.quantile.normalised",
                          "X5.MC.ChIP1.2AZ.4c.200.20250804_S5_L001_R1_001.quantile.normalised")

# Filtering for low sinal (based ont sum of all)
df$row_sum <- rowSums(df)
df <- setNames(df, c("P7", "P10", "P12", "NSN", "SN", "MII", "8cell", "Morula", "Blastocyst", "Zygote", "2cell", "4cell", "row_sums"))
df <- df %>% filter(row_sums > 10)
df <- subset(df, select = -row_sums)

# PCA calculation
pca_result <- prcomp(t(df))
summary(pca_result)

#using the other package
scatter3D(pca_result$x[,1], pca_result$x[,2], pca_result$x[,3], phi = 25, bty = "g",  type = "h", 
          ticktype = "detailed", pch = 19, cex = 1)


text3D(pca_result$x[,1], pca_result$x[,2], pca_result$x[,3],  labels = colnames(df),
       add = TRUE, colkey = FALSE, cex = 0.5)
