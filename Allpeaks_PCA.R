

GeneClustering2 <- read.csv("H2A.Z Peaks mouse Merged all samples.txt", sep = "\t")

df <- GeneClustering2 %>% select("X6.p7.Az.250.16_S14_R1_001.quantile.normalised", 
                          "X4.p10B.AZ.250.14_S10_R1_001.quantile.normalised", 
                          "Pooled.1.p12.AZ.200_S1_L001_R1_001.quantile.normalised", 
                          "Pooled.3.NSN.AZ.200_S3_L001_R1_001.quantile.normalised", 
                          "Pooled.4.SN.AZ.200_S4_L001_R1_001.quantile.normalised", 
                          "Pooled.5.MII.AZ.300_S5_L001_R1_001.quantile.normalised", 
                          "X8.8c.AZ.460_S8_L001_R1_001.quantile.normalised", 
                          "X11.mor.AZ460_S11_L001_R1_001.quantile.normalised", 
                          "X13.bla.AZ460_S13_L001_R1_001.quantile.normalised")


names(GeneClustering2)


df <- setNames(df, c("P7", "P10", "P12", "NSN", "SN", "MII", "8cell", "Morula", "Blastocyst"))
numerical_data <- apply(df, 2, function(x) as.numeric(gsub(",", ".", x)))


pcoutput <- pca(t(numerical_data))
dat <- as.data.frame(scores(pcoutput)[,1:2])
dat$celltype <- rownames(dat)


p <- ggplot(dat, aes(PC1, PC2, colour=celltype, label=celltype)) + 
  geom_point(size=3) + 
  geom_text_repel() + 
  theme_minimal() +
  xlab("PC1 43.6%") + 
  ylab("PC2 22.8%")
p



