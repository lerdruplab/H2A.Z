# correlation matrices
library(tidyverse)
library(corrplot)

Madeleinedf <- read.csv("Madeleine Homogenization of mm10 Chromosomes.txt", header = TRUE, sep= '\t')
Liudf <- read.csv("Liu Homogenization of mm10 Chromosomes.txt", header = TRUE, sep= '\t')
XieXudf <- read.csv("Xie_Xu Homogenization of mm10 Chromosomes pooled samples.txt", header = TRUE, sep= '\t')
InoueMeidf <- read.csv("Inoue_Mei Homogenization of mm10 Chromosomes pooled samples.txt", header = TRUE, sep= '\t')

# Renaming
names(Madeleinedf) <- c("Chromosome", "Start", "End", "Strand", "P10", "P7","8-cell", "Morula", "Blastocyst","P10-pooled", "P12", "NSN", "SN","MII","Zygote", "2-cell", "4-cell")
names(Liudf) <- c("Chromosome", "Start", "End", "Strand","2-cellinput","4-cellinput", "8-cellinput", "Blastocystinput","MIIinput","Morulainput","Zygoteinput", "2-cell", "4-cell", "8-cell", "Blastocyst", "MII", "Morula", "Zygote")
names(InoueMeidf) <- c("Chromosome", "Start", "End", "Strand","FGO", "p5", "p10", "p15")
names(XieXudf) <- c("Chromosome", "Start", "End", "Strand","FGO", "NSN", "p7", "SN")

# Correlation plots
# OBS: Corrplot() takes y-axis first, then x

y <- Liudf %>% select(MII, Zygote, "2-cell", "4-cell", "8-cell", Morula, Blastocyst)
x <- Madeleinedf %>% select(MII, Zygote, "2-cell", "4-cell", "8-cell", Morula, Blastocyst)
corrplot(cor(y, x, method = "spearman"), method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)  

y <- InoueMeidf %>% select(p5, p10, p15, FGO)
x <- Madeleinedf %>% select(P7, P10, P12, NSN)
corrplot(cor(y, x, method = "spearman"), method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)  

y <- XieXudf %>% select(p7, FGO, NSN, SN)
x <- Madeleinedf %>% select(P7, P12, NSN, SN)
corrplot(cor(y, x, method = "spearman"), method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)  

y <- InoueMeidf %>% select(p5, FGO)
x <- XieXudf %>% select(p7, FGO)
corrplot(cor(y, x, method = "spearman"), method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)  

# Self correlation plots
MAD <- Madeleinedf %>% select(P7, P10, P12, NSN, SN, MII, Zygote, "2-cell", "4-cell", "8-cell", Morula, Blastocyst)
corrplot(cor(MAD, method = "spearman"), type = "upper", method = "color",addCoef.col = "black", tl.col = "black", tl.srt = 45)

LIU <- Liudf %>% select(MII, Zygote, "2-cell", "4-cell", "8-cell", Morula, Blastocyst)
corrplot(cor(LIU, method = "spearman"), type = "upper", method = "color",addCoef.col = "black", tl.col = "black", tl.srt = 45)  

MEI <- InoueMeidf %>% select(p5, p10, p15, FGO)
corrplot(cor(MEI, method = "spearman"), type = "upper", method = "color", tl.col = "black", addCoef.col = "black", tl.srt = 45)  

XU <- XieXudf %>% select(p7, FGO, NSN, SN)
corrplot(cor(XU, method = "spearman"), type = "upper", method = "color", tl.col = "black", addCoef.col = "black", tl.srt = 45)  
