library(rstatix)
library(ggplot2)

# load in data
confocal <- read.table("confocal.txt", skip = 3, sep = '\t', fill = TRUE)
confocal <- setNames(confocal, c("p7", "p10", "p12"))

confocal <- confocal %>%
  pivot_longer(
    cols = c("p7", "p10", "p12"),
    names_to = "stage",
    values_to = "value"
  )


# boxplot (we got 8 NA values, therefore warnings)
my_comparisons <- list( c("p7", "p10"), c("p10", "p12"), c("p7", "p12") )

ggplot(confocal, aes(x = factor(stage, levels = c("p7", "p10", "p12")), y = value)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) + 
  geom_jitter(width = 0.15, size = 2) +
  theme_bw()

# kruskal-wallis test
res <- kruskal.test(value ~ stage, data = confocal)
res

# dunn test with benjamini hockberg correction
dunn_res <- confocal %>% dunn_test(value ~ stage, p.adjust.method = "BH")
dunn_res