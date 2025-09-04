# this script takes H2AZ peaks as df and performs
# GO analysis and outputs results for each cluster

# dependencies
library(gprofiler2)
library(tidyverse)
library(ggpubr)
library(gridExtra)

# set wd
setwd("~/work/postdoc/projects/h2az/peaks/")

# import data
df <- read.csv("./H2A_peaks_annotated.txt", sep = "\t")

# clusters are Order.based.on.stages (rearranged by Madeleine)

# optional filtering based on distance to tss/end
df <- filter(df, Distance.to.tss <= 5000)

# function to perform go on a specific cluster
cluster_go <- function(go_df, clust) {

  # filter data
  go_df <- go_df %>% filter(Order.based.on.stages == clust)
  clust_genes <- go_df$Gene.Name

  # GO analysis
  clust_go <- gost(query = clust_genes, organism = "mmusculus",
                   significant = FALSE)

  # get go result df
  clust_go <- as.data.frame(clust_go[1])
  clust_go <- clust_go[, -14]

  # -log10 conversion of padj and get top 10 for plot
  clust_go <- clust_go %>% 
      mutate(log10_padj = -log10(result.p_value)) %>%
      filter(result.source == "GO:BP") %>%
      arrange(log10_padj) #%>%
#      slice_tail(n = 10)

  # export the go result df
  write.table(clust_go, file = paste0("cluster", clust, "_bp_go.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# this visualize the results as barplot. we decided to go for heatmap
# visualization for better clarity and it was ommitted
#  go_plot <- ggpubr::ggbarplot(clust_go,
#                  x = "result.term_name",
#                  y = "log10_padj",
#                  fill = "darkgray",
#                  xlab = "GO Term",
#                  ylab = "p-adjusted(-log10)",
#                  size = 0.5,
#                  palette = "jco",
#                  title = paste("Cluster", clust),
#                  lab.size = 5,
#                  lab.vjust = 0.5,
#                  lab.hjust = 1.2,
#                  sort.by.groups = FALSE,
#                  rotate = TRUE,
#                  ggtheme = ggpubr::theme_pubr(base_size = 18)
#  )

  return(clust_go)

}

# below part is to visualize each peak cluster go results as a barplot
# this was not used for the manuscript
#c1 <- cluster_go(df, 1)
#c2 <- cluster_go(df, 2)
#c3 <- cluster_go(df, 3)
#c4 <- cluster_go(df, 4)
#c5 <- cluster_go(df, 5)
#c6 <- cluster_go(df, 6)
#c7 <- cluster_go(df, 7)
#c8 <- cluster_go(df, 8)
#c9 <- cluster_go(df, 9)
#c10 <- cluster_go(df, 10)

## export go plots as pdf
#pdf(file = "peaks_clusters_go.pdf",
#    width = 32, height = 20)

# arrange the plots layout
#gridExtra::grid.arrange(c1, c2, c3, c4, c5,
#                        c6, c7, c8, c9, c10,
#                        nrow = 4)
#dev.off()

