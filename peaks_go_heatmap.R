# dependencies
library(tidyverse)
library(ComplexHeatmap)

# set wd
setwd("~/work/postdoc/projects/h2az/peaks/go/")

# load in the data
c1 <- read.table("./cluster1_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c2 <- read.table("./cluster2_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c3 <- read.table("./cluster3_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c4 <- read.table("./cluster4_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c5 <- read.table("./cluster5_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c6 <- read.table("./cluster6_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c7 <- read.table("./cluster7_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c8 <- read.table("./cluster8_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c9 <- read.table("./cluster9_bp_go.tsv", sep = "\t", header = TRUE, quote = "")
c10 <- read.table("./cluster10_bp_go.tsv", sep = "\t", header = TRUE, quote = "")

# create a list of dfs
dfs <- list(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)

# select columns of interest in all dfs
dfs <- lapply(dfs, function(df) select(df, result.term_name, log10_padj))

# merge go tables by term
go_table <- reduce(dfs, full_join, by = "result.term_name")

# write.table(go_table, file = "test_dfs.tsv",
#             row.names = FALSE, quote = FALSE, sep = "\t")

# move term_id to rownames
go_table <- remove_rownames(go_table)
go_table <- column_to_rownames(go_table, var = "result.term_name")

# replace na with 0
go_table <- go_table  %>% replace(is.na(.), 0)

# rename columns
colnames(go_table) <- c("c1", "c2", "c3", "c4",
                        "c5", "c6", "c7", "c8", "c9", "c10")

#colnames(go_table) <- c("c1_-log10_padj", "c2_-log10_padj", "c3_-log10_padj", "c4_-log10_padj",
#                        "c5_-log10_padj", "c6_-log10_padj", "c7_-log10_padj", "c8_-log10_padj", "c9_-log10_padj", "c10_-log10_padj")

# filter out rows if it is 0 in all columns
go_table_filtered <- go_table[rowSums(go_table) > 150, ]
#go_table_filtered <- rownames_to_column(go_table_filtered, var = "go_term:bp")


#write.table(go_table_filtered, "./go_data_for_supp.tsv", quote = FALSE, sep = "\t", 
#            row.names = FALSE)

# export go plots as pdf
pdf(file = "peaks_heatmap_go.pdf",
    width = 6, height = 8)

# create a heatmap
Heatmap(go_table_filtered, cluster_columns = FALSE, column_title = "BP GO",
        width = unit(50, "mm"), height = unit(100, "mm"), show_row_names = FALSE,
        heatmap_legend_param = list(title = "-log10 padj"), col = colorRampPalette(c("white", "black"))(100))

dev.off()

