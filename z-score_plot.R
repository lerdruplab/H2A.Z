# this script reads in the H2AZ values quantified and z-score normalized on TSSs
# and plots the z-score values as a line plot

# dependencies
library(tidyverse)
library(gprofiler2)
#library(ggpubr)

# set wd
setwd("~/work/postdoc/projects/h2az/downscaling/")

# function to bin the data
bin_df <- function(df, size) {

  # split df into groups of size
  df_grouped <- split(df, rep(1:(nrow(df)/size), each = size))

  # take mean of each group
  df_binned <- do.call(rbind, lapply(df_grouped, function(x) colMeans(x, na.rm = TRUE)))

  # convert to df
  df_binned <- as.data.frame(df_binned)

  return(df_binned)

}

# import data
z_scores <- read.csv("./uTSS mm10 quant z_score.txt",
                     sep = "\t", header = TRUE)

# remove crappy columns
z_scores <- z_scores[, 5:11]

# rename columns
cols <- c("Mad_500", "Mad_1k", "Ku", "Hu", "Bla", "Bla_Liu_1", "Bla_Liu_2")
colnames(z_scores) <- cols

# get es_df and bla_df
es_df <- z_scores[, 1:4]
bla_df <- z_scores[, 5:7]

# sort based on reference z_score
es_df <- es_df[order(-es_df$Ku), ]
bla_df <- bla_df[order(-bla_df$Bla_Liu_1), ]

# randomly remove 166 rows to be able to bin df to 320
rows_to_keep <- sample(1:nrow(es_df), nrow(es_df) - 166)
es_df <- es_df[rows_to_keep, ]
bla_df <- bla_df[rows_to_keep, ]

# bin the data
es_df_binned <- bin_df(es_df, 320)
bla_df_binned <- bin_df(bla_df, 320)

# order again
es_df_binned <- es_df_binned[order(-es_df_binned$Ku), ]
bla_df_binned <- bla_df_binned[order(-bla_df_binned$Bla_Liu_1), ]

# add ranks as a column to dfs
es_df_binned <- mutate(es_df_binned, rank = 1:nrow(es_df_binned))
bla_df_binned <- mutate(bla_df_binned, rank = 1:nrow(bla_df_binned))

# pivot longer
es_df_binned <- pivot_longer(es_df_binned, cols = c(Ku, Hu, Mad_500, Mad_1k),
                                names_to = "variable", values_to = "value")

bla_df_binned <- pivot_longer(bla_df_binned, cols = c(Bla, Bla_Liu_1, Bla_Liu_2),
                                names_to = "variable", values_to = "value")

# plot the sorted data for each sample
pdf(file = "./es_downscaling_plot.pdf")
ggplot(es_df_binned, aes(x = rank, y = value, color = variable)) + geom_line() +
    theme_classic() + labs(title= "ES cells", x = "Rank", y = "z-score") + ylim(0,3)
dev.off()

pdf(file = "./bla_downscaling_plot.pdf")
ggplot(bla_df_binned, aes(x = rank, y = value, color = variable)) + geom_line() +
    theme_classic() + labs(title= "Blastocyst", x = "Rank", y = "z-score") + ylim (0,3)
dev.off()

