# Major waves of H2A.Z incorporation during mouse oogenesis and preimplantation embryo development
A collection of scripts used for manuscript: Major waves of H2A.Z incorporation during mouse oogenesis and preimplantation embryo development. \
bioRXiv link: https://www.biorxiv.org/content/10.1101/2025.06.14.659461v1

## List of R packages required to run the scripts:
tidyverse \
readxl \
rstatix \
plot3D \
ggrepel \
gprofiler2 \
ggpubr \
gridExtra \
ComplexHeatmap \
ChIPseeker \
TxDb.Mmusculus.UCSC.mm10.knownGene \
GenomicRanges \
corrplot

## Easy installation of all packages by creating a conda environment:
Executing the command below in terminal will create a conda environment with all the required R packages.
```bash
conda create -n h2az_lerdrup r-base=4.3.3 r-tidyverse r-lintr r-languageserver r-devtools nodejs r-readxl r-rstatix r-plot3d r-ggrepel r-gprofiler2 r-ggpubr r-gridextra bioconductor-complexheatmap bioconductor-chipseeker bioconductor-txdb.mmusculus.ucsc.mm10.knowngene bioconda::bioconductor-genomicranges
```

## Table of contents
| File  | Figure |
| ------------- | ------------- |
| fluorescenceBoxplots.R| 1e |
| peaks_PCA.R  | 1g |
| plotAnnoPie.R  | 1f |
| H2AZbubbleplots.Rmd  | 2c |
| LADbubbleplots.Rmd  | 3c |
| plotAnnoClust.R  | 3f |
| peaks_go.R  | 3g |
| peaks_go_heatmap.R  | 3g |
| TEbubbleplots.Rmd  | 4b & S5a |
| z-score_plot.R  | S1c |
| CorrelationMatrices.R | S1f, S2a, S2b & S2d |

## Contact
If you encounter any difficulties running the scripts, please create a new issue in the repository’s Issues section.
