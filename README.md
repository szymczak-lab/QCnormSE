
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QCnormSE

<!-- badges: start -->

<!-- badges: end -->

The QCnormSE package provides implementations of several quality control
procedures and processing steps for normalized gene expression data sets
from microarray and RNA-seq experiments in a common framework.

## Installation

First install the dependencies from CRAN:

``` r
## install CRAN dependencies
cran.dep = c("aplpack",
             "caret",
             "GGally",
             "ggplot2",
             "ggpubr",
             "Hmisc",
             "tidyr",
             "viridis")

cran.dep.to.install = setdiff(cran.dep,
                     installed.packages()[, "Package"])

if (length(cran.dep.to.install) > 0) {
    install.packages(cran.dep.to.install)
}
```

Next install the dependencies from Bioconductor:

``` r
## install Bioconductor dependencies
bioc.dep = c("affyio",
             "BiocGenerics",
             "ComplexHeatmap",
             "edgeR",
             "doppelgangR",
             "GEOquery",
             "illuminaio",
             "recount",
             "S4Vectors",
             "SummarizedExperiment")

bioc.dep.to.install = setdiff(bioc.dep,
                              installed.packages()[, "Package"])

if (length(bioc.dep.to.install) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install(bioc.dep.to.install)
}
```

Then install the most recent version of QCnormSE from GitHub:

``` r
library(devtools)
install_github("szymczak-lab/QCnormSE")
```

## Usage

A detailed user guide is available as vignette.
