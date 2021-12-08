<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

# muscatWrapper

<!-- badges: start -->

[![R build
status](https://github.com/saeyslab/muscatWrapper/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/saeyslab/muscatWrapper/actions)
[![Coverage
Status](https://codecov.io/gh/saeyslab/muscatWrapper/branch/master/graph/badge.svg?token=0X627I4TM7)](https://codecov.io/gh/saeyslab/muscatWrapper)
<!-- badges: end -->

**muscatWrapper: the R package containing wrapper functions for easier
muscat analysis and downstream visualization (= DE analysis on
single-cell transcriptomics data with complex multi-sample,
multi-condition designs).**

At the basis of muscatWrapper is the differential state analysis as
implemented in the muscat R package
(<https://doi.org/10.1038/s41467-020-19894-4>,
<https://bioconductor.org/packages/release/bioc/html/muscat.html>).
Muscat provides a framework for differential state/expression analyses
based on aggregated “pseudobulk” data . We use this muscat framework to
make inferences on the sample-level (as wanted in a multi-sample,
multi-condition setting) and not the classic cell-level differential
expression analysis of Seurat (Seurat::FindMarkers), because muscat
allows us to overcome some of the limitations of cell-level analyses for
differential state analyses. Some of these limitations include: a bias
towards samples with more cells of cell type, a lack of flexibility to
work with complex study designs, and a too optimistic estimation of the
statistical power since the analysis is done at the cell-level and not
at the sample level. In this package, we provide wrapper functions
around muscat for ease-of-use and flexibility in tweaking parameters.
Moreover, we provide some more advanced downstream visualizations of the
results.

## Installation of muscatWrapper

Installation typically takes a few minutes, depending on the number of
dependencies that has already been installed on your pc.

You can install muscatWrapper (and required dependencies) from github
with:

    # install.packages("devtools")
    devtools::install_github("saeyslab/muscatWrapper")

muscatWrapper is tested via Github Actions version control on Windows,
Linux (Ubuntu) and Mac (most recently tested R version: R 4.1.0).

## Learning to use muscatWrapper

In the following vignettes, you can find how to do a multi-sample
multi-condition DE analysis:

-   [Multi-sample Multi-condition Differential Expression Analysis via
    Muscat: HNSCC application](vignettes/basic_analysis.md):
    `vignette("basic_analysis", package="muscatWrapper")`

## References

Crowell, H.L., Soneson, C., Germain, PL. et al. muscat detects
subpopulation-specific state transitions from multi-sample
multi-condition single-cell transcriptomics data. Nat Commun 11, 6077
(2020). <https://doi.org/10.1038/s41467-020-19894-4>
