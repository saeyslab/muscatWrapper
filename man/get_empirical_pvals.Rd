% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_wrappers.R
\name{get_empirical_pvals}
\alias{get_empirical_pvals}
\title{get_empirical_pvals}
\usage{
get_empirical_pvals(de_output_tidy)
}
\arguments{
\item{de_output_tidy}{Differential expression analysis output for the sender cell types. `de_output_tidy` slot of the output of `perform_muscat_de_analysis`.}
}
\value{
`de_output_tidy`, but now 2 columns added with the empirical pvalues (normal and adjusted for multiple testing); Histogram plot of the empirical p-values is also returned.
}
\description{
\code{get_empirical_pvals} Calculate empirical p-values based on a DE output. Show p-value distribution histograms. Under the hood, the following functions are used: `add_empirical_pval_fdr` and `get_FDR_empirical_plots_all`
}
\examples{
\dontrun{
library(dplyr)
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network \%>\% dplyr::rename(ligand = from, receptor = to) \%>\% dplyr::distinct(ligand, receptor)
sample_id = "tumor"
group_id = "pEMT"
celltype_id = "celltype"
covariates = NA
contrasts_oi = c("'High-Low','Low-High'")
DE_info = get_DE_info(
   sce = sce,
   sample_id = sample_id,
   celltype_id = celltype_id,
   group_id = group_id,
   covariates = covariates,
   contrasts = contrasts_oi)
DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
}

}
