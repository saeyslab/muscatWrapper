% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{make_DEgene_dotplot_pseudobulk_covariate}
\alias{make_DEgene_dotplot_pseudobulk_covariate}
\title{make_DEgene_dotplot_pseudobulk_covariate}
\usage{
make_DEgene_dotplot_pseudobulk_covariate(genes_oi, celltype_info, abundance_data, celltype_oi, covariate_oi, groups_oi = NA)
}
\arguments{
\item{genes_oi}{Character vector with names of genes to visualize}

\item{celltype_info}{`celltype_info` slot of the output of the `muscat_analysis()` function}

\item{abundance_data}{`abundance_data` slot of the output of the `get_abundance_info()` function}

\item{celltype_oi}{Character vector with names of celltype of interest}

\item{covariate_oi}{Name of the covariate/batch that needs to be visualized for each sample}

\item{groups_oi}{Which groups to show? Default: NA -- will show all groups.}
}
\value{
Gene expression dotplot list: pseudobulk version and single-cell version
}
\description{
\code{make_DEgene_dotplot_pseudobulk_covariate}  Visualize the scaled pseudobulk expression of DE genes per sample, and compare the different groups. Genes in rows, samples in columns. Add for each sample to which covariate/batch it belongs.
}
\examples{
\dontrun{
library(dplyr)
sample_id = "tumor"
group_id = "pEMT"
celltype_id = "celltype"
covariates = "batch"
contrasts_oi = c("'High-Low','Low-High'")
contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
min_cells = 10
abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells, covariates = NA)
muscat_output = muscat_analysis(
     sce = sce,
     celltype_id = celltype_id,
     sample_id = sample_id,
     group_id = group_id,
     covariates = covariates,
     contrasts_oi = contrasts_oi,
     contrast_tbl = contrast_tbl)
celltype_oi = "Malignant"
DE_table = filter(muscat_output$celltype_de$celltype_de$de_output_tidy,
                 cluster_id == celltype_oi &
                   p_adj <= 0.05 & logFC >= 1)
DE_genes = unique(DE_table$gene)
dotplots = make_DEgene_dotplot_pseudobulk_covariate(genes_oi = DE_genes, celltype_info = muscat_output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi, covariate_oi = covariates)
}

}
