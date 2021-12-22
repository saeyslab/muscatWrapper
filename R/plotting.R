scaling_zscore = function (x)
{
  if (typeof(x) == "double") {
    if (sd(x, na.rm = TRUE) > 0) {
      return((x - mean(x, na.rm = TRUE ))/sd(x, na.rm = TRUE))
    }
    else {
      return((x - mean(x, na.rm = TRUE)))
    }
  }
  else {
    return(x)
  }
}
#' @title make_DEgene_dotplot_pseudobulk
#'
#' @description \code{make_DEgene_dotplot_pseudobulk}  Visualize the scaled pseudobulk expression of DE genes per sample, and compare the different groups. Genes in rows, samples in columns
#' @usage make_DEgene_dotplot_pseudobulk(genes_oi, celltype_info, abundance_data, celltype_oi, groups_oi = NA)
#'
#' @param genes_oi Character vector with names of genes to visualize
#' @param celltype_info `celltype_info` slot of the output of the `muscat_analysis()` function
#' @param abundance_data `abundance_data` slot of the output of the `get_abundance_info()` function
#' @param celltype_oi Character vector with names of celltype of interest
#' @param groups_oi Which groups to show? Default: NA -- will show all groups.
#'
#' @return Gene expression dotplot list: pseudobulk version and single-cell version
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' min_cells = 10
#' abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells, covariates = NA)
#' muscat_output = muscat_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl)
#' celltype_oi = "Malignant"
#' DE_table = filter(muscat_output$celltype_de$celltype_de$de_output_tidy,
#'                  cluster_id == celltype_oi &
#'                    p_adj <= 0.05 & logFC >= 1)
#' DE_genes = unique(DE_table$gene)
#' dotplots = make_DEgene_dotplot_pseudobulk(genes_oi = DE_genes, celltype_info = muscat_output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi)
#'}
#'
#' @export
#'
make_DEgene_dotplot_pseudobulk = function(genes_oi, celltype_info, abundance_data, celltype_oi, groups_oi = NA){

  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####

  keep_sender_receiver_values = c(1, 4)
  names(keep_sender_receiver_values) = c(FALSE, TRUE)

  abundance_data = abundance_data %>% dplyr::rename(sample = sample_id, celltype = celltype_id, group = group_id) %>% select(-n)
  plot_data = celltype_info$pb_df %>% inner_join(abundance_data, by = c("sample", "celltype"))
  plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = scaling_zscore(pb_sample)) %>% dplyr::ungroup()
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)

  plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)

  if(!is.na(groups_oi)){
    plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)

  }

  p1 = plot_data %>%
    ggplot(aes(sample, gene, color = scaled_gene_exprs, size = keep)) +
    geom_point() +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled pseudobulk\nexpression", size= "Celltype present") + xlab("") +ylab("") +
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p1 = p1 + custom_scale_fill

  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####

  frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
  plot_data = plot_data %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype"))
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)

  p2 = plot_data %>%
    # ggplot(aes(gene, sample , fill = scaled_gene_exprs)) +
    # geom_tile(color = "white") +
    ggplot(aes(sample, gene , color = scaled_gene_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    )  +
    # labs(color = "Scaled gene\navg expression")
    labs(color = "Scaled pseudobulk\nexpression", size= "Fraction of\nexpressing cells") + xlab("") +ylab("")
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p2 = p2 + custom_scale_fill

  return(list(pseudobulk_plot = p1, singlecell_plot = p2))

}
#' @title make_DEgene_dotplot_pseudobulk_reversed
#'
#' @description \code{make_DEgene_dotplot_pseudobulk_reversed}  Visualize the scaled pseudobulk expression of DE genes per sample, and compare the different groups. Genes and sample positions are reversed compared to `make_DEgene_dotplot_pseudobulk`: genes in columns, samples in rows.
#' @usage make_DEgene_dotplot_pseudobulk_reversed(genes_oi, celltype_info, abundance_data, celltype_oi, groups_oi = NA)
#'
#' @inheritParams make_DEgene_dotplot_pseudobulk
#' @return Gene expression dotplot list: pseudobulk version and single-cell version
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' min_cells = 10
#' abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells, covariates = NA)
#' muscat_output = muscat_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl)
#' celltype_oi = "Malignant"
#' DE_table = filter(muscat_output$celltype_de$celltype_de$de_output_tidy,
#'                  cluster_id == celltype_oi &
#'                    p_adj <= 0.05 & logFC >= 1)
#' DE_genes = unique(DE_table$gene)
#' dotplots = make_DEgene_dotplot_pseudobulk_reversed(genes_oi = DE_genes, celltype_info = muscat_output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi)
#'}
#'
#' @export
#'
make_DEgene_dotplot_pseudobulk_reversed = function(genes_oi, celltype_info, abundance_data, celltype_oi, groups_oi = NA){

  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####

  keep_sender_receiver_values = c(1, 4)
  names(keep_sender_receiver_values) = c(FALSE, TRUE)

  abundance_data = abundance_data %>% dplyr::rename(sample = sample_id, celltype = celltype_id, group = group_id) %>% select(-n)
  plot_data = celltype_info$pb_df %>% inner_join(abundance_data, by = c("sample", "celltype"))
  plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = scaling_zscore(pb_sample)) %>% dplyr::ungroup()
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)

  plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)

  if(!is.na(groups_oi)){
    plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)

  }

  p1 = plot_data %>%
    ggplot(aes(gene, sample, color = scaled_gene_exprs, size = keep)) +
    geom_point() +
    facet_grid(group ~ ., scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.x = element_text(face = "bold.italic", size = 9, angle = 90,hjust = 0),
      axis.text.y = element_text(size = 9),
      strip.text.y.left = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(1.5, "lines"),
      panel.spacing.x = unit(0.25, "lines"),
      strip.text.y = element_text(size = 11, color = "black", face = "bold", angle = 0),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled pseudobulk\nexpression", size= "Celltype present") + xlab("") +ylab("") +
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p1 = p1 + custom_scale_fill

  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####

  frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
  plot_data = plot_data %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype"))
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)

  p2 = plot_data %>%
    # ggplot(aes(gene, sample , fill = scaled_gene_exprs)) +
    # geom_tile(color = "white") +
    ggplot(aes(gene, sample , color = scaled_gene_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(group ~ ., scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.x = element_text(face = "bold.italic", size = 9, angle = 90,hjust = 0),
      axis.text.y = element_text(size = 9),
      strip.text.y.left = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(1.5, "lines"),
      panel.spacing.x = unit(0.25, "lines"),
      strip.text.y = element_text(size = 11, color = "black", face = "bold", angle = 0),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    )  +
    # labs(color = "Scaled gene\navg expression")
    labs(color = "Scaled pseudobulk\nexpression", size= "Fraction of\nexpressing cells") + xlab("") +ylab("")
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p2 = p2 + custom_scale_fill

  return(list(pseudobulk_plot = p1, singlecell_plot = p2))

}
#' @title make_DEgene_violin_plot
#'
#' @description \code{make_DEgene_violin_plot}  Violin plot of a gene gene of interest: per sample, and samples are grouped per group
#' @usage make_DEgene_violin_plot(sce, gene_oi, celltype_oi, group_id, sample_id, celltype_id, covariate_oi = NA, groups_oi = NA)
#'
#' @param gene_oi Name of the gene of interest
#' @param celltype_oi  Character vector with the names of the cell type of interest
#' @param covariate_oi Name of a covariate of interest based on which the visualization will be split. Default: NA: no covariate.
#' @inheritParams make_DEgene_dotplot_pseudobulk
#' @inheritParams muscat_analysis
#'
#' @return ggplot object: Violin plot of a gene gene of interest: per sample, and samples are grouped per group
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom generics setdiff
#' @importFrom muscat prepSCE
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' min_cells = 10
#' abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells, covariates = NA)
#' muscat_output = muscat_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl)
#' celltype_oi = "Malignant"
#' group_oi = "High"
#' gene_oi = "RAB31"
#' p_violin_gene = make_DEgene_violin_plot(sce = sce, gene_oi = gene_oi, celltype_oi = celltype_oi, group_id = group_id, sample_id, celltype_id = celltype_id)
#' }
#' @export
#'
make_DEgene_violin_plot = function(sce, gene_oi, celltype_oi, group_id, sample_id, celltype_id, covariate_oi = NA, groups_oi = NA){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% celltype_oi]

  if(!is.na(groups_oi)){
    sce_subset =  sce[, SummarizedExperiment::colData(sce)[,group_id] %in% groups_oi]

  }

  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #

  coldata_df = SummarizedExperiment::colData(sce_subset)
  if(! "cell" %in% colnames(coldata_df)){
    coldata_df = coldata_df %>% data.frame() %>% tibble::rownames_to_column("cell")
  }

  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[gene_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[gene_oi,])) %>% dplyr::inner_join(coldata_df %>% tibble::as_tibble() )

  if(is.na(covariate_oi)){
    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE, groupOnX = TRUE)  +
      facet_grid(.~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the gene ",gene_oi, " in cell type ", celltype_oi, sep = ""))
  } else {

    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(covariate_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","covariate_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)

    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE, groupOnX = TRUE)  +
      facet_grid(covariate_oi~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the gene ",gene_oi, " in cell type ", celltype_oi, sep = ""))
  }

  return(p_violin)

}
#' @title make_DEgene_dotplot_pseudobulk_covariate
#'
#' @description \code{make_DEgene_dotplot_pseudobulk_covariate}  Visualize the scaled pseudobulk expression of DE genes per sample, and compare the different groups. Genes in rows, samples in columns. Add for each sample to which covariate/batch it belongs.
#' @usage make_DEgene_dotplot_pseudobulk_covariate(genes_oi, celltype_info, abundance_data, celltype_oi, covariate_oi, groups_oi = NA)
#'
#' @inheritParams make_DEgene_dotplot_pseudobulk
#' @param covariate_oi Name of the covariate/batch that needs to be visualized for each sample
#'
#' @return Gene expression dotplot list: pseudobulk version and single-cell version
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = "batch"
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' min_cells = 10
#' abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells, covariates = NA)
#' muscat_output = muscat_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl)
#' celltype_oi = "Malignant"
#' DE_table = filter(muscat_output$celltype_de$celltype_de$de_output_tidy,
#'                  cluster_id == celltype_oi &
#'                    p_adj <= 0.05 & logFC >= 1)
#' DE_genes = unique(DE_table$gene)
#' dotplots = make_DEgene_dotplot_pseudobulk_covariate(genes_oi = DE_genes, celltype_info = muscat_output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi, covariate_oi = covariates)
#'}
#'
#' @export
#'
make_DEgene_dotplot_pseudobulk_covariate = function(genes_oi, celltype_info, abundance_data, celltype_oi, covariate_oi, groups_oi = NA){

  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####

  keep_sender_receiver_values = c(1, 4)
  names(keep_sender_receiver_values) = c(FALSE, TRUE)

  abundance_data = abundance_data %>% dplyr::rename(sample = sample_id, celltype = celltype_id, group = group_id) %>% select(-n)
  plot_data = celltype_info$pb_df %>% inner_join(abundance_data, by = c("sample", "celltype"))
  plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = scaling_zscore(pb_sample)) %>% dplyr::ungroup()
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)

  plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)

  if(!is.na(groups_oi)){
    plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)

  }

  p1 = plot_data %>%
    ggplot(aes(sample, gene, color = scaled_gene_exprs, size = keep)) +
    geom_point() +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled pseudobulk\nexpression", size= "Celltype present") + xlab("") +ylab("") +
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p1 = p1 + custom_scale_fill

  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####

  frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi)
  plot_data = plot_data %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype"))
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)

  p2 = plot_data %>%
    # ggplot(aes(gene, sample , fill = scaled_gene_exprs)) +
    # geom_tile(color = "white") +
    ggplot(aes(sample, gene , color = scaled_gene_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    )  +
    # labs(color = "Scaled gene\navg expression")
    labs(color = "Scaled pseudobulk\nexpression", size= "Fraction of\nexpressing cells") + xlab("") +ylab("")
  max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p2 = p2 + custom_scale_fill

  ### covariate plot
  grouping_tbl = plot_data %>% distinct(sample, group, covariate_oi) %>% rename(covariate = covariate_oi)
  grouping_tbl_plot = grouping_tbl %>% mutate(covariate_ = paste0(" ",covariate_oi," "), mock = "")
  p_covariate = grouping_tbl_plot %>%
    ggplot(aes(sample, mock, fill = covariate)) +
    geom_tile(color = "black") +
    facet_grid(.~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1.50, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + scale_fill_brewer(palette = "Set2") + xlab("") + labs(fill=covariate_oi)

  p1 = patchwork::wrap_plots(
    p_covariate, p1,
    guides = "collect", nrow = 2,
    heights = c(1, plot_data$gene %>% unique() %>% length())
  )
  p2 = patchwork::wrap_plots(
    p_covariate, p2,
    guides = "collect", nrow = 2,
    heights = c(1, plot_data$gene %>% unique() %>% length())
  )

  return(list(pseudobulk_plot = p1, singlecell_plot = p2))

}

