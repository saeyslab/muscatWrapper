#' @title make_sample_lr_prod_activity_covariate_plots
#'
#' @description \code{make_sample_lr_prod_activity_covariate_plots}  Visualize the scaled product of Ligand-Receptor (pseudobulk) expression per sample, and compare the different groups. In addition, show the NicheNet ligand activities in each celltype-celltype combination. On top of this summary plot, a heatmap indicates the covariate value for each displayed sample.
#' @usage make_sample_lr_prod_activity_covariate_plots(prioritization_tables, prioritized_tbl_oi, grouping_tbl, covariate_oi, widths = NULL, heights = NULL)
#'
#' @param prioritization_tables Output of `generate_prioritization_tables` or sublist in the output of `multi_nichenet_analysis`
#' @param prioritized_tbl_oi Subset of `prioritization_tables$group_prioritization_tbl`: the ligand-receptor interactions shown in this subset will be visualized: recommended to consider the top n LR interactions of a group of interest, based on the prioritization_score (eg n = 50; see vignettes for examples).
#' @param widths Vector of 3 elements: Width of the LR exprs product panel,  width of the scaled ligand activity panel, width of the ligand activity panel. Default NULL: automatically defined based on nr of samples vs nr of group-celltype combinations. If manual change: example format: c(5,1,1)
#' @param grouping_tbl Data frame linking the sample_id, group_id and covariate_oi
#' @param covariate_oi Name of the covariate/batch that needs to be visualized for each sample
#' @param heights Vector of 2 elements: Height of the covariate panel and height of the ligand-receptor prod+activity panel. Default NULL: automatically defined based on the nr of Ligand-Receptor pairs. If manual change: example format: c(1,5)
#'
#' @return Ligand-Receptor Expression Product & Ligand Activities Dotplot/Heatmap, complemented with a heatmap indicating the covariate/batch of interest
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      lr_network = lr_network,
#'      ligand_target_matrix = ligand_target_matrix,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl,
#'      sender_celltype_separate = FALSE
#'      )
#' group_oi = "High"
#' covariate_oi = "batch"
#' prioritized_tbl_oi = output$prioritization_tables$group_prioritization_tbl %>% filter(fraction_expressing_ligand_receptor > 0) %>% filter(group == group_oi) %>% top_n(50, prioritization_score)
#' plot_oi = make_sample_lr_prod_activity_covariate_plots(output$prioritization_tables, prioritized_tbl_oi, output$grouping_tbl, covariate_oi = covariate_oi)
#' plot_oi
#' }
#'
#' @export
#'
make_sample_lr_prod_activity_covariate_plots = function(prioritization_tables, prioritized_tbl_oi, grouping_tbl, covariate_oi, widths = NULL, heights = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  sample_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_celltype = paste(sender, celltype, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))   %>%  dplyr::arrange(celltype) %>% dplyr::group_by(celltype) %>%  dplyr::arrange(sender, .by_group = TRUE)
  sample_data = sample_data %>% dplyr::mutate(sender_celltype = factor(sender_celltype, levels = sample_data$sender_celltype %>% unique()))

  group_data = prioritization_tables$group_prioritization_tbl %>% dplyr::mutate(sender_celltype = paste(sender, celltype, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, celltype, sender_celltype, lr_interaction, group, ligand_receptor_lfc_avg, activity, activity_scaled, fraction_ligand_group, prioritization_score, scaled_avg_exprs_ligand) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(celltype) %>% dplyr::group_by(celltype) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data = group_data %>% dplyr::mutate(sender_celltype = factor(sender_celltype, levels = group_data$sender_celltype %>% unique()))

  p1 = sample_data %>%
    ggplot(aes(sample, lr_interaction, color = scaled_LR_pb_prod, size = scaled_LR_frac)) +
    geom_point() +
    facet_grid(sender_celltype~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand-Receptor expression in samples\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled L-R\npseudobulk exprs product", size= "Scaled L-R\navg exprs fraction product")
  max_lfc = abs(sample_data$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p1 = p1 + custom_scale_fill


  p2 = group_data %>%
    # ggplot(aes(celltype, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(celltype, lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_celltype~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in celltype cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
      # axis.text.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity in celltype")
  max_activity = abs(group_data$activity_scaled) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.40, 0.50, 0.60, 0.70, 0.825, 1),  limits = c(-1*max_activity, max_activity))

  p2 = p2 + custom_scale_fill

  p3 = group_data %>%
    ggplot(aes(celltype, lr_interaction, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_celltype~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in celltype cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Ligand\nActivity in celltype")
  max_activity = (group_data$activity) %>% max()
  min_activity = (group_data$activity) %>% min()
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges") %>% .[-7]),values = c(0, 0.250, 0.5550, 0.675, 0.80, 0.925, 1),  limits = c(min_activity-0.01, max_activity))

  p3 = p3 + custom_scale_fill

  grouping_tbl_plot = grouping_tbl %>% mutate(covariate_ = paste0(" ",covariate_oi," "), mock = "", covariate = grouping_tbl[[covariate_oi]])
  p_covariate = grouping_tbl_plot %>%
    ggplot(aes(sample, mock, fill = covariate)) +
    geom_tile(color = "black") +
    facet_grid(covariate_~group, scales = "free", space = "free", switch = "y") +
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
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + scale_fill_brewer(palette = "Set2") + xlab("")

  if(!is.null(widths)){
    design <- "D##
               ABC"
    p = patchwork::wrap_plots(
      A = p1, B = p2, C= p3, D = p_covariate,
      guides = "collect", design = design,
      widths = widths,
      heights = heights
    ) + patchwork::plot_layout(design = design)
  } else {
    design <- "D##
               ABC"
    p = patchwork::wrap_plots(
      A = p1, B = p2, C= p3, D = p_covariate,
      guides = "collect", design = design,
      widths = c(sample_data$sample %>% unique() %>% length(), sample_data$celltype %>% unique() %>% length(), sample_data$celltype %>% unique() %>% length()),
      heights = c(1, group_data$lr_interaction %>% unique() %>% length())
    ) + patchwork::plot_layout(design = design)
  }


  return(p)

}
#' @title make_sample_target_plots
#'
#' @description \code{make_sample_target_plots}  Heatmap/Dotplot of scaled target gene expression per sample, compared between the groups of interest. (samples in columns, genes in rows)
#' @usage make_sample_target_plots(celltype_info, targets_oi, celltype_oi, grouping_tbl)
#'
#' @param celltype_info `celltype_info` or `celltype_info` slots from the output of `multi_nichenet_analysis`. Also: part of the output of `get_abundance_expression_info_separate`
#' @param targets_oi Character vector of genes of which the expression should be visualized
#' @param celltype_oi Name of the celltype cell type of interest for which expression of genes should be visualized
#' @param grouping_tbl Data frame linking the sample_id, group_id and covariate_oi
#'
#' @return Heatmap/Dotplot of scaled target gene expression (samples in columns, genes in rows)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nichenetr scaling_zscore
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      lr_network = lr_network,
#'      ligand_target_matrix = ligand_target_matrix,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl,
#'      sender_celltype_separate = FALSE
#'      )
#' group_oi = "High"
#' celltype_oi = "Malignant"
#' targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(celltype == celltype_oi) %>% pull(gene) %>% unique()
#' p_target = make_sample_target_plots(celltype_info = output$celltype_info, targets_oi, celltype_oi, output$grouping_tbl)
#' }
#'
#' @export
#'
make_sample_target_plots = function(celltype_info, targets_oi, celltype_oi, grouping_tbl){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  pb_df =  celltype_info$pb_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% celltype_oi)
  frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% celltype_oi)

  filtered_data = pb_df %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype")) %>% dplyr::inner_join(grouping_tbl, by = "sample")

  filtered_data = filtered_data %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled_target_exprs = nichenetr::scaling_zscore(pb_sample), scaled_target_frac = nichenetr::scaling_zscore(fraction_sample)) %>% dplyr::ungroup()
  filtered_data$gene = factor(filtered_data$gene, levels=targets_oi)

  p1 = filtered_data %>%
    ggplot(aes(sample, gene, color = scaled_target_exprs, size = fraction_sample)) +
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
      panel.spacing.x = unit(2.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 11, color = "black", face = "bold"),
      strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled target\navg expression", size= "Target\n exprs fraction")
  max_lfc = abs(filtered_data$scaled_target_exprs) %>% max()
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  p1 = p1 + custom_scale_fill
  return(p1)
}

#' @title make_sample_target_plots_reversed
#'
#' @description \code{make_sample_target_plots_reversed}  Heatmap/Dotplot of scaled target gene expression per sample, compared between the groups of interest. (genes in columns, samples in rows)
#' @usage make_sample_target_plots_reversed(celltype_info, targets_oi, celltype_oi, grouping_tbl)
#'
#' @inheritParams make_sample_target_plots
#'
#' @return Heatmap/Dotplot of scaled target gene expression (samples in rows, genes in columns)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nichenetr scaling_zscore
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      lr_network = lr_network,
#'      ligand_target_matrix = ligand_target_matrix,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl,
#'      sender_celltype_separate = FALSE
#'      )
#' group_oi = "High"
#' celltype_oi = "Malignant"
#' targets_oi = output$ligand_activities_targets_DEgenes$de_genes_df %>% inner_join(contrast_tbl) %>% filter(group == group_oi) %>% arrange(p_val) %>% filter(celltype == celltype_oi) %>% pull(gene) %>% unique()
#' p_target = make_sample_target_plots_reversed(celltype_info = output$celltype_info, targets_oi, celltype_oi, output$grouping_tbl)
#' }
#'
#' @export
#'
make_sample_target_plots_reversed = function(celltype_info, targets_oi, celltype_oi, grouping_tbl){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  pb_df =  celltype_info$pb_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% celltype_oi)
  frq_df =  celltype_info$frq_df %>% dplyr::filter(gene %in% targets_oi & celltype %in% celltype_oi)

  filtered_data = pb_df %>% dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype")) %>% dplyr::inner_join(grouping_tbl, by = "sample")

  filtered_data = filtered_data %>% dplyr::group_by(gene) %>% dplyr::mutate(scaled_target_exprs = nichenetr::scaling_zscore(pb_sample), scaled_target_frac = nichenetr::scaling_zscore(fraction_sample)) %>% dplyr::ungroup()

  filtered_data$gene = factor(filtered_data$gene, levels=targets_oi)

  p1 = filtered_data %>%
    # ggplot(aes(gene, sample , fill = scaled_target_exprs)) +
    # geom_tile(color = "white") +
    ggplot(aes(gene, sample , color = scaled_target_exprs, size = fraction_sample)) +
    geom_point() +
    facet_grid(group~., scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0),
      axis.text.x = element_text(face = "italic", size = 9, angle = 90,hjust = 0),
      axis.text.y = element_text(size = 9),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(2.5, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 9, color = "black", face = "bold"),
      strip.text.y = element_text(size = 11, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) +
    # labs(color = "Scaled target\navg expression")
    labs(color = "Scaled target\nPseudobulk expression", size= "Target\n exprs fraction")
  max_lfc = abs(filtered_data$scaled_target_exprs) %>% max()
  # custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  plot = p1 + custom_scale_fill


  return( plot  )
}

##' @title make_featureplot
#'
#' @description \code{make_featureplot}  Make a panel of UMAP plots showing the expression of a gene of interest in the group of interest vs other groups
#' @usage make_featureplot(sce_subset_oi, sce_subset_bg, title_umap, gene_oi, group_oi, background_groups, group_id)
#'
#' @param sce_subset_oi SingleCellExperiment object containing the cells of the celltype and group of interest
#' @param sce_subset_bg SingleCellExperiment object containing the background cells (cells not of interest)
#' @param title_umap Character vector: title of the umap
#' @param gene_oi Character vector: name of the gene of interest
#' @param group_oi Character vector: name of the group of interest
#' @param background_groups Character vector: names of the background group(s)
#' @param group_id Name of the colData(sce) column in which the id of the group is defined
#'
#' @return UMAP plots showing the expression of gene of interest, compared between a group of interest and background groups
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom scater plotReducedDim
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' p = make_featureplot(sce, sce, title_umap = "title", gene_oi = "TNF", group_oi = "High", background_groups = "Low", group_id = "pEMT")
#' }
#'
#' @export
#'
make_featureplot = function(sce_subset_oi, sce_subset_bg, title_umap, gene_oi, group_oi, background_groups, group_id){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  p_dim = scater::plotReducedDim(sce_subset_oi,  "UMAP", colour_by = group_id) + ggtitle(title_umap) + theme(title = element_text(face = "bold"))

  p_oi = scater::plotReducedDim(sce_subset_oi,  "UMAP", colour_by = gene_oi)
  p_bg = scater::plotReducedDim(sce_subset_bg,  "UMAP", colour_by = gene_oi)

  p_oi = p_oi + ggtitle(paste(gene_oi, group_oi, sep = " in ")) + theme(title = element_text(face = "bold"))
  p_bg = p_bg + ggtitle(paste(gene_oi, background_groups %>% paste0(collapse = " & "), sep = " in ")) + theme(title = element_text(face = "bold"))

  wrapped_plots = patchwork::wrap_plots(p_dim,
                                        p_oi,
                                        p_bg,
                                        ncol = 3,guides = "collect")
}

#' @title make_target_violin_plot
#'
#' @description \code{make_target_violin_plot}  Violin plot of a target gene of interest: per sample, and samples are grouped per group
#' @usage make_target_violin_plot(sce, target_oi, celltype_oi, group_oi, group_id, sample_id, celltype_id, covariate_oi = NA, background_groups = NULL)
#'
#' @param target_oi Name of the gene of interest
#' @param sample_id Name of the colData(sce) column in which the id of the sample is defined
#' @param celltype_oi  Character vector with the names of the celltype cell type of interest
#' @param covariate_oi Name of a covariate of interest based on which the visualization will be split. Default: NA: no covariate.
#' @param group_oi Character vector of name of the group of interest
#' @param celltypes_oi Character vector with the names of the celltype cell types of interest
#' @param background_groups Default NULL: all groups in the group_id metadata column will be chosen. But user can fill in a character vector with the names of all gruops of interest.
#' @inheritParams muscat_analysis
#'
#' @return ggplot object: Violin plot of a target gene of interest: per sample, and samples are grouped per group
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
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      lr_network = lr_network,
#'      ligand_target_matrix = ligand_target_matrix,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl,
#'      sender_celltype_separate = FALSE
#'      )
#' celltype_oi = "Malignant"
#' group_oi = "High"
#' target_oi = "RAB31"
#' p_violin_target = make_target_violin_plot(sce = sce, target_oi = target_oi, celltype_oi = celltype_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id = celltype_id)
#' }
#' @export
#'
make_target_violin_plot = function(sce, target_oi, celltype_oi, group_oi, group_id, sample_id, celltype_id, covariate_oi = NA, background_groups = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }

  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% celltype_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]

  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #

  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[target_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[target_oi,])) %>% dplyr::inner_join(SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble())

  if(is.na(covariate_oi)){
    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  +
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
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the target ",target_oi, " in celltype cell type ", celltype_oi, sep = ""))
  } else {

    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(covariate_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","covariate_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)

    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  +
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
      ) + scale_color_brewer(palette = "Set2") + ggtitle(paste("Expression of the target ",target_oi, " in celltype cell type ", celltype_oi, sep = ""))
  }

  return(p_violin)

}

#' @title make_target_feature_plot
#'
#' @description \code{make_target_feature_plot}  UMAP showing target gene expression in groups and celltypes of interest vs other groups and cell types.
#' @usage make_target_feature_plot(sce, target_oi, group_oi, group_id, celltype_id, celltypes_oi, background_groups = NULL)
#'
#' @param target_oi Character vector of the name of the target gene of interest
#' @param celltype_oi  Character vector with the names of the celltype cell type of interest
#' @param covariate_oi Name of a covariate of interest based on which the visualization will be split. Default: NA: no covariate.
#' @param group_oi Character vector of name of the group of interest
#' @param celltypes_oi Character vector with the names of the celltype cell types of interest
#' @param background_groups Default NULL: all groups in the group_id metadata column will be chosen. But user can fill in a character vector with the names of all gruops of interest.
#' @inheritParams muscat_analysis
#'
#' @return UMAP showing target gene expression in groups and celltypes of interest vs other groups and cell types.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom SummarizedExperiment colData

#' @importFrom generics setdiff
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' covariates = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#' output = multi_nichenet_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      lr_network = lr_network,
#'      ligand_target_matrix = ligand_target_matrix,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl,
#'      sender_celltype_separate = FALSE
#'      )
#' celltype_oi = "Malignant"
#' group_oi = "High"
#' target_oi = "RAB31"
#' make_target_feature_plot(sce = sce, target_oi = target_oi, group_oi = group_oi, group_id = group_id, celltype_id = celltype_id, celltypes_oi = c("Malignant","myofibroblast","CAF"))
#' }
#'
#' @export
#'
make_target_feature_plot = function(sce, target_oi, group_oi, group_id, celltype_id, celltypes_oi, background_groups = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }

  # subset celltype - Nebulosa
  sce_subset_oi = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% group_oi]
  sce_subset_bg = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% background_groups]

  sce_subset_oi =  sce_subset_oi[, SummarizedExperiment::colData(sce_subset_oi)[,celltype_id] %in% celltypes_oi]
  sce_subset_bg =  sce_subset_bg[, SummarizedExperiment::colData(sce_subset_bg)[,celltype_id] %in% celltypes_oi]

  celltype_plots_feature = make_featureplot(sce_subset_oi, sce_subset_bg, "celltype UMAP", target_oi, group_oi, background_groups, group_id)

  p_feature = celltype_plots_feature

  return(p_feature)
}
