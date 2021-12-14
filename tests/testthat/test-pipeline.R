context("muscat wrapper pipeline")
test_that("Pipeline for all-vs-all analysis works & plotting functions work", {
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  covariates = NA
  contrasts_oi = c("'High-Low','Low-High'")
  contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
  abundance_output = get_abundance_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = 10, covariates = NA)
  expect_type(abundance_output,"list")
  output = muscat_analysis(
       sce = sce,
       celltype_id = celltype_id,
       sample_id = sample_id,
       group_id = group_id,
       covariates = covariates,
       contrasts_oi = contrasts_oi,
       contrast_tbl = contrast_tbl)
  expect_type(output,"list")

  # test plotting functions
  group_oi = "High"
  celltype_oi = "Malignant"

  DE_table = filter(output$celltype_de$celltype_de$de_output_tidy,
                    cluster_id == celltype_oi &
                      p_adj <= 0.05 & logFC >= 1 & contrast == "High-Low")
  DE_genes = unique(DE_table$gene)

  dotplots = make_DEgene_dotplot_pseudobulk(genes_oi = DE_genes, celltype_info = output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi)
  expect_true("ggplot" %in% class(dotplots$pseudobulk_plot))
  expect_true("ggplot" %in% class(dotplots$singlecell_plot))

  dotplots_reversed = make_DEgene_dotplot_pseudobulk_reversed(genes_oi = DE_genes, celltype_info = output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi)
  expect_true("ggplot" %in% class(dotplots_reversed$pseudobulk_plot))
  expect_true("ggplot" %in% class(dotplots_reversed$singlecell_plot))

  gene_oi = DE_genes[1]
  violin_plot = make_DEgene_violin_plot(sce = sce, gene_oi = gene_oi, celltype_oi = celltype_oi, group_id = group_id, sample_id = sample_id, celltype_id = celltype_id)
  expect_true("ggplot" %in% class(violin_plot))

  # for coming calculations: reduce running time by having only one contrast of interest
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))

  # test whether input checks are stringent enough: celltype_id, sample_id, group_id
  SummarizedExperiment::colData(sce)$`Test Celltype` = SummarizedExperiment::colData(sce)$celltype
  SummarizedExperiment::colData(sce)$`Test Group` = SummarizedExperiment::colData(sce)$pEMT
  SummarizedExperiment::colData(sce)$`Test Sample` = SummarizedExperiment::colData(sce)$tumor
  SummarizedExperiment::colData(sce)$false_celltype = SummarizedExperiment::colData(sce)$seurat_clusters %>% as.double()

  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = "Test Celltype",
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = "Test Sample",
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = "Test Group",
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))

  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = "Valencia",
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id =  "Valencia",
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id =  "Valencia",
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = "false_celltype",
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))
  # test whether input checks are stringent enough: covariates
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = "dataset",
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  ))

  # test whether input checks are stringent enough: contrast definition and contrast tbl
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = c("'High-Medium'"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = c("High-Low","Low-High"),
    contrast_tbl = contrast_tbl
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Medium","Low-High"), group = c("High","Low"))
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("Medium","Low"))
  ))
  expect_warning(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = c("'High-Low','High-(Low+Low+Low)/3'"),
    contrast_tbl = tibble(contrast = c("High-Low","High-(Low+Low+Low)/3"), group = c("High","High"))
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    assay_oi_sce = "Spatial"
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    assay_oi_pb = "Spatial"
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    de_method_oi = "EdgePython"
  ))
  expect_error(muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl,
    min_cells = "0"
  ))


})
test_that("Pipeline with wrapper function works - while correcting for batch effects", {
  sample_id = "tumor"
  group_id = "pEMT"
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  covariates = "batch"
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))

  abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells = 10, covariates = covariates)
  expect_type(abundance_output,"list")

  output = muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  )
  expect_type(output,"list")

  # test batch effect indicated plots
  gene_oi = "PTHLH"
  celltype_oi = "Malignant"

  violin_plot = make_DEgene_violin_plot(sce = sce, gene_oi = gene_oi, celltype_oi = celltype_oi, group_id = group_id, sample_id = sample_id, celltype_id = celltype_id, covariate_oi = covariates)
  expect_true("ggplot" %in% class(violin_plot))

  DE_table = filter(output$celltype_de$celltype_de$de_output_tidy,
                    cluster_id == celltype_oi &
                      p_adj <= 0.05 & logFC >= 1 & contrast == "High-Low")
  DE_genes = unique(DE_table$gene)

  dotplots = make_DEgene_dotplot_pseudobulk_covariate(genes_oi = DE_genes, celltype_info = output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi, covariate_oi = covariates)
  expect_true("ggplot" %in% class(dotplots$pseudobulk_plot))
  expect_true("ggplot" %in% class(dotplots$singlecell_plot))

})

test_that("Empirical null procedure works", {
  sample_id = "tumor"
  group_id = "pEMT"
  sample_id = "tumor"
  group_id = "pEMT"
  celltype_id = "celltype"
  covariates = "batch"
  contrasts_oi = c("'High-Low'")
  contrast_tbl = tibble(contrast = c("High-Low"), group = c("High"))

  abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells = 10, covariates = covariates)
  expect_type(abundance_output,"list")

  output = muscat_analysis(
    sce = sce,
    celltype_id = celltype_id,
    sample_id = sample_id,
    group_id = group_id,
    covariates = covariates,
    contrasts_oi = contrasts_oi,
    contrast_tbl = contrast_tbl
  )
  expect_type(output,"list")

  empirical_pval = TRUE
  if(empirical_pval == TRUE){
    DE_info_emp = get_empirical_pvals(output$celltype_de$celltype_de$de_output_tidy)
  }

  # test batch effect indicated plots
  gene_oi = "PTHLH"
  celltype_oi = "Malignant"

  violin_plot = make_DEgene_violin_plot(sce = sce, gene_oi = gene_oi, celltype_oi = celltype_oi, group_id = group_id, sample_id = sample_id, celltype_id = celltype_id, covariate_oi = covariates)
  expect_true("ggplot" %in% class(violin_plot))

  DE_table = filter(DE_info_emp$de_output_tidy_emp,
                    cluster_id == celltype_oi &
                    p_emp <= 0.05 & logFC >= 1 & contrast == "High-Low")
  DE_genes = unique(DE_table$gene)

  dotplots = make_DEgene_dotplot_pseudobulk_covariate(genes_oi = DE_genes, celltype_info = output$celltype_info, abundance_data = abundance_output$abundance_data, celltype_oi = celltype_oi, covariate_oi = covariates)
  expect_true("ggplot" %in% class(dotplots$pseudobulk_plot))
  expect_true("ggplot" %in% class(dotplots$singlecell_plot))



})
