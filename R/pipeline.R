#' @title muscat_analysis
#'
#' @description \code{muscat_analysis}  Perform a multi-sample multi-condition DE analysis with the pseudobulk approach implemented in muscat.
#' @usage muscat_analysis(
#' sce, celltype_id, sample_id, group_id, covariates, contrasts_oi, contrast_tbl, assay_oi_pb ="counts", fun_oi_pb = "sum", de_method_oi = "edgeR", min_cells = 10, verbose = FALSE
#' )
#'
#' @param sample_id Name of the meta data column that indicates from which sample/patient a cell comes from (in sce)
#' @param group_id Name of the meta data column that indicates from which group/condition a cell comes from (in sce)
#' @param covariates NA if no covariates should be corrected for. If there should be corrected for covariates, this argument should be the name(s) of the columns in the meta data that indicate the covariate(s).
#' @param contrasts_oi String indicating the contrasts of interest (= which groups/conditions will be compared) for the differential expression and MultiNicheNet analysis.
#' We will demonstrate here a few examples to indicate how to write this. Check the limma package manuals for more information about defining design matrices and contrasts for differential expression analysis.\cr
#' If wanting to compare group A vs B: `contrasts_oi = c("'A-B'")` \cr
#' If wanting to compare group A vs B & B vs A: `contrasts_oi = c("'A-B','B-A'")` \cr
#' If wanting to compare group A vs B & A vs C & A vs D: `contrasts_oi = c("'A-B','A-C', 'A-D'")` \cr
#' If wanting to compare group A vs B and C: `contrasts_oi = c("'A-(B+C)/2'")` \cr
#' If wanting to compare group A vs B, C and D: `contrasts_oi = c("'A-(B+C+D)/3'")` \cr
#' If wanting to compare group A vs B, C and D & B vs A,C,D: `contrasts_oi = c("'A-(B+C+D)/3', 'B-(A+C+D)/3'")` \cr
#' Note that the groups A, B, ... should be present in the meta data column 'group_id'.
#' @param contrast_tbl Data frame providing names for each of the contrasts in contrasts_oi in the 'contrast' column, and the corresponding group of interest in the 'group' column. Entries in the 'group' column should thus be present in the group_id column in the metadata.
#' Example for `contrasts_oi = c("'A-(B+C+D)/3', 'B-(A+C+D)/3'")`:
#' `contrast_tbl = tibble(contrast = c("A-(B+C+D)/3","B-(A+C+D)/3"), group = c("A","B"))`
#' @param sce SingleCellExperiment object of the scRNAseq data of interest.
#' @param celltype_id Name of the column in the meta data of sce that indicates the cell type of a cell.
#' @param assay_oi_pb Indicates which information of the assay of interest should be used (counts, scaled data,...). Default: "counts". See `muscat::aggregateData`.
#' @param fun_oi_pb Indicates way of doing the pseudobulking. Default: "sum". See `muscat::aggregateData`.
#' @param de_method_oi Indicates the DE method that will be used after pseudobulking. Default: "edgeR". See `muscat::pbDS`.
#' @param min_cells Indicates the minimal number of cells that a sample should have to be considered in the DE analysis. Default: 10. See `muscat::pbDS`.
#' @param verbose Indicate which different steps of the pipeline are running or not. Default: FALSE.
#'
#' @return List containing information and output of the Muscat analysis.\cr
#' celltype_info: contains average expression value and fraction of each cell type - sample combination,\cr
#' celltype_de: contains output of the differential expression analysis,\cr
#' grouping_tbl: data frame showing the group per sample \cr
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom generics setdiff intersect union
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom tidyr spread
#' @importFrom SummarizedExperiment colData
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
#' output = muscat_analysis(
#'      sce = sce,
#'      celltype_id = celltype_id,
#'      sample_id = sample_id,
#'      group_id = group_id,
#'      covariates = covariates,
#'      contrasts_oi = contrasts_oi,
#'      contrast_tbl = contrast_tbl)
#' }
#'
#' @export
#'
muscat_analysis = function(sce, celltype_id,
                                            sample_id,
                                            group_id,
                                            covariates,
                                            contrasts_oi,
                                            contrast_tbl,
                                            assay_oi_pb ="counts",
                                            fun_oi_pb = "sum",
                                            de_method_oi = "edgeR",
                                            min_cells = 10,
                                            verbose = FALSE){


  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  # input checks
  if (class(sce) != "SingleCellExperiment") {
    stop("sce should be a SingleCellExperiment object")
  }

  if (!celltype_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("celltype_id should be a column name in the metadata dataframe of sce")
  }
  if (celltype_id != make.names(celltype_id)) {
    stop("celltype_id should be a syntactically valid R name - check make.names")
  }
  if (!sample_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("sample_id should be a column name in the metadata dataframe of sce")
  }
  if (sample_id != make.names(sample_id)) {
    stop("sample_id should be a syntactically valid R name - check make.names")
  }
  if (!group_id %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("group_id should be a column name in the metadata dataframe of sce")
  }
  if (group_id != make.names(group_id)) {
    stop("group_id should be a syntactically valid R name - check make.names")
  }

  if(is.double(SummarizedExperiment::colData(sce)[,celltype_id])){
    stop("SummarizedExperiment::colData(sce)[,celltype_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce)[,group_id])){
    stop("SummarizedExperiment::colData(sce)[,group_id] should be a character vector or a factor")
  }
  if(is.double(SummarizedExperiment::colData(sce)[,sample_id])){
    stop("SummarizedExperiment::colData(sce)[,sample_id] should be a character vector or a factor")
  }

  # if some of these are factors, and not all levels have syntactically valid names - prompt to change this
  if(is.factor(SummarizedExperiment::colData(sce)[,celltype_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,celltype_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,celltype_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,celltype_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,celltype_id] should be a syntactically valid R names - see make.names")
    }
  }

  if(is.factor(SummarizedExperiment::colData(sce)[,group_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,group_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,group_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,group_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,group_id] should be a syntactically valid R names - see make.names")
    }
  }

  if(is.factor(SummarizedExperiment::colData(sce)[,sample_id])){
    is_make_names = levels(SummarizedExperiment::colData(sce)[,sample_id]) == make.names(levels(SummarizedExperiment::colData(sce)[,sample_id]))
    if(sum(is_make_names) != length(levels(SummarizedExperiment::colData(sce)[,sample_id]))){
      stop("The levels of the factor SummarizedExperiment::colData(sce)[,sample_id] should be a syntactically valid R names - see make.names")
    }
  }


  if(!is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector")
  }
  if(!is.data.frame(contrast_tbl)){
    stop("contrast_tbl should be a data frame / tibble")
  }
  # conditions of interest in the contrast should be present in the in the group column of the metadata
  groups_oi = SummarizedExperiment::colData(sce)[,group_id] %>% unique()
  conditions_oi = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    # stringr::str_split("[:digit:]") %>% unlist() %>% unique() %>%
    stringr::str_split("\\)") %>% unlist() %>% unique() %>%
    stringr::str_split("\\(") %>% unlist() %>% unique() %>%
    stringr::str_split("-") %>% unlist() %>% unique() %>%
    stringr::str_split("\\+") %>% unlist() %>% unique() %>%
    stringr::str_split("\\*") %>% unlist() %>% unique() %>%
    stringr::str_split("\\/") %>% unlist() %>% unique() %>% generics::setdiff(c("",","," ,", ", ")) %>% unlist() %>% unique()
  conditions_oi = conditions_oi[is.na(suppressWarnings(as.numeric(conditions_oi)))]

  if(length(contrasts_oi) != 1 | !is.character(contrasts_oi)){
    stop("contrasts_oi should be a character vector of length 1. See the documentation of the function for having an idea of the right format of setting your contrasts.")
  }

  # conditions of interest in the contrast should be present in the in the contrast_tbl
  contrasts_oi_simplified = stringr::str_split(contrasts_oi, "'") %>% unlist() %>% unique() %>%
    stringr::str_split(",") %>% unlist() %>% unique() %>% generics::setdiff(c("",",")) %>% unlist() %>% unique()

  if (sum(conditions_oi %in% groups_oi) != length(conditions_oi)) {
    stop("conditions written in contrasts_oi should be in the condition-indicating column! This is not the case, which can lead to errors downstream.")
  }
  if (sum(contrasts_oi_simplified %in% unique(contrast_tbl$contrast)) != length(contrasts_oi_simplified)) {
    stop("conditions written in contrasts_oi should be in the contrast column of contrast_tbl column! This is not the case, which can lead to errors downstream.")
  }

  #
  groups_oi_contrast_tbl = contrast_tbl$group %>% unique()
  if(sum(groups_oi_contrast_tbl %in% groups_oi) != length(groups_oi_contrast_tbl)){
    stop("You have defined some groups in contrast_tbl$group that are not present SummarizedExperiment::colData(sce)[,group_id]. This will result in lack of information downstream. We recommend to change your metadata or this contrast_tbl appropriately.")
  }

  if(length(groups_oi_contrast_tbl) != length(contrast_tbl$group)){
    warning("According to your contrast_tbl, some of your contrasts will be assigned to the same group. This should not be a problem if this was intended, but be aware not to make mistakes in the further interpretation and plotting of the results.")
  }

  if(!is.na(covariates)){
    if (sum(covariates %in% colnames(SummarizedExperiment::colData(sce))) != length(covariates) ) {
      stop("covariates should be NA or all present as column name(s) in the metadata dataframe of sce")
    }
  }

  if(!is.character(assay_oi_pb)){
    stop("assay_oi_pb should be a character vector")
  } else {
    if(assay_oi_pb != "counts"){
      warning("are you sure you don't want to use the counts assay?")
    }
  }
  if(!is.character(fun_oi_pb)){
    stop("fun_oi_pb should be a character vector")
  }
  if(!is.character(de_method_oi)){
    stop("de_method_oi should be a character vector")
  }

  if(!is.double(min_cells)){
    stop("min_cells should be numeric")
  } else {
    if(min_cells <= 0) {
      warning("min_cells is now 0 or smaller. We recommend having a positive, non-zero value for this parameter")
    }
  }
  if(!is.logical(verbose)){
    stop("verbose should be TRUE or FALSE")
  }

  ### celltype abundance plots + Calculate expression information
  if(verbose == TRUE){
    print("Make diagnostic abundance plots + Calculate expression information")
  }

  expression_info = get_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, covariates = covariates)

  ### Perform the DE analysis ----------------------------------------------------------------

  if(verbose == TRUE){
    print("Calculate differential expression for all cell types")
  }

  DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells,
                          assay_oi_pb = assay_oi_pb,
                          fun_oi_pb = fun_oi_pb,
                          de_method_oi = de_method_oi)

  ### Remove types of information that we don't need anymore:

  metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

  if(!is.na(covariates)){
    grouping_tbl = metadata_combined[,c(sample_id, group_id, covariates)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group",covariates)
  } else {
    grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
    colnames(grouping_tbl) = c("sample","group")
  }

  muscat_output = list(
    celltype_info = expression_info,
    celltype_de = list(hist_pvals = DE_info$hist_pvals, celltype_de = DE_info$celltype_de)
  )

  return(muscat_output)
}

