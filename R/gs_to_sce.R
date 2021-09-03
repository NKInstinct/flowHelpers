#' Convert GatingSet to SCE with option to run UMAP as well
#' 
#' This function is very similar to gs_to_Seurat, but rather than shoehorning
#' all the flow data into Seurat's structures, it falls back to using just 
#' (tidy)SingleCellExperiment packages and constructors.
#' 
#' Note that if you want to be able to interact with the sce in tidy format, you
#' should attach tidySingleCellExperiment yourself!
#' 
#' @param gs The GatingSet to draw the data from
#' @param node A string specifying which gate to take data from. Defaults to
#'   "root" (all acquired events)
#' @param condition A vector of strings specifying which group each sample belongs to.
#'   Must be of the same length as the number of samples in gs, and must be in
#'   the same order as the samples in gs (can check this with sampleNames(gs)).
#'   Defaults to just the sampleNames themselves, which is not ideal so you
#'   should probably specify something here!
#' @param do_UMAP A boolean specifying whether to run UMAP dimensional reduction
#'   or not. This is the whole point of turning flow data into an SCE, but it
#'   can be computationally intense for even moderately sized experiments so you
#'   need to opt-in.
#' @importFrom magrittr "%$%"
#'   
#' @export
gs_to_sce <- function(gs, 
                      node = "root", 
                      condition = Biobase::sampleNames(gs),
                      do_UMAP = FALSE){
  fs <- flowWorkspace::gs_pop_get_data(gs, 
                                       y = node, 
                                       inverse.transform = TRUE) |>
    flowWorkspace::cytoset_to_flowSet()
  
  md <- tibble::tibble(file_name = Biobase::sampleNames(gs),
                       condition = condition)
  
  sce <- CATALYST::prepData(fs, 
                            md = md, 
                            cofactor = 150, 
                            md_cols = list(file = "file_name", 
                                           id = "file_name", 
                                           factors = "condition"),
                            FACS = TRUE)
  
  sce <- sce |>
    scater::runPCA(exprs_values = "exprs")
  
  if(do_UMAP == TRUE){
    SingleCellExperiment::colLabels(sce) <- sce |>
      scran::buildSNNGraph(use.dimred="PCA") |>
      igraph::cluster_walktrap() %$%
      membership |>
      as.factor()
    
    sce <- sce |>
      scater::runUMAP(ncomponents = 2, exprs_values = "exprs")
  }
  
  return(sce)
}