#' Convert GatingSet to Seurat and add metadata.
#' 
#' This function outputs a Seurat object based on a provided GatingSet and
#' (optionally) a node within the GatingSet to draw data from. The output is a
#' Seurat object, and the function prompts the next steps in completing the
#' Seurat analysis.
#' 
#' Note that as part of this call, gs_to_seurat will prompt for user input to
#' add metadata to the created Seurat object. Right-cliking on a cell in the
#' created shiny window will allow you to add columns to the datatable, which
#' can then be given colnames and values interactively. When finished,
#' gs_to_seurat will apply all of the created metadata columns to the Seurat
#' object as a separate metadata entry with the colname as a name - this can be
#' useful for downstream regression and labeling.
#' 
#' @param gs The GatingSet to take data from
#' @param node Character vector specifying a (nonredundant) node in the
#'   GatingSet to use for data extraction. Defaults to "root" (all events in
#'   GatingSet).
#' 
#' @export
gs_to_seurat <- function(gs,
                        node = "root"){
# Prepare SCE Inputs -----------------------------------------------------------
  message("Converting GatingSet to Seurat")
  
  fs <- flowWorkspace::gs_pop_get_data(gs, y = node, inverse.transform = TRUE)
  panel <- dplyr::left_join(
    tibble::enframe(flowCore::colnames(fs), 
                    name = NULL, 
                    value = "fcs_colname"),
    tibble::enframe(flowCore::markernames(fs), 
                    name = "fcs_colname", 
                    value = "antigen"),
    by = "fcs_colname") |>
    dplyr::mutate(marker_class = "none")
  metadata <- tibble::tibble(file_name = flowCore::sampleNames(fs),
                             sample_id = flowCore::sampleNames(fs),
                             condition = flowCore::sampleNames(fs))
  
  # Metadata is just a placeholder so CATALYST doesn't lose it. Whatever you put
  # here gets lost on conversion to Seurat and has to be restored manually later
  # on, so no use specifying anything useful at this stage.
  
  inputs <- list(fs = fs, panel = panel, metadata = metadata)
  
  sce <- CATALYST::prepData(inputs$fs,
                            inputs$panel,
                            inputs$metadata,
                            cofactor = 150,
                            md_cols = list(file = "file_name",
                                           sample = "sample_id",
                                           factors = "condition"),
                            FACS = TRUE)
  # Note - the above centers and scales the data, which is why we don't in the
  # subsequent Seurat steps (and why we invert the gs transforms above,
  # incidentally).

# Convert SCE to Seurat --------------------------------------------------------
  seurat <- Seurat::as.Seurat(sce, counts = "counts", data = "exprs")
  
  message("Requesting metadata")
  
# Fix Metadata -----------------------------------------------------------------
  meta <- seurat@meta.data$sample_id |>
    tibble::enframe(name = NULL, value = "sample_id")
  # Need to make a second df with the sample id and the added metadata and then
  # left_join them I think
  
  metaDisplay <- dplyr::distinct(meta) |>
    as.data.frame()
  
  metaEdit <- DataEditR::data_edit(x = metaDisplay, 
                                      col_readonly = "sample_id", 
                                      col_edit = TRUE, 
                                      row_edit = FALSE)

  meta <- dplyr::left_join(meta, metaEdit, by = "sample_id")
  
  # Should return a dataframe with any metadata columns added in
  if(length(meta) > 1){
    for(i in 2:length(meta)){
      seurat <- SeuratObject::AddMetaData(seurat,
                                          meta[[i]],
                                          col.name = colnames(meta[i]))
    }
  }
  
  message("Seurat Produced - please proceed with scaling (ScaleData, with 
          optional variable regression), PCA (runPCA), Elbow Plot to determine 
          optimal dimensions (ElbowPlot), and then FindNeighbors (specify dims), 
          FindClusters, and RunUMAP (specify dims)")
  
  return(seurat)

# Okay so something about Seurat (possibly the way it implements Futures???) is
# making it near impossible to build into a pipeline. All the commented-out code
# (both futures- and non-futures attempts at solving this) work PERFECTLY when
# run in debug mode but silently fail when run normally. Googling around gives
# nothing for R, but for other languages, this pattern (debug works but normal
# doesn't) indicates a multi-threading issue, where threads are resolving out of
# order and the program isn't written to catch them. I suspect that running the
# Seurat command interactively "gets around" this problem but calling it in a
# function (even in a wrapped function) doesn't work.
  
# For now, I think the best solution is to call this function "gs_to_seurat" and
# then have it print a helpful message at the end saying what to do next.
# Minimum is Scale, PCA, ElbowPlot, then Dimensional Reduction. Maybe as I'm
# getting those worked out, I'll learn that the problem is only with ScaleData,
# and then I can proceed with writing some more helper functions below.
  

  # # Normalize and dimension reduce -----------------------------------------------
  #   if(!is.null(regress)){
  #     if(!regress %in% colnames(meta)){
  #       stop("Factors to regress must be provided in the metadata table. Please
  #            check your spelling and try again.")
  #     }
  #   }
  # 
  #   
  #   if(is.null(cruftFeatures)){
  #     feat <- rownames(seurat)
  #   } else{
  #     feat <- rownames(seurat)[!str_detect(rownames(seurat),
  #                                          cruftFeatures)]
  #   }

  # future::plan(strategy = "sequential")
  # scaleAndPCA <- function(obj, feat, regress){
  #   y <- Seurat::ScaleData(obj,
  #                          features = feat,
  #                          do.scale = FALSE,
  #                          do.center = FALSE,
  #                          verbose = FALSE,
  #                          vars.to.regress = regress)
  #   y <- Seurat::RunPCA(y,
  #                       features = feat,
  #                       npcs = 20,
  #                       verbose = FALSE)
  #   
  #   return(y)
  # }
  # 
  # seurat <- scaleAndPCA(seurat, feat, regress)
  # 
  # options(future.globals.maxSize = Inf)
  # future::plan(strategy = "sequential")
  # seurat <- future::future({Seurat::ScaleData(seurat,
  #                             features = feat,
  #                             do.scale = TRUE,
  #                             do.center = TRUE,
  #                             verbose = FALSE,
  #                             vars.to.regress = regress)},
  #                          globals = c("seurat", "feat", "regress"))
  # 
  # seurat <- future::value(seurat)
  # 
  # seurat <- future::future({Seurat::RunPCA(seurat, features = feat, npcs = 20, verbose = FALSE)},
  #                          globals = c("seurat", "feat"))
  # 
  # seurat <- future::value(seurat)
  # 
  # return(seurat)

# 
# # Set Elbowplot ----------------------------------------------------------------
#   # Write shiny app here to show the Elbow plot and then ask for the number
#   # of dimensions
# 
# # Run cluster and UMAP ---------------------------------------------------------
# 
# # Name clusters ----------------------------------------------------------------
#   # Write shiny app here to show UMAP and top 4 markers 
#   # and then have a box for each cluster to name it
# 
# # Save Seurat ------------------------------------------------------------------
}