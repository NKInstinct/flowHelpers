#' Bunch of small helper functions for the different exported functions here.
#'
#' Might move these to individual files later on if this gets out of hand, but
#' for now, all the helper functions can just live here without too many issues.
#' This way too, I can work on unit testing these first in the same file and
#' then move on to the harder ones later.
#'
#' @noRd

readFS <- function(FCSDirectory,
                   pattern,
                   recursive,
                   ncdf){
  if(ncdf == TRUE){
    fs <- ncdfFlow::read.ncdfFlowSet(files = list.files(FCSDirectory,
                                                        pattern = pattern,
                                                        full.names = TRUE,
                                                        recursive = recursive))
  } else{
    fs <- flowCore::read.flowSet(files = list.files(FCSDirectory,
                                                    pattern = pattern,
                                                    full.names = TRUE,
                                                    recursive = recursive))
  }

  return(fs)

}

applyComp <- function(fs, comp){

  if(comp == "acquisition"){
    compObj <- flowCore::spillover(fs[[1]])
    compObj <- compObj[[1]]
  } else if(stringr::str_detect(comp, ".csv$")){
    comp.mat <- readr::read_csv(comp)
    comp.mat <- tibble::column_to_rownames(comp.mat, var = "X1")
    comp.mat <- as.matrix(comp.mat)
    compObj <- flowCore::compensation(comp.mat)
  } else{
    stop("If not FALSE, 'comp' must equal 'acquisition' or be a path to
           a csv file that defines a comp matrix.")
  }
  return(compObj)
}

getLinearChannels <- function(gs, linearChannels){
  if(length(linearChannels) == 1 & is.numeric(linearChannels)){
    linear <- linearChannels + 1
    chnls <- flowWorkspace::colnames(gs)[linear:(length(flowWorkspace::colnames(gs))-1)]
  } else if(is.character(linearChannels)){
    chnls <- flowWorkspace::colnames(gs)[!flowWorkspace::colnames(gs) %in% linearChannels]
  } else{
    stop("'linearChannels' must be an integer specifying the number of
         linear channels (FSC and SSC) active in the fcs (don't count TIME
         it is already accounted for here) or a character vector naming
         each channel you don't want transformed.")
  }

  return(chnls)
}

applyTransform <- function(gs, linearChannels, transform, arcsinh_coeff){

  chnls <- getLinearChannels(gs, linearChannels)

  if(transform == "biex"){
    transObj <- flowWorkspace::flowjo_biexp_trans()
  } else if(transform == "arcsinh") {
    transObj <- ggcyto::flowCore_asinht_trans(a=arcsinh_coeff)
  } else{
    stop("Currently, 'transform' must be either 'biex' or 'arcsinh'.")
  }

  transList <- flowWorkspace::transformerList(chnls, transObj)

  return(transList)
}

exportHandler <- function(export_opts,
                          gating_strategy,
                          parent_directory,
                          pattern){

  if(is.null(export_opts$exported_gate)){
    export_opts$exported_gate <- as.character(gating_strategy$filterId[length(gating_strategy$filterId)])
  }

  if(is.null(export_opts$fs_dir)){
    export_opts$fs_dir <- parent_directory
  }

  if(is.null(export_opts$fs_names)){
    export_opts$fs_names <- pattern
  }

  if(!is.null(export_opts$gs)){
    if(is.null(export_opts$gs_dir)){
      export_opts$gs_dir <- parent_directory
    }

    if(is.null(export_opts$gs_names)){
      export_opts$gs_names <- pattern
    }
  }
  return(export_opts)
}

exportFS <- function(gs, export_opts){

  if(class(gs) == "GatingSet"){
    fs <- flowWorkspace::gs_pop_get_data(gs,
                                         y = export_opts$exported_gate,
                                         inverse.transform = FALSE)

    flowCore::write.flowSet(fs,
                            outdir = paste(export_opts$fs_dir,
                                           export_opts$fs_names,
                                           "/", sep = ""))

    if(!is.null(export_opts$gs)){
      flowWorkspace::save_gs(gs,
                             paste(export_opts$gs_dir,
                                   export_opts$gs_names,
                                   ".gs", sep = ""))
    }

  } else if(is.list(gs) & class(gs[[1]]) == "GatingSet"){
    fs <- purrr::map(gs,
                     ~flowWorkspace::gs_pop_get_data(..1,
                                                     y = export_opts$exported_gate,
                                                     inverse.transform = FALSE))
    purrr::walk2(fs,
                 export_opts$fs_names,
                 ~flowCore::write.flowSet(..1,
                                          outdir = paste(export_opts$fs_dir,
                                                         ..2,
                                                         "/",
                                                         sep = "")))
    if(!is.null(export_opts$gs)){
      purrr::walk2(gs,
                   export_opts$gs_names,
                   ~flowWorkspace::save_gs(..1,
                                           paste(export_opts$gs_dir,
                                                 ..2,
                                                 ".gs",
                                                 sep = "")))
    }

  } else{
    stop("fs construction failed: gs is not a GatingSet or List of GatingSets")
  }
}

# Add metadata interactively ---------------------------------------------------

addMetadata <- function(df){
  data <- dplyr::distinct(df) |>
    as.data.frame()
  
  # Very weird error here where it works fine in debug mode but not in normal. Something about calling stopApp halts the whole function, not just the shiny app. Might not be solvable, so might need to split this whole thing into two pieces...
  result <- DataEditR::data_edit(x = data, 
                                 col_readonly = "sample_id", 
                                 col_edit = TRUE, 
                                 row_edit = FALSE,
                                 hide = FALSE)
    dplyr::right_join(df, by = "sample_id")
  
  return(result)
}
