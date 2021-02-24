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
