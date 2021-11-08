#' Create a GatingSet from a directory of fcs files
#'
#' This function takes a directory (or directories recursively) containing some
#' fcs files and reads them into the GatingSet format from flowWorkspace so that
#' they can be subsequentially gated on. The function can optionally read the
#' data in as ncdf files to save memory, can apply a compensation matrix from an
#' acquisition-defined slot or a provided matrix in .csv format, and can apply a
#' biexponential or arcsinh transform.
#'
#' @param FCSDirectory The parent directory containing the fcs files.
#' @param pattern a string/regexp that selects the files to read in. Defaults to
#'   all fcs files with ".fcs$".
#' @param autoqc Boolean specifying whether to run flowAI's auto QC pipeline on
#'   the generated fs before converting to gs. Defaults to FALSE.
#' @param recursive Boolean specifying whether the fcs files should be searched
#'   for within FCSDirectory recursively (i.e. are there subfolders?).
#' @param ncdf Boolean specifying whether the fcs files should be read in as
#'   ncdf files.
#' @param comp Either "acquisition" or the filepath to a comp matrix saved as a
#'   .csv file. Specify FALSE if you want no comp applied.
#' @param transform Either "biex" or "arcsinh" to apply the appropriate
#'   transform function to all channels not defined as linear (see below). Keep
#'   FALSE to apply no transformation (default, since I'm starting to move to
#'   scale_._transform instead of raw data transforms)
#' @param arcsinh_coeff passed to the arcsinh function - defaults to 150, which
#'   is considered appropriate for fluorescent cytometry data (use coefficient 1
#'   if you have cytof data instead).
#' @param linearChannels Either a number specifying how many FSC and SSC
#'   channels you have (don't include TIME as linear even though it is - this is
#'   taken care of already), or a character vector naming all of the channels
#'   you *don't* want transformed (i.e. the linear channels).
#'
#' @return a GatingSet object with the comp and transformation applied as
#'   specified.
#'
#' @examples
#'
#'path_to_fcs <- system.file("extdata", package = "flowHelpers")
#' gs <- gs_create(FCSDirectory = path_to_fcs,
#' pattern = ".fcs$",
#' recursive = FALSE,
#' ncdf = FALSE,
#' comp = "acquisition",
#' transform = "biex",
#' linearChannels = 6)
#'
#'
#' @export
gs_create <-  function(FCSDirectory,
                       pattern = ".fcs$",
                       autoqc = FALSE,
                       recursive = FALSE,
                       ncdf = FALSE,
                       comp = "acquisition",
                       transform = FALSE,
                       arcsinh_coeff = 150,
                       linearChannels = 2){

  # Read files from directory ---------------------------------------------

  fs <- readFS(FCSDirectory, pattern, recursive, ncdf)
  
  # Apply flowAI QC -------------------------------------------------------
  if(autoqc == TRUE){
    fs <- flowAI::flow_auto_qc(fs)
  }

  gs <- flowWorkspace::GatingSet(fs)


  # Apply comp from acq or matrix -----------------------------------------
  if(comp != FALSE){
    compObj <- applyComp(fs, comp)
    gs <- flowWorkspace::compensate(gs, compObj)
  }

  # Transform fluorescent channels ----------------------------------------

  if(transform != FALSE){
    transList <- applyTransform(gs, linearChannels, transform, arcsinh_coeff)
    gs <- flowWorkspace::transform(gs, transList)
  }

  return(gs)
}
