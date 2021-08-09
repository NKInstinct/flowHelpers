#' Pregate a gs for infinityFlow
#' 
#' Pass this function a gatingset or the instructions to build one (still in
#' development - might not work very well) and it will draw gates and export the
#' resulting gates as .fcs files for use in infinityFlow.
#' 
#' Note that a previous version of this function would allow lists of gs to be
#' passed directly into it. Unfortunately, this behaviour seemed to break every
#' time flowWorkspace was updated? Or something else, but in any case it would
#' regularly stop working and was a bit jank anyway. The better way to do this
#' is to explicitly pass the list of gatingsets into this function with map (or
#' more likely, imap and then use the names as the fs_names argument in the
#' export options).
#'
#' @param gs A GatingSet. If NULL, must specify a
#'   parent_directory containing .fcs files and a pattern to find them so
#'   gs_create can build a gs or list of gs.
#' @param gating_strategy A gating strategy tibble as required by
#'   flowGate::gs_apply_gating_strategy. See flowGate docs for more details.
#'   Note that for this function specially, there is currently no way to supply
#'   invariant but non-default params to gs_gate_iteractive, so if you have
#'   something you want done non-default, it must be specified in the
#'   gating_strategy even if it is invariant.
#' @param parent_directory A string passed to gs_create if gs is not specified
#'   to build a gs from .fcs files contained within.
#' @param pattern A string or list of strings containing patterns to pass to
#'   gs_create and build a gs or list of gs.
#' @param ... Additional params to pass to gs_create if desired.
#' @param export_opts A list allowing fine control of how the resulting gated
#'   data are exported. Default behaviour will export only an fcs file per
#'   sample in the parent directory (which must be specified in this case),
#'   named based on the pattern (which also must be specified in this case), and
#'   will use the last gate specified in gating_strategy as the export gate.
#'   Override this behaviour by passing a list here with the following named
#'   items:
#'       exported_gate: a string specifying the filterId of the gate to export.
#'       fs_dir: the directory to save the resulting fcs files.
#'       fs_names: the strategy to name the subfolders containing each exported
#'           fcs (specify one name per GatingSet used in the function call).
#'       gs: set to TRUE if you want to export the GatingSets as well.
#'       gs_dir: where to save the GatingSets
#'       gs_names: how to name the GatingSets (like fs_names)
#'
#' @return Nothing---instead, saves the facs data contained in the export gate
#'   as .fcs files according to the export options, and optionally the GatingSet
#'   object for further use.
#' @examples
#' if(interactive()){
#'     path_to_fcs <- system.file("extdata", package = "flowHelpers")
#'     gating_strategy <- readRDS(system.file("extdata",
#'                                            "gating_strategy.Rds",
#'                                            package = "flowHelpers"))
#'
#'     gs_pregate(gating_strategy = gating_strategy,
#'                parent_directory = path_to_fcs,
#'                pattern = "Vitro",
#'                export_opts = list(fs_dir = "../data/processed/examplePregate/",
#'                                   fs_names = "Example Vitro"))
#' }
#'
#' @export
gs_pregate <- function(gs = NULL,
                       gating_strategy,
                       parent_directory = NULL,
                       pattern = ".fcs$",
                       ...,
                       export_opts = list()){

  # Build gs from pattern if not supplied -------------------------

  if(is.null(gs)){
    if(length(pattern) == 1){
      gs <- gs_create(parent_directory, pattern = pattern, ...)
    } else if(length(pattern) > 1){
      gs <- purrr::map(pattern,
                       ~gs_create(FCSDirectory = parent_directory,
                                  pattern = ..1,
                                  ...))
    } else{
      stop("'pattern' and 'parent_directory' must be supplied if gs is not")
    }

  }

  # Apply gating strategy ------------------------------------------------------

  # note for this funnction, since I don't know how to handle multiple ... yet
  # and suspect it makes the args super messy anyway, I'm going to insist that
  # the gating_strategy simply contain all args you want non-default rather than
  # passing them to purrr later.

  flowGate::gs_apply_gating_strategy(gs, gating_strategy)


  # Setup Export Opts -------------------------------------------------------

  export_opts <- exportHandler(export_opts, gating_strategy, parent_directory, pattern)


  # Convert to fs & export --------------------------------------------------


  exportFS(gs, export_opts)

  return(NULL)
}
