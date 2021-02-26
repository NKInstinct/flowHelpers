# This script will allow you to quickly gate through a gating set according to
# any gatingStrat object you make, according to the following template.

# gatingStrat <- tribble(
#   ~filterId, ~dims, ~subset, ~coords,
#   "Lymphocytes", list("FSC-A", "SSC-A"), "root", lincoords,
#   "Single cells", list("FSC-A", "FSC-W"), "Lymphocytes", lincoords,
#   "Live Leukocytes", list("Viability", "CD45"), "Single cells", logcoords)


gs_apply_gating_strategy <- function(gs = NULL,
                                     gating_strategy,
                                     parent_directory = NULL,
                                     pattern = ".fcs$",
                                     mode = "individual",
                                     reference_sample = 1,
                                     ...){
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

  if(class(gs) == "GatingSet"){
    # Apply each step of the gating strategy to the gs.
    purrr::pwalk(gating_strategy, flowGate::gs_gate_interactive, gs = gs)
  } else if(is.list(gs)){
    if(mode == "individual"){
      purrr::walk(gs,
                 ~purrr::pwalk(gating_strategy,
                              flowGate::gs_gate_interactive,
                              gs = ..1))
    } else if(mode == "reference"){
      gsRef <- gs[[reference_sample]]
      gsNon <- gs[-reference_sample]

      purrr::pwalk(gating_strategy, flowGate::gs_gate_interactive, gs = gsRef)

      for(i in seq_along(gating_strategy$filterId)){
        purrr::walk(gsNon,
                   ~flowWorkspace::gs_pop_add(..1,
                                              gate = flowWorkspace::gh_pop_get_gate(gsRef[[1]],
                                                                                    gating_strategy$filterId[[i]]),
                                              parent = gating_strategy$subset[[i]]))
      }

      purrr::walk(gs, flowWorkspace::recompute)
    } else{
      stop("'mode' must be 'individual' or 'reference'")
    }
  } else{
    stop("'gs' must be a gatingset, a list of gatingsets, or a directory and pattern suitable to build one of the above")
  }

  return(gs)

}

# purrr::pmap(gating_strategy, flowGate::gs_gate_interactive, gs = gs)
#
# # Export should probably be a separate function so that it can be called apart from gating the data - might be the easiest way to solve this interpreting problem where
#
#     # fs <- flowWorkspace::gs_pop_get_data(gs,
#     #                                      y = exported_gate,
#     #                                      inverse.transform = TRUE)
#
#   # purrr::pmap(gating_strategy, flowGate::gs_gate_interactive, gs = gs)
#
#   if(is.null(export_opts$fs_dir)){
#     outDir <- parent_directory
#   } else{
#     outDir <- export_opts$fs_dir
#   }
#
#   if(is.null(export_opts$fs_names)){
#     outName <- pattern
#   } else{
#     outName <- export_opts$fs_names
#   }
#
#   flowCore::write.flowSet(fs, outdir = paste(outDir, outName, "/", sep = ""))
#
#   if(!is.null(export_opts$gs)){
#     if(is.null(export_opts$gs_dir)){
#       gsOutDir <- parent_directory
#     } else{
#       gsOutDir <- export_opts$gs_dir
#     }
#
#     if(is.nulll(export_opts$gs_names)){
#       gsOutName <- pattern
#     } else{
#       gsOutName <- export_opts$gs_names
#     }
#
#     save_gs(gs, paste(gsOutDir, gsOutName, ".gs", sep = ""))
#
#   }
#
#   return(NULL)
# }
