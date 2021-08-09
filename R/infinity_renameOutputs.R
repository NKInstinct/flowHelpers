#' Rename the outputs of infinityFlow to something more helpful than concat.fcs
#' 
#' Annoyingly, all infinityFlow concatenated files are output as 'concat.fcs',
#' making it very difficult to use them in downstream analysis software like
#' Seurat or CATALYST. This function renames them the name of their
#' infinityOutput directory so that adding them to some downstream process
#' doesn't throw errors later.
#' 
#' @param parentDir a filepath to the directory containing both the
#'   infinityOutput and pregated directories.
#' @param pattern a regular expression that will be used to find the output
#'   files. The default of "\\.fcs" works if you have only output a single
#'   concatenated file, but you'll need to be more specific if your
#'   infinityOutput also includes non-concatenated .fcs files (or any other .fcs
#'   for that matter, since the directories must be searched recursively).
#' @param concat_tag_to_remove a regular expression showing what string to
#'   remove from the filename. The default will work if your infinityOutput only
#'   has a concatenated file in it.
#' @param prefix a string to add in front of the filename for each renamed file.
#'   Pasted directly on, so if you want a specific separater, add it here. The
#'   default (NULL) pastes no string. Adding a prefix is useful if you are going
#'   to be combining several batches of samples with similar naming strategies
#'   (WT1, WT2, etc).
#' @return Nothing - this modifies the filenames in place, so be careful! It's
#'   hard to undo this!
#' @importFrom stringr str_detect str_remove_all
#' 
#' 
#' @export
infinity_renameOutputs <- function(parentDir, 
                                   pattern = "\\.fcs$", 
                                   concat_tag_to_remove = "FCS/concatenated",
                                   prefix = NULL){
  
  if(!stringr::str_detect(parentDir, "/$")){
    parentDir <- paste(parentDir, "/", sep = "")
  }
  
  old_files <- list.files(paste(parentDir,
                                "infinityOutput",
                                sep = ""),
                          recursive = TRUE,
                          pattern = pattern,
                          full.names = TRUE)
  
  parent_paths <- dirname(old_files) |>
    stringr::str_remove_all(concat_tag_to_remove)
  
  sample_names <- list.files(paste(parentDir,
                                   "infinityOutput",
                                   sep = ""),
                             full.names = FALSE) |>
    paste(".fcs", sep = "")
  
  if(!is.null(prefix)){
    sample_names <- paste(prefix, sample_names, sep = "")
  }
  
  new_files <- paste(parent_paths, sample_names, sep = "")
  
  file.rename(old_files, new_files)
}