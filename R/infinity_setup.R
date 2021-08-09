#' Prepare the infinity setup Rds
#' 
#' Run this command after you have pregated the infinity files to make an
#' infinity setup tibble that you can then pmap into infinityFlow so that each
#' sample gets correctly imputed. Note that for now, this version is insisting
#' on an intermediary results folder - although the bug where not specifying
#' intermediaries appended the files has been resolved, it's safest to just run
#' this with the intermediaries saved and then delete them later if you don't
#' need them.
#'
#' @param pregated_dir A path to the directory that contains your pregated .fcs
#'   files for infinity imputation.
#' @param infinity_targets A character vector with the identity of each infinity
#'   marker, file-by-file. Note that this needs to be in the order that R will
#'   read it into the system, which is often well A10, A11, A12, A1, A2... when
#'   acquiring fortessa data. If you aren't sure what the order will be, call
#'   (list.files) on one of the individual pregated sample folders and see what
#'   order they come out in.
#' @param infinity_isotypes A character vector specifying the name of each
#'   isotype target for each infinity target. The isotype string must be found
#'   in one of the entries in infinity_targets. Often, just pass "FMO" to this
#'   and all background correction will be done against the FMO sample.
#' @param create_output_dir Boolean. Should an infinityOutput directory be
#'   created and populated with the same subfolder structure as the pregated
#'   directory?
#' @param create_intermed_dir Boolean. Should an infinityIntermed directory be
#'   created and populated with the same subfolder structure as the pregated
#'   directory?
#'
#' @return an infinitySetup tibble that should be saved to a .Rds and then fed
#'   into an imputation pipeline as a set of args through pwalk.
#' @examples
#' # Note - commented because R-CMD-BUILD can't create the output folders in a
#' # safe way.
#' # infinitySetup <- infinity_setup(pregated_dir = system.file("extdata",
#' #                                                     "pregated",
#' #                                                     package = "flowHelpers"),
#' #                                 infinity_targets = c("Marker1",
#' #                                                      "Marker2"),
#' #                                 infinity_isotypes = "FMO",
#' #                                 create_output_dir = FALSE,
#' #                                 create_intermed_dir = FALSE)
#'
#' @export
infinity_setup <- function(pregated_dir,
                           infinity_targets,
                           infinity_isotypes,
                           create_output_dir = TRUE,
                           create_intermed_dir = TRUE){

  sampleList <- list.files(pregated_dir, full.names = TRUE)

  descList <- purrr::map(sampleList,
                  ~tibble::tibble(file = list.files(..1, pattern = ".fcs"),
                          Infinity_target = infinity_targets,
                          Infinity_isotype = infinity_isotypes))

  names(descList) <- list.files(pregated_dir, full.names = FALSE)

  targets <- purrr::map(descList, ~as.character(..1$Infinity_target))
  for(i in seq_along(targets)){
    names(targets[[i]]) <- descList[[i]]$file
  }

  isotypes <- purrr::map(descList, ~as.character(..1$Infinity_isotype))
  for(i in seq_along(isotypes)){
    names(isotypes[[i]]) <- descList[[i]]$file
  }

  if(create_output_dir == TRUE){
    dir.create(paste(dirname(pregated_dir), "/infinityOutput", sep = ""))

    purrr::walk(list.files(pregated_dir,
                    full.names = FALSE),
         ~dir.create(paste(dirname(pregated_dir), "/infinityOutput/",
                           ..1,
                           sep = "")))
  }
  if(create_intermed_dir == TRUE){
    dir.create(paste(dirname(pregated_dir), "/infinityIntermed", sep = ""))
    
    purrr::walk(list.files(pregated_dir,
                           full.names = FALSE),
                ~dir.create(paste(dirname(pregated_dir), "/infinityIntermed/",
                                  ..1,
                                  sep = "")))
  }

  infinitySetup <- tibble::tibble(path_to_fcs = list.files(pregated_dir,
                                                   full.names = TRUE),
                          path_to_output = list.files(paste(dirname(pregated_dir),
                                                            "/infinityOutput",
                                                            sep = ""),
                                                      full.names = TRUE),
                          path_to_intermediary_results = list.files(paste(dirname(pregated_dir),
                                                                          "/infinityIntermed",
                                                                          sep = ""), 
                                                                    full.names = TRUE),
                          annotation = targets,
                          isotype = isotypes)

  message("Don't forget to create a backbone annotation file by calling infinity_backbone()")

  return(infinitySetup)
}

# This one needs some work - probably pop it out into its own function to
# document it better.

infinity_backbone <- function(fcs_directory, pattern, recursive){
  fcsFiles <- list.files(fcs_directory,
                         pattern = pattern,
                         recursive = recursive,
                         full.names = TRUE)

  fs <- flowCore::read.flowSet(fcsFiles)

  backbone <- tibble::tibble(name = colnames(fs)) |>
    dplyr::left_join(tibble::as_tibble(flowCore::markernames(fs), rownames = "name")) |>
    dplyr::rename(desc = "value")

  message("IMPORTANT: \n
          Ensure that you edit the backbone file produced by this function
          to add a 'type' column that indicates whether a given marker is to
          be of type 'discard', 'backbone', or 'exploratory'. Once you have
          done this, save the resulting backbone file in the processed dir
          as 'infinity backbone.csv")

  return(backbone)

}
