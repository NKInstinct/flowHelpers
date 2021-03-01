#' Prepare the infinity setup Rds
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
#'
#' @return an infinitySetup tibble that should be saved to a .Rds and then fed
#'   into an imputation pipeline as a set of args through pwalk.
#' @examples
#' infinitySetup <- infinity_setup(pregated_dir = system.file("extdata",
#'                                                            "pregated",
#'                                                            package = "flowHelpers"),
#'                                 infinity_targets = c("Marker1",
#'                                                      "Marker2"),
#'                                 infinity_isotypes = "FMO",
#'                                 create_output_dir = FALSE)
#'
#'
#' @export
infinity_setup <- function(pregated_dir,
                           infinity_targets,
                           infinity_isotypes,
                           create_output_dir = TRUE){

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

  infinitySetup <- tibble::tibble(path_to_fcs = list.files(pregated_dir,
                                                   full.names = TRUE),
                          path_to_output = list.files(paste(dirname(pregated_dir),
                                                            "/infinityOutput",
                                                            sep = ""),
                                                      full.names = TRUE),
                          annotation = targets,
                          isotype = isotypes)

  message("Don't forget to create a backbone annotation file by calling infinity_backbone()")

  return(infinitySetup)
}

infinity_backbone <- function(fcs_directory, pattern, recursive){
  fcsFiles <- list.files(fcs_directory,
                         pattern = pattern,
                         recursive = recursive,
                         full.names = TRUE)

  fs <- flowCore::read.flowSet(fcsFiles)

  backbone <- tibble::tibble(name = colnames(fs)) %>%
    dplyr::left_join(as_tibble(markernames(fs), rownames = "name")) %>%
    dplyr::rename(desc = "value")

  message("IMPORTANT: \n
          Ensure that you edit the backbone file produced by this function
          to add a 'type' column that indicates whether a given marker is to
          be of type 'discard', 'backbone', or 'exploratory'. Once you have
          done this, save the resulting backbone file in the processed dir
          as 'infinity backbone.csv")

  return(backbone)

}
