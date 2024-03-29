% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs_to_sce.R
\name{gs_to_sce}
\alias{gs_to_sce}
\title{Convert GatingSet to SCE with option to run UMAP as well}
\usage{
gs_to_sce(
  gs,
  node = "root",
  condition = Biobase::sampleNames(gs),
  do_UMAP = FALSE,
  feature_exclusion = NULL
)
}
\arguments{
\item{gs}{The GatingSet to draw the data from}

\item{node}{A string specifying which gate to take data from. Defaults to
"root" (all acquired events)}

\item{condition}{A vector of strings specifying which group each sample belongs to.
Must be of the same length as the number of samples in gs, and must be in
the same order as the samples in gs (can check this with sampleNames(gs)).
Defaults to just the sampleNames themselves, which is not ideal so you
should probably specify something here!}

\item{do_UMAP}{A boolean specifying whether to run UMAP dimensional reduction
or not. This is the whole point of turning flow data into an SCE, but it
can be computationally intense for even moderately sized experiments so you
need to opt-in.}

\item{feature_exclusion}{A regular expression specifying which rownames to
exclude from building the UMAP. Defaults to NULL (no exclusions), but you
likely want to exclude Time, FSC/SSC, and Viability at least since these
should be pre-gated.}
}
\description{
This function is very similar to gs_to_Seurat, but rather than shoehorning
all the flow data into Seurat's structures, it falls back to using just
(tidy)SingleCellExperiment packages and constructors.
}
\details{
Note that if you want to be able to interact with the sce in tidy format, you
should attach tidySingleCellExperiment yourself!
}
