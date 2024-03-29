% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infinity_renameOutputs.R
\name{infinity_renameOutputs}
\alias{infinity_renameOutputs}
\title{Rename the outputs of infinityFlow to something more helpful than concat.fcs}
\usage{
infinity_renameOutputs(
  parentDir,
  pattern = "\\\\.fcs$",
  concat_tag_to_remove = "FCS/concatenated",
  prefix = NULL,
  trialRun = TRUE
)
}
\arguments{
\item{parentDir}{a filepath to the directory containing both the
infinityOutput and pregated directories.}

\item{pattern}{a regular expression that will be used to find the output
files. The default of "\\.fcs" works if you have only output a single
concatenated file, but you'll need to be more specific if your
infinityOutput also includes non-concatenated .fcs files (or any other .fcs
for that matter, since the directories must be searched recursively).}

\item{concat_tag_to_remove}{a regular expression showing what string to
remove from the filename. The default will work if your infinityOutput only
has a concatenated file in it.}

\item{prefix}{a string to add in front of the filename for each renamed file.
Pasted directly on, so if you want a specific separater, add it here. The
default (NULL) pastes no string. Adding a prefix is useful if you are going
to be combining several batches of samples with similar naming strategies
(WT1, WT2, etc).}

\item{trialRun}{Boolean. If true, this function will return a list of old and
new filepaths but won't actually change them. You can either then run
file.rename yourself using the old and new paths, or rerun with trialRun =
FALSE to actually rename the files.}
}
\value{
Nothing - this modifies the filenames in place, so be careful! It's
hard to undo this, hence the trialRun option to make sure it's working as
intended before hitting go.
}
\description{
Annoyingly, all infinityFlow concatenated files are output as 'concat.fcs',
making it very difficult to use them in downstream analysis software like
Seurat or CATALYST. This function renames them the name of their
infinityOutput directory so that adding them to some downstream process
doesn't throw errors later.
}
