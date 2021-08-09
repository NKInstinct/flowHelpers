# Goal here is a shiny-enabled pipeline for creating the whole seurat object, up
# to and including annotating the clusters if possible. This is going to be a
# big project but will hopefully showcase the kind of coding I can do now. Basic
# plan is to go through all the steps of building the Seurat object, from
# (renamed) infinity .fcs to a single Seurat object. At the various interjection
# points, I'll add in calls to shiny apps that will interactively add metadata,
# pick dimensions to use for neighbor search, add cluster labels, etc.