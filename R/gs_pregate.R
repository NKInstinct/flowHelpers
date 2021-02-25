# Pregate Infinity

# This script will allow you to quickly gate through a gating set according to
# any gatingStrat object you make, according to the following template.

# gatingStrat <- tribble(
#   ~filterId, ~dims, ~subset, ~coords,
#   "Lymphocytes", list("FSC-A", "SSC-A"), "root", lincoords,
#   "Single cells", list("FSC-A", "FSC-W"), "Lymphocytes", lincoords,
#   "Live Leukocytes", list("Viability", "CD45"), "Single cells", logcoords)

# parentDir = directory to search for the fcs files in
# pattern = regexp to find the fcs files you want
# strat = gatingStrat object you made
# expGate = which gate to export to the fcs file
# pregateDir = parent directory to save the pregated data to. Note that each
#     flowSet will get its own directory from `pattern`, so this should be a
#     directory chain ending in "pregated"
# recursive = should the parentDir be searched recursively?
# linearChannels = number of linearChannels or list of their names/positions, as
#     per gs_create
# ncdf = should netCDF be used for the fcs read-in? defaults to false since we
#     aren't saving these gs objects, only their exports
# bins = how fine to resolve the data when gating.

gs_pregate <- function(parentDir,
                            pattern,
                            strat,
                            expGate,
                            pregateDir,
                            outDirName = NULL,
                            recursive = TRUE,
                            comp = "acquisition",
                            transform = "biex",
                            arcsinh_coeff = 150,
                            linearChannels = 6,
                            ncdf = FALSE,
                            bins = 512){
  # Turn this into something that takes ... instead of each param to gs_create
  # for a cleaner function call
  gs <- gs_create(parentDir,
                  recursive = recursive,
                  pattern = pattern,
                  ncdf = ncdf,
                  linearChannels = linearChannels,
                  comp = comp,
                  transform = transform,
                  arcsinh_coeff = arcsinh_coeff)

  # Three kinds of pwalk needed here. If gs is a single gatingset, do this
  # exactly as writ. If it is a list, either take one reference gs and apply the
  # whole strat to it, then apply all of those gates to all other gs, or
  # alternatively gate each things individually.
  pwalk(strat,
        gs_gate_interactive,
        gs = gs,
        bins = bins)

  if(is.null(outDirName)){
    outName <- pattern
  } else{
    outName <- outDirName
  }

  # Also include an option to write the gatingsets to an outdir as well
  fs <- gs_pop_get_data(gs,
                        y = expGate,
                        inverse.transform = TRUE) %>%
    flowCore::write.flowSet(outdir = paste(pregateDir,
                                           outName,
                                           "/",
                                           sep = ""))

  return(NULL)
}
