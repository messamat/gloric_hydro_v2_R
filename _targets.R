source("R/cartoCE_packages.R")
source("R/cartoCE_functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

tar_option_set(format = "qs")

#--- Data repositories --------------------------------------------------------
statsgdb = file.path(resdir, 'env_stats.gdb')

#--- Parameters ---------------------------------------------------------------
excluded_deps = c(75, 92, 93, 94)

############################# Define targets plan ##############################
list(
  #------------------------------- Read files ----------------------------------
  #Read in metadata
  tar_target(ddt_metadata_path, file.path(datdir, 'metadonnes_cartographie_cours_deau_20231106.xlsx')), #format='file'
  tar_target(ddt_nets_colnas_path, file.path(resdir, 'cartos_loi_eau_colNAs.csv'), format='file') ,

)