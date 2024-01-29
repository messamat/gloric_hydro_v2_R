source("R/gloric_packages.R")
source("R/gloric_functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

tar_option_set(format = "qs",
               controller= crew_controller_local(workers = 4))

#--- Data repositories --------------------------------------------------------


#--- Parameters ---------------------------------------------------------------


############################# Define targets plan ##############################
list(
  #------------------------------- File paths ----------------------------------
  tar_target(path_GRDC_qdat_dir, file.path(datdir, "gauges", "grdc", "q_day"), format='file') ,
  tar_target(path_GRDC_metadata, file.path(datdir, "gauges", "grdc", "GRDC_Stations.xlsx"), format='file'),
  tar_target(path_GIRES_metadata, file.path(datdir, 'gauges', 'gires', "high_qual_daily_stations.csv"), format='file'),
  tar_target(path_gaugep, file.path(resdir, 'stations_preprocess.gdb', 'grdc_p_o20y_cleanjoin')),
  
  tar_target(path_gauges_anthropo_stats, file.path(resdir, 'stations_anthropo_stats_upst.csv'), format='file')
  ,  
  
  #------------------------------ Read files -----------------------------------
  tar_target(GRDC_metadata,
             read_xlsx(path_GRDC_metadata, sheet = 'station_catalogue') %>%
               as.data.table
  )
  ,
  
  tar_target(g_anthropo_stats,
             fread(path_gauges_anthropo_stats)
             )
  ,

  tar_target(GRDCgauged_filenames,
             read_GRDCgauged_paths(
               inp_GRDC_qdat_dir = path_GRDC_qdat_dir,
               in_GRDC_metadata = GRDC_metadata,
               inp_GIRES_metadata = path_GIRES_metadata,
               inp_gaugep = path_gaugep
               ))
  #,
  # inp_GRDC_qdat_dir = tar_read(path_GRDC_qdat_dir)
  # in_GRDC_metadata = tar_read(path_GRDC_metadata)
  # inp_GIRES_metadata = tar_read(path_GIRES_metadata)
  # inp_gaugep = tar_read(path_gaugep)
             
)

