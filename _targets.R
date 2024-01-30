source("R/gloric_packages.R")
source("R/gloric_functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

tar_option_set(format = "qs",
               controller= crew_controller_local(workers = 4))

#--- Data repositories --------------------------------------------------------


#--- Parameters ---------------------------------------------------------------
min_yrs = 20

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
  
  tar_target(gaugep_dt,
             vect(dirname(path_gaugep), layer=basename(path_gaugep)) %>%
               as.data.table
  ),
  
  tar_target(g_anthropo_stats,
             read_format_anthropo_stats(inp = path_gauges_anthropo_stats)
             )
  ,
  tar_target(GRDC_filenames,
             read_GRDCgauged_paths(
               inp_GRDC_qdat_dir = path_GRDC_qdat_dir,
               in_GRDC_metadata = GRDC_metadata,
               inp_GIRES_metadata = path_GIRES_metadata,
               in_gaugep_dt = gaugep_dt
               ))
  ,
  tar_target(gmeta_formatted,
             format_gauges_metadata(
               in_GRDC_filenames = GRDC_filenames,
               in_GRDC_metadata = GRDC_metadata,
               in_gaugep_dt = gaugep_dt,
               in_g_anthropo_stats = g_anthropo_stats,
               min_yrs = min_yrs
             ))
  ,
  
  tar_target(anthropo_plot,
             plot_anthropo_stats (in_gmeta_formatted = gmeta_formatted)
  ),

  tar_target(ref_gauges,
             filter_reference_gauges(in_gmeta_formatted = gmeta_formatted)
  ),
  
             
)

