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
  tar_target(path_GRDCgaugedir, file.path(datdir, "gauges", "grdc", "q_day"), format='file') ,
  tar_target(path_GRDCmetadata, file.path(datdir, "gauges", "grdc", "GRDC_Stations.xlsx")),
  tar_target(path_GIRESgaugemetadata, file.path(datdir, 'gauges', 'gires', "high_qual_daily_stations.csv"), format='file'),
  tar_target(path_gaugep, file.path(resdir, 'stations_preprocess.gdb', 'grdc_p_o20y_cleanjoin'))
  ,
  
  
  #------------------------------ Read files -----------------------------------
  tar_target(GRDC_metadata,
             read_xlsx(path_GRDCmetadata, sheet = 'station_catalogue') %>%
               as.data.table
  )
  # ,
  # inp_GRDCgaugedir = tar_read(path_GRDCgaugedir)
  # inp_GIRESgaugemetadata = tar_read(path_GIRESgaugemetadata)
  # inp_gaugep = tar_read(path_gaugep)
  # in_GRDC_metadata = tar_read(GRDC_metadata)
  #
  # tar_target(GRDCgauged_filenames,
  #            read_GRDCgauged_paths(
  #              inp_GRDCgaugedir = path_GRDCgaugedir,
  #              inp_GRDCmetadata = path_GRDCmetadata,
  #              inp_GIRESgaugemetadata = path_GIRESgaugemetadata,
  #              inp_gaugep = path_gaugep
  #              ))
  # )
  # ,
  # 
  # GRDCgaugestats = future_map(GRDCgauged_filenames, #To map
  #                             comp_GRDCdurfreq, #Function to run on each file name
  #                             maxgap = 20,
  #                             in_gaugep = gaugep,
  #                             windowsize = 20,
  #                             fullwindow = FALSE,
  #                             monthsel = NULL, #Other arguments in function
  #                             mdurthresh = 1,
  #                             verbose = FALSE,
  #                             .progress = TRUE),
)
