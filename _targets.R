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
  tar_target(path_GRDC_qdat_dir, file.path(datdir, "gauges", "grdc", "q_day")), # format='file') ,
  tar_target(path_GRDC_metadata, file.path(datdir, "gauges", "grdc", "GRDC_Stations.xlsx")), # format='file'),
  tar_target(path_GIRES_metadata, file.path(datdir, 'gauges', 'gires', "high_qual_daily_stations.csv")), # format='file'),
  tar_target(path_gaugep, file.path(resdir, 'stations_preprocess.gdb', paste0('grdc_p_o', min_yrs, 'y_cleanjoin'))),
  tar_target(url_riggs2023, "https://zenodo.org/records/7150168/files/zenodo.zip", format='url'),
  
  tar_target(path_gauge_netdist, file.path(resdir, 'grdc_p_o20y_cleanjoin_netdist_tab.csv')), # format='file')
  tar_target(path_gauge_geodist, file.path(resdir, 'grdc_p_o20y_cleanjoin_geodist_tab.csv')), # format='file')
  
  tar_target(path_gauges_anthropo_stats, file.path(resdir, 'stations_anthropo_stats_upst.csv')), # format='file')
  
  tar_target(path_gauges_tmin, file.path(resdir, 'grdc_p_o20y_cleanjoin_terra_tmin.csv')), # format='file')
  tar_target(path_gauges_tmax, file.path(resdir, 'grdc_p_o20y_cleanjoin_terra_tmax.csv')), # format='file')
  tar_target(path_rivice_elv, file.path(resdir, 'landsat_grdc_p_o20y_clean_sub_elvstats.csv')), # format='file')
  tar_target(path_rivice_ts, file.path(resdir, 'global_river_ice_dataset_sub.csv')), # format='file')
  tar_target(path_gauges_rivicetiles_join, file.path(resdir, 'grdc_p_o20y_clean_landsatjoin.csv'))
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
  
  tar_target(riggs2023_dirpath,
             download_riggs2023(in_url=path_riggs2023,
                                out_file=file.path(datdir, 'gauges', 
                                                   'riggs2023', 'zenodo.zip'))
  )
  ,
  
  tar_target(riggs2023_dt,
          format_riggs2023(riggs2023_dirpath)
  )
  ,
  
  tar_target(rivice_dt,
             read_format_rivice(inp_elv=path_rivice_elv,
                                inp_ts=path_rivice_ts,
                                inp_gtiles_join=path_gauges_rivicetiles_join)
  )
  ,
  
  tar_target(gauges_tmin,
             fread(path_gauges_tmin,
                   select= c('grdc_no', 'time', 'tmin'),
                   colClasses= c('character', 'date', 'integer', 'numeric')
             )
  )
  ,
  
  tar_target(gauges_tmax,
             fread(path_gauges_tmax,
                   select= c('grdc_no', 'time', 'tmax'),
                   colClasses= c('character', 'date', 'integer', 'numeric')
             )
  )
  ,
  
  tar_target(netdist,
          fread(path_gauge_netdist) %>%
            .[, (c('grdc_no', 'grdc_no_destination')) := 
                tstrsplit(Name, ' - ')]  %>%
            .[, c('grdc_no', 'grdc_no_destination', 'Total_Leng', 'travel_mod'),
              with=F] %>%
            setnames(c('Total_Leng', 'travel_mod'), c('dist_net', 'dist_net_dir'))
  )
  ,
  
  tar_target(geodist,
             fread(path_gauge_geodist) %>%
               .[, c('grdc_no_origin', 'grdc_no_destination', 
                     'NEAR_DIST', 'NEAR_RANK')] %>%
               setnames(c('grdc_no_origin', 'NEAR_DIST'),
                        c('grdc_no', 'dist_geo'))
  )
  ,
  
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
  
  #------------------------------ Format data ----------------------------------
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
             plot_anthropo_stats(in_gmeta_formatted = gmeta_formatted)
  ),
  
  tar_target(ref_gauges,
             filter_reference_gauges(in_gmeta_formatted = gmeta_formatted,
                                     max_dor = 2, #2% of mean annual flow regulated by a dam
                                     max_crop = 25, #max 25% of upstream area cultivated
                                     max_pop = 100, #max 100 people/km2 upstream
                                     max_built = 1 #max 1% of upstream area is "built"
             )  
  )
  # ,
  # tar_target(data_for_qc,
  #            prepare_QC_data(in_ref_gauges = ref_gauges,
  #                            in_gmeta_formatted = gmeta_formatted,
  #                            in_geodist = geodist,
  #                            in_netdist = netdist,
  #                            in_tmax = gauges_tmax,
  #                            #in_pdsi_dir = pdsi_dir,
  #                            in_rivice = rivice_dt,
  #                            in_riggs2023 = riggs2023_dt
  #            )  
  # )
  
  
)
