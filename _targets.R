source("R/gloric_packages.R")
source("R/gloric_functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')
temp_qs_dir = file.path(resdir, 'temp_qs')
if (!dir.exists(temp_qs_dir)) {
  dir.create(temp_qs_dir)
}

nthreads <- round(parallel::detectCores()*0.75)
plan(multisession, workers=nthreads)
tar_option_set(format = "qs",
               controller= crew_controller_local(workers = nthreads))

#--- Parameters ---------------------------------------------------------------
min_yrs = 20
overwrite=F

############################# Define targets plan ##############################
list(
  #------------------------------- File paths ----------------------------------
  tar_target(path_GRDC_qdat_dir, file.path(datdir, "gauges", "grdc", "q_day")), # format='file') ,
  tar_target(path_GRDC_metadata, file.path(datdir, "gauges", "grdc", "GRDC_Stations.xlsx")), # format='file'),
  tar_target(path_GIRES_metadata, file.path(datdir, 'gauges', 'gires', "high_qual_daily_stations.csv")), # format='file'),
  tar_target(path_gaugep, file.path(resdir, 'stations_preprocess.gdb', paste0('grdc_p_o', min_yrs, 'y_cleanjoin'))),
  tar_target(url_riggs2023, "https://zenodo.org/records/7150168/files/zenodo.zip", format='url'),
  
  tar_target(path_gauge_netdist, file.path(resdir, 'grdc_p_o20y_cleanjoin_netdist_tab.csv')), # format='file'),
  tar_target(path_gauge_geodist, file.path(resdir, 'grdc_p_o20y_cleanjoin_geodist_tab.csv')), # format='file')
  
  tar_target(path_gauges_anthropo_stats, file.path(resdir, 'stations_anthropo_stats_upst.csv')), # format='file')
  
  tar_target(path_gauges_pdsi, file.path(resdir, 'terra', 'PDSI.csv')), #format='file')
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
  
  tar_target(gauges_satQ,
          format_riggs2023(riggs2023_dirpath)
  )
  ,
  
  tar_target(gauges_rivice,
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
  
  tar_target(gauges_pdsi,
              fread(path_gauges_pdsi) %>%
               melt(id.vars = 'grdc_no',
                    measure.vars = grep('PDSI_.*', names(.), value=T),
                    value.name = 'PDSI'
               ) %>%
               .[, `:=`(date=as.character(as.Date(gsub('PDSI_', '', variable),
                                                  '%Y%m%d')),
                        grdc_no=as.character(grdc_no))] %>%
               .[, -'variable', with=F]
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
  ,
  #Prepare data for QCing
  tar_target(data_for_qc,
             future_lapply(
               unique(ref_gauges$dt$grdc_no), 
               function(in_no) {
                 print(in_no)
                 out_qs <- file.path(temp_qs_dir, 
                                     paste0('data_for_qc_', in_no, '.qs')
                 )
                 
                 out_list <- prepare_QC_data_util(
                   in_grdc_no=in_no,
                   in_ref_gauges=ref_gauges,
                   in_gmeta_formatted=gmeta_formatted,
                   #in_geodist,
                   in_netdist=netdist,
                   in_tmax=gauges_tmax,
                   in_pdsi=gauges_pdsi,
                   in_rivice=gauges_rivice,
                   in_riggs2023=gauges_satQ) 
                 
                 if (!file.exists(out_qs) | overwrite) {
                   qsave(out_list$q_dt_attri, out_qs)
                 }
                 
                 return(data.table(
                   grdc_no = in_no,
                   qs_path = out_qs,
                   potential_npr = out_list$potential_npr,
                   integer_perc = out_list$integer_perc,
                   near_gcols_sel = out_list$near_gcols_sel)
                 )
               }
               ) %>% rbindlist(., fill=T)
  ),
  
  tar_target(q_outliers_flags,
             lapply(
               data_for_qc[potential_npr==T & integer_perc<0.95, grdc_no], 
               function(in_no) {
                 print(in_no)
                 row <- data_for_qc[grdc_no == in_no,]
                 qs_path <- row$qs_path[[1]]
                 nearg_cols_sel <- row$near_gcols_sel[[1]]
                 
                 out_qs <- file.path(temp_qs_dir, 
                                     gsub('data_for_qc_', 
                                          'q_outliers_flags_',
                                          basename(qs_path))
                                     )
                 
                 overwrite=F
                 if (!file.exists(out_qs) | overwrite) {
                   out_flags <- detect_outliers_ts(
                     in_data=qread(qs_path),
                     in_nearg_cols=nearg_cols_sel)
                   
                   qsave(out_flags$outliers_dt, out_qs)

                 return(data.table(
                   grdc_no = in_no,
                   out_qs = out_qs,
                   arima_model = out_flags$arima_model,
                   p_fit = out_flags$p_fit,
                   p_seasonal = out_flags$p_fit_seasonal,
                   p_fit_forplotly = out_flags$p_fit_forplotly)
                 )
                 }                 
               }
             ) %>% rbindlist(., fill=T)
  ),
  
  tar_target(
    metastats_dt,
    compute_metastatistics_wrapper(q_outliers_flags)
  ),
  
  tar_target(
    metastats_analyzed,
    analyze_metastats(in_metastats_dt=metastats_dt)
  ),
  
  tar_target(
    noflow_hydrostats,
    compute_noflow_hydrostats_wrapper(in_metastats_analyzed=metastats_analyzed,
                                      in_metastats_dt=metastats_dt,
                                      in_gaugep_dt=gaugep_dt,
                                      max_interp_sel = 10,
                                      max_miss_sel = 0,
                                      min_nyears = 15,
                                      q_thresh = 0.001) 
  )
)




# tar_target(outlier_plot,
#            plotGRDCtimeseries(GRDCgaugestats_record=data_for_qc$q_dt_attri,
#                               outpath=NULL,
#                               maxgap = 366,
#                               showmissing = T)
#            )
# ,
