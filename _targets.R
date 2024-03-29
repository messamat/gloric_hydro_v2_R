source("R/gloric_packages.R")
source("R/gloric_functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')
figdir = file.path(resdir, 'figures')
if (!dir.exists(figdir)) {
  dir.create(figdir)
}

temp_qs_dir = file.path(resdir, 'temp_qs')
if (!dir.exists(temp_qs_dir)) {
  dir.create(temp_qs_dir)
}

nthreads <- round(parallel::detectCores()*0.75)
plan(multisession, workers=nthreads)
tar_option_set(format = "qs")
#,
#              controller= crew_controller_local(workers = nthreads))

#------ There is an issue with large fread in the pipeline so need to run this---
gires_qs_path <- file.path(resdir, 'gires_net_sub.qs')
if (!file.exists(gires_qs_path)) {
  gires_net_dt <- fread(file.path(datdir, 'gires', 'GIRES_v10_rivers.csv'),
                        select = c("HYRIV_ID","LENGTH_KM", 
                                   "predprob1", "predprob30", 
                                   predvars$varcode))
  qsave(gires_net_dt,
        gires_qs_path)
}

#--- Parameters ---------------------------------------------------------------
min_nyrs = 15
overwrite=F
manual_class_order <- c(4, 1, 3, 2, 5, 9, 6, 8, 7) #tree rotation while keeping relationship intact
class_colors <- c("#62001D", "#E60000", "#FF5500", "#FFC100",
                  "#00A884", "#00FFC5", "#897044", "#D69DBC", "#CDAA66")
# class_colors <- c("#A20101", "#E69F00", "#DECF03","#0072B2", "#03C18E",
#                   "#CC79A7", "#686868", "#7570b3","#A8BF7A")

hydrostats_sel <- c('f0', 'medianN', 'sdN',
                    'medianD', 'sdD', 
                    'Ic', 'bfi', 'medianDr',
                    'Fper', 'FperM10',
                    'PDSIdiff', 'P90PDSI',
                    'r_cos_theta', 'r_sin_theta')

############################# Define targets plan ##############################
list(
  #------------------------------- File paths ----------------------------------
  tar_target(path_GRDC_qdat_dir, file.path(datdir, "gauges", "grdc", "q_day")), # format='file') ,
  tar_target(path_GRDC_metadata, file.path(datdir, "gauges", "grdc", "GRDC_Stations.xlsx")), # format='file'),
  tar_target(path_GIRES_metadata, file.path(datdir, 'gauges', 'gires', "high_qual_daily_stations.csv")), # format='file'),
  tar_target(path_gaugep, file.path(resdir, 'stations_preprocess.gdb', paste0('grdc_p_o20y_cleanjoin'))), ####EDIT down the line
  tar_target(url_riggs2023, "https://zenodo.org/records/7150168/files/zenodo.zip", format='url'),

  tar_target(path_riveratlas_meta, file.path(datdir, "hydroatlas", "HydroATLAS_metadata_MLMv11.xlsx")),
  tar_target(path_riveratlas_legends, file.path(datdir, "hydroatlas", "HydroATLAS_v10_Legends.xlsx")),


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

  tar_target(
    predvars,
    selectformat_predvars(inp_riveratlas_meta = path_riveratlas_meta)
  ),

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
               min_yrs = min_nyrs
             ))
  ,

  tar_target(anthropo_plot,
             plot_anthropo_stats(in_gmeta_formatted = gmeta_formatted,
                                 vlines = list(dor=2, crop=25, pop=100, built=1),
                                 export = T,
                                 fig_outdir = figdir)
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

  tar_target(
    gauges_for_qc,
    {terra::vect(dirname(path_gaugep), layer=basename(path_gaugep)) %>%
        terra::merge(data_for_qc[potential_npr==T & integer_perc<0.95,],
              by='grdc_no') %>%
        terra::writeVector(filename=file.path(resdir, 'gaugesp_for_qc.shp'),
                           overwrite=T)
    }
  ),

  tar_target(q_outliers_flags,
             future_lapply(
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

                 overwrite=T
                 if (!file.exists(out_qs) | overwrite) {
                   out_flags <- detect_outliers_ts(
                     in_data=qread(qs_path),
                     in_nearg_cols=nearg_cols_sel,
                     run_arima = F)

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
    compute_metastatistics_wrapper(q_outliers_flags,
                                   in_gaugep_dt=gaugep_dt)
  ),

  tar_target(
    metastats_analyzed,
    analyze_metastats(in_metastats_dt=metastats_dt,
                      export = T,
                      fig_outdir = figdir)
  ),

  tar_target(
    noflow_hydrostats,
    compute_noflow_hydrostats_wrapper(in_metastats_analyzed=metastats_analyzed,
                                      in_metastats_dt=metastats_dt,
                                      max_interp_sel = 5,
                                      max_miss_sel = 0,
                                      min_nyears = 15,
                                      q_thresh = 0.001)
  ),

  tar_target(
    noflow_hydrostats_preformatted,
    preformat_hydrostats(in_hydrostats = noflow_hydrostats)
  ),

  tar_target(
    noflow_clusters,
    cluster_noflow_gauges_full(in_hydrostats_preformatted=
                                 noflow_hydrostats_preformatted,
                               stats_sel = hydrostats_sel,
                               in_colors=class_colors)
  ),

  tar_target(
    sel_cluster_postanalysis,
    noflow_hydrostats_preformatted$mat %>%
      .[, -which(colnames(.) 
                 %in% c('r_cos_theta', 'r_sin_theta'))] %>%
    cuttree_and_visualize(
      in_mat =  .,
      in_hclust = noflow_clusters$hclust_reslist_all$ward.D2$hclust,
      in_kclass = noflow_clusters$kclass,
      in_colors = class_colors, 
      in_meltdt = noflow_hydrostats_preformatted$dt,
      id_col = 'grdc_no',
      classnames= data.table(
        gclass = seq(1, 9),
        classnames= manual_class_order
      ),
      colorder = manual_class_order
    )
  ),
  
  tar_target(
    out_gaugesp_path,
    export_gauges_classes(in_cluster_postanalysis=sel_cluster_postanalysis,
                          in_path_gaugep =path_gaugep,
                          out_shp_root=file.path(resdir, 'gaugep_classstats_ward')
    )
  ),

  tar_target(
    unit_hydrographs,
    plot_class_hydrograph_wrapper(in_cluster_postanalysis=sel_cluster_postanalysis,
                                  in_metastats_dt = metastats_dt,
                                  max_interp_sel = 5,
                                  max_miss_sel = 0,
                                  noflow_qthresh = 0.001
    )
  ),

  tar_target(
    cluster_sensitivity,
    analyze_cluster_sensitivity(
      in_noflow_clusters = noflow_clusters,
      in_cluster_postanalysis = sel_cluster_postanalysis,
      in_hydrostats_preformatted = noflow_hydrostats_preformatted,
      stats_sel = hydrostats_sel,
      classnames = manual_class_order,
      export = T,
      fig_outdir = figdir
      )
  ),

  tar_target(
    gires_dt,
    qread(gires_qs_path) %>%
      comp_derivedvar
  )
  ,

  tar_target(
    clusters_env_analysis,
    analyze_environmental_correlates(
      in_cluster_postanalysis = sel_cluster_postanalysis,
      in_gaugep_dt = gaugep_dt,
      in_predvars = predvars,
      in_gires_dt = gires_dt,
      in_colors = class_colors,
      export = T,
      fig_outdir = figdir
    )
  ),

  tar_target(
    gauge_representativeness,
    analyze_gauge_representativeness(in_noflow_clusters = noflow_clusters,
                                     in_gaugep_dt = gaugep_dt,
                                     in_predvars = predvars,
                                     in_gires_dt = gires_dt,
                                     export = T,
                                     fig_outdir = figdir
    )
  ),

  tar_target(
    export_plots,
    {
      ggsave(file.path(figdir, paste0('p_varscor', '_',
                                      format(Sys.Date(), '%Y%m%d'), '.png')),
             noflow_clusters$p_varscor)
      
      pname_suffix_pdf <- paste0(gsub('[.]','_', noflow_clusters$chosen_hclust), '_',
                             format(Sys.Date(), '%Y%m%d'), '.pdf')
      pname_suffix_png <- paste0(gsub('[.]','_', noflow_clusters$chosen_hclust), '_',
                                 format(Sys.Date(), '%Y%m%d'), '.png')
      
      lapply(list(pname_suffix_pdf, pname_suffix_png), function(pname_suffix) {
        ggsave(file.path(figdir, paste0('p_scree_', pname_suffix)),
               noflow_clusters$hclust_reslist_all[[noflow_clusters$chosen_hclust]]$p_scree)
        
        ggsave(file.path(figdir, paste0('p_dendo_', pname_suffix)),
               sel_cluster_postanalysis$p_dendo,
               width = 25, height = 10, units = "cm")
        
        ggsave(file.path(figdir, paste0('p_boxplot_small_', pname_suffix)),
               sel_cluster_postanalysis$p_boxplot,
               width = 20, height = 20, units = "cm")
        
        ggsave(file.path(figdir, paste0('p_hydrographs_',  pname_suffix)),
               unit_hydrographs,
               width = 20, height = 20, units = "cm")
      }
        )

    }
  ),
  
  tar_target(
    export_tables,
    {
    #Hydrological metrics table
    lapply(c(0,1,2), function(digits) {
      stats_table <- sel_cluster_postanalysis$class_dt %>%
        rbind(., copy(.)[, gclass := 'All']) %>%
        .[variable=='f0', value := 365.25*value] %>%
        .[, list(mean = mean(value, na.rm=T), 
                 q10 = quantile(value, 0.1, na.rm=T),
                 q90 = quantile(value, 0.9, na.rm=T),
                 classn = length(unique(grdc_no))),
          by=.(variable, gclass)] %>%
        digitform(cols=c('mean', 'q10', 'q90'),
                  extradigit = digits, inplace=F) %>%
        .[, stat_format := paste0(mean, ' (', q10, '-', q90, ')')] %>%
        dcast(gclass+classn~variable, value.var = 'stat_format')
      
      
      fwrite(stats_table, 
             file.path(figdir,paste0('hydrostats_', digits, 'digits_',
                                     format(Sys.Date(), '%Y%m%d'), '.csv')
             ))
    })
    
    fwrite(gauge_representativeness$KLdiv_marginal,
           file.path(figdir,paste0('global_KLdiv_',
                                   format(Sys.Date(), '%Y%m%d'), '.csv'))
    )
    }
  )
)