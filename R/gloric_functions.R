############################ UTILITY FUNCTIONS ##################################
#------ readformatGRDC -----------------
#' Read and pre-format GRDC data
#'
#' Reads text file of daily discharge data for a single GRDC station.
#' Creates columns for year, month, and date of last non-zero flow day +
#' computes yearly number of days of missing data
#'
#' @param path (character) path to the text file of daily discharge data in
#'   standard GRDC format.
#'
#' @return \link[data.table]{data.table} of daily discharge data with additional columns
#'
#' @export
readformatGRDC<- function(path) {
  #extract GRDC unique ID by formatting path
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  
  #Read GRDC text data
  gaugetab <- cbind(fread(path, header = T, skip = 40, sep=";",
                          colClasses = c('character', 'character', 'numeric',
                                         'numeric', 'integer')),
                    GRDC_NO = gaugeno)%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)
  
  #Format data
  gaugetab[, `:=`(year = as.numeric(substr(dates, 1, 4)), #Create year column
                  month = as.numeric(substr(dates, 6, 7)), #create month column
                  integervalue = fifelse(Original == round(Original), 1, 0) #Flag integer discharge values]
  )]
  
  #For each record, compute date of last non-zero flow day
  gaugetab[, prevflowdate := gaugetab[zero_lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA] #If non-zero flow, set prevflowdate to NA
  
  #Compute number of missing days per year, excluding NoData values
  gaugetab[!(Original %in% c(-999, -99, -9999, 99, 999, 9999) | is.na(Original)),
           `:=`(missingdays = diny(year)-.N,
                datadays = .N),
           by= 'year']
  
  return(gaugetab)
}


#------ format_for_ridges -----------------------------------------------------
format_for_ridges <- function(in_dt, fieldroot, id_var='grdc_no') {
  dt_formatted <- melt(in_dt,
                       id.vars=id_var,
                       measure.vars = grep(paste0(fieldroot, '_[0-9]{4}'), 
                                           names(in_dt),
                                           value=T)) %>%
    .[, year:= as.integer(tstrsplit(variable, '_', keep=2)[[1]])] %>%
    merge(in_dt[, list(year = seq(d_start, d_end)), by=id_var],
          by=c(id_var, 'year'), all.x=F)
  
  dt_formatted[, n := .N, by=year]
  return(dt_formatted)
}
############################ ANALYSIS FUNCTIONS ##################################
#------ read_GRDCgauged_paths -----------------
#' Read file paths to streamflow data from GRDC gauging stations
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @param inp_GRDC_qdat_dir path to directory containing streamflow data GRDC standard files.
#' @param in_gaugep table containing column named \code{GRDC_NO} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables, assuming
#' that files are called "GRDC_NO.txt", GRDC_NO being replaced with a 7-digit integer.
#'
#' @export


read_GRDCgauged_paths <- function(inp_GRDC_qdat_dir, in_GRDC_metadata,
                                  inp_GIRES_metadata, in_gaugep_dt) { #, gaugeid = 'GRDC_NO' down the line
  
  #Get data paths of daily records for gauge stations
  file_names_dt <- in_gaugep_dt[!is.na(in_gaugep_dt$grdc_no), 
                                list(filename = file.path(inp_GRDC_qdat_dir,
                                                          paste(grdc_no,"_Q_Day.Cmd.txt", sep=""))
                                ), 
                                by=grdc_no]
  
  #Check those that are missing
  file_names_dt[, exists := file.exists(filename)]
  
  gires_metada <- fread(inp_GIRES_metadata) #get metadata for stations from Messager et al. 2021
  #the grdc_no of some stations have changed
  
  #Make sure that there aren't missing 2024 gauges
  missing_gauges_grdc <- merge(file_names_dt[!(exists),],
                               in_GRDC_metadata,
                               by='grdc_no')
  
  #Check which gauges from 2014 that were used in 2021 et al. are not in the 2024 GRDC dataset
  missing_gauges_gires <- merge(file_names_dt[!(exists),],
                                gires_metada[, .(GRDC_NO, RIVER, STATION, AREA, LONG_ORG, LAT_ORG)],
                                by.x='grdc_no',
                                by.y='GRDC_NO')
  
  #Get those that match by coordinates
  missing_gauges_gires_rejoin_coord <-  merge(missing_gauges_gires,
                                              in_GRDC_metadata,
                                              by.x=c('LONG_ORG', 'LAT_ORG'),
                                              by.y=c('long', 'lat'),
                                              suffix=c('_old', '_new'))
  
  #Get those that strictly match by name
  setnames(missing_gauges_gires, c('RIVER', 'STATION'), c('river', 'station'))
  
  
  missing_gauges_gires_rejoin_name <-  merge(missing_gauges_gires,
                                             in_GRDC_metadata,
                                             by=c('river', 'station'),
                                             suffix=c('_old', '_new')) %>%
    .[!(grdc_no_new %in% missing_gauges_gires_rejoin_coord$grdc_no_new),]
  
  
  #from the rest, get those that match by name through partial match
  nojoin <- missing_gauges_gires[
    !(grdc_no %in% c(missing_gauges_gires_rejoin_coord$grdc_no_old,
                     missing_gauges_gires_rejoin_name$grdc_no_old)),
  ]
  
  missing_gauges_gires_fuzzyjoin <- stringdist_join(nojoin, 
                                                    in_GRDC_metadata, 
                                                    mode = "inner",
                                                    by = c("river", "station"),
                                                    max_dist = 2) %>%
    as.data.table %>%
    .[river.x != 'KURA',] %>%
    setnames(c('grdc_no.x', 'grdc_no.y'),
             c('grdc_no_old', 'grdc_no_new'))
  
  
  #Output a clean dataset
  file_names_dt_clean <- do.call(
    rbind,
    list(missing_gauges_gires_rejoin_coord[, .(grdc_no_new, grdc_no_old)],
         missing_gauges_gires_rejoin_name[, .(grdc_no_new, grdc_no_old)],
         missing_gauges_gires_fuzzyjoin[, .(grdc_no_new, grdc_no_old)]
    )
  ) %>%
    merge(file_names_dt, ., by.x='grdc_no', by.y='grdc_no_old', all.x=T) %>%
    .[is.na(grdc_no_new), grdc_no_new := grdc_no] %>%
    .[, filename := file.path(inp_GRDC_qdat_dir,
                              paste(grdc_no_new,"_Q_Day.Cmd.txt", sep=""))
    ] %>%
    setnames(c('grdc_no'), c('grdc_no_old')) 
  
  setnames(file_names_dt_clean, c('grdc_no_new'), c('grdc_no'))
  file_names_dt_clean[, exists := file.exists(filename)]
  file_names_dt_clean[, .N, by=exists]
  
  
  return(file_names_dt_clean)
}

#------ read_format_anthropo_stats --------------------------------------------
#inp <- tar_read(path_gauges_anthropo_stats)
read_format_anthropo_stats <- function(inp) {
  stats_raw <- fread(inp)
  
  #Create new column names
  col_dt <- data.table(
    pre=c('GHS_POP_E', 'GHS_BUILT_S_E', 'G_AEI_', 'cropland', 'gdw_dor_res_'),
    repl=c('pop', 'built', 'aei', 'crop', 'dor')) %>%
    .[,
      list(old = grep(pre, names(stats_raw), value=T, perl=T)),
      by=c('pre', 'repl')] %>%
    .[,           
      new := paste(repl,
                   str_extract(old, paste0("(?<=", pre, ")[0-9]{4}")),
                   sep='_'
      )]
  
  setnames(stats_raw, 
           c('grdc_p_o20y_cleanjoin', col_dt$old), 
           c('grdc_no',col_dt$new))
  
  #Check for NAs
  stats_raw[, sapply(.SD, function(x) sum(is.na(x)))]
  dor_cols = col_dt[repl=='dor',new]
  stats_raw[, 
            (col_dt[repl=='dor',new]) := sapply(
              .SD, function(x) fifelse(is.na(x), 0, x),
              simplify=F),
            .SDcols = dor_cols,
            by='grdc_no']
  
  #The other NAs are from Mauritius and Hawaii because of the HydroSHEDS landmask
  
}


#------ format_gauges_metadata  ---------------------------------------
# in_g_anthropo_stats = tar_read(g_anthropo_stats)
# in_gaugep_dt = tar_read(gaugep_dt)
# in_GRDC_metadata = tar_read(GRDC_metadata)
# in_GRDC_filenames = tar_read(GRDC_filenames)

format_gauges_metadata <- function(
    in_g_anthropo_stats,
    in_gaugep_dt,
    in_GRDC_metadata,
    in_GRDC_filenames,
    min_yrs) {
  
  gstats_merge <- merge(in_GRDC_filenames, 
                        in_GRDC_metadata, 
                        by='grdc_no') %>%
    merge(in_gaugep_dt, 
          by.x='grdc_no_old', by.y='grdc_no') %>%
    merge(in_g_anthropo_stats, 
          by.x='grdc_no_old', by.y='grdc_no') %>%
    .[, -c('exists', 'wmo_reg', 'lat', 'long', 'Join_Count', 
           'station_river_distance', 'TARGET_FID', 'MAIN_RIV', 'LENGTH_KM',
           'NEXT_DOWN', 'DIST_DN_KM', 'DIST_UP_KM', 'CATCH_SKM',
           'ORD_CLAS', 'ORD_FLOW', "OBJECTID", 'X', 'Y')] %>%
    .[, grdc_match := as.numeric(grdc_no_old != grdc_no)]
  
  gstats_merge_nodupli <- sort(gstats_merge, 'grdc_match') %>%
    .[!duplicated(grdc_no),] %>% #remove duplicates
    .[!is.na(d_yrs) | (d_yrs < min_yrs),]  #Remove those with insufficient daily data (anymore. used to)

  #Compute pop density
  pop_cols <- grep('pop_[0-9]{4}', names(gstats_merge_nodupli), value=T)  
  gstats_merge_nodupli[
    , (pop_cols) := sapply(
      .SD, function(x) x/up_area_skm_15s, simplify=F),
    .SDcols=pop_cols]
  
  #Re-scale DOR
  dor_cols <- grep('dor_[0-9]{4}', names(gstats_merge_nodupli), value=T)  
  gstats_merge_nodupli[
    , (dor_cols) := sapply(
      .SD, function(x) x/10, simplify=F),
    .SDcols=dor_cols]
  
  #Re-scale built extent
  built_cols <- grep('built_[0-9]{4}', names(gstats_merge_nodupli), value=T)  
  gstats_merge_nodupli[
    , (built_cols) := sapply(
      .SD, function(x) x/100, simplify=F),
    .SDcols=built_cols]
  
  
  return(gstats_merge_nodupli[, -c('grdc_match'), with=F])
}

#------ analyze_anthropo_stats -------------------------------------------------
#in_gmeta_formatted = tar_read(gmeta_formatted)

plot_anthropo_stats <- function(in_gmeta_formatted) {
  #Degree of regulation plot ---------------------------------------------------
  dt_format_dor <- format_for_ridges(in_dt=in_gmeta_formatted,
                                     fieldroot = 'dor') 
  
  p_dor <- ggplot(dt_format_dor, aes(x=value, y=year, fill = n, 
                                     group=year, height = ..count..)) +
    geom_density_ridges(stat='density', scale=15, alpha=1/2,
                        rel_min_height = 0.0005)  +
    scale_x_continuous(
      name = 'Degree of regulation (%)',
      trans=scales::pseudo_log_trans(base = 10),
      breaks=c(0, 1, 2, 5, 10, 100, 1000)) + 
    scale_y_continuous(name='Year',
                       breaks=seq(1900, 2020, 20)) +
    scale_fill_viridis(name=str_wrap('Number of active gauges', 20),
                       direction=-1,
                       breaks=c(200,1500,3000,4500,6000)) +
    coord_cartesian(ylim=c(1900,2020), expand=F, clip='off') +
    theme_ridges() + 
    theme(plot.margin=margin(3, 0.5, 0.1, 0.1, "cm")
    )
  #, scale=20
  
  #Crop plot -----------------------------------------------------
  dt_format_crop <- format_for_ridges(in_dt=in_gmeta_formatted,
                                      fieldroot = 'crop') %>%
    .[!is.na(value),]
  
  p_crop <- ggplot(dt_format_crop, aes(x=value, y=year, fill = n, 
                                       group=year, height = ..count..)) +
    geom_density_ridges(stat='density', scale=8, alpha=1/2,
                        rel_min_height = 0.001)  +
    scale_x_continuous(
      name = 'crop extent upstream (%)',
      breaks=c(0, 10, 25, 50, 75, 100)) + 
    scale_y_continuous(name='Year',
                       breaks=seq(1900, 2020, 20)) +
    scale_fill_viridis(name=str_wrap('Number of active gauges', 20),
                       direction=-1,
                       breaks=c(200,1500,3000,4500,6000)) +
    coord_cartesian(ylim=c(1900,2020), expand=F, clip='off') +
    theme_ridges() + 
    theme(
      legend.position='none'
      , plot.margin=margin(1, 0.5, 0.1, 0.1, "cm")
    )

  #Population density plot -----------------------------------------------------
  dt_format_pop <- format_for_ridges(in_dt=in_gmeta_formatted,
                                     fieldroot = 'pop') %>%
    .[!is.na(value),]

  p_pop <- ggplot(dt_format_pop, aes(x=value, y=year, fill = n, 
                                     group=year, height = ..count..)) +
    geom_density_ridges(stat='density', scale=5, alpha=1/2,
                        rel_min_height = 0.001)  +
    scale_x_continuous(
      name = expression('Population density upstream'~(km^-2)),
      trans=scales::pseudo_log_trans(base = 10),
      breaks=c(0, 1, 2, 5, 10, 100, 1000)) + 
    scale_y_continuous(name='Year',
                       breaks=seq(1900, 2020, 20)) +
    scale_fill_viridis(name=str_wrap('Number of active gauges', 20),
                       direction=-1,
                       option='F',
                       breaks=c(200,1500,3000,4500,6000)) +
    coord_cartesian(expand=F, clip='off') +
    theme_ridges() + 
    theme(
      legend.position = 'none'
      , plot.margin=margin(1, 0.5, 0.1, 0.1, "cm")
    )

  #Built plot -----------------------------------------------------
  dt_format_built <- format_for_ridges(in_dt=in_gmeta_formatted,
                                     fieldroot = 'built') %>%
    .[!is.na(value),]

  p_built <- ggplot(dt_format_built, aes(x=value, y=year, fill = n, 
                                     group=year, height = ..count..)) +
    geom_density_ridges(stat='density', scale=10, alpha=1/2,
                        rel_min_height = 0.001)  +
    scale_x_continuous(
      name = 'Built extent upstream (%)',
      trans=scales::pseudo_log_trans(base = 10),
      breaks=c(0, 1, 2, 5, 10, 25)) + 
    scale_y_continuous(name='Year',
                       breaks=seq(1900, 2020, 20)) +
    scale_fill_viridis(name=str_wrap('Number of active gauges', 20),
                       direction=-1,
                       option='F',
                       breaks=c(200,1500,3000,4500,6000)) +
    coord_cartesian(expand=F, clip='off') +
    theme_ridges() + 
    theme(plot.margin=margin(1, 0.5, 0.1, 0.1, "cm")
    )

  p_patchwork <- (p_crop + p_dor+ plot_layout(axes='collect'))/
    (p_pop + p_built + plot_layout(axes='collect'))
  
  return(p_patchwork)
}

#in_gmeta_formatted = tar_read(gmeta_formatted)