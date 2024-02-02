############################ UTILITY FUNCTIONS ##################################
#------ mergDTlist -----------------------------------------
#Function to merge a list of data.tables (dt), adding a suffix to all columns
#by dt based on the name of the input data.table
mergeDTlist <- function(dt_list, by = NULL, all = TRUE, sort = FALSE,
                        set_suffix=TRUE) {
  
  if (set_suffix) {
    dt_list <-  Map(function(dt_name, dt) {
      dt_copy <- copy(dt)
      cols_to_rename <- names(dt)[!(names(dt) %in% by)]
      setnames(dt_copy,
               old=cols_to_rename,
               new=paste(cols_to_rename, dt_name, sep='_'))
      return(dt_copy)
    },
    names(dt_list), dt_list
    )
  }
  
  out_merge <- Reduce(
    function(...) {
      merge(..., by = by, all = all, sort = sort)
    }, dt_list)
  
  return(out_merge)
}
#------ fread_cols -----------------
#' fread columns
#'
#' Fast data.table-based reading of a subset of columns from a table
#'
#' @param file_name path of the table to be read
#' @param cols_tokeep character vector, names of columns to read
#'
#' @return a data.table with the \code{cols_tokeep} that are found in the table
#'
#' @examples
#' fread_cols(iris, c('Sepal.Length', 'Sepal.Width')
#'
#' @export
fread_cols <- function(file_name, cols_tokeep) {
  #Only read the first row from the file
  header <- fread(file_name, nrows = 1, header = FALSE)
  #Check which columns are in the table
  keptcols <- cols_tokeep[cols_tokeep %chin% unlist(header)]
  missingcols <- cols_tokeep[!(cols_tokeep %chin% unlist(header))]
  paste('Importing', file_name, 'with ', length(keptcols),
        'columns out of ', length(cols_tokeep), 'supplied column names')
  #Reading in table
  paste(missingcols, 'columns are not in the file...')
  
  dt <- fread(input=file_name, header=TRUE,
              select=keptcols, verbose=TRUE, ...)
  return(dt)
}

#------ diny -----------------
#' Number of days in the year
#'
#' Computes the number of days in a given year, taking in account leap years
#'
#' @param year Numeric or integer vector of the year (format: 4-digit year, %Y).
#'
#' @return Number of days in that year.
#'
#' @examples
#' diny(1999)
#' diny(2000)
#' diny(2004)
#' diny(2100)
#' diny(1600)
#'
#' @export
diny <- function(year) {
  365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
}

#------ zero_lomf -----------------
#' Last \[non-zero\] Observation Moved Forward (lomf)
#'
#' Finds the index, for each row, of the previous row with a non-zero value
#'
#' @param x Numeric vector.
#' @param first (logical) Whether to consider first value as a non-zero value
#'   whose index is moved forward even if it is zero. This prevents having NAs
#'   in the results and somewhat assumes that, for a time series, the day prior
#'   to the first value is non-zero.
#'
#' @return Numeric vector of the indices of the previous non-zero for each
#'   element of the input vector.
#'
#' @examples
#' test1 <- c(1,1,1,0,0,0,0,1,1)
#' zero_lomf(test1)
#' test2 <- c(0,0,0,0,0,1,1,0,1)
#' zero_lomf(test2, first=FALSE)
#' zero_lomf(test2, first=TRUE)
#'
#' @export
zero_lomf <- function(x, first=TRUE) {
  if (length(x) > 0) {
    non.zero.idx <- which(x != 0)
    if(first==T & x[1]==0)
      non.zero.idx=c(1,non.zero.idx)
    #Repeat index of previous row with non-zero as many times gap until next non-zero values
    rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
  }
}

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
readformatGRDC<- function(path, newformat=T) {
  #Read GRDC text data
  if (!newformat) {
    #extract GRDC unique ID by formatting path
    gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
    
    gaugetab <- cbind(fread(path, header = T, skip = 40, sep=";",
                            colClasses = c('character', 'character', 'numeric',
                                           'numeric', 'integer')),
                      grdc_no = gaugeno)%>%
      setnames('YYYY-MM-DD', 'dates') %>%
      setorder(grdc_no, dates)
  } else {
    #extract GRDC unique ID by formatting path
    gaugeno <- strsplit(basename(path), '_')[[1]][1]
    
    gaugetab <- cbind(fread(path, header = T, skip = 36, sep=";",
                            colClasses = c('character', 'character', 'numeric')),
                      grdc_no = gaugeno) %>%
      setnames(c('YYYY-MM-DD','Value'),
               c('date', 'Qobs')) %>%
      .[, -c('hh:mm'), with=F] %>%
      setorder(grdc_no, date)
  }

  #Format data
  gaugetab[, `:=`(year = as.numeric(substr(date, 1, 4)), #Create year column
                  month = as.numeric(substr(date, 6, 7)), #create month column
                  integervalue = fifelse(Qobs == round(Qobs), 1, 0) #Flag integer discharge values]
  )]
  
  #For each record, compute date of last non-zero flow day
  gaugetab[, prevflowdate := gaugetab[zero_lomf(Qobs),'date', with=F]] %>% #Get previous date with non-zero flow
    .[Qobs != 0, prevflowdate:=NA] #If non-zero flow, set prevflowdate to NA
  
  #Compute number of missing days per year, excluding NoData values
  gaugetab[!(Qobs %in% c(-999, -99, -9999, 99, 999, 9999) | is.na(Qobs)),
           `:=`(missingdays = diny(year)-.N,
                datadays = .N),
           by= 'year']
  
  gaugetab[(Qobs %in% c(-999, -99, -9999, 99, 999, 9999)),
           Qobs := NA]

  return(gaugetab)
}


#------ melt_anthropo_stats -----------------------------------------------------
melt_anthropo_stats <- function(in_dt, fieldroot, 
                                id_var='grdc_no', keep_yrs = F) {
  dt_formatted <- melt(in_dt,
                       id.vars=id_var,
                       measure.vars = grep(paste0(fieldroot, '_[0-9]{4}'), 
                                           names(in_dt),
                                           value=T)) %>%
    .[, year:= as.integer(tstrsplit(variable, '_', keep=2)[[1]])] %>%
    merge(in_dt[, list(year = seq(d_start, d_end)), by=id_var], .,
          by=c(id_var, 'year'), all.x=keep_yrs)
  
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
  out_dt <- stats_raw[, 
                      (col_dt[repl=='dor',new]) := sapply(
                        .SD, function(x) fifelse(is.na(x), 0, x),
                        simplify=F),
                      .SDcols = dor_cols,
                      by='grdc_no']
  #The other NAs are from Mauritius and Hawaii because of the HydroSHEDS landmask
  
  out_dt[, grdc_no := as.character(grdc_no)]
  
  return(out_dt)
}


#------ read_format_river_ice ------------------------------------------------------
# inp_elv = tar_read(path_rivice_elv)
# inp_ts = tar_read(path_rivice_ts)
# inp_gtiles_join = tar_read(path_gauges_rivicetiles_join)
read_format_rivice <- function(inp_elv, inp_ts, inp_gtiles_join) {
  ts <- fread(inp_ts)
  elv <- fread(inp_elv, select = c('PATH', 'ROW', 'mean_elv'))
  gtiles_join <- fread(inp_gtiles_join, select = c('grdc_no', 'PATH', 'ROW'))
  g_riverice_ts <- merge(gtiles_join, ts, 
                         by=c('PATH', 'ROW'), all.x=T, allow.cartesian=TRUE) %>%
    .[!is.na(river_ice_fraction),] %>%
    merge(elv, by=c('PATH', "ROW")) %>%
    .[, c('grdc_no', 'date', 'river_ice_fraction', 'N_river_pixel', 'mean_elv'),
      with=F]
  
  g_riverice_ts[, grdc_no := as.character(grdc_no)]
  
  return(g_riverice_ts)
}

#------ download_riggs2023 -----------------------------------------------------
download_riggs2023 <- function(in_url, out_file) {
  if (!(dir.exists(dirname(out_file)))) {
    dir.create(dirname(out_file))
  }
  
  #Download and unzip file
  if (!file.exists(out_file)) {
    options(timeout = max(300, getOption("timeout")))
    #https://zenodo.org/records/7150168
    download.file(url = 'https://zenodo.org/records/7150168/files/zenodo.zip',
                  destfil = out_file,
                  method='libcurl')
    unzip(out_file)
  }
  
  return(tools::file_path_sans_ext(out_file))
}

#------ format_riggs2023--------------------------------------------------------
#in_dir <- tar_read(riggs2023_dirpath)
format_riggs2023 <- function(in_dir) {
  dat_dir <- file.path(in_dir, 'supplemented')
  file_list <- list.files(dat_dir,
                          pattern='^[0-9]{7}_GRDC.csv$',
                          full.names=T)
  
  q_bind <- lapply(file_list, function(path) {
    #print(path)
    fread(path, 
          select= c('ID', 'Date', 'RC'))
  }) %>% rbindlist 
  
  setnames(q_bind, c('ID', 'Date', 'RC'), c('grdc_no', 'date', 'Qmod'))
  
  q_bind[, grdc_no := as.character(grdc_no)]
  
  return(q_bind)
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
    .[, `:=`(grdc_no = as.character(grdc_no),
             grdc_no_old = as.character(grdc_no_old))] %>%
    merge(in_g_anthropo_stats, 
          by.x='grdc_no_old', by.y='grdc_no') %>%
    .[, -c('exists', 'wmo_reg', 'lat', 'long', 'Join_Count', 
           'station_river_distance', 'TARGET_FID', 'MAIN_RIV', 'LENGTH_KM',
           'NEXT_DOWN', 'DIST_DN_KM', 'DIST_UP_KM', 'CATCH_SKM',
           'ORD_CLAS', 'ORD_FLOW', "OBJECTID", 'X', 'Y')] %>%
    .[, grdc_match := as.numeric(grdc_no_old != grdc_no)]
  
  gstats_merge_nodupli <- sort(gstats_merge, 'grdc_match') %>%
    .[!duplicated(grdc_no),] %>% #remove duplicates
    .[!is.na(d_yrs) & (d_yrs > min_yrs),]  #Remove those with insufficient daily data (anymore. used to)

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
  
  gstats_merge_nodupli[, grdc_no := as.character(grdc_no)]
  
  return(gstats_merge_nodupli[, -c('grdc_match'), with=F])
}

#------ analyze_anthropo_stats -------------------------------------------------
#in_gmeta_formatted = tar_read(gmeta_formatted)

plot_anthropo_stats <- function(in_gmeta_formatted) {
  #Degree of regulation plot ---------------------------------------------------
  dt_format_dor <- melt_anthropo_stats(in_dt=in_gmeta_formatted,
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
  dt_format_crop <- melt_anthropo_stats(in_dt=in_gmeta_formatted,
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
  dt_format_pop <- melt_anthropo_stats(in_dt=in_gmeta_formatted,
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
  dt_format_built <- melt_anthropo_stats(in_dt=in_gmeta_formatted,
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

#------ filter_reference_gauges ------------------------------------------------
#in_gmeta_formatted = tar_read(gmeta_formatted)
# max_dor = 2
# max_crop = 25
# max_pop = 100
# max_built = 1

filter_reference_gauges <- function(in_gmeta_formatted,
                                    max_dor, max_crop, max_pop, max_built) {
  
  nafill_anthropo <- function(in_dt, fieldroot) {
    fn <- paste0(fieldroot, '_interp')
    dor_dat <- melt_anthropo_stats(in_gmeta_formatted, 
                                   fieldroot = fieldroot, keep_yrs = T) %>%
      .[order(grdc_no, year)] %>%
      .[, (fn) := nafill(value, type='nocb'), by=grdc_no] %>%
      .[, (fn) := nafill(get(fn), type='locf'), by=grdc_no] 
    return(dor_dat[, c('grdc_no', 'year',  fn), with=F])
  }

  anthropo_interpolated <- lapply(
    c('dor', 'pop', 'built', 'crop'), 
    function(in_fieldroot) {
      out_l <- nafill_anthropo(in_gmeta_formatted, in_fieldroot) 
      return(out_l)
    }) %>%
    mergeDTlist(dt_list=., by=c('grdc_no', 'year'), set_suffix = F)
  
  
  gauges_sub <- anthropo_interpolated[
    (dor_interp<=max_dor) &
      (pop_interp<=max_pop) &
      (crop_interp<=max_crop) & 
      (built_interp<=max_built),]
  
  gauges_sub[, `:=`(n_years = .N,
                    grdc_no = as.character(grdc_no)),
             by=grdc_no]
  
  ngauges <- gauges_sub[n_years >= 20, .N, by=year]  
  ref_gauges_ts <- ggplot(ngauges, aes(x=year, y=N)) +
    geom_line(linewidth=1.5) +
    scale_y_log10(breaks=c(1, 10, 100, 1000, 2000, 3000)) +
    theme_ridges()
  
  return(list(
    dt = gauges_sub,
    plot = ref_gauges_ts)
  )
}
#------ prepare_QC_data --------------------------------------------------------
# in_ref_gauges = tar_read(ref_gauges)
# in_gmeta_formatted = tar_read(gmeta_formatted)
# in_geodist = tar_read(geodist)
# in_netdist = tar_read(netdist)
# in_tmax = tar_read(gauges_tmax)
# #in_pdsi_dir = tar_read(pdsi_dir)
# in_rivice = tar_read(rivice_dt)
# in_riggs2023 = tar_read(riggs2023_dt)
# 
# 
# prepare_QC_data_util <- function(
#     
#   )
# 
# prepare_QC_data <- function(in_ref_gauges,
#                             in_gmeta_formatted,
#                             in_geodist,
#                             in_netdist,
#                             in_tmax,
#                             in_rivice,
#                             in_riggs2023) {
#   gauges_anthropo <- in_ref_gauges$dt %>%
#     .[, grdc_no := as.character(grdc_no)]
#     
#   
#   in_grdc_no <- 1634700
#   
#   
#   filepath <- in_gmeta_formatted[grdc_no == in_grdc_no, filename]
#   q_dt <- readformatGRDC(filepath)
#   
#   q_dt_attri <- merge(q_dt, #merge anthropo data
#                       gauges_anthropo[, -c('n_years'), with=F],
#                       by=c('grdc_no','year'),
#                       all.x=T
#                       ) %>%
#     merge(in_tmax[grdc_no == in_grdc_no],
#           by.x=c('grdc_no', 'date'), by.y=c('grdc_no', 'time'))
#   
#   t_max_g <- 
#   
#   
#   in_netdist[grdc_no ==in_grdc_no,]
#   in_geodist[grdc_no ==in_grdc_no,]
#   
#   q_dt
#   
#   #get terra climate data
#   
#   #identify gauges within 100-m from each other
#   
#   #initial value flag
#   
#   #train ARIMA model
#   
#   
# }