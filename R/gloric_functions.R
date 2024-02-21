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
      setnames('YYYY-MM-DD', 'date') %>%
      setorder(grdc_no, date)
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
#------ flagGRDCoutliers ------
#' Flag GRDC outliers
#'
#' Flag potential outliers in daily discharge records for a given GRDC gauging
#' station following the criteria developed for GSIM by
#' [Gudmundsson et al. (2018)](https://essd.copernicus.org/articles/10/787/2018/).
#'
#' @param in_gaugetab \link[data.table]{data.table} containing formatted daily
#' discharge record from GRDC gauging station (as formatted by \code{\link{readformatGRDC}}.
#'
#' @details Criteria to flag a daily discharge value (Qt) as a potential outlier include:
#' \itemize{
#'   \item Negative values (Qt < 0)
#'   \item At least ten identical consecutive discharge values (for Qt > 0)
#'   \item |log(Qt + 0.01) - mean| are larger than the mean values of log(Q + 0.01)
#'   plus or minus 6 times the standard deviation of log(Q + 0.01) computed for
#'   that calendar day for the entire length of the series. The mean and SD are
#'   computed for a 5-day window centred on the calendar day to ensure that a
#'   sufficient amount of data is considered. The log-transformation is used to
#'   account for the skewness of the distribution of daily streamflow values.
#'   \item Qt for which Qobs != Calculated discharge in GRDC record
#' }
#'
#' @return \link[data.table]{data.table} of daily discharge records with additional
#' columns for outlier flags
#'
#' @source Gudmundsson, L., Do, H. X., Leonard, M., & Westra, S. (2018). The Global
#'   Streamflow Indices and Metadata Archive (GSIM) – Part 2: Quality control,
#'   time-series indices and homogeneity assessment. Earth System Science Data,
#'   10(2), 787–804. https://doi.org/10.5194/essd-10-787-2018
#'
#' @export
flagGRDCoutliers <- function(in_gaugetab) {
  in_gaugetab[, `:=`(jday = format(as.Date(date), '%j'),#Julian day
                     q_rleid = rleid(Qobs),#Identify each group of consecutive values
                     flag_mathis = 0)] #Create flag field)
  
  #Flag negative values
  in_gaugetab[Qobs < 0, flag_mathis := flag_mathis + 1]
  
  #Flag when more than 10 identical values in a row or when a single zero-flow
  in_gaugetab[, flag_mathis := flag_mathis +
                ((Qobs > 0) & (.N > 10)) +
                ((Qobs == 0) & (.N == 1)),
              by=q_rleid]
  
  #Flag |log(Q + 0.01) - mean| > 6SD for julian day mean and SD of 5d mean of log(Q + 0.01)
  in_gaugetab[, logmean5d := frollapply(log(Qobs + 0.01), n = 5, align='center',
                                        FUN=mean, na.rm = T)] %>% #Compute 5-day mean of log(Q+0.01)
    .[, `:=`(jdaymean = mean(logmean5d, na.rm = T),
             jdaysd = sd(logmean5d, na.rm = T)),
      by = jday] %>% #Compute mean and SD of 5-day mean of log(Q + 0.01) by Julian day
    .[abs(log(Qobs + 0.01) - jdaymean) > (6 * jdaysd),
      flag_mathis := flag_mathis + 1]
  
  return(in_gaugetab)
}

#------ plotGRDCtimeseries ----------------------
#' Plot a GRDC time series
#'
#' Creates a plot of daily discharge ({m^3}/s) for a GRDC gauging station,
#' with flags for 0-flow values and potential outliers.
#' Save plot to png if path is provided.
#'
#' @param GRDCgaugestats_record \link[data.table]{data.table} of formatted daily
#' discharge records for a single GRDC station. In this project, e.g. the output
#' from \code{\link{comp_GRDCdurfreq}}. Must contain at least five columns:
#' \code{GRDC_NO, date, Qobs, flag_mathis, missingdays} \cr
#' Alternatively, a named list or vector with
#' a column called "path" towards a standard GRDC text file containing daily discharge
#' records for a single gauging station.
#' @param outpath (character) path for writing output png plot (default is no plotting).
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param showmissing (logical) whether to show records in years with number of missing daily records beyond \code{maxgap}.
#'
#' @details the output graphs show the time series of daily streamflow values
#' for the station. For the flagging criteria, see documentation for \code{\link{flagGRDCoutliers}}.
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Individual points show daily discharge values (in {m^3}/s).
#'   \item blue lines link daily values (which may result in unusual patterns due to missing years).
#'   \item red points are zero-flow flow values.
#'   \item green points are non-zero flow daily values statistically flagged as potential outliers .
#'   \item black points are zero-flow values flagged as potential outliers.
#' }
#'
#' @return plot
#'
#' @export
plotGRDCtimeseries <- function(GRDCgaugestats_record,
                               outpath=NULL, maxgap = 366,  showmissing = FALSE) {
  #Read and format discharge records
  if (GRDCgaugestats_record[,.N>1] &
      ('Qobs' %in% names(GRDCgaugestats_record))) {
    gaugetab <- GRDCgaugestats_record
  } else {
    gaugetab <- readformatGRDC(GRDCgaugestats_record$path) %>%
      flagGRDCoutliers %>%
      .[, date := as.Date(date)] %>%
      .[!is.na(Qobs), missingdays := diny(year)-.N, by= 'year']
  }
  
  #Plot time series
  qtiles <- union(gaugetab[, min(Qobs, na.rm=T)],
                  gaugetab[, quantile(Qobs, probs=seq(0, 1, 0.1), na.rm=T)])
  
  subgaugetab <- gaugetab[missingdays < maxgap,]
  
  rawplot <- ggplot(subgaugetab[missingdays < maxgap,],
                    aes(x=date, y=Qobs)) +
    geom_line(color='#045a8d', size=1, alpha=1/5) +
    geom_point(data = subgaugetab[flag_mathis == 0 & Qobs > 0,],
               color='#045a8d', size=1, alpha=1/3) +
    geom_point(data = subgaugetab[flag_mathis > 0 & Qobs > 0,],
               color='green') +
    geom_point(data = subgaugetab[flag_mathis == 0 & Qobs == 0,],
               color='red') +
    geom_point(data = subgaugetab[flag_mathis > 0 & Qobs == 0,],
               color='black') +
    scale_y_sqrt(breaks=qtiles, labels=qtiles) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(y='Discharge (m3/s)',
         title=paste0('GRDC: ', GRDCgaugestats_record$grdc_no)) +
    coord_cartesian(expand=0, clip='off')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text.y = element_text())
  
  if (showmissing) {
    rawplot <- rawplot +
      geom_point(data=gaugetab[missingdays >= maxgap,], color='black', alpha=1/10)
  }
  
  if (!is.null(outpath)) {
    if (!(file.exists(outpath))) {
      ggsave(filename = outpath, plot = rawplot, device = 'png',
             width = 10, height = 10, units='in', dpi = 300)
    }
  } else {
    return(rawplot)
  }
}


############################ ANALYSIS FUNCTIONS ##################################
#------ read_GRDCgauged_paths -----------------
#' Read file paths to streamflow data from GRDC gauging stations
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @param inp_GRDC_qdat_dir path to directory containing streamflow data GRDC standard files.
#' @param in_gaugep table containing column named \code{grdc_no} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables, assuming
#' that files are called "grdc_no.txt", GRDC_NO being replaced with a 7-digit integer.
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
  
  g_riverice_ts[, `:=`(grdc_no = as.character(grdc_no),
                       date = as.character(date)
  )]
  
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
  
  q_bind[, `:=`(grdc_no =as.character(grdc_no),
                date=as.character(date)
  )]
  
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

#------ plot_anthropo_stats -------------------------------------------------
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
#------ prepare_QC_data_util ---------------------------------------------------
# in_grdc_no=3651650
# in_ref_gauges = tar_read(ref_gauges)
# in_gmeta_formatted = tar_read(gmeta_formatted)
# in_geodist = tar_read(geodist)
# in_netdist = tar_read(netdist)
# in_tmax = tar_read(gauges_tmax)
# in_pdsi = tar_read(gauges_pdsi)
# in_rivice = tar_read(rivice_dt)
# in_riggs2023 = tar_read(riggs2023_dt)

prepare_QC_data_util <- function(in_grdc_no,
                                 in_ref_gauges,
                                 in_gmeta_formatted,
                                 #in_geodist,
                                 in_netdist,
                                 in_tmax,
                                 in_pdsi,
                                 in_rivice,
                                 in_riggs2023) {
  #Pre-format anthropo data
  gauges_anthropo <- in_ref_gauges$dt %>%
    .[, grdc_no := as.character(grdc_no)]
  
  #Get basic metadata and discharge observations
  filepath <- in_gmeta_formatted[grdc_no == in_grdc_no, filename]
  q_dt <- readformatGRDC(filepath)
  
  #Merge all data
  q_dt_attri <- merge(q_dt, #merge anthropo data, keeping only reference years
                      gauges_anthropo[, -c('n_years'), with=F],
                      by=c('grdc_no','year'),
                      all.x=F
  ) %>%
    merge(in_tmax[grdc_no == in_grdc_no], #merge tmax
          by.x=c('grdc_no', 'date'), by.y=c('grdc_no', 'time'),
          all.x=T) %>%
    merge(in_pdsi[grdc_no == in_grdc_no,], #merge upstream average PDSI
          by=c('grdc_no', 'date'),
          all.x=T
    ) %>%
    merge(in_riggs2023[grdc_no == in_grdc_no,], #merge Riggs satellite-base discharge
          by=c('grdc_no', 'date'),
          all.x=T
    ) %>%
    merge(in_rivice[grdc_no == in_grdc_no, -'mean_elv',with=F],
          by=c('grdc_no', 'date'),
          all.x=T
    )
  
  #Get discharge  data from upstream and downstream gauges
  netnear_gauges <- in_netdist[ #Make sure these are also ref_gauges
    (grdc_no ==in_grdc_no) & 
      (grdc_no_destination %in% in_ref_gauges$dt$grdc_no),]
  
  #Get the corresponding data
  get_nearg_qobs <- function(in_grdc_no, in_direction){
    path_up <- in_gmeta_formatted[grdc_no == in_grdc_no, filename]
    q_up <- readformatGRDC(path_up) %>%
      .[, c('date', 'Qobs'), with=F] 
    return(q_up)
  }
  
  if (nrow(netnear_gauges)>0) {
    #Create a cast data.table of discharge for those gauges
    near_gaugesq <- netnear_gauges[
      , get_nearg_qobs(in_grdc_no=grdc_no_destination, 
                       in_direction=dist_net_dir)
      , by=c('grdc_no_destination', 'dist_net_dir')] %>%
      .[, colname := paste0('Qobs_', dist_net_dir, '_', grdc_no_destination)] %>%
      dcast(date~colname,
            value.var='Qobs'
      )
    
    #Merge them to the main table
    q_dt_attri <-  q_dt_attri %>% 
      merge(near_gaugesq, by='date', all.x=T) %>%
      .[, date := as.Date(date)]
    
    #Sort near gauges based on overlap and similarity in median Q----------
    nearg_cols <- grep('Qobs_(down|up)stream_travel_.*',
                       names(q_dt_attri),
                       value=T)
    
    #Check number of years where both the target ts has less than 30 missing days
    # and number of years where the ref ts has less than 30 missing days
    near_gcols_stats <- lapply(nearg_cols, function(gcol) {
      #print(gcol)
      ntotal_yrs <-  q_dt_attri[, .SD[(sum(is.na(Qobs))<30),], by=year] %>%
        .[, length(unique(year))]
      common_yrs <- q_dt_attri[, .SD[(sum(is.na(get(gcol)))<30) &
                                       (sum(is.na(Qobs))<30),], by=year]
      ncommon_yrs <- common_yrs[, length(unique(year))]
      
      if (ncommon_yrs > 5) {
        cc_all <- common_yrs[, ccf(ts(log(get(gcol)+0.01)), ts(log(Qobs+0.01)), 
                                   na.action=na.pass,
                                   i = 10,
                                   plot=F)]
        max_cc <- max(cc_all$acf)
      } else {
        max_cc <-NA
      }
      
      
      return(data.table(
        col=gcol,
        ncommon_yrs=ncommon_yrs,
        ncommon_ratio=ncommon_yrs/ntotal_yrs,
        max_cc=max_cc
      ))
    }) %>% rbindlist
    
    near_gcols_sel <- near_gcols_stats[ncommon_yrs>5 & ncommon_ratio>0.75,] %>%
      .[order(-max_cc),col]
  } else {
    near_gcols_sel <- NULL
  }
  
  #Format a few columns
  q_dt_attri[date > as.Date('1958-01-01'), 
             `:=`(
               tmax = nafill(tmax, type="locf"),
               PDSI = nafill(PDSI, type="nocb")
             )]
  
  #Check whether there are 0
  potential_npr <- q_dt_attri[, .SD[Qobs<0.01,.N]>0]
  
  #Computer the percentage of integer values
  integer_perc <- q_dt_attri[!is.na(Qobs), 
                             sum(fifelse(Qobs == round(Qobs), 1, 0))/.N]
  
  
  return(list(
    q_dt_attri=q_dt_attri,
    near_gcols_sel = near_gcols_sel,
    potential_npr = potential_npr,
    integer_perc = integer_perc
  ))
}

#------ train_qARIMA ------------------------------------------------------------
#Simplification of check residuals from forecast package
LBtest_util <- function(in_mod, lag) {
  moddf <-  sum(arimaorder(in_mod)[c("p","q","P","Q")], na.rm = TRUE)
  modres <- residuals(in_mod)
  freq <- frequency(residuals)
  
  if (missing(lag)) {
    lag <- ifelse(freq > 1, 2 * freq, 10)
    lag <- min(lag, round(length(modres)/5))
    lag <- max(moddf + 3, lag)
  }
  
  LBtest <- Box.test(zoo::na.approx(modres), fitdf = moddf, 
                     lag = lag, type = "Ljung")
  LBtest$method <- "Ljung-Box test"
  LBtest$data.name <- "Residuals"
  names(LBtest$statistic) <- "Q*"
  
  return(LBtest)
}


run_qARIMA <- function(in_q_dt, in_qcol, nearg_cols, obs_bclambda, max_K) {
  setDT(in_q_dt)
  q_dt <- copy(in_q_dt)
  q_ts <- ts(q_dt[, get(in_qcol)], frequency=365.25)
  
  #Identify outliers with ARIMAX -----------------------------------------------
  #If there is a nearby gauge
  # by default, set this to FALSE to not throw an error if there is no nearby gauge
  if ((length(nearg_cols) >= 1) & !any(is.na(nearg_cols))) {
    nearg <- T
    test_without_nearg <- F
  } else {
    test_without_nearg <- T
    nearg <- F
  }
  
  if (!test_without_nearg) {
    #Get ts object of transformed data for nearby gauge
    q_near <- q_dt[, get(nearg_cols)+0.01] %>%
      forecast::BoxCox(lambda='auto') %>%
      ts(frequency=365.25)
    nearg_bclambda <- attr(q_near, 'lambda')
    
    #Determine lag time with highest cross-correlation between the two gauges
    maxccf_lag_nearg <- ccf(q_ts, q_near, 
                            na.action=na.pass, plot=F)$acf %>%
      as.data.table %>%
      .[, shift := seq_along(.I)-median(.I)] %>%
      .[which.max(value), -shift]
    
    #Find best ARIMA model
    bestfit_nearg <- list(aicc=Inf)
    break_loop <- F
    for(i in 1:max_K) #Select the number of Fourier series by minimizing AICc
    {
      tryCatch(
        fit <- auto.arima(q_ts, 
                          seasonal=FALSE, 
                          xreg=cbind(fourier(q_ts, K=i),
                                     lag(q_near, maxccf_lag_nearg)
                          )
        ),
        error = function(e) { #If no suitable ARIMA model
          break_loop <- T
        }
      )
      if (break_loop==T) {
        break
      }
      if(fit$aicc < bestfit_nearg$aicc) {
        print(i)
        bestfit_nearg <- fit #Keep model with lowest AIC
      } else {break}
    }
    
    #Test whether the residuals are stationary and uncorrelated
    restests_nearg <- LBtest_util(bestfit_nearg)
    if (restests_nearg$p.value > 0.05) {
      test_without_nearg <- F
    } 
  }
  
  if ((!nearg) | (test_without_nearg)) {
    #Find best ARIMA model (without other gauge)
    bestfit_nonearg <- list(aicc=Inf)
    break_loop <- F
    for(i in 1:max_K) #Select the number of Fourier series by minimizing AIC
    {
      tryCatch(
        fit <- auto.arima(q_ts, 
                          seasonal=FALSE, 
                          xreg=cbind(fourier(q_ts, K=i)
                          )
        ),
        error = function(e) { #If no suitable ARIMA model
          break_loop <- T
        }
      )
      if (break_loop==T) {
        break
      }
      if(fit$aicc < bestfit_nonearg$aicc){
        print(i)
        bestfit_nonearg <- fit #Keep model with lowest AIC
      }else {break};
    }
    
    restests_nonearg <- LBtest_util(in_mod=bestfit_nonearg)
  }
  
  #Get the model with whitenoise residuals or lowest AICc
  if (!test_without_nearg) { #If there was no testing without a nearby gauge
    bestfit <- bestfit_nearg
  } else {
    if (nearg) { #If there was testing with and without a nearby gauge
      if (restests_nearg$p.value < 0.05 
          & restests_nonearg$p.value > 0.05) { #If the residuals without a nearby gauge are whitenoise and those with a nearby gauge are not
        bestfit <- bestfit_nonearg
      } else { #Otherwise, pick the model with the lowest AICc
        if (arimaorder(bestfit_nearg)['d'] == arimaorder(bestfit_nonearg)['d']) { #can only compare AICc if same number of differencing
          bestfit <- list(bestfit_nearg, bestfit_nonearg) %>%
            .[which.min(lapply(., function(x) x$aicc))] %>%
            .[[1]]
        } else {
          bestfit <- bestfit_nearg
        }
      }
    } else {
      bestfit <- bestfit_nonearg
    }
  }
  LBtest <- checkresiduals(bestfit, plot=F)
  
  #Get one-step ahead forecast with confidence interval
  sigma2 <- bestfit$sigma2
  q_dt$fit <- as.numeric(fitted(bestfit, h=1))
  q_dt$fit_l95 <- InvBoxCox(q_dt$fit - sqrt(sigma2)*1.96, lambda=obs_bclambda)-0.01
  q_dt$fit_u95 <- InvBoxCox(q_dt$fit + sqrt(sigma2)*1.96, lambda=obs_bclambda)-0.01
  q_dt$fit_l99 <- InvBoxCox(q_dt$fit - sqrt(sigma2)*2.68, lambda=obs_bclambda)-0.01
  q_dt$fit_u99 <- InvBoxCox(q_dt$fit + sqrt(sigma2)*2.68, lambda=obs_bclambda)-0.01
  q_dt$fit <- InvBoxCox(q_dt$fit, lambda=obs_bclambda)-0.01
  q_dt$obsfitdiff <- q_dt[, InvBoxCox(get(in_qcol), lambda=obs_bclambda)] - q_dt$fit
  
  #Make sure the predictions are bounded to 0 (because of the 0.01)
  q_dt[fit<0, fit:=0]
  q_dt[fit_l95<0, fit_l95:=0]
  q_dt[fit_l99<0, fit_l99:=0]
  
  return(list(
    fit_dt=q_dt[, .(date, fit, fit_l95, fit_u95, fit_l99, fit_u99, obsfitdiff)],
    mod=bestfit,
    LBtest=LBtest,
    sigma2=sigma2)
  )
}


#------ detect_outliers_ts ----------------------------------------------------
# data_for_qc<- tar_read(data_for_qc)
# in_row <- data_for_qc[potential_npr==T & integer_perc<0.95,][77,]
# in_data <- qread(in_row[, qs_path][[1]])
# in_nearg_cols <- in_row[, near_gcols_sel][[1]]
# arima_split=F
# plot_fit = F
# 
# detect_outliers_ts(in_data, in_nearg_cols, arima_split=F, plot_fit=F) 

detect_outliers_ts <- function(in_data, in_nearg_cols, arima_split=F, plot_fit=F) {
  #Remove negative values
  ts_cols <- c('date', 'Qobs', 'year', 'missingdays','PDSI')
  in_data[Qobs < 0, Qobs := NA]
  first_nona_ix <- in_data[,min(which(!is.na(Qobs)))]
  if (first_nona_ix > 1) {
    in_data <- in_data[!seq(1,first_nona_ix),]
  }
  
  #Transform discharge and create time series object
  obs_bclambda  <- BoxCox.lambda(ts(in_data$Qobs+0.01, 
                                    frequency=365.25)) #Determine BoxCox transformation lambda
  in_data[, Qobs_trans := BoxCox(Qobs+0.01, lambda=obs_bclambda)]
  
  if (length(in_nearg_cols) > 1) {
    in_nearg_cols <- in_nearg_cols[[1]]
  }
  
  #Run ARIMA model to detect potential outliers---------------------------------
  #start <- Sys.time()
  max_K <- ifelse((nrow(in_data)<(365.25*80)), 10, 3) 
  qARIMA <- run_qARIMA(in_q_dt=in_data, in_qcol='Qobs_trans', 
                       nearg_cols=in_nearg_cols, obs_bclambda=obs_bclambda,
                       max_K=max_K)
  #}
  #print(Sys.time()-start)
  #ARIMA is good at detecting sudden peaks, but classifies them all as peaks, 
  #STL decomposition-based outlier detecting seems better at detecting 
  #greater peakiness by season
  
  #Join ARIMA to dt 
  in_data <- merge(in_data, qARIMA$fit_dt, by='date', all.x=T)
  in_data[, date := as.Date(date)] %>%
    .[, jday := as.numeric(format(date, '%j'))]
  
  #Find differences from forecast outside of the julian day 90% interval
  obsfitdiff_roll_dt <- lapply(seq(366), function(i) {
    dt_processed <- in_data[
      jday==i, 
      list(i,
           obsfitdiff_roll_l = .SD[(jday >= i-7) & (jday <= i +7),
                                   quantile(obsfitdiff, 0.025, na.rm=T)],
           obsfitdiff_roll_u = .SD[(jday >= i-7) & (jday <= i +7),
                                   quantile(obsfitdiff, 0.975, na.rm=T)]
      )
    ]
    return(dt_processed)
  }) %>% rbindlist
  
  in_data <- merge(in_data, obsfitdiff_roll_dt, by.x='jday', by.y='i', all.x=T) %>%
    .[order(date),]
  in_data[, arima_outlier_rolldiff := fifelse(
    ((Qobs < fit_l99) & (obsfitdiff<0) & (obsfitdiff<(obsfitdiff_roll_l*2))) | 
      ((Qobs > fit_u99) & (obsfitdiff>0) & (obsfitdiff>(obsfitdiff_roll_u*2))), 1, 0)]
  
  #Detect outliers through periodic STL decomposition --------------------------
  in_data_sub <- in_data[missingdays < 90,]
  if (nrow(in_data_sub) > 1) {
    stl_outliers <- tsoutliers(ts(in_data_sub$Qobs_trans, frequency=365.25),
                               iterate = 1)
    in_data[date %in% in_data_sub[stl_outliers$index, date], stl_outlier := 1]
  } else {
    in_data[, stl_outlier := 0]
  }
  
  #Detect potential outliers through hard rules --------------------------------
  flagGRDCoutliers(in_data)
  setnames(in_data, 'flag_mathis', 'auto_flag')
  
  #Detect abnormally smooth stretches ------------------------------------------
  #Identify periods of at least 7 days whose CV in second order difference
  # is below the 10th percentile for their respective calendar day
  in_data[, Qobs_diff2 := c(NA, NA, diff(Qobs_trans, differences=2))] %>%
    .[, Qdiff2_rollcv := frollapply(Qobs_diff2, n=10, sd)/abs(frollmean(Qobs_diff2, n=10))] 
  in_data[, Qdiff2_rollcv_jdayq90 := quantile(Qdiff2_rollcv, 1/10, na.rm=T), by=jday]
  in_data[, smooth_flag := (.N>7)&(Qdiff2_rollcv< Qdiff2_rollcv_jdayq90), 
          by=rleid(Qdiff2_rollcv< Qdiff2_rollcv_jdayq90)]
  
  # ggplotly(ggplot(in_data, aes(x=date, y=Qobs_trans)) +
  #   geom_line() +
  #   geom_point(data=in_data[smooth_flag,], 
  #              aes(color= smooth_flag))
  # )
  
  #Mark potential outliers -----------------------------------------------------
  # in_data[, `:=`(
  #   arima_outlier_95 = fifelse((Qobs < fit_l95) | (Qobs > fit_u95), 1, 0),
  #   arima_outlier_99 = fifelse((Qobs < fit_l99) | (Qobs > fit_u99), 1, 0)
  # )]
  
  in_data[, all_flags := fcase(
    stl_outlier==1, 'Stl flag',
    auto_flag>0, 'Auto flag',
    smooth_flag==T, 'Smooth flag',
    arima_outlier_rolldiff==1, 'Arima flag',
    default='No flag')]
  
  #Plot outliers ---------------------------------------------------------------
  if (plot_fit) {
    #Format outliers for plotting
    all_flags <- melt(
      in_data[(stl_outlier==1)|(arima_outlier_99==1)|(auto_flag>0)|smooth_flag
              |(arima_outlier_rolldiff==1),],
      id.vars=c('date', 'jday','year', 'Qobs'),
      measure.vars = c('stl_outlier', 'arima_outlier_rolldiff',
                       'auto_flag', 'smooth_flag'))  %>%
      .[!(value %in% c(NA, 0)),]
    
    p_rect_dat <- in_data[!is.na(PDSI), .(date, PDSI)] %>%
      .[, trimester := lubridate::round_date(date, unit='3 months')] %>%
      .[!duplicated(trimester),] %>%
      .[order(trimester), end_date := .SD[.I+1, trimester]-1]
    
    p_ymin <- in_data[, min(c(min(Qobs,na.rm=T), 
                              min(fit_l99, na.rm=T)))]
    
    p_fit <- ggplot(data=in_data) +
      geom_rect(data=p_rect_dat,
                aes(xmin=trimester, xmax=end_date, 
                    ymin=p_ymin, ymax=Inf, 
                    fill = PDSI), alpha=1/2) +
      geom_ribbon(aes(x=date, ymin=fit_l99, ymax=fit_u99), color='darkgrey') + 
      geom_line(aes(x=date, y=Qobs), size=1) +
      #geom_line(aes(x=date, y=fit),  color='orange') +
      geom_point(data=all_flags, aes(x=date, y=Qobs, color=variable)) +
      scale_fill_distiller(name='PDSI', palette='RdBu', direction=1) +
      scale_y_sqrt() + 
      theme_classic()
    
    in_data[, grp_int := interaction(format(date, '%Y%m'), 
                                     as.numeric(factor(PDSI, exclude = 999)))]
    
    p_seasonal <-  ggplot(in_data[!is.na(Qobs)], 
                          aes(x=as.numeric(jday), y=Qobs)) + 
      # geom_ribbon(aes(ymin=fit_l95, ymax=fit_u95, 
      #                 fill=PDSI, group=grp_int), 
      #             alpha=1/5) + 
      #geom_point(aes(color=year), alpha=1/4) +
      geom_line(aes(color=PDSI, group=grp_int), linewidth=1) +
      #geom_line(aes(x=jday, y=get(in_nearg_cols[[1]]), group=year), color='blue') +
      scale_color_distiller(palette='Spectral', direction=1) +
      scale_fill_distiller(palette='Spectral', direction=1) +
      ggnewscale::new_scale_color() +
      geom_point(data=all_flags, aes(color=variable), alpha=1/2, size=2) +
      scale_color_brewer(palette='Dark2') +
      scale_y_sqrt() +
      theme_classic()
    
    p_fit_forplotly <- ggplot(data=in_data) +
      geom_point(aes(x=date, y=Qobs), alpha=1/3) +
      geom_line(aes(x=date, y=Qobs), linewidth=1.2) +
      #geom_line(aes(x=date, y=fit),  color='orange') +
      geom_point(data=all_flags, aes(x=date, y=Qobs, color=variable)) +
      scale_x_date(date_breaks="1 year") +
      scale_y_sqrt() + 
      theme_classic()
  } else {
    p_fit = NULL
    p_seasonal = NULL
    p_fit_forplotly = NULL
  }
  
  #return statement ------------------------------------------------------------
  return(list(
    outliers_dt = in_data[,
                          .(grdc_no, date, jday, Qobs, year, integervalue, 
                            missingdays, dor_interp, pop_interp, built_interp,
                            crop_interp, tmax, PDSI, Qmod, river_ice_fraction,
                            jdaymean, jdaysd, all_flags)], 
    arima_model = qARIMA$mod,
    p_fit = p_fit,
    p_seasonal = p_seasonal,
    p_fit_forplotly = p_fit_forplotly
  ))
}


#------ na_interp_dt_custom ----------------------------------------------------
na_interp_dt_custom <- function(in_dt, in_var, in_freq=365.25, in_floor=0) {
  #Convert negative values to NA
  in_dt[get(in_var) < in_floor, (in_var) := NA]
  
  #For interpolation, get rid of leading NAs
  first_nona_ix <- in_dt[,min(which(!is.na(get(in_var))))]-1
  trans_fn <- paste0(in_var, '_trans')
  interp_fn <- paste0(in_var, '_interp')
  
  #Interpolate with robust STL decomposition on pre-transformed data
  obs_bclambda  <- BoxCox.lambda(ts(in_dt[!seq(0,first_nona_ix),][[in_var]]+0.01, 
                                    frequency=in_freq)) #Determine BoxCox transformation lambda
  in_dt[, (trans_fn) := BoxCox(get(in_var)+0.01, lambda=obs_bclambda)]
  
  in_dt[!seq(0,first_nona_ix),  
        (interp_fn) := as.vector(
          InvBoxCox(
            forecast::na.interp(ts(get(trans_fn), 
                                   frequency=in_freq)
            ),
            lambda = obs_bclambda)-0.01
        )]
  in_dt[get(interp_fn)<in_floor, (interp_fn):=in_floor]
  
  #Compute length of each interpolated period
  in_dt[, NAperiod := rleid(is.na(get(in_var)))] %>%
    .[is.na(get(in_var)), NAperiod_n := .N, by=NAperiod]
  in_dt[, NAperiod := NULL]
  
}


#------ compute metastatistics_util --------------------------------------------
# in_outliers_output_dt <- tar_read(q_outliers_flags)
# lapply(unique(in_outliers_output_dt, by=c('grdc_no'))[,grdc_no], function(in_no) {
#   print(in_no)
#   in_outliers_path<- in_outliers_output_dt[grdc_no==in_no, out_qs][[1]]
#   n_freq=365.25
#   in_floor=0
#   compute_metastatistics_util(in_outliers_path, in_no)
# })

compute_metastatistics_util <- function(in_outliers_path, in_no) {
  in_dt_edit_path <- paste0(
    tools::file_path_sans_ext(in_outliers_path),
    '_edit.csv')
  in_dt_clean <- fread(in_dt_edit_path, select=c('date', 'Qobs', 'year'))
  
  if (in_dt_clean[, sum(!is.na(Qobs) & (Qobs > 0))] > 1) {
    na_interp_dt_custom(in_dt=in_dt_clean,
                        in_var='Qobs')
    
    # ggplot(in_dt_clean, aes(x=date, y=Qobs_interp, color=NAperiod_n, group=NAperiod)) +
    #   geom_line(color='grey') +
    #   geom_line(data=in_dt_clean[is.na(Qobs)], size=1.1) +
    #   theme_bw()
    
    #Compute number of missing days depending on the maximum interpolated gap
    #as well as the number of "no-flow" records depending on flow threshold
    metastats_dt <- in_dt_clean[, list(
      missingdays_edit = sum(is.na(Qobs)),
      missingdays_edit_interp5 = .SD[is.na(Qobs) & NAperiod_n > 5, .N],
      missingdays_edit_interp7 = .SD[is.na(Qobs) & NAperiod_n > 7, .N],
      missingdays_edit_interp10 = .SD[is.na(Qobs) & NAperiod_n > 10, .N],
      ndays_Qu5 = sum(Qobs < 0.005, na.rm=T),
      ndays_Que1 = sum(Qobs <= 0.001, na.rm=T),
      ndays_zeroflow = sum(Qobs==0, na.rm=T)
    ), by=year] %>%
      .[,  `:=`(
        grdc_no = in_no,
        edited_data_path = in_dt_edit_path,
        Q50 = in_dt_clean[
          year %in% .SD[missingdays_edit_interp10 < 15,year],
          quantile(Qobs_interp, 0.5, na.rm=T)]
      )] 

    return(metastats_dt)
  }
}

#------ compute metastatistics wrapper -----------------------------------------
#in_outliers_output_dt <- tar_read(q_outliers_flags)

compute_metastatistics_wrapper <- function(in_outliers_output_dt) {
    out_dt <- unique(in_outliers_output_dt, by=c('grdc_no'))[,
      compute_metastatistics_util(
        in_outliers_path=out_qs,
        in_no=grdc_no),
      by=grdc_no
    ]
  return(out_dt)
}
  
#------ plot_metastats -----------------------------------------------
#in_metastats_dt <- tar_read(metastats_dt)

analyze_metastats <- function(in_metastats_dt) {
  #Check which ones are non-perennial (npr) 
  # have at least one day < 5L/s if Q50>1 or <= 1L/s if Q50<1
  in_metastats_dt[, unique(grdc_no)]
  
  metastats_npr <- in_metastats_dt[, list(
    npr=max(ndays_Que1)>0), by='grdc_no'] %>%
    .[npr==T, unique(grdc_no)]
  
   # metastats_npr <- in_metastats_dt[, list(npr = fifelse(
   #  Q50>1, 
   #  max(ndays_Qu5)>0, 
   #  max(ndays_Que1)>0)),
   #  by=grdc_no] %>%
   #  .[npr==T, unique(grdc_no)]
  
  metastats_pformat <- melt(
    in_metastats_dt[grdc_no %in% metastats_npr,],
    id.vars=c('grdc_no', 'year'),
    measure.vars = grep('missingdays_', names(in_metastats_dt), value=T)
  ) %>%
    .[, variable := as.numeric(gsub('missingdays_edit(_interp)*', '', 
                                    variable))] %>%
    setnames('variable', 'max_interp') %>%
    .[is.na(max_interp), max_interp:=0]
  
  metastats_nyears <- lapply(seq(0,30, 5), function(max_miss) {
    out_stats <- metastats_pformat[value<=max_miss, 
                            list(nyears=.N),
                            by=c('grdc_no', 'max_interp')]
    out_stats[, max_miss := max_miss]
    return(out_stats)
  }) %>% rbindlist
  
  metastats_ngauges <- lapply(seq(0,30, 5), function(min_nyears) {
    out_stats <- metastats_nyears[nyears>=min_nyears, 
                           list(ngauges=.N),
                           by=c('max_miss','max_interp')]
    out_stats[, min_nyears := min_nyears]
    return(out_stats)
  }) %>% rbindlist

  plot_ngauges <- ggplot(metastats_ngauges, 
                         aes(x=min_nyears, y=max_miss, color=ngauges)) +
    geom_text(aes(label=ngauges)) +
    facet_wrap(~max_interp) +
    scale_color_distiller(palette='RdPu') +
    theme_bw()
  
  # ggplot(metastats_ngauges, 
  #        aes(x=min_nyears, y=ngauges, color=max_miss, group=max_miss)
  #        )+
  #   geom_line() +
  #   scale_color_distiller(palette='PuBu') +
  #   facet_wrap(~max_interp) +
  #   theme_bw()
  
  metastats_meanyrs <- metastats_nyears[nyears>15, 
                                        .SD[, mean(nyears)],
                                        by=c('max_miss', 'max_interp')]
  return(list(
    plot_ngauges = plot_ngauges,
    metastats_meanyrs = metastats_meanyrs,
    metastats_nyears = metastats_nyears
  ))
}

#------ compute_hydrostats_utils -------------------------------------
in_metastats_dt <- tar_read(metastats_dt)
in_metastats_analyzed <- tar_read(metastats_analyzed)
in_gaugep_dt <- tar_read(gaugep_dt)

max_interp_sel = 10
max_miss_sel = 0
min_nyears = 15
q_thresh = 0.001

interp_fn <- paste0('missingdays_edit',
                    fifelse(max_interp_sel >0,
                            paste0('_interp', max_interp_sel),
                            '')
)

gauges_sel_no <- in_metastats_analyzed$metastats_nyears[
  max_interp==max_interp_sel & max_miss==max_miss_sel & nyears>=min_nyears,
  grdc_no]

in_no <- gauges_sel_no[1]
#For given gauge, get years with less than the maximum number of missing days
#given maximum interpolation period
meta_sub_yrs <- in_metastats_dt[grdc_no==in_no
                                & get(interp_fn)<=max_miss_sel,]

#Read data, keeping only years with number of missing days under threshold
q_dt <- fread(meta_sub_yrs$edited_data_path[[1]]) %>%
  .[year %in% meta_sub_yrs$year, .(grdc_no, date, jday, Qobs, year, tmax, PDSI)] %>%#subset columns for speed
  .[min(which(!is.na(Qobs))):max(which(!is.na(Qobs))),] #remove leading and trailing NAs

#Interpolate discharge, keeping only interpolation periods under threshold
na_interp_dt_custom(in_dt=q_dt,
                    in_var='Qobs')
q_dt[NAperiod_n > max_interp_sel, Qobs_interp := NA]

in_lat <- in_gaugep_dt[grdc_no==in_no]$y_geo

#Prepare data for computing statistics -----------------------------------------
#Uniquely identify no-flow period
q_dt[, noflow_period := rleid(Qobs_interp <= q_thresh)] %>%
  .[Qobs_interp > q_thresh, noflow_period := NA] %>%
  .[!is.na(noflow_period), noflow_period_dur := .N, by = noflow_period]

#Identify continuous blocks of 5 years
whole_5yrblock <- q_dt[order(date) & !duplicated(year),] %>%
  .[, year_rleid :=  cumsum(c(1, diff(year) != 1))] %>%
  .[, year_5blockid := ceiling(seq_along(year)/5), by='year_rleid'] %>%
  .[, year_5blockidn := .N, by=c('year_rleid', 'year_5blockid')] %>%
  .[year_5blockidn==5, blockid := 10*year_rleid+year_5blockid]
q_dt <- merge(q_dt, whole_5yrblock[, .(year, blockid)], by='year')

#Convert each day i with no flow into an angular (ti) and represent it
#by a unit vector with rectangular coordinates (cos(ti); sin(ti))
q_dt[!is.na(noflow_period), `:=`(cos_t=cos(2*pi*(jday-1)/(diny(year)-1)),
                                 sin_t=sin(2*pi*(jday-1)/(diny(year)-1))
                                 )]

#Identify contiguous six months with the most zero-flow days for computing Sd6
identify_drywet6mo <- function(in_dt, q_col='Qobs_interp', 
                               q_thresh=0.001, jday_col = 'jday') {
  
  #Add 3 months before and after record to avoid having NAs on the edges
  fill_dt <- data.table(rep(NA,91)) %>% setnames(q_col)
  
  in_dt <- rbind(fill_dt,
                rbind(in_dt, fill_dt, 
                      fill=T),
                fill=T)
  
  #Computer total number of no-flow days in the 6-month period centered around
  #each day in the record
  in_dt <- in_dt[, noflow_6morollsum := frollsum(
    get(q_col) <= q_thresh,
    n=183, align='center', hasNA=T, na.rm=T)] %>%
    .[ 92:(.N-91), ]
  
  #Find the Julian day with the most number of no-flow days on interannual average
  driest_6mocenter <- in_dt[
    , as.Date(
      unique(get(jday_col))[
        which.max(.SD[, mean(noflow_6morollsum, na.rm=T),by=jday_col]$V1)],
      origin=as.Date("1970-01-01"))
  ]
  
  #Identify the 6-month period centered on that julian day as the dry period
  driest_6moperiod <- data.table(
    dry_6mo = T,
    jday = as.numeric(format(
      seq(driest_6mocenter-91, driest_6mocenter+91, by='day'), '%j'))
  )
  
  #Identify the other julian days as the wet period
  in_dt <- merge(in_dt, driest_6moperiod, 
                 by.x=jday_col, by.y='jday', all.x=T) %>%
    .[is.na(dry_6mo), dry_6mo := F]
  
  in_dt[, noflow_6morollsum := NULL]
  
  return(in_dt)
}
q_dt <- identify_drywet6mo(q_dt)

#Compute statistics ------------------------------------------------------------
q_stats <- q_dt[, list(
  #Intermittence
  f0 = sum(Qobs_interp <= 0.001)/.N,
  
  #Duration
  meanD = mean(
    .SD[!is.na(noflow_period) & !duplicated(noflow_period), noflow_period_dur]),
  medianD = median(
    .SD[!is.na(noflow_period) & !duplicated(noflow_period), noflow_period_dur]),
  sD = sd(
    .SD[!is.na(noflow_period) & !duplicated(noflow_period), noflow_period_dur]),
  #d80 is computed, excluding all years not within continuous blocks of 5 years 
  #in chronological order (i.e., no re-ordering within continuous periods of 
  #record to minimize or maximize d80)
  d80 = .SD[!is.na(blockid), 
           fifelse(sum(!is.na(noflow_period_dur)) > 0,
                   max(noflow_period_dur, na.rm=T), 0)
           , by=blockid][, ceiling(quantile(V1, 0.8))],
  
  #Frequency
  meanN = mean(
    .SD[, uniqueN(noflow_period, na.rm=T), by=year]$V1),
  medianN = median(
    .SD[, uniqueN(noflow_period, na.rm=T), by=year]$V1),
  sdN = sd(
    .SD[, uniqueN(noflow_period, na.rm=T), by=year]$V1),
  
  #Timing
  theta = atan2(mean(cos_t, na.rm=T), 
                mean(sin_t, na.rm=T)),
  r = sqrt(mean(cos_t, na.rm=T)^2 
           + mean(sin_t, na.rm=T)^2)
)]

#Adjust theta for stations in southern hemisphere. To avoid discontinuities, 
#shift the seasons progressively up to 23.5 degrees. 
#After 23.5, shift by a full half year 
#i.e. Meteorological summer starts December 1st instead of June 1st
q_stats[, theta := fifelse(in_lat >= 0,
                           theta,
                           (theta - pi*max(-1, (in_lat/23.5)))%%(2*pi)
                           )]

#Compute seasonal predictability of no-flow events (Sd6)
q_stats$Sd6 <- q_dt[, ym := format(date, '%Y%m')] %>%
  .[, any(!is.na(noflow_period)), by=.(year, ym, dry_6mo)] %>%
  .[, sum(V1, na.rm=T), by=.(year, dry_6mo)] %>% #Compute number of months with zero-flows for each wet and dry period and year
  .[, 1-(.SD[!dry_6mo,mean(V1)]/.SD[dry_6mo,mean(V1)])]



# Rate of change	Drec	Seasonal recession time scale (Catalogne, 2012) (day)
# Ic	Concavity index derived from the flow duration curve, introduced by Sauquet and Catalogne (2011) (*) (dimensionless)
# BFI	Baseflow index computed with the smoothed minima method introduced by the Institute of Hydrology (1980) (dimensionless)
# medianDr	Median duration of runoff event (*) (day)








# dat_qs_list <- tar_read(data_for_qc)[potential_npr==T & integer_perc < 0.80, qs_path]
# as.data.table(qread(dat_qs_list[[1]]) )
# 
# outlier_qs_list <- file.path(temp_qs_dir,
#                              gsub('data_for_qc_', 'q_outliers_flags_',
#                                   basename(dat_qs_list)))
# 
# for (i in seq(1, 20)) {
#   in_dat_path <- dat_qs_list[[i]]
#   in_qs <- file.path(temp_qs_dir,
#                      gsub('data_for_qc_',
#                           'q_outliers_flags_',
#                           basename(in_dat_path))
#   )
#   in_dat <- qread(in_dat_path)
#   integer_perc <- in_dat[!is.na(Qobs),
#                          sum(fifelse(Qobs == round(Qobs), 1, 0))/.N]
#   print(integer_perc)
# }
# 
# in_dat_path <- dat_qs_list[[5]]
# in_qs <- file.path(temp_qs_dir,
#                    gsub('data_for_qc_',
#                         'q_outliers_flags_',
#                         basename(in_dat_path))
# )
# in_dat <- qread(in_dat_path)$q_dt_attri
# in_outliers_qs <- qread(in_qs)
# in_outliers <- in_outliers_qs$outliers_dt

# in_outliers_qs$p_fit_forplotly$data



################### EXTRA STUFF ################################################
# ggplot(q_dt, aes(x=jday, y=Qobs, group=year)) + 
#   geom_ribbon(aes(ymin=fit_l95, ymax=fit_u95, fill=year), alpha=1/5) + 
#   #geom_point(aes(color=year), alpha=1/4) +
#   geom_line(aes(color=year)) +
#   scale_color_distiller(palette='Spectral') +
#   scale_fill_distiller(palette='Spectral') +
#   ggnewscale::new_scale_color() +
#   geom_point(data=all_flags, aes(color=variable), alpha=1/2) +
#   scale_y_sqrt() +
#   theme_classic()

# kr <- KalmanSmooth(q_ts, bestfit$model) #Impute missing values with Kalman Smoother (from https://stats.stackexchange.com/questions/104565/how-to-use-auto-arima-to-impute-missing-values)
# id.na <- which(is.na(q_ts)) #Limit to gaps < 9 months
# pred <- q_ts
# for (i in id.na)
#   pred[i] <- bestfit$model$Z %*% kr$smooth[i,]
# q_dt[id.na,'Qinterp'] <- pred[id.na] #Replace NA values in the time series by predicted values
# #precip1KA41sub[precip1KA41sub$Flow<0.05 & !is.na(precip1KA41sub$Flow),'Flow'] <- 0 #For values < 0 and improbably low values, replace with 0
# ggplot(q_dt, aes(x=date, y=Qobs)) + geom_point() + #Plot result
#   geom_point(data=q_dt[id.na,], aes(y=Qinterp), color='red')


#Extra fable/forecast functions
# autoplot(q_ts, log(Qobs+0.1))   #Standard ts plot
# gg_season(q_ts, log(Qobs+0.1))  #Seasonal plot
# GGally::ggpairs(as.data.table(q_ts),  #Correlation matrix
#                 columns=grep('Qobs.*', names(q_ts)))
# ACF(q_ts, Qobs, lag_max=30) %>% autoplot() #correlogram for the past month
# ACF(q_ts, Qobs, lag_max=365*5) %>% autoplot() #correlogram for the past 5 years (trend)
# PACF(q_ts, Qobs, lag_max=30) %>% autoplot() #correlogram for the past month
# PACF(q_ts, Qobs, lag_max=365) %>% autoplot() #correlogram for the past month

# check <- tslm(ts(q_ts$Qobs, frequency=365.25) ~ trend + season, lambda='auto') 
# plot(forecast(check, h=365.25))
#check<- features(q_ts,Qobs, feature_set(pkgs = "feasts")) #Get all features from feasts
#check <- tsoutliers(ts(log(q_ts$Qobs+0.01), 365.25)) #look for outliers

