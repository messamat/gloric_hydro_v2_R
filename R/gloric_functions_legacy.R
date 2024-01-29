############################ UTILITY FUNCTIONS ##################################
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
#'   \item Qt for which Original != Calculated discharge in GRDC record
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
  in_gaugetab %>%
    .[Original %in% c(-999, -99, -9999, 999, 9999), Original := NA] %>%
    .[Calculated %in% c(-999, -99, -9999, 999, 9999), Calculated := NA] %>%
    .[, `:=`(jday = format(as.Date(dates), '%j'),#Julian day
             q_rleid = rleid(Original),#Identify each group of consecutive values
             flag_mathis = 0)] #Create flag field)
  
  #Flag negative values
  in_gaugetab[Original < 0, flag_mathis := flag_mathis + 1]
  
  #Flag when more than 10 identical values in a row or when a single zero-flow
  in_gaugetab[, flag_mathis := flag_mathis +
                ((Original > 0) & (.N > 10)) +
                ((Original == 0) & (.N == 1)),
              by=q_rleid]
  
  #Flag |log(Q + 0.01) - mean| > 6SD for julian day mean and SD of 5d mean of log(Q + 0.01)
  in_gaugetab[, logmean5d := frollapply(log(Original + 0.01), n = 5, align='center',
                                        FUN=mean, na.rm = T)] %>% #Compute 5-day mean of log(Q+0.01)
    .[, `:=`(jdaymean = mean(logmean5d, na.rm = T),
             jdaysd = sd(logmean5d, na.rm = T)),
      by = jday] %>% #Compute mean and SD of 5-day mean of log(Q + 0.01) by Julian day
    .[abs(log(Original + 0.01) - jdaymean) > (6 * jdaysd),
      flag_mathis := flag_mathis + 1]
  
  #Flag values where Original != Calculated discharge
  in_gaugetab[Original != Calculated, flag_mathis := flag_mathis + 1]
  
  return(in_gaugetab)
}

#------ remove_constant -----------------------
#Function to remove spurious constant values
remove_constant <- function(gage_data, gap_n=20) {
  gage_data[1,'Flag2'] <- 0
  delist <- list()
  for (i in 2:nrow(gage_data)){
    #Accumulate number of days with constant values
    if (gage_data[i,'Flow'] == gage_data[i-1,'Flow']){
      gage_data[i,'Flag2'] <- gage_data[i-1,'Flag2']+1 
    } else {
      gage_data[i,'Flag2'] <- 0
    }
    #Flag period with at least gap_n days of constant values
    if (gage_data[i,'Flag2'] == gap_n) {
      delist <- c(delist, gage_data[(i-19):i,'Date'])
    } 
    if (gage_data[i,'Flag2'] > gap_n) {
      delist <- c(delist, gage_data[i,'Date'])
    }
  }
  return(delist)
}

#------ na.lomf ----------------------------------------------------------------
#Function to find the index, for each row, of the previous row with a non-NA value
#Takes in a univariate time series with NAs
na.lomf <- function(x) {
  if (length(x) > 0) {
    non.na.idx <- which(!is.na(x))
    if(is.na(x[1]))
      non.na.idx=c(1,non.na.idx)
    rep.int(non.na.idx, diff(c(non.na.idx, length(x) + 1))) #Repeat index of previous row with non-NA as many times gap until next non-NA values
  }
}

#------ na.lomb ----------------------------------------------------------------
#Function to find the index, for each row, of the next row with a non-NA value
#Takes in a univariate time series with NAs
na.lomb <- function(x) {
  if (length(x) > 0) {
    non.na.idx <- which(!is.na(x))
    if(is.na(x[length(x)]))
      non.na.idx=c(non.na.idx,length(x))
    y<-rep(NA,length(x))
    y[non.na.idx] <- non.na.idx
    zoo::na.locf(y, fromLast = TRUE)
  }
}
#------ get_legend -------------------------------------------------------------
#Function to extract legend from graph as a grob to be re-inserted later
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
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
#' \code{GRDC_NO, dates, Original, flag_mathis, missingdays} \cr
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
      ('Original' %in% names(GRDCgaugestats_record))) {
    gaugetab <- GRDCgaugestats_record
  } else {
    gaugetab <- readformatGRDC(GRDCgaugestats_record$path) %>%
      flagGRDCoutliers %>%
      .[, dates := as.Date(dates)] %>%
      .[!is.na(Original), missingdays := diny(year)-.N, by= 'year']
  }
  
  #Plot time series
  qtiles <- union(gaugetab[, min(Original, na.rm=T)],
                  gaugetab[, quantile(Original, probs=seq(0, 1, 0.1), na.rm=T)])
  
  subgaugetab <- gaugetab[missingdays < maxgap,]
  
  rawplot <- ggplot(subgaugetab[missingdays < maxgap,],
                    aes(x=dates, y=Original)) +
    geom_line(color='#045a8d', size=1, alpha=1/5) +
    geom_point(data = subgaugetab[flag_mathis == 0 & Original > 0,],
               color='#045a8d', size=1, alpha=1/3) +
    geom_point(data = subgaugetab[flag_mathis > 0 & Original > 0,],
               color='green') +
    geom_point(data = subgaugetab[flag_mathis == 0 & Original == 0,],
               color='red') +
    geom_point(data = subgaugetab[flag_mathis > 0 & Original == 0,],
               color='black') +
    scale_y_sqrt(breaks=qtiles, labels=qtiles) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(y='Discharge (m3/s)',
         title=paste0('GRDC: ', GRDCgaugestats_record$GRDC_NO)) +
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

#------ checkGRDCzeroes --------
#' Check zero-flow values in GRDC discharge record
#'
#' Plot and subset a given period (in days) on each side of every zero-flow
#' period in record for a given GRDC station. The output is a plot with a panel
#' showing discharge for each zero-flow period.
#'
#' @param GRDCstatsdt formatted data.table including a "path" column to access
#' GRDC standard daily discharge record text file for the gauging station of interest.
#' @param in_GRDC_NO GRDC unique identifier of the station to be investigated.
#' @param period (integer) number of days on each side of zero-flow values to subset.
#' @param yearthresh (integer) minimum year from which to analyze discharge record
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param in_scales (character) should panel scales be fixed ("fixed", the default),
#' free ("free"), or free in one dimension ("free_x", "free_y")?
#' @param labelvals (logical) whether to label individual discharge values
#'
#' @return \link[data.table]{data.table} subset of GRDCstatsdt, including only
#' records for zero-flow days +- period
#'
#' @export
checkGRDCzeroes <- function(GRDCstatsdt, in_GRDC_NO, period=15, yearthresh,
                            maxgap, in_scales='free_x', labelvals) {
  check <- readformatGRDC(GRDCstatsdt[GRDC_NO ==  in_GRDC_NO, path]) %>%
    flagGRDCoutliers %>%
    .[, dates := as.Date(dates)] %>%
    .[!is.na(Original), missingdays := diny(year)-.N, by= 'year'] %>%
    .[missingdays < maxgap,] %>%
    merge(check[, list(dates=seq(min(dates), max(dates), by='day'))],
          by='dates', all.x=T, all.y=T)
  
  checkdates <- unique(do.call("c",
                               lapply(
                                 as.Date(check[Original==0, dates]),
                                 function(i) seq.Date(i-period, i+period, by='day'))))
  
  zeroes <- check[dates %in% checkdates,]
  zeroes[, zerogrp := rleid(round(difftime(dates, lag(dates))))] %>%
    .[, grpN := .N, by=zerogrp]
  
  p <- plotGRDCtimeseries(zeroes[grpN > 1 & year > yearthresh,],
                          outpath=NULL, maxgap=maxgap, showmissing=T) +
    scale_y_sqrt() +
    scale_x_date(breaks='1 month') +
    facet_wrap(~zerogrp, scales=in_scales)
  
  if (labelvals == T) {
    p <- p + geom_text(data=zeroes[grpN > 1 & Original >0,],
                       aes(label=Original), vjust=1, alpha=1/2)
  }
  
  print(p)
  return(zeroes)
}

#------ plot_winterir -------
#' Plot winter-only intermittent rivers
#'
#' Utility function: plot discharge time series for GRDC or GSIM gauging stations
#' that have been deemed to monitor winter-only intermittent rivers.
#'
#' @param dt formatted data.table (generated in \code{\link{analyzemerge_gaugeir}}).
#' @param dbname (character) 'GRDC' or 'GSIM". Will determine whether to use \code{\link{plotGRDCtimeseries}}
#' or \code{\link{plotGSIMtimeseries}}
#' @param inp_resdir (character) path to directory where to save plots. A directory will be
#' created called " 'winterir_rawplots_o(yearthresh)_(date)'
#' @param yearthresh (integer) minimum year from which to plot discharge record
#' @param plotseries (logical) whether to generate plots
#'
#' @details “Winter-only” non-perennial gauging stations were defined as those
#' whose stream record contained less than one zero-flow day per year on average
#' during months with long-term mean air temperature over 10°C (averaged across
#' the local catchment immediately draining to the river reach, according to
#' WorldClim 2). In other words, “winter-only” non-perennial gauging stations
#' were those which would not have qualified as non-perennial according to our
#' criterion if only non-winter months were taken into account.
#'
#' @return subset of input dt only including winter-only intermittent irvers
#'
#' @export
plot_winterir <- function(dt, dbname, inp_resdir, yearthresh, plotseries = TRUE) {
  #Get data subset
  wintergauges <- dt[get(paste0('winteronlyir_o', yearthresh)) == 1 &
                       get(paste0('totalYears_kept_o', yearthresh)) >= 10,]
  
  #Create output directory
  resdir_winterirplots <- file.path(
    inp_resdir,
    paste0(dbname, 'winterir_rawplots_o', yearthresh, '_',
           format(Sys.Date(), '%Y%m%d')
    )
  )
  
  if (!(dir.exists(resdir_winterirplots))) {
    print(paste0('Creating ', resdir_winterirplots ))
    dir.create(resdir_winterirplots )
  }
  
  #Generate plots depending on database
  if (plotseries) {
    if (str_to_upper(dbname) == 'GRDC') {
      lapply(wintergauges$GRDC_NO, function(gauge) {
        print(gauge)
        plotGRDCtimeseries(GRDCgaugestats_record =  wintergauges[GRDC_NO == gauge,],
                           outpath = file.path(resdir_winterirplots,
                                               paste0('GRDC', gauge, '.png')))
      })
    } else if (str_to_upper(dbname) == 'GSIM') {
      lapply(wintergauges$gsim_no, function(gauge) {
        print(gauge)
        plotGSIMtimeseries(GSIMgaugestats_record = wintergauges[gsim_no == gauge,],
                           outpath = file.path(resdir_winterirplots,
                                               paste0('GSIM', gauge, '.png')))
      })
    }
  }
  
  return(wintergauges)
}


#------ plot_coastalir -------
#' Plot coastal intermittent rivers
#'
#' Utility function: plot discharge time series for GRDC or GSIM gauging stations
#' that have been deemed to be near marine waters.
#'
#' @param in_gaugep \link[sf]{sf} object of gauges with column called class99_19_rsp9_buf3k1
#' indicating whether the gauge is located within 3 km of marine waters.
#' @param dt formatted data.table (generated in \code{\link{analyzemerge_gaugeir}}).
#' @param yearthresh (integer) minimum year from which to plot discharge record
#' @param dbname (character) 'GRDC' or 'GSIM". Will determine whether to use \code{\link{plotGRDCtimeseries}}
#' or \code{\link{plotGSIMtimeseries}}
#' @param inp_resdir (character) path to directory where to save plots. A directory will be
#' created called " 'winterir_rawplots_o(yearthresh)_(date)'
#' @param plotseries (logical) whether to generate plots
#'
#' @details so-called marine stations were defined as those within 3 km of a coastline.
#'
#' @return subset of input dt only including only intermittent rivers within 3 km of a coastline
#'
#' @export
plot_coastalir <- function(in_gaugep = in_gaugep, dt = GRDCstatsdt, yearthresh,
                           dbname = 'grdc', inp_resdir = inp_resdir,
                           plotseries = TRUE) {
  #Get column name of uniquer identifiers
  idno <- fifelse(str_to_upper(dbname) == 'GRDC', 'GRDC_NO', 'gsim_no')
  
  #Select gauges within 3 km of a coastline of the given database (GRDC or GSIM)
  #and with at least 10 valid years of data
  coastalgauges <- in_gaugep[!is.na(in_gaugep$class99_19_rsp9_buf3k1) &
                               !is.na(as.data.frame(in_gaugep)[,idno]),] %>%setDT
  
  
  coastalir <- dt[get(paste0('totalYears_kept_o', yearthresh)) >= 10 &
                    get(idno) %in% coastalgauges[,get(idno)],]
  
  #Create output directory for plots
  resdir_coastalirplots <- file.path(inp_resdir,
                                     paste0(str_to_upper(dbname),
                                            'coastalir_rawplots_o',
                                            yearthresh, '_',
                                            format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_coastalirplots))) {
    print(paste0('Creating ', resdir_coastalirplots ))
    dir.create(resdir_coastalirplots )
  }
  
  #Generate plots depending on database
  if (plotseries) {
    if (str_to_upper(dbname) == 'GRDC') {
      lapply(coastalir$GRDC_NO, function(gauge) {
        print(gauge)
        plotGRDCtimeseries(GRDCgaugestats_record = coastalir[GRDC_NO == gauge,],
                           outpath = file.path(resdir_coastalirplots, paste0('GRDC', gauge, '.png')))
      })
    } else if (str_to_upper(dbname) == 'GSIM') {
      lapply(coastalir$gsim_no, function(gauge) {
        print(gauge)
        plotGSIMtimeseries(GSIMgaugestats_record = coastalir[gsim_no == gauge,],
                           outpath = file.path(resdir_coastalirplots, paste0('GSIM', gauge, '.png')))
      })
    }
  }
  
  
  return(coastalir)
}

#------ comp_ymean -------------
#' Compute annual mean
#'
#' Compute annual mean of a variable based on columns of monthly means.
#'
#' @param in_dt data table, contains columns of monthly values.
#' @param fieldex character, example column name containing monthly value.
#' It must include the month number in %m format (two digits e.g. 02 for February)
#' @param mstart integer; position of the first digit of the month in the column
#' of monthly values.
#' @param outcol character; name of new column to be created containing monthly mean
#'
#' @return Modifies in_dt in place. Add column of weighted mean across months
#' (based on number of days in months)
#'
#' @export
comp_ymean<- function(in_dt, fieldex, mstart, outcol) {
  refd <- data.table(month=c('01', '02', '03', '04', '05', '06',
                             '07', '08', '09', '10', '11', '12'),
                     dimonth = c(31, 28.25, 31, 30, 31, 30,
                                 31, 31, 30, 31, 30, 31))
  refd[, mcols := paste0(substr(fieldex, 1, mstart-1),
                         month)]
  in_dt2 <- copy(in_dt[, refd$mcols, with=F])
  in_dt2[, comp_ymeanID := .I]
  in_dt[, outcol] <- melt(in_dt2,
                          id.vars='comp_ymeanID',
                          variable.name = 'mcols') %>%
    .[refd, on='mcols'] %>%
    .[, weighted.mean(value, dimonth), by=comp_ymeanID] %>%
    .[,'V1',with=F]
}

#------ comp_irstats -----------
#' Compute flow intermittence statistics
#'
#' Format and compute a set of summary statistics based on a yearly streamflow 
#' gauging station time series. Used in \code{\link{comp_GRDCdurfreq}}. 
#' 
#'
#' @param tabyearly data.table 
#' @param maxgap (integer) maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics.
#' @param mdurthresh (numeric) threshold of mean annual number of zero-flow days beyond
#' which to classify gauge as intermittent.
#' @param yearthresh (integer) minimum year from which to analyze discharge record.
#' @param windowsize (integer) window size to check for zero-flow days. 
#' @param fullwindow (logical) whether years for which the window is truncated 
#' (e.g., beginning and end of time series) are taken in account in moving window analysis.
#'
#' 
#' @return link[data.table]{data.table} with 1 row and 10 columns:
#' \itemize{
#'   \item firstYear_kept: the first year in the time series with a number of missing daily values <= \code{maxgap} 
#'   (i.e., the first year in the subset of the time series kept for further analysis)
#'   \item lastYear_kept: the last year in the time series with a number of missing daily values <= \code{maxgap} 
#'   (i.e., the first year in the subset of the time series kept for further analysis)
#'   \item totalYears_kept: total number of years in the time series with a number of missing daily values <= \code{maxgap} 
#'   (i.e., the number of years in the subset of the time series kept for further analysis)
#'   \item totaldays: total number of days with data in years with a number of missing daily values <= \code{maxgap} 
#'   (i.e., total number of days of data kept for further analysis)
#'   \item integerperc: proportion of daily values which are in integer format. 
#'   Only including years with <= maxgap missing daily discharge values.
#'   \item sumDur: total number of daily zero-flow values. 
#'   Only including years with <= maxgap missing daily discharge values.
#'   \item mDur: mean annual number of daily zero-flow values. 
#'   Only including years with <= maxgap missing daily discharge values.
#'   \item mFreq: mean annual frequency of zero flow events. A zero-flow event is defined 
#'   as one or more consecutive days with a recorded daily discharge of zero. 
#'   Only including years with <= maxgap missing daily discharge values.
#'    \item intermittent: binary flow intermittence class. 1: non-perennial (if mDur >= 1, i.e., 
#'    if gauging station recorded zero-flow for at least one day per year on average); 0: perennial.
#'    \item movinginter: whether there is at least one zero-flow day in every \code{windowsize}-year (e.g., 20-year) moving window across the record.
#' } 
#' 
#'
#' @export
comp_irstats <- function(tabyearly, maxgap, mdurthresh, yearthresh,
                         windowsize, fullwindow) {
  
  checkpos <- function(x) {any(x>0)}
  
  #
  if ((windowsize %% 2)==0) {
    windowsize <- windowsize + 1
  }
  
  if (!fullwindow) {
    movinginter <- all(
      tabyearly[missingdays <= maxgap & year >= yearthresh,
                (frollapply(dur, n=windowsize, FUN=checkpos, align="center") >= mdurthresh)],
      na.rm=T)
  } else {
    movinginter <- all(
      tabyearly[missingdays <= maxgap & year >= yearthresh,
                c((frollapply(dur, n=round(windowsize/2), FUN=checkpos, align="left") >= mdurthresh),
                  (frollapply(dur, n=round(windowsize/2), FUN=checkpos, align="right") >= mdurthresh))],
      na.rm=T)
  }
  
  
  irstats <- tabyearly[missingdays <= maxgap & year >= yearthresh,
                       .(firstYear_kept=min(year),
                         lastYear_kept=max(year),
                         totalYears_kept=length(unique(year)),
                         totaldays = sum(datadays),
                         integerperc = sum(integerperc)/(sum(datadays)+sum(missingdays)),
                         sumDur = sum(dur),
                         mDur = mean(dur),
                         mFreq = mean(freq),
                         intermittent =
                           factor(fifelse(mean(dur)>=mdurthresh, 1, 0),
                                  levels=c('0','1')),
                         movinginter = movinginter
                       )]
  setnames(irstats, new = paste0(names(irstats), '_o', yearthresh))
  return(irstats)
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
              select=keptcols, verbose=TRUE)
  return(dt)
}

#------ formatscales ------------
#' Format plot scales
#'
#' Utility function to format plot scales for density distribution plots of environmental variables.
#' 
#' @param in_df data.frame with all records for environmental variables. 
#' Used to determine the appropriate range of values for each variable. In this case, 
#' a data.table of the river network hydro-environmental attributes.
#' @param varstopplot vector of variable names that will be plots and for which 
#' to return a list of scales
#'  
#' @return list of x and y scale objects + cartesian coordinates for ggplot
#' 
#' 
#' @export
formatscales <- function(in_df, varstoplot) {
  scales_x <- list(
    ari_ix_uav = scale_x_continuous(expand=c(0,0)),
    bio12_mm_uav  = scale_x_sqrt(expand=c(0,0),
                                 labels=c(0, 1000, 2000, 5000, 10000)),
    bio14_mm_uav  = scale_x_sqrt(expand=c(0,0),
                                 breaks = c(0, 50, 100, 200, 500),
                                 labels=c(0, 50, 100, 200, 500)),
    cly_pc_uav = scale_x_continuous(labels=percent_format(scale=1), expand=c(0,0)),
    cmi_ix_uyr = scale_x_continuous(),
    dis_m3_pyr = scale_x_log10(breaks=c(1, 10^2,
                                        10^(0:log10(max(in_df$dis_m3_pyr)))),
                               labels=c(0, 10^2,
                                        10^(0:log10(max(in_df$dis_m3_pyr)))),
                               expand=c(0,0)),
    dor_pc_pva = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    for_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    gla_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    kar_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    lka_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    pet_mm_uyr = scale_x_continuous(expand=c(0,0)),
    sdis_ms_uyr = scale_x_continuous(expand=c(0,0)),
    snw_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    run_mm_cyr = scale_x_continuous(expand=c(0,0)),
    swc_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    tmp_dc_uyr = scale_x_continuous(expand=c(0,0)),
    hdi_ix_cav = scale_x_continuous(expand=c(0,0)),
    hft_ix_c93 = scale_x_continuous(expand=c(0,0)),
    ORD_STRA = scale_x_continuous(expand=c(0,0)),
    UPLAND_SKM = scale_x_log10(breaks=c(1, 10^2,
                                        10^(0:log10(max(in_df$UPLAND_SKM)))),
                               labels=c(1, 10^2,
                                        10^(0:log10(max(in_df$UPLAND_SKM)))),
                               expand=c(0,0)),
    gwt_m_cav = scale_x_sqrt(expand=c(0,0)),
    ire_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0))
  ) %>%
    .[(names(.) %in% names(in_df)) & names(.) %in% varstoplot]
  #Only keep those variables that are actually in df and that we want to plot
  
  scales_y <- unlist(rep(list(scale_y_continuous(expand=c(0,0))),
                         labels = scientific_format(),
                         length(scales_x)),
                     recursive=F) %>%
    setNames(names(scales_x))
  
  scales_y[['dis_m3_pmn']] <- scale_y_sqrt(expand=c(0,0))
  scales_y[['glc_pc_u16']] <- scale_y_continuous(trans='log1p',
                                                 breaks=c(10, 1000, 100000, 10000000))
  
  coordcart <- lapply(varstoplot, function(var) {
    coord_cartesian(xlim=as.data.table(in_df)[, c(min(get(var), na.rm=T),
                                                  max(get(var), na.rm=T))])
  }) %>%
    setNames(varstoplot)
  
  coordcart[['clz_cl_cmj']] <-  coord_cartesian(
    xlim=c(1,max(in_df$clz_cl_cmj)))
  coordcart[['kar_pc_use']] <-  coord_cartesian(
    xlim=c(0, 100))
  coordcart[['pet_mm_uyr']] <-  coord_cartesian(
    xlim=c(0, max(in_df$pet_mm_uyr)))
  coordcart[['ORD_STRA']] <-  coord_cartesian(
    xlim=c(1, 10))
  coordcart[['ari_ix_uav']] <-  coord_cartesian(
    xlim=c(0, 100))
  
  return(list(scales_x=scales_x, scales_y=scales_y, coordcart=coordcart))
}

#------ ggenvhist -------------
#' Plot of environmental histogram
#'
#' Utility function to create an individual density plot of the distribution of a 
#' given environmental variables across gauges and the whole global river network.
#' 
#' @param vartoplot (column) variable for which to produce a density plot.
#' @param in_gaugedt data.table of gauging stations' environmental attributes.
#' @param in_rivdt data.table of global river network's environmental attributes 
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' @param scalesenvhist list of scale objects to format plot. From \link{formatscales}.
#' 
#' @return ggplot with two density distributions of the environmental variable,
#' one for the gauging stations and one for the global river network. 
#' 
#' 
#' @export
ggenvhist <- function(vartoplot, in_gaugedt, in_rivdt, in_predvars,
                      scalesenvhist) {
  print(vartoplot)
  
  varname <- in_predvars[varcode==vartoplot, Attribute]
  #paste0(Attribute, ' ',Keyscale,Keystat,' (',unit,')')]
  
  if (vartoplot == "clz_cl_cmj") {
    rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
                       by=as.factor(clz_cl_cmj)]
    gclz <- in_gaugedt[,.N/in_gaugedt[,.N],by=as.factor(clz_cl_cmj)]
    bindclz <- rbind(rivclz, gclz, idcol='source')%>%
      setnames(c( 'source', vartoplot, 'density'))
    
    penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
      geom_bar(aes(fill=as.factor(source)), stat='identity',
               position = 'dodge', alpha=1/2, width=.6) +
      scale_fill_manual(values=c('#2b8cbe', '#dd3497'))
    
  } else if (vartoplot == "glc_pc_u16") {
    rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
                       by=glc_pc_u16]
    gclz <- in_gaugedt[,.N/in_gaugedt[,.N],by=glc_pc_u16]
    bindclz <- rbind(rivclz, gclz, idcol='source')%>%
      setnames(c( 'source', vartoplot, 'density'))
    
    penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
      geom_bar(aes(fill=as.factor(source)), stat='identity',
               position = 'identity', alpha=1/2, width=.6) +
      scale_fill_manual(values=c('#2b8cbe', '#dd3497'))
    #
    #     penvhist <- ggplot(in_gaugedt, aes_string(x=vartoplot)) +
    #       geom_histogram(data=in_rivdt, aes(weight = LENGTH_KM),
    #                      fill='#2b8cbe', alpha=0.5, bins=101) +
    #       geom_histogram(fill='#dd3497', alpha=0.5, bins=101)
    
  } else {
    penvhist <- ggplot(in_gaugedt, aes_string(x=vartoplot)) +
      geom_density(data=in_rivdt, aes(weight = LENGTH_KM),
                   fill='#2b8cbe', alpha=0.5) +
      geom_density(fill='#dd3497', alpha=0.5) +
      ylab('Density')
  }
  
  penvhist <- penvhist +
    scalesenvhist$scales_x[[vartoplot]] +
    #scalesenvhist$scales_y[[vartoplot]] +
    scalesenvhist$coordcart[[vartoplot]] +
    xlab(varname) +
    theme_classic() +
    theme(strip.background=element_rect(colour="white", fill='lightgray'),
          legend.position = 'none',
          axis.title.y = element_blank(),
          axis.title = element_text(size=12))
  
  # if (which(vartoplot %in% varstoplot_hist)!=length(varstoplot_hist)) {
  #   penvhist <- penvhist +
  #     theme(legend.position='none')
  # }
  
  return(ggplotGrob(penvhist))
}

#------ layout_ggenvhist --------------------------
#' Layout plots of environmental histograms
#'
#' Run plotting functions across predictor variables and arrange plots
#' 
#' @param in_rivernetwork data.table of global river network's environmental attributes. Here, output from \link{rformat_network}.
#' @param in_gaugepred selected gauging stations' environmental attributes. Here, output from \link{write_gaugepreds}.
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' 
#' @details function used to produce Extended Data Fig. 8  c-p in Messager et al. 2021.
#' 
#' @return ggplots with two density distributions of the environmental variable,
#' one for the gauging stations and one for the global river network. 
#' 
#' 
#' @export
layout_ggenvhist <- function(in_rivernetwork, in_gaugepred, in_predvars) {
  varstoplot_hist <- c(
    "bio1_dc_uav", "bio7_dc_uav", "bio12_mm_uav", "bio14_mm_uav", "clz_cl_cmj",
    "ari_ix_uav", "dis_m3_pyr", "sdis_ms_uyr", "gwt_m_cav", "UPLAND_SKM",
    "lka_pc_use", "snw_pc_uyr", "kar_pc_use", "for_pc_use") #, "glc_pc_u16")
  
  if ("dis_m3_pyr" %in% varstoplot_hist) {
    setDT(in_rivernetwork)[, dis_m3_pyr := dis_m3_pyr + 1]
    setDT(in_gaugepred)[, dis_m3_pyr := dis_m3_pyr + 1]
  }
  
  #Get legend
  pleg <- ggplot(in_gaugepred, aes(x=dis_m3_pyr, fill=factor(IRpredcat_full))) +
    geom_density(alpha=1/2) +
    scale_fill_manual(values=c('#2b8cbe', '#dd3497'),
                      name = 'Dataset',
                      labels=c('Global river network',
                               'Training gauges')) +
    theme(text=element_text(size=14))
  
  tmp <- ggplot_gtable(ggplot_build(pleg))
  leg <- tmp$grobs[[
    which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
  ]]
  
  #Get scales
  scalesenvhist <- formatscales(in_df=in_rivernetwork, varstoplot=varstoplot_hist)
  
  #Plot each facet
  penvhist_grobs <- lapply(varstoplot_hist, ggenvhist,
                           in_gaugedt = in_gaugepred,
                           in_rivdt = in_rivernetwork,
                           in_predvars = in_predvars,
                           scalesenvhist = scalesenvhist)
  #Add legend
  penvhist_grobs[[length(penvhist_grobs) + 1]] <- leg
  
  #Plot
  grid.newpage()
  do.call("grid.arrange", list(grobs=penvhist_grobs, nrow=5))
}

#


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
  
  Reduce(
    function(...) {
      merge(..., by = by, all = all, sort = sort)
    }, dt_list)
}

#------ prettydend --------------------------------------
#Make a nice looking dendogram based on a clustering output
prettydend <- function(hclus_out, colorder=NULL, colors=NULL, labels=NULL,
                       kclass=7, classnames = NULL) {
  classr <- dendextend::cutree(hclus_out, k=kclass, 
                               order_clusters_as_data = FALSE)
  classr_df <- data.frame(ID=names(classr), gclass=classr) 
  
  if (!is.null(classnames)) {
    classr_df <- merge(classr_df, classnames, by='gclass')
    grouplabels <- classnames$classnames
  } else {
    grouplabels <- TRUE
  }
  
  hclus_out_name <- hclus_out
  if (!is.null(labels)) {
    hclus_out_name$labels <- labels
  }
  
  #,paste(RGS_No,"-",RGS_Loc," R. at ", RGS_Name,sep=""))
  
  dendname <- as.dendrogram(hclus_out_name)
  
  if (is.null(colorder)) colorder = 1:kclass
  
  par(cex=0.7, mar=c(2.5, 1, 0, 9)) #bottom left top right
  dendname %>% set("branches_lwd", 3) %>% 
    color_branches(k=kclass, col=colors[colorder]) %>% 
    #color_branches(clusters=as.numeric(temp_col), col=levels(temp_col), groupLabels=as.character(as.numeric(temp_col))) %>% 
    color_labels(k=kclass, col=colors[colorder]) %>%
    plot(horiz=TRUE,xlab="Gower's distance", ylab="",mgp=c(1.5,0.5,0))
  title(ylab="Department", line=0, cex=1)
  
  return(list(classr_df, dendname))
}


#------ bin_dt -----------------
#' Bin data.table
#'
#' Bins a data.table over a numeric column.
#'
#' @param in_dt \link[data.table]{data.table} to bin.
#' @param binvar (character) column that will be used to define bins.
#' @param binfunc (character) binning approach. One of 'manual', 'equal_length', 'equal_freq'.
#' @param binarg (numeric) binning argument, depends on binning approach (\code{binfunc}).
#' @param bintrans (character or numeric) transformation of \code{binvar}, default is NULL.
#' @param ndigits (integer) number of decimals to keep for displaying formatted bin limits
#' @param na.rm (logical) whether to include NAs.
#' @param valuevar (character) if na.rm = FALSE, column to use to detect NAs and remove records.
#'
#' @return input \link[data.table]{data.table} with four new columns:
#' \itemize{
#'   \item bin - bin number (1 being the lowest value bin)
#'   \item bin_lmin - bin lower limit
#'   \item bin_lmax - bin higher limit
#'   \item bin_lformat - formatted character of bin limits
#'   (format: \code{round(bin_lmin, ndigits) - round(bin_lmax, ndigits))
#' }
#'
#' @details inspired from [rbin package](https://github.com/rsquaredacademy/rbin).
#' Differences include that it concentrates all binning approaches within a single
#' function and works on a data.table.
#'
#' binfunc: \cr
#' \itemize{
#'   \item 'manual' - bin continuous data manually. \code{binarg} sets the inner bin limits,
#'  such that the final table will have \code{length(binarg) + 1} bins. The lower end of the
#'  first bin is automatically set to be the minimum value in \code{binvar} and the upper end of
#'  the last bin is set to be the maximum value in \code{binvar}
#'
#'   \item 'equal_length' - Bin continuous data such that each bin has the same \code{binvar} interval length.
#'   If \code{bintrans} is not \code{NULL}, then interval length is computed on transformed scale.
#'   \code{binarg} (automatically rounded to the nearest integer) sets the number of bins.
#'
#'   \item 'equal_freq' - Bin continuous data such that each bin has the same number of records.
#'   \code{binarg} (automatically rounded to the nearest integer) sets the number of bins.
#' }
#'
#' bintrans: can either be 'log' (for natural log) or a numeric exponent to transform
#' according to x^bintrans.
#'
#' @seealso for examples, see applications in \code{\link{bin_misclass}},
#' \code{\link{eval_watergap}}, \code{\link{tabulate_globalsummary}},
#' \code{\link{formathistab}}, \code{\link{compare_us}}
#'
#' @export
bin_dt <- function(in_dt, binvar, binfunc, binarg,
                   bintrans=NULL, ndigits=2,
                   na.rm=FALSE, valuevar=NULL) {
  #Inspired from rbin, adapted to dt and simplified
  in_dt <- copy(in_dt)
  
  el_freq <- function(byd, bins) {
    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    append(min(byd, na.rm = TRUE), min(byd, na.rm = TRUE) + (bin_length * seq_len(bins)))[1:bins]
    
  }
  
  eu_freq <- function(byd, bins) {
    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    ufreq      <- min(byd, na.rm = TRUE) + (bin_length * seq_len(bins))
    n          <- length(ufreq)
    ufreq[n]   <- max(byd, na.rm = TRUE) + 1
    return(ufreq)
    
  }
  
  binvar_orig <- copy(binvar)
  
  #Remove NAs
  if (na.rm) {
    in_dt <- in_dt[!is.na(get(eval(valuevar))),]
  }
  
  #Transform data if trans
  if (!is.null(bintrans)) {
    transvar <- paste0(binvar, '_bintrans')
    if (bintrans == 'log') {
      nneg <- in_dt[get(eval(binvar)) <= 0, .N]
      warning(paste0('There are ', nneg, ' records with', binvar, ' <= 0...',
                     'removing them for log transformation'))
      in_dt[, eval(transvar) :=
              log(get(eval(binvar)))]
      
    } else if (is.numeric(bintrans)) {
      in_dt[, eval(transvar) :=
              get(eval(binvar))^eval(bintrans)]
      
    }
    binvar = transvar
  }
  
  byd <- in_dt[, get(eval(binvar))]
  
  if (binfunc == 'manual') {
    l_freq    <- append(min(byd), binarg)
    u_freq    <- c(binarg, (max(byd, na.rm = TRUE) + 1))
    bins      <- length(binarg) + 1
  }
  
  if (binfunc == 'equal_length') {
    bins = round(binarg)
    l_freq    <- el_freq(byd, bins)
    u_freq    <- eu_freq(byd, bins)
  }
  
  if (binfunc == 'equal_freq') {
    bins = round(binarg)
    bin_prop     <- 1 / bins
    bin_length   <- in_dt[, round(.N/bins)]
    first_bins   <- (bins - 1) * bin_length
    residual     <- in_dt[, .N - first_bins]
    bin_rep      <- c(rep(seq_len((bins - 1)), each = bin_length),
                      rep(residual, residual))
    l_freq        <- c(1, (bin_length * seq_len((bins - 1)) + 1))
    u_freq       <- c(bin_length * seq_len((bins - 1)), in_dt[,.N])
    setorderv(in_dt, cols= eval(binvar))
    in_dt[, binid := .I]
    binvar = 'binid'
  }
  
  for (i in seq_len(bins)) {
    in_dt[get(eval(binvar)) >= l_freq[i] & get(eval(binvar)) < u_freq[i],
          bin := i]
    in_dt[bin == i, `:=`(bin_lmin = min(get(eval(binvar_orig)), na.rm=T),
                         bin_lmax = max(get(eval(binvar_orig)), na.rm=T))] %>%
      .[bin == i, bin_lformat := paste(round(bin_lmin, ndigits),
                                       round(bin_lmax, ndigits),
                                       sep='-')]
    
    if (i == bins) {
      in_dt[get(eval(binvar)) == u_freq[i],  bin := i]
    }
  }
  
  if (binfunc == 'equal_freq') {in_dt[, binid := NULL]}
  
  return(in_dt)
}

#------ label_manualbins ------------
#' Label manual bins
#'
#' Utility function: label bin limits for \code{\link{formathistab}}.
#'
#' @param binarg (character vector) Arguments for bin_dt manual
#' @param minval (numeric) value to set for lower limit of first bin.
#'
#' @return vector of labels
#'
#' @export
label_manualbins <- function(binarg, minval) {
  minlabel <- paste(minval, binarg[1], sep=" - ")
  otherlabels <- mapply(function(x, y) {paste(x, y-1, sep=" - ")},
                        binarg[1:(length(binarg)-1)], binarg[2:length(binarg)])
  return(c(minlabel, otherlabels))
}


#------ formathistab -----------------
#' Format histogram table
#'
#' Creates a binary frequency histogram table by computing the proportion of records
#' for which a selected column has a given value
#' (in this study, the percentage length of rivers that are deemed intermittent)
#' after binning the table by a continuous variable (e.g. river discharge).
#'
#' @inheritParams bin_dt
#' @param castvar (character) olumn that will be used to define bins —
#' the equivalent of \code{binvar} in \code{bin_dt}.
#' @param valuevar (character) variable to summarize (e.g. intermittency class).
#' @param valuevarsub (character) value of \code{valuevar} to summarize (e.g. '1' for intermittent rivers)
#' @param weightvar (character) variable to weigh proportion by. (e.g. river reach length).
#' @param binlabels (character vector) formatted labels for bins.
#' @param datname (character) optional, adds a column to output table describing what the data describes (e.g. France)
#'
#' @return \link[data.table]{data.table} with five columns:
#' \itemize{
#'   \item bin - bin number
#'   \item perc - percentage of records (or of \code{weightvar}) that meet \code{valuevar == valuevarsub}.
#'   \item binsumlength - number of records (if default \code{weightvar}) or sum of \code{weightvar} by bin (e.g. total river length).
#'   \item dat - \code{datname} argument
#'   \item binformat - label for each bin provided through \code{binlabels} argument
#' }
#'
#' @seealso used in \code{\link{compare_fr}}, \code{\link{compare_us}},
#' and \code{\link{compare_au}}
#'
#' @export
formathistab <- function(in_dt, castvar, valuevar, valuevarsub,
                         weightvar=1, binfunc, binarg, binlabels,
                         datname=NULL) {
  rivbin <- bin_dt(in_dt = as.data.table(in_dt),
                   binvar = castvar,
                   valuevar = valuevar,
                   binfunc = binfunc,
                   binarg = binarg)
  
  netstat <- as.data.table(rivbin)[,sum(get(weightvar)),
                                   by=c(eval(valuevar), "bin")]
  tidyperc_riv <- netstat[, list(perc = 100*.SD[get(valuevar)==valuevarsub,
                                                sum(V1)]/sum(V1),
                                 binsumlength = sum(V1)),
                          by=c('bin')] %>%
    setorder(bin) %>%
    .[, `:=`(dat=datname,
             binformat = binlabels[bin])]
  return(tidyperc_riv)
}

#------ scientific_10 ------------
#' Format number to nearest log10 bin in scientific format
#'
#' Format a number to nearest order of magnitude (log10 integer) and gives out a
#' scientific format expression.
#'
#' @param x (numeric)
#'
#' @return (expression)
#'
#' @example
#' scientific_10(45435) #expression(10^+04)
#' scientific_10(0.0000431) #expression(10^-05)
#'
#' @export
scientific_10 <- function(x) {
  parse(text=gsub(".*e", "10^", scales::scientific_format()(x)))
}

#------ allHITcomp -------------------------------------------------------------
allHITcomp <- function(dfhydro, dfenv, gageID, templateID='1KA9',hstats="all", floodquantile=0.95) {
####Get template
dailyQClean <- validate_data(dfhydro[dfhydro$ID==templateID,c("Date", "Flow")], yearType="water")
#Calculate all hit stats
HITall_template <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=10, pref="mean",
                               drainArea=dfenv[dfenv$RGS_No==templateID,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
colnames(HITall_template)[2] <- templateID
HITall <- data.frame(indice=HITall_template$indice) 
####Compute metrics for all gages
for (gage in unique(dfhydro[,gageID])) {
  print(gage)
  try({
    #Check data for completeness
    dailyQClean <- validate_data(dfhydro[dfhydro$ID==gage,c("Date", "Flow")], yearType="water")
    #Calculate all hit stats
    calc_allHITout <- calc_allHIT(dailyQClean, yearType="water", stats=hstats, digits=10, pref="mean",
                                  drainArea=dfenv[dfenv$RGS_No==gage,'WsArea'], floodThreshold = quantile(dailyQClean$discharge, floodquantile))
    colnames(calc_allHITout)[2] <- gage
    HITall <- merge(HITall, calc_allHITout, by='indice')
  })
}
HITall_formatmelt <-melt(setDT(HITall), id.vars = "indice",variable.name = gageID) 
HITall_formatmelt[is.infinite(HITall_formatmelt$value),'value'] <- NA
HITall_formatmelt[is.nan(HITall_formatmelt$value),'value'] <- NA
return(HITall_formatmelt)
}
#------ HITboxplot -------------------------------------------------------------
HITboxplot <- function(HITdf, plotname) {
  HITallbox<- HITdf
  HITallbox$group1 <- as.factor(substr(HITdf$indice,1,1)) #Subset metric name into main category
  HITallbox$group2 <- as.factor(substr(HITdf$indice,2,2)) #and l, h, and a
  HITallbox$indice_sub <- substr(HITdf$indice,3,5) #and metric number
  HITallbox$indice_sub <- factor(HITallbox$indice_sub, levels = unique(HITallbox$indice_sub[order(as.numeric(as.character(HITallbox$indice_sub)))]))
  
  #Format panel orders and names
  HITallbox$group1 <- factor(HITallbox$group1, levels = c('m','f','d','t','r'), 
                             labels = c("Magnitude-m", "Frequency-f", "Duration-d",'Timing-t',"Rate of change-r")) 
  HITallbox$group2 <- factor(HITallbox$group2, levels = c('h','a','l'), labels=c("High flow-h",'Average flow-a',"Low flow-l"))
  
  #Format x-labels
  metric_breaks <- unique(HITallbox$indice_sub[order(HITallbox$indice_sub)])
  metric_labels <- as.character(metric_breaks)
  metric_labels[seq(2, length(unique(HITallbox$indice_sub)), by=2)] <- ""  #Show only every other metric number
  
  HITallboxplot <-ggplot(HITallbox, aes(x=indice_sub, y=value, color=group1)) + 
    scale_y_log10(name='Metric value') +
    geom_boxplot(notch=F) +
    facet_grid(group2~group1, scales = "free", space="free_x") + 
    scale_x_discrete(name='Hydrologic metric', 
                     breaks= metric_breaks,
                     labels= metric_labels)+ 
    theme_classic() +
    theme(axis.title = element_text(size=20),
          axis.text.y = element_text(size=16),
          axis.text.x = element_text(size=16),
          strip.text = element_text(size = 15),
          legend.position='none')
  png(file.path(outdir,plotname),width=20, height=12,units='in',res=600)
  print(HITallboxplot)
  dev.off()
}

#------ HITdist ---------------------------------------------------------------
#Format hydrologic metrics to use in classification and compute Gower's distance
HITdist <- function(HITdf, logmetrics) { 
  if (logmetrics==TRUE){
    HITdf$Value <- log(HITdf$Value+1) #log-transform metric
  }
  HITdf_format <- dcast(HITdf, ID ~ indice)
  HITdf_format <- merge(HITdf_format, gagesenvrec[,c('RGS_No','WsArea')], by.x='ID', by.y='RGS_No')
  dimindices <- c('ma1','ma2',paste('ma',seq(12,23),sep=''),paste('ml',seq(1,12),sep=''),paste('mh',seq(1,12),sep=''), 
                  paste('dl',seq(1,5),sep=''),paste('dh',seq(1,5),sep=''),'ra1','ra3','ra6','ra7') #List of dimensional indices 
  dimindices <- dimindices[dimindices %in% colnames(HITdf_format)] #Make sure they are all in the dataset
  HITdf_format <- as.data.frame(setDT(HITdf_format)[,(dimindices) := lapply(.SD, function(x) round(x/WsArea, digits=10)), .SDcols=dimindices]) #Standardize dimensional indices by drainage area
  row.names(HITdf_format) <- HITdf_format$ID 
  HITdf_format <- HITdf_format[,-which(colnames(HITdf_format) %in% c('ID','WsArea'))] #Get rid of non-indices columns
  HITdf_stand <- data.stand(HITdf_format[,2:(ncol(HITdf_format))],method='standardize',margin='column',plot=F) #z-standardize columnwise 
  gauge_gow<- gowdis(HITdf_stand, w=rep(1,ncol(HITdf_stand)), asym.bin = NULL) #Compute Gower's distance so that missing values will not be taken in account
  return(gauge_gow)
}

#------ cluster_diagnostic -----------------------------------------------------
cluster_diagnostic <- function(clusterres, clusname, gowdis, 
                               ylabel="Gower's distance", format='pdf') {
  #hclus.table(clusterres)
  print(paste0('Agglomerative coefficient: ', coef.hclust(clusterres))) #Compute agglomerative coefficient
  print(paste0('Cophenetic correlation coefficient: ',cor(gowdis, cophenetic(clusterres)))) #Compute cophenetic coefficient
  #Plot cophenetic relationship 
  png(file.path(outdir, paste(clusname,'r6r_cophe','.png',sep="")), width=8, height=8, units='in',res=600)
  hclus.cophenetic(gowdis, clusterres) 
  dev.off()
  #Scree plot
  if (format=='png'){
    png(file.path(outdir, paste(clusname,'r6r_scree','.png',sep="")), width=8, height=8, units='in',res=600)
  } 
  if (format=='pdf'){
    pdf(file.path(outdir, paste(clusname,'r6r_scree','.pdf',sep="")), width=8, height=8)
  }
  hclus.scree(clusterres, ylabel=ylabel, frame.plot=FALSE,cex=1, xlim=c(0,length(clusterres$height)+5), ylim=c(0,max(clusterres$height)+0.1), xaxs='i', yaxs='i') 
  dev.off()
  #Plot dendogram
  png(file.path(outdir, paste(clusname,'r6r_dendogram','.png',sep="")), width=8, height=8, units='in',res=600)
  plot(clusterres, main=paste(clusname, "gauge cluster dendogram",sep=" "), xlab='Gauge ID', ylab="Gower's distance", hang=-1)   
  rect.hclust(clusterres, k=4) #Draw rectangle around k classes
  rect.hclust(clusterres, k=5) 
  rect.hclust(clusterres, k=6) 
  rect.hclust(clusterres, k=7) 
  rect.hclust(clusterres, k=8) 
  dev.off()
}

#------ hydrographplots --------------------------------------------------------
hydrographplots <- function(hydrodat, classtab, dir, kclass, classnames=NULL) {
  outdirclass <- file.path(outdir,dir)
  hydrodat_class_join <- merge(hydrodat, classtab, by="ID")
  
  if (!is.null(classnames)) {
    hydrodat_class_join <- merge(hydrodat_class_join, classnames, by='gclass') 
  } else {
    hydrodat_class_join[, classnames := gclass]
  }
  
  write.csv(hydrodat_class_join, file.path(outdirclass,'rufidat_class_join.csv'), row.names=F)
  setDT(hydrodat_class_join)[,yrmean:=mean(Flow),.(ID,hyear)] #Compute average daily flow for each station and year
  #Compute statistics on long-term daily flow (average, min, max, Q10, Q25, Q75, Q90) across all stations and years for each class
  classflowstats <- setDT(hydrodat_class_join)[,list(classmeanfull=mean(Flow, na.rm=T), classmean= mean(Flow/yrmean,na.rm=T),classQ75= quantile(Flow/yrmean, .25,na.rm=T),
                                                     classQ25=quantile(Flow/yrmean, .75,na.rm=T),classQ90=quantile(Flow/yrmean, .10,na.rm=T),
                                                     classQ10=quantile(Flow/yrmean, .90,na.rm=T),classmax=max(Flow/yrmean,na.rm=T),
                                                     classmin=min(Flow/yrmean,na.rm=T),classsd=sd(Flow/yrmean,na.rm=T), 
                                                     cal_hdoy=format(as.Date(hdoy, origin='2015-10-01'), "%Y-%m-%d")),
                                               .(gclass, classnames,hdoy)] 
  
  #Superimposed non-standardized average yearly hydrograph for each class
  classhydro_allfull <- ggplot(as.data.frame(classflowstats), aes(x=as.Date(cal_hdoy), y=classmeanfull, color=factor(classnames))) + 
    geom_line(size=1, alpha=0.8) + 
    scale_color_manual(name='Hydrologic class',values=classcol[1:kclass]) +
    scale_y_continuous(name=expression('Daily mean discharge'~(m^{3}%.%s^{-1})),expand=c(0,0),limits=c(0,NA)) + 
    scale_x_date(name='Date',date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) + 
    theme_classic() + 
    theme(legend.position='none',
          text=element_text(size=18)) +
    labs(subtitle = "(a)")
  
  #Superimposed standardized average yearly hydrograph for each class
  classhydro_all <- ggplot(as.data.frame(classflowstats), aes(x=as.Date(cal_hdoy), y=classmean, color=factor(classnames))) + 
    geom_line(size=1, alpha=0.8) + 
    scale_color_manual(name='Hydrologic class',values=classcol[1:kclass]) +
    scale_fill_manual(name='Hydrologic class',values=classcol[1:kclass]) +
    scale_y_continuous(name='Daily mean discharge/Mean daily discharge',expand=c(0,0),limits=c(0,NA)) + 
    scale_x_date(name='Date',date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) + 
    theme_classic() + 
    theme(legend.position='none',
          text=element_text(size=18))
  
  #Facetted standardized average yearly hydrograph for each class + Q90-Q10 ribbon
  classhydro_facet <-ggplot(as.data.frame(classflowstats), aes(x=as.Date(cal_hdoy))) + 
    #geom_ribbon(aes(ymin=ifelse(classmean-2*classsd>=0,classmean-2*classsd,0), ymax=classmean+2*classsd,
    #                fill=factor(classnames)),alpha=0.3) +
    geom_ribbon(aes(ymin=classQ90, ymax=classQ10,
                    fill=factor(classnames)),alpha=0.3) +
    geom_line(aes(y=classmean, color=factor(classnames)),size=1.2) + 
    facet_grid(classnames~.,scale='free_y') +
    scale_color_manual(name='Hydrologic class',values=classcol[1:kclass]) +
    scale_fill_manual(name='Hydrologic class',values=classcol[1:kclass]) +
    scale_y_continuous(name='Standardized daily mean discharge',expand=c(0,0),limits=c(0,NA)) + 
    scale_x_date(name='Date',date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) + 
    annotate("segment", x=as.Date('2015-10-01'), xend=as.Date('2016-09-30'), y=0, yend=0)+ 
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank(),
          text=element_text(size=18)) +
    labs(subtitle = "(b)")
  
  p1 <- ggplot_gtable(ggplot_build(classhydro_allfull))
  p2 <- ggplot_gtable(ggplot_build(classhydro_facet))
  lay= t(c(1,1,2,2))
  png(file.path(outdirclass,paste0(kclass,'class_hydrograph.png')),width = 16, height=9,units='in',res=600)
  print(grid.arrange(p1,p2, ncol=2, layout_matrix = lay))
  dev.off()
}


#------ format_HIT_table -------------------------------------------------------
#Format (transpose, truncate and order col and row names) and export table to HTML format with default formatting
#original name: tableformat
format_HIT_table <- function(df,tabname) { 
  df <- digitform(df, 2:(ncol(df)))
  df[is.na(df)] <- '–' 
  n <- df[,1]
  ycountdf <-predsmelt[!duplicated(predsmelt$ID),c('ID','ycount_full')]#Add number of years on record
  colnames(ycountdf) <- c('ID', 'Record Length (years)')
  df <- merge(ycountdf, df, by='ID') 
  df <- as.data.frame(t(df[,-1]))
  stations <- data.frame(basin=substr(n,3,3),ID=as.numeric(str_extract(as.character(substr(n,3,6)), "\\-*\\d+\\.*\\d*"))) #Only keep most basic ID subset
  colnames(df) <- with(stations, paste(basin,ID,sep="")) #Remove 1K and trailing letter
  metrics <-data.frame(type=substr(rownames(df),1,2), num=as.numeric(substr(rownames(df),3,5))) 
  HITo15y_format <- df[with(metrics, order(type, num)), #Order rows by metrics number
                       with(stations, order(basin,ID))] #Order columns first by basin then by numeric station ID
  #Re-place record length as first row
  rec <- HITo15y_format[rownames(HITo15y_format)=='Record Length (years)',]
  HITo15y_format <- HITo15y_format[-which(rownames(HITo15y_format)=='Record Length (years)'),]
  HITo15y_format <- rbind(rec, HITo15y_format)
  HITo15y_format[which(rownames(HITo15y_format) %in% HITo15ysub3$indice),'selected'] <- 'Yes'
  
  kable(HITo15y_format) %>%
    kable_styling(bootstrap_options = "striped", font_size = 9) %>%
    save_kable(tabname, self_contained=T)
}

#------ format_HIT_classable ---------------------------------------------------
classtableformat <- function(df, KWtab, tabname) {
  classHIT_stats<- setDT(df)[,`:=`(classmean=mean(value, na.rm=T),
                                   classsd=sd(value,na.rm=T)),
                             .(indice, classnames)] #Get mean and SD of hydrologic metric for each class
  classHIT_stats <- classHIT_stats[!duplicated(classHIT_stats[,c('indice','classnames')]),]
  classHIT_statsmean <- as.data.frame(dcast(classHIT_stats, classnames~indice, value.var='classmean'))
  classHIT_statssd <- as.data.frame(dcast(classHIT_stats, classnames~indice, value.var='classsd'))
  #Format digits for mean and sd
  classHIT_stats_meanform <- melt(setDT(digitform(classHIT_statsmean,
                                                  cols=2:(ncol(classHIT_statsmean)), extradigit=1)),id.vars='classnames',variable.name='Metric', value.name='classmean')
  classHIT_stats_sdform <- melt(setDT(digitform(classHIT_statssd,
                                                cols=2:(ncol(classHIT_statssd)), extradigit=1)),id.vars='classnames',variable.name='Metric', value.name='classsd')
  classHIT_format <- merge(classHIT_stats_meanform, classHIT_stats_sdform, by=c('classnames','Metric'))
  classHIT_format$tabcol <- with(classHIT_format, paste0(classmean,' (',classsd,')'))
  
  classHIT_format[is.na(classHIT_format$classsd),'tabcol']  <- '–'
  table <- as.data.frame(dcast(classHIT_format, Metric~classnames, value.var='tabcol'))
  metrics <-data.frame(type=substr(table$Metric,1,2), num=as.numeric(substr(table$Metric,3,5)))
  table <- table[with(metrics, order(type, num)),] #Order rows by metrics number
  table <- merge(table, metricKW[,c('Metric','Significance')], by='Metric')
  
  kable(table) %>%
    kable_styling(bootstrap_options = "striped", font_size = 10) %>%
    row_spec(which(table$Metric %in% HITo15ysub3$indice), bold=T) %>%
    column_spec(2,color=classcol[1]) %>%
    column_spec(3,color=classcol[2]) %>%
    column_spec(4,color=classcol[3]) %>%
    column_spec(5,color=classcol[4]) %>%
    column_spec(6,color=classcol[5]) %>%
    column_spec(7,color=classcol[6]) %>%
    column_spec(8,color=classcol[7]) %>%
    save_kable(tabname, self_contained=T)
}
############################ ANALYSIS FUNCTIONS ################################
#------ read_gaugep -----------------
#' Read gauge points
#'
#' Import streamgauging station spatial data and attributes for GRDC and GSIM.
#' Add new and updated environmental predictors from RiverATLAS v1.0.9 and, optionally,
#' estimates of monthly discharge from WaterGAP v2.2.
#'
#' @param inp_GRDCgaugep path to point data for GRDC gauging stations.
#' @param inp_riveratlas2 path to attribute table of RiverATLAS v1.0.9

#'
#' @return object of class \link[sf]{sf}
#'
#' @export
read_gaugep <- function(inp_GRDCgaugep,
                        inp_riveratlas2) {
  #Import gauge stations
  GRDCgaugep <- vect(dsn = dirname(inp_GRDCgaugep),
                        layer = basename(inp_GRDCgaugep))
  


  #Row bind GSIM and GRDC stations
  GRDCgaugep$area_correct <- GRDCgaugep$GRDC_AREA
  GRDCgaugep$gsim_no <- NA
  GSIMgaugep$GRDC_NO <- NA
  GRDCgaugep[, c("gsim_no", "reference_db", "reference_no", "grdb_merge",
                 "grdb_no", "paired_db", "paired_db_no")] <- NA
  keepcols_bind <- intersect(names(GSIMgaugep), names(GRDCgaugep))
  
  gaugep <- rbind(GRDCgaugep[, keepcols_bind],
                  GSIMgaugep[, keepcols_bind])
  
  #Get new and updated environmental predictors from RiverATLAS v1.0.9
  riveratlas2 <- fread(inp_riveratlas2)
  setnames(riveratlas2,
           names(riveratlas2),
           gsub('_11$', '', names(riveratlas2)))
  
  #Replace variables from RiverATLAS v1.0 by those that have been re-calculated
  #in v1.0.9
  keepcols <- names(gaugep)[!(names(gaugep) %in% names(riveratlas2))]
  gaugep_attriall <- merge(gaugep[,keepcols, with=F], riveratlas2,
                           by.x = 'HYRIV_ID', by.y = 'REACH_ID',
                           all.x=TRUE, all.y=FALSE)
  
  #Merge with WaterGAP downscaled monthly naturalized discharge
  if (!is.null(in_monthlydischarge)) {
    gaugep_outformat <- merge(gaugep_attriall, in_monthlydischarge,
                              by.x = 'HYRIV_ID', by.y = 'REACH_ID',
                              all.x=TRUE, all.y=FALSE)
  } else {
    gaugep_outformat <- gaugep_attriall
  }
  
  return(gaugep_outformat)
}

#------ read_GRDCgauged_paths -----------------
#' Read file paths to streamflow data from GRDC gauging stations
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @param inp_GRDCgaugedir path to directory containing streamflow data GRDC standard files.
#' @param in_gaugep table containing column named \code{GRDC_NO} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables, assuming
#' that files are called "GRDC_NO.txt", GRDC_NO being replaced with a 7-digit integer.
#'
#' @export
read_GRDCgauged_paths <- function(inp_GRDCgaugedir, in_gaugep) { #, gaugeid = 'GRDC_NO' down the line
  #Get data paths of daily records for gauge stations
  fileNames <- file.path(inp_GRDCgaugedir,
                         paste(
                           in_gaugep[!is.na(in_gaugep$GRDC_NO),]$GRDC_NO,
                           ".txt", sep=""))
  #Check whether any GRDC record does not exist
  print(paste(length(which(!do.call(rbind, lapply(fileNames, file.exists)))),
              'GRDC records do not exist...'))
  return(fileNames)
}

#------ comp_GRDCdurfreq -------------------------------
#' Compute intermittency statistics for GRDC gauging stations
#'
#' Determine general characteristics of the whole time series and of the subset
#' of years that have less than a given threshold of missing data as well as
#' intermittency statistics. The intermittency statistics can be computed for a
#' subset of months of the year (e.g. only winter months)
#'
#' @param path (character) file path to a GRDC-formatted streamflow time series table
#' @param in_gaugep (table; data.frame or data.table) of gauges' hydro-environmental attributes (including mean monthly temperature data)
#' @param maxgap (integer) maximum number of days with missing data beyond which a year is
#' not used in the computation of statistics
#' @param mdurthresh (numeric) threshold of mean annual number of zero-flow days beyond
#' which to classify gauge as intermittent.
#' @param windowsize (integer) window size to check for zero-flow days. 
#' @param fullwindow (logical) whether years for which the window is truncated 
#' (e.g., beginning and end of time series) are taken in account in moving window analysis.
#' @param monthsel (integer vector) selected months to compute the statistics over
#' @param verbose whether to print input path
#'
#' @return One row data.table with 110 columns: \cr
#' \describe{
#' \itemize{
#'   \item{GRDC_NO} - (char) unique identifier for the gauge
#'   \item{firstYear} - (num) first year on full record
#'   \item{lastYear} - (num) last year on full record
#'   \item{totalYears} - (int) total number of years on full record \cr
#'   For three subsets of the time series post-1800, post-1961, and post-1971 (e.g., suffix "mDur_o1800"):
#'   \itemize{
#'     \item{firstYear_kept} - (num) first year on record with < maxgap missing days
#'     \item{lastYear_kept} - (num) first year on record with < maxgap missing days
#'     \item{totalYears_kept} - (int) total number of years with < maxgap missing days
#'     \item{totaldays} - (num) total number of days with discharge data
#'     \item{integerperc} - (num) proportion of daily values which are in integer format. 
#'     Only including years with <= maxgap missing days.
#'     \item{sumDur} - (int) total number of days with discharge = 0
#'     \item{mDur} - (num) mean number of days/year with discharge = 0
#'     \item{mFreq} - (num) mean number of periods with discharge = 0
#'     \item{intermittent} - (factor): binary flow intermittence class. 1: non-perennial (if mDur >= 1, i.e., 
#'    if gauging station recorded zero-flow for at least one day per year on average); 0: perennial.
#'    `\item{movinginter} - (logical): whether there is at least one zero-flow day 
#'    in every \code{windowsize}-year (e.g., 20-year) moving window across the record. (with at least one day of flow between periods)
#'    `\item{monthly statistics} - (numeric) long-term mean frequency of zero-flow events 
#'    and number of zero-flow days for each month (column name example: 'Jan_mdur_o1800')
#'     \item{winteronlyir} - Indicates whether gauge is only non-perennial during winter months. 
#'     if intermittent == 1 AND the average annual number of zero-flow days during warm months < 1. 
#'     Warm months are those with mean monthly catchment air temperature >= 10 (WorldClim v2; Fick and Hijmans 2017).
#'     0: either perennial, or non-perennial outside of winter months.
#'     Only including years with <= maxgap missing daily discharge values.
#' }
#' }
#' }
#'
#' @export
comp_GRDCdurfreq <- function(path, in_gaugep, maxgap, mdurthresh = 1,
                             windowsize = 100, fullwindow = FALSE,
                             monthsel = NULL, verbose = FALSE) {
  if (verbose) {
    print(path)
  }
  
  #Read and format discharge records
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  gaugetab <- readformatGRDC(path)
  
  #Function to compute mean zero-flow duration and event frequency by month
  #(average number of drying events in each month over record length)
  #and check whether intermittency only happens in cold months
  comp_monthlyirtemp <- function(gaugetab, gaugeno, in_gaugep, maxgap,
                                 mdurthresh, tempthresh, yearthresh) {
    if (gaugetab[(missingdays < maxgap) & (year >= yearthresh), .N>0]) { #Make sure that there are years with sufficient data
      monthlyfreq <- gaugetab[
        (missingdays < maxgap)  & (year >= yearthresh),
        .(GRDC_NO = unique(GRDC_NO),
          monthrelfreq = length(unique(na.omit(prevflowdate)))/length(unique(year)),
          monthmdur = length(na.omit(prevflowdate))/length(unique(year))),
        by='month'] #Count proportion of years with zero flow occurrence per month
    } else {
      monthlyfreq <- data.table(GRDC_NO = gaugetab[, unique(GRDC_NO)],
                                month=1:12,
                                monthrelfreq=rep(NA,12),
                                monthmdur=rep(NA,12))
    }
    abbrev_months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    monthlyfreq_format <- monthlyfreq %>%
      dcast(GRDC_NO~month, value.var=c('monthrelfreq', 'monthmdur')) %>%
      setnames(c('GRDC_NO',
                 paste0(rep(abbrev_months, times=2), '_',
                        rep(c('mfreq', 'mdur'), each=12)
                 )
      ))
    
    #Check whether average yearly number of zero-flow days is only above threshold
    # if includes zero-flow days during months with average temperature < 5C
    gaugeirtemp <- as.data.table(in_gaugep)[GRDC_NO == gaugeno,
                                            c('GRDC_NO',
                                              grep('^tmp.*_c[0-9]{2}',
                                                   colnames(in_gaugep),
                                                   value = T)),
                                            with=F] %>%
      melt(id.vars='GRDC_NO') %>%
      .[, month := as.numeric(gsub('^tmp_dc_c', '', variable))] %>%
      merge(monthlyfreq, ., by='month')
    
    mdur_otempthresh <- gaugeirtemp[value >= 10*tempthresh, sum(monthmdur)]
    mdur_utempthresh <- gaugeirtemp[value < 10*tempthresh, sum(monthmdur)]
    
    monthlyfreq_format[, winteronlyir := as.numeric(
      mdur_otempthresh < mdurthresh &
        mdur_utempthresh >= mdurthresh)] %>%
      .[, GRDC_NO := NULL]
    
    setnames(monthlyfreq_format,
             paste0(names(monthlyfreq_format), '_o', yearthresh))
    
    return(monthlyfreq_format)
  }
  
  
  monthlyirtemp_all <- comp_monthlyirtemp(gaugetab, gaugeno, in_gaugep,
                                          maxgap, mdurthresh,
                                          yearthresh=1800, tempthresh=10)
  monthlyirtemp_o1961 <- comp_monthlyirtemp(gaugetab, gaugeno,in_gaugep,
                                            maxgap, mdurthresh,
                                            yearthresh=1961, tempthresh=10)
  monthlyirtemp_o1971 <- comp_monthlyirtemp(gaugetab, gaugeno,in_gaugep,
                                            maxgap, mdurthresh,
                                            yearthresh=1971, tempthresh=10)
  
  #If analysis is only performed on a subset of months
  if (!is.null(monthsel)) {
    gaugetab <- gaugetab[month %in% monthsel, ]
  }
  
  #Compute number of days of zero flow for years with number of gap days under threshold
  gaugetab_yearly <- merge(gaugetab[, .(missingdays=max(missingdays, na.rm=T),
                                        datadays=max(datadays, na.rm=T),
                                        integerperc=sum(integervalue,na.rm=T)
  ), by='year'],
  gaugetab[Original == 0, .(dur=.N,
                            freq=length(unique(prevflowdate))), by='year'],
  by = 'year', all.x = T) %>%
    .[!is.finite(missingdays), missingdays := diny(year)] %>%
    .[!is.finite(datadays), datadays := diny(year)] %>%
    .[dur==diny(year), freq:=1] %>%
    .[is.na(dur), `:=`(dur=0, freq=0)]
  
  
  #Combine all statistics (and determine which stations are labeled as
  #intermittent based on mdurthresh)
  gaugetab_all <- gaugetab_yearly[, .(GRDC_NO = gaugeno,
                                      firstYear=min(year),
                                      lastYear=max(year),
                                      totalYears=length(unique(year))
  )]
  
  irstats_all <- comp_irstats(tabyearly = gaugetab_yearly, maxgap=maxgap,
                              mdurthresh = mdurthresh,
                              yearthresh = 1800,
                              windowsize = windowsize,
                              fullwindow = fullwindow)
  irstats_1961 <- comp_irstats(tabyearly = gaugetab_yearly, maxgap=maxgap,
                               mdurthresh = mdurthresh,
                               yearthresh = 1961,
                               windowsize = windowsize,
                               fullwindow = fullwindow)
  irstats_1971 <- comp_irstats(tabyearly = gaugetab_yearly, maxgap=maxgap,
                               mdurthresh = mdurthresh,
                               yearthresh = 1971,
                               windowsize = windowsize,
                               fullwindow = fullwindow)
  
  statsout <- cbind(gaugetab_all,
                    irstats_all, irstats_1961, irstats_1971,
                    monthlyirtemp_all, monthlyirtemp_o1961, monthlyirtemp_o1971)
  
  #Include local path to discharge records in table
  statsout[, path := path]
  
  return(statsout)
}

#------ plot_GRDCflags ------
#' Plot all GRDC time series with data quality flags
#'
#' Creates pngs of streamflow time series plots of daily discharge ({m^3}/s; with flags 
#' for 0-flow values and potential outlier) for all GRDC gauges with at least 
#' 10 years of data, excluding years with more than 20 days of missing data. 
#'
#' @param in_GRDCgaugestats data.table (or list of data.tables) with time series 
#' and intermittency statistics for all GRDC gauging stations. 
#' @param yearthresh  (integer) minimum year from which to plot discharge record.
#' @param inp_resdir (character) path to the results directory in which to create folder and write output plots
#' @param maxgap (integer) threshold number of missing daily records to consider a calendar year unfit for analysis.
#' @param showmissing (logical) whether to show records in years with number of missing daily records beyond \code{maxgap}.
#'
#' @details the output graphs are written into two separate newly created directories in \code{inp_resdir} called
#' GRDCir_rawplots_\code{yearthresh}_\code{YYYYMMDD} and GRDCper_rawplots_\code{yearthresh}_\code{YYYYMMDD}
#' (e.g., GRDCir_rawplots_1800_20200512 and GRDCper_rawplots_1800_20200512). The first contains plots for non-perennial
#' gauging stations and the second contains plots for perennial gauging stations. \cr
#' \cr
#' Each plot is generated by \code{\link{plotGRDCtimeseries}} and shows
#' the time series of daily streamflow values for a station. 
#' For the flagging criteria, see documentation for \code{\link{flagGRDCoutliers}}.
#' \itemize{
#'   \item The y-axis is square-root transformed.
#'   \item Individual points show daily discharge values (in {m^3}/s).
#'   \item blue lines link daily values (which may result in unusual patterns due to missing years).
#'   \item red points are zero-flow flow values.
#'   \item green points are non-zero flow daily values statistically flagged as potential outliers .
#'   \item black points are zero-flow values flagged as potential outliers.
#' }
#'
#' @return nothing (empty data.table)
#'
#' @export
plot_GRDCflags <- function(in_GRDCgaugestats, yearthresh,
                           inp_resdir, maxgap, showmissing = FALSE) {
  if (inherits(in_GRDCgaugestats, 'data.table')) {
    GRDCstatsdt <- in_GRDCgaugestats
  }  else if (is.list(in_GRDCgaugestats)) {
    GRDCstatsdt <- rbindlist(in_GRDCgaugestats)
  }
  
  #Create output directory for IRs
  resdir_GRDCirplots <- file.path(inp_resdir,
                                  paste0('GRDCir_rawplots_', yearthresh, '_',
                                         format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_GRDCirplots))) {
    print(paste0('Creating ', resdir_GRDCirplots ))
    dir.create(resdir_GRDCirplots )
  }
  
  #Plot
  GRDCstatsdt[(get(paste0('totalYears_kept_o', yearthresh)) >= 10) &
                (get(paste0('intermittent_o', yearthresh)) == 1),
              plotGRDCtimeseries(.SD,
                                 outpath = file.path(resdir_GRDCirplots,
                                                     paste0(GRDC_NO, '.png')),
                                 maxgap=maxgap,
                                 showmissing = showmissing
              ), by=GRDC_NO]
  
  #Create output directory for non IRs
  resdir_GRDCperplots <- file.path(inp_resdir,
                                   paste0('GRDCper_rawplots_',yearthresh, '_',
                                          format(Sys.Date(), '%Y%m%d')))
  if (!(dir.exists(resdir_GRDCperplots))) {
    print(paste0('Creating ', resdir_GRDCperplots ))
    dir.create(resdir_GRDCperplots )
  }
  
  GRDCstatsdt[(get(paste0('totalYears_kept_o', yearthresh)) >= 10) &
                (get(paste0('intermittent_o', yearthresh)) == 0),
              plotGRDCtimeseries(.SD,
                                 outpath = file.path(resdir_GRDCperplots,
                                                     paste0(GRDC_NO, '.png')),
                                 maxgap=maxgap,
                                 showmissing = showmissing
              ), by=GRDC_NO]
  
  
}
#------ analyze_gaugeir ----------------------------
#' Quality assurance/quality checking of streamflow time series
#'
#' Record of quality-checking procedure which aimed to ensure the validity of 
#' zero-flow readings and the flow intermittence class assigned to each gauge 
#' (i.e., perennial or non-perennial).
#' 
#'
#' @param in_GRDCgaugestats list of data.tables of intermittency statistics, output from 
#' \code{\link{comp_GRDCdurfreq}} applied to all GRDC gauging stations.
#' @param in_GSIMgaugestats list of data.tables of intermittency statistics, output from 
#' \code{\link{comp_GSIMdurfreq}} applied to all GSIM gauging stations.
#' @param yearthresh (integer) minimum year from which to analyze discharge record.
#' @param in_gaugep \link[sf]{sf} object of gauging stations.
#' @param inp_resdir (character) path to the results directory in which to create folder and write output plots
#' @param plotseries (logical) whether to create plots of streamflow time series for gauging stations 
#' that are non-perennial only in winter or within 3 km from the coast.
#'
#' @details For details on the QA/QC procedure, see Supplementary Information at 
#' \link{https://www.nature.com/articles/s41586-021-03565-5}. The reason for 
#' station exclusion is also presented in an interactive online map at  
#' \link{https://messamat.github.io/globalIRmap/} in the Data and methods/Reference streamflow gauging stations tab.
#'
#' @return list of two data.tables. 
#' \itemize{
#'   \item A data.table listed as 'data' which contains the collated outputs from \code{\link{comp_GRDCdurfreq}}
#'   and\code{\link{comp_GSIMdurfreq}} for the stations that were not excluded through this QA/QC process.
#'   \item A data.table listed as flags, listing all stations excluded from further analysis and the reason their exclusion.
#' }
#'
#' @export
analyzemerge_gaugeir <- function(in_GRDCgaugestats, in_GSIMgaugestats, yearthresh,
                                 in_gaugep, inp_resdir, plotseries = FALSE) {


### Analyze GRDC data ########################################################################
GRDCstatsdt <- rbindlist(in_GRDCgaugestats)

#Remove all gauges with 0 values that have at least 99% of integer values as not reliable (see GRDC_NO 6140700 as example)
GRDCtoremove_allinteger <- data.table(
  GRDC_NO = GRDCstatsdt[integerperc_o1800 >= 0.95 &
                          intermittent_o1800 == 1, GRDC_NO],
  flag = 'removed',
  comment = 'All integer discharge values'
)

#Remove those which have at least one day per year of zero-flow day but instances
#of no zero-flow day within a 20-year window — except for three gauges that have a slight shift in values but are really IRES
GRDCtoremove_unstableIR <- data.table(
  GRDC_NO = GRDCstatsdt[(mDur_o1800 >= 1) & (!movinginter_o1800) &
                          !(GRDC_NO %in% c(1160115, 1160245, 4146400)), GRDC_NO],
  flag = 'removed',
  comment = 'automatic filtering: at least one no-flow day/year on average but no zero-flow event during >= 20 years'
)

#Outliers from examining plots of ir time series (those that were commented out were initially considered)
GRDCtoremove_irartifacts <- list(
  c(1104800, 'inspected', 'low-quality gauging but confirmed seasonally intermittent flow'),
  c(1134300, 'removed', 'changed flow permanence from perennial to non-perennial, large data gaps'),
  c(1134500, 'removed', 'only 1 occurrence of 0 flow values'),
  c(1159110, 'inspected', 'regulated but drying was confirmed upstream of reservoir'),
  c(1159120, 'inspected', 'at weird or dam; but drying was confirmed upstream of reservoir'),
  c(1159132, 'inspected', 'regulated but drying was confirmed upstream of dam'),
  c(1159302, 'removed', 'abrupt decreases to 0 flow values'),
  c(1159303, 'removed', 'unreliable record, isolated 0s, sudden jumps and capped at 77'),
  c(1159320, 'inspected', "0 values for the first 14 years but still apparently originally IRES"),
  c(1159325, 'removed', "0 values for most record. probably episodic and due to series of agricultural ponds"),
  c(1159510, 'removed', "0 values for most record, on same segment as 1159511 but seems unreliable"),
  c(1159520, 'inspected', "values seem capped after 1968, otherwise seem fine. Could just be rating curve"),
  c(1159830, 'removed', 'only one occurence of 0 flow values'),
  c(1160101, 'removed', 'abrupt decreases to 0 flow values'),
  c(1160210, 'inspected', 'lower plateaus in the beginning of record are likely 0-flow values'),
  c(1160245, 'removed', 'regulated, at reservoir'),
  c(1160301, 'inspected', 'lower plateaus in the beginning of record are likely 0-flow values'),
  c(1160340, 'removed', 'abrupt decreases to 0 flow values'),
  c(1160378, 'inspected', 'appears IRES before regulation, reservoir fully dry on satellite imagery, so keep as intermittent even when not regulated'),
  c(1160420, 'removed', 'decrease to 0 appears a bit abrupt'),
  c(1160435, 'removed', 'unreliable record, abrupt decreases to 0, capped'),
  c(1160470, 'removed', 'unreliable record, probably change of rating curve in 1947, mostly missing data until 1980 but truly intermittent based on imagery'),
  c(1160540, 'inspected', 'only 0 - nodata for first 15 years. seemingly good data post 1979 and IRES'),
  c(1160635, 'removed', 'valid 0 values only during early 80s'),
  c(1160670, 'removed', 'regulated'),
  c(1160675, 'removed', 'some outliers but otherwise most 0 values seem believable'),
  c(1160780, 'removed', 'unreliable record, abrupt decreases to 0 flow values, large data gaps'),
  c(1160775, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1160785, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1160793, 'removed', 'changed flow permanence from perennial to non-perennial, unreliable record'),
  c(1160795, 'removed', 'abrupt decreases to 0 flow values'),
  c(1160800, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1160840, 'removed', 'only 2 zero flow values are believable, others are outliers'),
  c(1160850, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1160880, 'removed', 'unreliable record. Tugela river, perennial'),
  c(1160881, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1160900, 'removed', 'most 0 values look like outliers, abrupt decreases'),
  c(1160911, 'removed', 'most 0 values look like outliers, abrupt decreases'),
  c(1160971, 'removed', 'most 0 values look like outliers, abrupt decreases'),
  c(1160975, 'removed', 'most 0 values look like outliers, abrupt decreases'),
  c(1196102, 'removed', 'unreliable record, large data gaps, hard to tell original flow permanence'),
  c(1196141, 'removed', "doesn't look reliable, hard to assess long term flow permanence"),
  c(1196160, 'removed', 'changed flow permanence, some outlying 0 flow values but most are good'),
  c(1197500, 'removed', 'only one flow intermittency event, abrupt decrease to 0'),
  c(1197540, 'removed', 'abrupt decreases to 0'),
  c(1197591, 'removed', 'abrupt decreases to 0'),
  c(1197700, 'removed', 'abrupt decreases to 0'),
  c(1197740, 'removed', 'some outlying 0 flow values but most are good'),
  c(1199100, 'removed', 'most 0 values look like outliers'),
  c(1199200, 'removed', 'abrupt decreases to 0'),
  c(1199410, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1234130, 'inspected', 'low quality record but confirmed intermittent by https://doi.org/10.3390/w11010156'),
  c(1259500, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1259800, 'removed', '0 values come from integer-based part of the record'),
  c(1286690, 'removed', 'changed flow permanence, record too short to determine original flow permanence'),
  c(1289230, 'removed', 'unreliable record, gaps, shifts'),
  c(1259800, 'removed', 'changed flow permanence, only one 0 flow value post 1963'),
  c(1428400, 'removed', '0 values come from integer-based part of the record'),
  c(1428500, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1434200, 'removed', 'almost all integers'),
  c(1434300, 'removed', 'almost all integers'),
  c(1434810, 'removed', '0 values come from integer-based part of the record'),
  c(1491790, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1491815, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1491870, 'removed', 'some outlying 0 flow values but most are good'),
  c(1494100, 'removed', 'abrupt decreases to 0'),
  c(1494100, 'inspected', 'regulated but naturally intermittent'),
  c(1495360, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1495700, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(1495720, 'removed', 'too many gaps, probably changed flow permanence from perennial to non-perennial'),
  c(1591110, 'removed', "doesn't look reliable, changed flow permanence from perennial to non-perennial"),
  c(1591730, 'removed', 'abrupt decreases to 0'),
  c(1733600, 'removed', '0 values come from integer-based part of the record and outliers'),
  c(1837410, 'removed', 'abrupt decreases to 0'),
  c(1837430, 'inspected', 'nearly same as 1837410. Naturally intermittent before dam'),
  c(1897550, 'removed', 'abrupt decreases to 0'),
  c(1898501, 'removed', 'abrupt decreases to 0'),
  c(1992840, 'removed', 'changed flow permanence'),
  c(1992400, 'removed', 'most 0 values look like outliers'),
  c(2181960, 'removed', 'series of small reservoirs upstream on both tributaries'),
  c(2588500, 'removed', 'abrupt decreases to 0'),
  c(2588551, 'removed', 'abrupt decreases to 0'),
  c(2588630, 'removed', 'abrupt decreases to 0'),
  c(2588640, 'removed', 'abrupt decreases to 0'),
  c(2588708, 'removed', 'abrupt decreases to 0'),
  c(2588820, 'removed', 'abrupt decreases to 0'),
  c(2589230, 'removed', 'abrupt decreases to 0'),
  c(2589370, 'removed', 'abrupt decreases to 0'),
  c(2591801, 'removed', 'abrupt decreases to 0'),
  c(2694450, 'removed', 'abrupt decreases to 0'),
  c(2969081, 'removed', 'abrupt decreases to 0'),
  c(2999920, 'removed', '0 values come from integer-based part of the record'),
  c(3650380, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(3650460, 'removed', 'large portion of integer values, probably changed flow permanence'),
  c(3650470, 'removed', '0 values come from integer-based part of the record'),
  c(3650475, 'removed', 'large portion of integer values, probably changed flow permanence'),
  c(3650610, 'removed', 'integers pre-1960s but still intermittent after'),
  c(3650640, 'removed', 'abrupt decreases to 0'),
  c(3650649, 'inspected', 'change of flow regime due to dam building but intermittent before'),
  c(3650690, 'inspected', 'abrupt decreases to 0'),
  c(3650860, 'removed', 'only one 0 flow event in first 28 years'),
  c(3650928, 'removed', 'most 0 values look like outliers'),
  c(3652050, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(3652135, 'removed', 'only one valid 0-flow event'),
  c(3652200, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(3844460, 'removed', 'abrupt decreases to 0'),
  c(3844460, 'removed', 'abrupt decreases to 0'),
  c(4101451, 'inspected', 'station downstream also has 0s'),
  c(4103700, 'removed', '0 values come from integer-based part of the record'),
  c(4150605, 'inspected', 'just downstream of lwesville dam in Dallas. previously intermittent as well but will be removed anyways as >50% dor'),
  c(4151513, 'inspected', 'looks regulated but will be removed as > 50% regulated'),
  c(4185051, 'removed', 'regulated'),
  c(4203820, 'removed', 'only 1 occurrence of 0 flow values'),
  c(4208043, 'removed', 'changed from perennial to non-perennial'),
  c(4208195, 'removed', 'unreliable record, 0 flow values stem from interpolation'),
  c(4208372, 'removed', 'abrupt decreases to 0 flow, probably data gaps'),
  c(4208585, 'removed', 'changed from perennial to non-perennial'),
  c(4208655, 'removed', 'insufficient data to tell flow permanence'),
  c(4208855, 'removed', 'insufficient data to tell flow permanence'),
  c(4208857, 'removed', 'only 1 occurrence of 0 flow values'),
  c(4213090, 'removed', 'only 1 occurrence of 0 flow values'),
  c(4213091, 'removed', 'only 2 occurrence of 0 flow values'),
  c(4213566, 'removed', 'only 1 occurrence of 0 flow values'),
  c(4213905, 'removed', "regulated, changed flow permanence"),
  c(4214075, 'removed', '0 flow values are data gaps'),
  c(4214200, 'removed', 'changed from perennial to non-perennial'),
  c(4214297, 'removed', 'only 1 occurrence of 0 flow values'),
  c(4214298, 'removed', 'only 1 occurrence of 0 flow values'),
  c(4234300, 'removed', "regulated, changed flow permanence"),
  c(4243610, 'removed', "regulated, abrupt decrease to 0 probably due to reservoir filling/construction"),
  c(4351710, 'removed', '0 values come from integer-based part of the record and outliers'),
  c(4355500, 'removed', "regulated, outlier 0 flow values"),
  c(4357510, 'removed', "single flow intermittency event, probably gap in data"),
  c(4769200, 'removed', 'only 1 occurrence of 0 flow values'),
  c(4773050, 'removed', 'abrupt decreases to 0'),
  c(5101020, 'removed', "single flow intermittency event, probably gap in data"),
  c(5101101, 'removed', "single flow intermittency event, probably gap in data"),
  c(5101130, 'removed', 'abrupt decreases to 0'),
  c(5101201, 'removed', 'unreliable record, large data gap as 0 flow values'),
  c(5101290, 'removed', '0 flow values before 2000 are outliers, changed flow permanence'),
  c(5101305, 'removed', 'most 0 values look like outliers'),
  c(5101380, 'removed', 'abrupt decreases to 0'),
  c(5109200, 'removed', 'unreliable record, interpolation, large data gap as 0 flow values'),
  c(5109230, 'removed', 'abrupt decreases to 0'),
  c(5202140, 'removed', 'abrupt decreases to 0'),
  c(5202145, 'removed', 'most 0 values look like outliers'),
  c(5202185, 'only one maybe erroneous 0 flow values in first 34 years'),
  c(5204140, 'only one 0 flow values in first 30 years'),
  c(5202228, 'removed', 'maybe regulated, unreliable record post 1983 accounts for 0 flow values'),
  c(5204170, 'removed', 'changed flow permanence'),
  c(5302251, 'removed', 'large data gap as 0 flow values, otherwise only one flow intermittency event'),
  c(5302261, 'removed', 'large data gap as 0 flow values'),
  c(5405095, 'removed', 'changed flow permanence from perennial to non-perennial'),
  c(5405105, 'only one 0 flow value'),
  c(5608100, 'removed', 'large data gap as 0 flow values'),
  c(5708200, 'removed', 'changed flow permanence'),
  c(5803160, 'removed', 'large data gap as 0 flow values'),
  c(5864500, 'removed', 'rounded to 10L/s'),
  c(5870100, 'removed', 'rounded to 100L/s'),
  c(6119100, 'removed', 'rounded to 10L/s'),
  c(6125680, 'removed', 'changed flow permanence'),
  c(6139140, 'removed', 'only one 0 flow occurrence in first 30 years'),
  c(6233410, 'removed', '0 flow values from tidal reversals'),
  c(6442300, 'removed', '0 values come from integer-based part of the record and outliers'),## perfect example of what an integer-based record involves
  c(6444250, 'removed', '0 values come from integer-based part of the record and outliers'),
  c(6444350, 'removed', '0 values come from integer-based part of the record and outliers'),
  c(6444400, 'removed', 'abrupt decreases to 0'),
  c(6935570, 'removed', 'rounded to 10L/s')
) %>%
  do.call(rbind, .) %>%
  as.data.table %>%
  setnames(c('GRDC_NO', 'flag', 'comment'))

#### Check intermittent record
# checkno <- 6444400 #GRDC_NO
# check <- checkGRDCzeroes( #Check area around 0 values
#   GRDCstatsdt, in_GRDC_NO=checkno, period=15, yearthresh=1800,
#   maxgap=20, in_scales='free', labelvals = F)
# checkno %in% GRDCtoremove_allinteger #Check whether all integers
# in_gaugep[in_gaugep$GRDC_NO==checkno & !is.na(in_gaugep$GRDC_NO), "dor_pc_pva"] #check DOR
# GRDCstatsdt[GRDC_NO == checkno, integerperc_o1800] #Check % integers

#Outliers from examining plots of perennial time series (those that were commented out were initially considered)
#Try to find those:
# whose low flow plateaus could be 0s
# whose perennial character is dam-driven or maybe irrigation driven (changed from IR to perennial but hard to find)
# whose missing data are actually 0s
# whose quality is too low to be reliable
GRDCtoremove_pereartifacts <- list(
  c(1159800, 'removed', 'regulated'),
  c(1160324, 'removed', 'regulated'),
  c(1160331, 'removed', 'low-flow plateaus are likely overestimated 0 values'),
  c(1160520, 'removed', 'regulated by diversions'),
  c(1160602, 'removed', 'regulated'),
  c(1160709, 'removed', 'regulated'),
  c(1160788, 'removed', 'low-flow plateaus may be overestimated 0 values'),
  c(1197310, 'removed', 'regulated'),
  c(1255100, 'inspected', 'now regulated, but naturally perennial, Cunene River'),
  c(1593100, 'inspected', 'bad quality but clearly not IRES'),
  c(1593751, 'inspected', 'missing values may contain intermittency, strange regular patterns, maybe interpolated'),
  c(2357750, 'inspected', 'not regulated, in Sri Lanka'),
  c(2588707, 'removed', 'regulated'),
  c(3628200, 'removed', 'appears to change flow permanence'),
  c(3650634, 'removed', 'regulated, probably changed flow permanence'),
  c(3652030, 'removed', '0s in missing years and low flows in other years may also be 0s'),
  c(4101200, 'removed', 'low-flow plateaus are likely overestimated 0 values'),
  c(4115225, 'removed', 'regulated, may have been IRES otherwise'),
  c(4118850, 'removed', 'low-flow plateaus may be overestimated 0 values'),
  c(4125903, 'removed', 'regulated, may have been IRES otherwise'),
  c(4126351, 'removed', 'regulated, may have been IRES otherwise'),
  c(4148850, 'removed', 'many 0 flow values in missing years'),
  c(4151801, 'removed', 'regulated, may have been IRES otherwise, Rio Grande'),
  c(4152651, 'removed', 'regulated by blue mesa reservoir, may have been IRES otherwise, Rio Grande'),
  c(4208610, 'removed', 'too much missing data but if not would be IRES'),
  c(4213055, 'removed', 'too much missing data but if not would be IRES'),
  c(4213802, 'removed', 'identical to 4213801'),
  c(4214320, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
  c(4362100, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
  c(5606090, 'removed', 'low-flow plateaus may be overestimated 0 values'),
  c(5606414, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
  c(6123630, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
  c(6335020, 'removed', 'identical to 6335060'),
  c(6335050, 'removed', 'identical to 6335060'),
  c(6337503, 'removed', 'regulated, cannot tell whether may have been intermittent before'),
  c(6442100, 'removed', 'identical to 6442600'),
  c(6935146, 'removed', 'identical to 6935145'),
  c(6335050, 'removed', 'identical to 6335060'),
  c(6935600, 'removed', 'identical to 6935145')
) %>%
  do.call(rbind, .) %>%
  as.data.table %>%
  setnames(c('GRDC_NO', 'flag', 'comment'))

#---------- Check flags in winter IR
plot_winterir(dt = GRDCstatsdt, dbname = 'grdc', inp_resdir = inp_resdir,
              yearthresh = 1800, plotseries = plotseries)
#Checked for seemingly anomalous 0s. Sudden decreases.
#Check for flags, check satellite imagery, station name, check for construction of reservoir

#------ Check time series of stations within 3 km of seawater
GRDCcoastalirall <- plot_coastalir(in_gaugep = in_gaugep, dt = GRDCstatsdt,
                                   dbname = 'grdc', inp_resdir = inp_resdir,
                                   yearthresh = 1800, plotseries = plotseries)
#GRDCcoastalirall[, unique(readformatGRDC(path)$Flag), by=GRDC_NO]
#Nothing obviously suspect beyond those that ad already been flagged

#Inspect statistics for 4208857, 4213531 as no flow days occurred only one year
# ID = '6976300'
# GRDCstatsdt[GRDC_NO == ID,]
# check <- readformatGRDC(GRDCstatsdt[GRDC_NO == ID,path])
# unique(check$Flag)
#
# plotGRDCtimeseries(GRDCstatsdt[GRDC_NO == ID,], outpath=NULL)


### Summarize removal ########################################################################

#Before cleaning
GRDCflags <- rbindlist(list(GRDCtoremove_allinteger,
                            GRDCtoremove_unstableIR,
                            GRDCtoremove_irartifacts,
                            GRDCtoremove_pereartifacts
))

GRDCtoremove_all <- GRDCflags[flag=='removed', GRDC_NO]

GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10, .N]
GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10 &
              !(GRDC_NO %in% GRDCtoremove_all), .N]

GRDCstatsdt_clean <- GRDCstatsdt[!(GRDC_NO %in% GRDCtoremove_all),]
}

#------ clean_record_dryver ----------------------------------------------------
#Use others from https://github.com/mahabbasi/europeanIRmap/blob/main/src/pre_processing_functions.R
clean_record_dryver <- function() {
  daily_q_records_to_remove = list(
    # c("station_id", "dataset_name","comments", "flag", "start_date", "end_date")
    c("UA_0000023", "GSIM","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("sm_1281", "smires","abrupt decrease to zero", 'months with zeros', "2009-01-01", "2013-12-31"),
    c("sm_1245", "smires","abrupt decrease to zero", 'months with zeros','1992-01-01', '1992-12-31'),
    c("sm_1181", "smires", "getting start with zero and then never been again", 'remove entire record', NA, NA),
    c("sm_1138", "smires", "getting start with zero and then never been again", 'remove entire record', NA, NA),
    c("sm_1124", "smires", "abrupt decrease to zero", 'months with zeros', NA, NA),
    c("sm_1118", "smires", "abrupt decrease to zero", 'months with zeros', '1987-04-01', '1987-04-30'),
    c("sm_1080", "smires", "abrupt decrease to zero", 'remove entire record', NA, NA),
    c("sm_1048", "smires", "abrupt decrease to zero", 'months with zeros', '1993-10-01', '1993-10-31'),
    c("sm_1038", "smires", "abrupt decrease to zero", 'months with zeros', '1984-09-01', '1984-09-30'),
    c("arpal_2004", "arpal","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("arpal_2002", "arpal","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6854103", "GRDC","abrupt decrease to zero & round records", 'remove entire record', NA, NA),
    c("6642100", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6503352", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6444500", "GRDC","abrupt decrease to zero & round records", 'remove entire record', NA, NA),
    c("6337530", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6233410", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6125310", 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c("6123760", 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c("6123641", "GRDC","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6123370", 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c("9263", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9231", "spanish_stations","start with zero values but never exprience again", 'remove entire record', NA, NA),
    c("9181", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9170", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9157", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9135", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9096", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("9090", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("9005", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("8132", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("8130", "spanish_stations","start with zero values but never exprience again", 'remove entire record', NA, NA),
    c("8071", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7102", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7062", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7050", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7006", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7004", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("7003", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("5084", "spanish_stations","abrupt decrease to zero", 'months with zeros', '2004-01-01', '2004-03-31'),
    c("5043", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("5039", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("5003", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("5002", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("4122", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3285", "spanish_stations","abrupt decrease to zero and rounded records to two digits", 'remove entire record', NA, NA),
    c("3283", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3270", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3263", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3252", "spanish_stations","abrupt decrease to zero and duplicated values in months", 'remove entire record', NA, NA),
    c("3155", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("3154", "spanish_stations","records is too incomplete and unreliable", 'remove entire record', NA, NA),
    c("3113", "spanish_stations","records is too incomplete and unreliable", 'remove entire record', NA, NA),
    c("2125", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("2064", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("1544", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("1464", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("1196", "spanish_stations","abrupt decrease to zero", 'months with zeros', NA, NA),
    c("6", "rbis","abrupt decrease to zero", 'months with zeros', NA, NA),
    c('6854510', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6683300', 'GRDC', "rounded to nearest m3", 'remove entire record', NA, NA),
    c('6683200', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6683010', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6682300', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6444400', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6444250', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('6373911', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6373460', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6373444', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6373224', 'GRDC', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c('6233740', 'GRDC', "abrupt decrease to zero", 'remove entire record', NA, NA),
    c('6221660', 'GRDC', "an anomalously jump in the period", 'remove entire record', NA, NA),
    c('6125300', 'GRDC', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('9292', 'spanish_stations', "abrupt decrease to zero. rounded to nearest m3", 'remove entire record', NA, NA),
    c('9145', 'spanish_stations', "start with zero values but later an anomalously jump in the period", 'remove entire record', NA, NA),
    c('9143', 'spanish_stations', "an anomalously jump in the period. record is unreliable", 'remove entire record', NA, NA),
    c('9105', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('7007', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3203', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3128', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3019', 'spanish_stations', "abrupt decrease to zero. rounded to two digits", 'remove entire record', NA, NA),
    c('3016', 'spanish_stations', "abrupt decrease to zero. records is too incomplete and  unreliable", 'remove entire record', NA, NA),
    c(5, 'rbis', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1264, 'spanish_stations','abrupt decrease to zero', 'months with zeros', '2007-01-01', '2007-01-31'),
    c(1295, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1398, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1431, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1438, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1446, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals',  'months with zeros', NA, NA),
    c(1485, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals. odd patterns probably due to regulation',  'months with zeros', NA, NA),
    c(1519, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals and rating curve discontinuity',  'months with zeros', NA, NA),
    c(1520, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(1542, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals and rating curve discontinuity',  'months with zeros', NA, NA),
    c(2000, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero','months with zeros', NA, '2001-01-01'),
    c(2003, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(2006, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2009, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '2001-01-01'),
    c(2015, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2016, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '1983-01-01'),
    c(2018, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2030, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2031, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '1985-01-01'),
    c(2033, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(2035, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', '2000-08-01', '2000-09-30'),
    c(2040, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2053, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '2001-01-01'),
    c(2054, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '1986-01-01'),
    c(2056, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals at the beginning of the time series','months with zeros', NA, '1986-01-01'),
    c(2066, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2068, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2081, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero','months with zeros', NA, '2008-01-01'),
    c(2083, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2093, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2097, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2102, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2109, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2129, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '2001-01-01'),
    c(2144, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2145, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(2149, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals',  'months with zeros', NA, NA),
    c(3001, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3009, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3015, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3030, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3065, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3067, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3127, 'spanish_stations', 'zero flows may be due to values being rounded to one decimals',  'remove entire record', NA, NA),
    c(3142, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3144, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros','2008-12-01', '2009-12-31'),
    c(3146, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'months with zeros', NA, '2010-01-01'),
    c(3152, 'spanish_stations', 'rounded to nearest m3 or first decimal', 'remove entire record', NA, NA),
    c(3159, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3162, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3175, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3180, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero','months with zeros', NA, '2010-01-01'),
    c(3189, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3190, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'remove entire record', NA, NA),
    c(3195, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'remove entire record', NA, NA),
    c(3196, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(3200, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', '2005-01-01', '2005-10-31'),
    c(3212, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '1987-01-01'),
    c(3262, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(3269, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(3940, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(4009, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(4013, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'remove entire record', NA, NA),
    c(4014, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4105, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4131, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4140, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4164, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4165, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, '2004-01-01'),
    c(4229, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4248, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(4283, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(4287, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(5008, 'spanish_stations', 'zero flows may be due to values being rounded to two decimals', 'remove entire record', NA, NA),
    c(5016, 'spanish_stations', 'record is too incomplete', 'months with zeros', NA, '2015-01-01'),
    c(5023, 'spanish_stations', 'record is too incomplete', 'months with zeros', NA, '2015-01-01'),
    c(5028, 'spanish_stations', 'record is too incomplete', 'remove entire record', NA, NA),
    c(5047, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, NA),
    c(5082, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '1985-01-01'),
    c(5083, 'spanish_stations', 'abrupt decrease and increase', 'months with zeros', NA, NA),
    c(5086, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(5101, 'spanish_stations', 'record is too incomplete', 'months with zeros', '2010-01-01', '2010-12-31'),
    c(5133, 'spanish_stations', 'rating curve discontinuity. abrupt decrease to zero', 'months with zeros', NA, '1983-01-01'),
    c(5138, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(7043, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(7121, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8028, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8032, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8042, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(8074, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '1988-01-01'),
    c(8148, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9014, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9018, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9035, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9039, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2010-01-01'),
    c(9049, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', '1986-09-01','1986-09-30'),
    c(9057, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2011-01-01'),
    c(9088, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, NA),
    c(9095, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, NA),
    c(9097, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9109, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(9115, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, NA),
    c(9118, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2010-01-01'),
    c(9124, 'spanish_stations', 'rating curve discontinuity', 'months with zeros', NA, '2010-01-01'),
    c(9138, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(9144, 'spanish_stations', 'rating curve discontinuity', 'remove entire record', NA, NA),
    c(9159, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', '1987-01-01', '1996-01-01'),
    c(9182, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9186, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(9201, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9257, 'spanish_stations', 'rating curve discontinuity','remove entire record', NA, NA),
    c(9259, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9287, 'spanish_stations', 'rating curve discontinuity','months with zeros', NA, NA),
    c(9321, 'spanish_stations', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6114500, 'GRDC',  'zero flows may be due to values being rounded to two decimals',  'months with zeros', NA, NA),
    c(6119100, 'GRDC', 'zero flows may be due to values being rounded to two decimals', 'months with zeros', NA, NA),
    c(6123461, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6123501, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124300, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124430, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124440, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6124502, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6125360, 'GRDC', 'abrupt decreases to zero at the end of the series', 'months with zeros', '1996-01-01', NA),
    c(6125440, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6128050, 'GRDC', 'zero flows may be due to values being rounded to two decimals. flow regulation started within record',  'months with zeros', '1989-01-01', NA),
    c(6128101, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6128220, 'GRDC', 'abrupt decrease to zeor after gap in series', 'months with zeros', '2007-04-01', '2007-05-01'),
    c(6139340, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6139501, 'GRDC', 'abrupt decrease to anomalously low values and zeros at the end of the series', 'months with zeros', '2009-01-01', NA),
    c(6139700, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6139850, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6144490, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6172028, 'GRDC', 'decrease to low discharge may be plausible but zero flows may be due to values being rounded to two decimals', 'months with zeros', NA, NA),
    c(6172031, 'GRDC', 'zero flows may be due to values being rounded to two decimals',  'remove entire record', NA, NA),
    c(6233128, 'GRDC', 'low-flows could be miscalibrated zeros', 'remove entire record', NA, NA),
    c(6233250, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6233414, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6242810, 'GRDC', 'abrupt decrease to zero. probably due to regulation for channel maintenance or other', 'months with zeros', NA, NA),
    c(6273752, 'GRDC', 'zero flows may be due to values being rounded to two decimal', 'months with zeros', NA, NA),
    c(6273900, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6274550, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6274655, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373040, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373050, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373200, 'GRDC', 'abrupt decrease to zero. unreliable record', 'remove entire record', NA, NA),
    c(6373215, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373217, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6373306, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6442300, 'GRDC', 'rounded to nearest m3 until 1992', 'months with zeros', NA, '1992-01-01'),
    c(6444350, 'GRDC', 'rounded to nearest m3 until 1992', 'months with zeros', NA, '1992-01-01'),
    c(6458420, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6574154, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6680400, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6680410, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6682500, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6729180, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6729225, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6729310, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6731556, 'GRDC', 'abrupt decrease to zero at the end of the series', 'months with zeros', '2015-01-01', NA),
    c(6731815, 'GRDC', 'abrupt decrease to zero at the end of the series', 'months with zeros', '2006-01-01', NA),
    c(6731940, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6854580, 'GRDC', 'rounded to one decimal until 1992', 'months with zeros', NA, '1990-01-01'),
    c(6854591, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6854660, 'GRDC', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c(6854708, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6854711, 'GRDC', 'rounded to nearest m3', 'remove entire record', NA, NA),
    c(6854722, 'GRDC', 'rounded to first decimal', 'remove entire record', NA, NA),
    c(6871100, 'GRDC', 'rounded to first decimal', 'remove entire record', NA, NA),
    c(6731940, 'GRDC', 'abrupt decrease to zero', 'remove entire record', NA, NA),
    c(6401120, 'GRDC', 'anomalously increased during the period', 'remove entire record', NA, NA),
    c('arpae_1013', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1016', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1017', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1018', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1020', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1024', 'emr', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c('arpae_1032', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1034', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1035', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1057', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1072', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1073', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpae_1079', 'emr', 'zero flows may be due to values being rounded to two decimal or to regulation', 'keep', NA, NA),
    c('arpal_2010', 'arpal', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('arpal_2028', 'arpal', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('arpas_3007', 'arpas', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('HS103', 'rbis', 'abrupt decrease to zero', 'months with zeros', NA, NA),
    c('PT_0000025', 'GSIM', 'zero flows may be due to values being rounded to two decimal or to regulation', 'months with zeros', NA, NA),
    c('sm_1003', 'smires', 'low-flows should probably be zero flows starting in 2016', 'change to zero', '2016-08-24', NA),
    c('sm_1171', 'smires', 'record is too incomplete', 'months with zeros', '2007-05-31', NA),
    c('sm_1174', 'smires', 'record is too incomplete', 'months with zeros', '2011-01-01', NA),
    c('sm_1180', 'smires', 'record is too incomplete', 'months with zeros', '2011-01-01', NA),
    c('sm_1098', 'smires', 'no-data for the reach- static predictors', 'remove entire record', NA, NA)
  )  %>% 
    do.call('rbind', .) %>% 
    as.data.table() %>% 
    `colnames<-`(c("station_id", "dataset_name","comments", "flag", "start_date", "end_date"))
  
  # Add the year and month columns to the data table to expand NA for the whole month
  gaugetab[, year_month := format(dates, "%Y-%m")]
  # Determine the different situations of the flags to exculd the gauges  
  gaugeid_removed <- daily_q_records_to_remove[flag == "remove entire record", station_id]
  
  gaugeid_zero_all <- daily_q_records_to_remove[flag == "months with zeros" &
                                                  is.na(start_date) & is.na(end_date), station_id]
  
  gaugeid_zero_part <- daily_q_records_to_remove[flag == "months with zeros" &
                                                   !is.na(start_date) | !is.na(end_date)] %>% 
    .[is.na(end_date), end_date := '2019-12-31'] %>% 
    .[is.na(start_date), start_date := '1981-01-01'] %>% 
    .[, c('start_date', 'end_date') := lapply(.SD, as.Date),
      .SDcols = c('start_date', 'end_date')]
  
  gaugeid_lowflow_zero <- daily_q_records_to_remove[flag == "change to zero" &
                                                      (!is.na(start_date) | !is.na(end_date))] %>% 
    .[is.na(end_date), end_date := '2019-12-31'] %>% 
    .[is.na(start_date), start_date := '1981-01-01'] %>% 
    .[, c('start_date', 'end_date') := lapply(.SD, as.Date),
      .SDcols = c('start_date', 'end_date')]
  
  # Implement the flags to the data table
  gaugetab[gaugeid %in% gaugeid_removed, value :=NA]
  
  gaugetab[gaugeid %in% gaugeid_zero_all & value == 0, value :=NA]
  
  for (iterator in 1:nrow(gaugeid_lowflow_zero)) {
    
    gaugeid_item <- gaugeid_lowflow_zero[iterator, station_id]
    start_date_item <- gaugeid_lowflow_zero[iterator, start_date]
    end_date_item <- gaugeid_lowflow_zero[iterator, end_date]
    gaugetab[gaugeid %in% gaugeid_item & 
               dates >= start_date_item & 
               dates <= end_date_item &
               value < 0.004,
             value:=0]
    
  }
  
  # Iterate over the data table with the gauges that modification requried in the middle
  # the period
  
  for (iterator in 1:nrow(gaugeid_zero_part)) {
    
    gaugeid_item <- gaugeid_zero_part[iterator, station_id]
    start_date_item <- gaugeid_zero_part[iterator, start_date]
    end_date_item <- gaugeid_zero_part[iterator, end_date]
    gaugetab[gaugeid %in% gaugeid_item & 
               dates >= start_date_item & 
               dates <= end_date_item &
               value == 0,
             value:=NA]
    
  }
  
  # Define the days of NA grouped by gauges and year_month then,
  # make all of the months NA with at least on one NA day
  gaugetab[, value := if (any(is.na(value))) NA else value, by = .(gaugeid, year_month)]
  
  
  
  # Select the non-NA records and return it
  output <- gaugetab[!is.na(value)] %>%
    .[,c('gaugeid', 'dates', 'value')]
  
  function_output <- list(output_ts = output,
                          gaugeids_removed = gaugeid_removed)
  
  return(function_output)
  
}

#------ format_gaugestats --------------------------------------------------------
#' Format gauge statistics
#'
#' Final selection of gauging stations and formatting of gauge statistics + 
#' hydro-environmental attributes.
#'
#' @param in_gaugestats data.table of summary time series statistics and intermittency statistics for 
#' all gauging stations (GRDC +  GSIM) not excluded in QA/QC process.
#' @param in_gaugep table of hydro-environmental attributes for all gauging stations.
#' @param yearthresh (integer) minimum year from which to analyze/consider discharge record.
#' 
#' @details An additional set of gauges are excluded in this step:
#' - All gauges with less than 10 years of data are excluded (considering only years with no more than 20 days of missing data)
#' - All gauges with a Degree of Regulation >= 50% are excluded (Lehner et al. 2011)
#' 
#' @return data.table of gauging stations included in all subsequent analysis, with
#' summary time series statistics and hydro-environmental attributes.
#'
#' @source Lehner, B., Liermann, C.R., Revenga, C., Vörösmarty, C., Fekete, B., 
#' Crouzet, P., Döll, P., Endejan, M., Frenken, K., Magome, J., Nilsson, C., 
#' Robertson, J.C., Rödel, R., Sindorf, N. and Wisser, D. (2011), 
#' High-resolution mapping of the world's reservoirs and dams for sustainable 
#' river-flow management. Frontiers in Ecology and the Environment, 9: 494-502. 
#' \link{https://doi.org/10.1890/100125}
#' 
#' @export

format_gaugestats <- function(in_gaugestats, in_gaugep, yearthresh) {
  #Join intermittency statistics to predictor variables and subset to only
  #include those gauges with at least 10 years of data
  gaugestats_join <- in_gaugestats[
    , GAUGE_NO := fifelse(is.na(gsim_no), GRDC_NO, gsim_no)] %>%
    .[!is.na(get(paste0('totalYears_kept_o', yearthresh))) &
        get(paste0('totalYears_kept_o', yearthresh))>=10,] %>%  # Only keep stations with at least 10 years of data pas yearthresh
    merge(as.data.table(in_gaugep)[, -c('GRDC_NO', 'gsim_no'), with=F],
          by='GAUGE_NO', all.x=T, all.y=F) %>%
    .[, c('X', 'Y') := as.data.table(sf::st_coordinates(geometry))] %>%
    .[, DApercdiff := (area_correct-UPLAND_SKM)/UPLAND_SKM]
  #Check for multiple stations on same HydroSHEDS segment
  dupliseg <- gaugestats_join[duplicated(HYRIV_ID) |
                                duplicated(HYRIV_ID, fromLast = T),] %>%
    .[, interdiff := length(unique(
      eval(paste0('intermittent_o', yearthresh))))-1, by=HYRIV_ID]
  
  
  #Only two stations have differing status. Inspect their record
  # dupliseg[interdiff==1,
  #          c('HYRIV_ID', 'GAUGE_NO', 'DApercdiff',
  #            'mDur_o1961', paste0('intermittent_o', yearthresh)), with=F]
  # plotGRDCtimeseries(dupliseg[GAUGE_NO == '1197560',])
  # plotGRDCtimeseries(dupliseg[GAUGE_NO == '1197591',]) #Remove because lots of erroneous data
  # plotGRDCtimeseries(dupliseg[GAUGE_NO == '4150605',])
  # plotGSIMtimeseries(dupliseg[GAUGE_NO == 'US_0006104',]) #Remove because identical record but without precise zero flow day count
  
  gaugestats_joinsel <- gaugestats_join[!GAUGE_NO %in% c('1197591', 'US_0006104')] %>%
    setorder(HYRIV_ID, DApercdiffabs) %>%
    .[!duplicated(HYRIV_ID),]
  
  print(paste0('Removing ', gaugestats_joinsel[dor_pc_pva >= 500, .N],
               ' stations with >=50% flow regulation'))
  gaugestats_derivedvar <- gaugestats_joinsel[dor_pc_pva < 500, ] %>% #Only keep stations that have less than 50% of their discharge regulated by reservoir
    comp_derivedvar #Compute derived variables, rescale some variables, remove -9999
  
  return(gaugestats_derivedvar)
}


########################### TO DEAL WITH #######################################
####TUN THE FOLLOWING COMMENTED-OUT BLOCK IN ORDER TO GET DIAGNOSTIC PLOTS USED IN DATA QA/QCing####
# for (gage in unique(rufidat_screenform$ID)) {
#   print(gage)
#   gname <- as.character(unique(rufidat[rufidat$Gage.ID==gage,'Station.Name']))
#   #Generate FlowScreen time series
#   gts<- create.ts(rufidat_screenform[rufidat_screenform$ID==gage,]) #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
#   #Compute and output flowScreen metrics and plots
#   try({
#     res <- metrics.all(gts)
#     ginfo <- data.frame(StationID=gage, StnName=gname, ProvState='Rufiji Basin',Country='Tanzania',Lat=0, Long=0, Area=0, RHN='RBWB')
#     png(file.path(outdir,paste(gage,'screenb.png',sep="_")),width = 20, height=12,units='in',res=300)
#     screen.summary(res, type="b", StnInfo=ginfo)
#     dev.off()
#     png(file.path(outdir,paste(gage,'screenl.png',sep="_")),width = 20, height=12,units='in',res=300)
#     screen.summary(res, type="l", StnInfo=ginfo)
#     dev.off()
#     png(file.path(outdir,paste(gage,'screenh.png',sep="_")),width = 20, height=12,units='in',res=300)
#     screen.summary(res, type="h", StnInfo=ginfo)
#     dev.off()
#   })
#   #Fit Savitzky-Golay 1st order derivative
#   # p = polynomial order w = window size (must be odd) m = m-th derivative (0 = smoothing) 
#   d1 <- as.data.frame(savitzkyGolay(gts$Flow, p = 3, w = 21, m = 1)) 
#   #A shorter period than 21 days would be used but it seems like granularity of stage measurements at low flow makes them appear very constant
#   gts[11:(nrow(gts)-10),'sg.d1'] <- d1
#   gts[gts$sg.d1<(10^-10) & gts$sg.d1>-(10^-10) & !is.na(gts$sg.d1),'Flag'] <- 'Y'
#   
#   #Make raw time series plot
#   rawsgplot <-ggplot(gts, aes(x=Date, y=Flow)) + 
#     geom_point(color='#045a8d', size=1) + 
#     geom_point(data=gts[gts$Flag=='Y',],aes(x=Date, y=Flow), color='#d01c8b', size=1.5) +
#     geom_point(data=gts[gts$Flow==0,],aes(x=Date, y=Flow), color='#e31a1c', size=1.5) +
#     scale_y_sqrt()+
#     scale_x_date(date_breaks = "2 years", date_labels = "%Y") + 
#     labs(y='Discharge (m3/s)', title=paste(gage, gname,sep=" - ")) +
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust=1))
#   png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 20, height=12,units='in',res=300)
#   print(rawsgplot)
#   dev.off()
# }


########################## Interpolate data with na.interp ########################
# CustomImpute_nainterp <- function(tscast, sn, maxgap,pplot=F) {
#   if (sn > 1) {
#     print(colnames(tscast)[sn])
#     #Only select records within 37 days of the first and last date of records for the gauge
#     mindate <- tscast[min(which(!is.na(tscast[,sn]))),'Date']-37
#     maxdate <- tscast[max(which(!is.na(tscast[,sn]))),'Date']+37
#     tscastsub <- tscast[tscast$Date>mindate & tscast$Date<maxdate,]
#     tscastsub$prevdate <- tscastsub[na.lomf(tscastsub[,sn]),'Date']
#     tscastsub$nextdate <- tscastsub[na.lomb(tscastsub[,sn]),'Date']
#     tscastsub$gap <- as.numeric(as.Date(tscastsub$nextdate)-as.Date(tscastsub$prevdate))
#     tscastsub <- as.data.frame(tscastsub)
#     
#     #Determine Box-Cox lamba for tsl decomposition
#     bc <- BoxCox.lambda(tscastsub[,sn]+10,method='loglik',lower=0)
#     #Interpolate data
#     pred_try <- data.frame(pred=na.interp(ts(tscastsub[,sn]+10,frequency=365), lambda=bc)-10)
#     #Bound interpolation to min and max if outside bound of record
#     pred_try[pred_try$pred<min(tscastsub[,sn],na.rm=T),'pred'] <- min(tscastsub[,sn],na.rm=T) #if interpolated value is outside bound, bound it
#     pred_try[pred_try$pred>max(tscastsub[,sn],na.rm=T),'pred'] <- max(tscastsub[,sn],na.rm=T)
#     #Round abnormally low discharge values to 0
#     pred_try[pred_try$pred<0.01,'pred'] <- 0
#     #Create output dataframe
#     out<-data.frame(tscastsub[,'Date'],pred_try$pred, tscastsub$gap)
#     gapcol <- paste('gap',colnames(tscastsub)[sn],sep='_') 
#     colnames(out) <- c("Date",colnames(tscastsub)[sn],gapcol)
#     
#     #Re-convert records for gaps longer than 'maxgap' to NA
#     out[out[[gapcol]]>=maxgap, colnames(tscastsub)[sn]] <- NA
#     
#     #Flag interpolated vs. original data
#     srccol <-  paste('source',colnames(tscastsub)[sn],sep='_') 
#     out[out[[gapcol]]>0,srccol] <- 'interpolated' 
#     out[out[[gapcol]]==0, srccol] <- 'original'
#     
#     #Plot original and interpolated data for subsequent QA/QC
#     if (pplot==T){
#       print(ggplot(out, aes(x=Date, y=get(colnames(tscastsub)[sn]), color=get(srccol))) +
#               geom_point()+
#               ggtitle(colnames(tscastsub)[sn]) +
#               scale_y_sqrt(name='Discharge (cms)')+
#               scale_x_date(limits=c(as.Date(min(rufidat_cast[!is.na(tscast[,sn]),'Date'])),
#                                     as.Date(max(rufidat_cast[!is.na(tscast[,sn]),'Date']))))) #Only plot the area with NAs
#     }
#     return(out)
#   }
# }
# 
# #Fill in every gauge
# impute_preds <- data.frame(Date=rufidat_cast[,"Date"])
# for (i in 2:(ncol(rufidat_cast))) {
#   try({
#     impute_preds<-merge(impute_preds, CustomImpute_nainterp(tscast=rufidat_cast,sn=i,maxgap=37, pplot=T), by='Date', all.x=T)
#   })
# }

# #Function to linearly interpolate for a set of dates.
# interpcor <- function(df, station, cordates) {
#   df[df$Date %in% cordates, station] <- NA
#   df[,station] <- na.interpolation(df[,station], option='linear')
#   gapcol <- paste('gap',station,sep='_') 
#   df[df[[gapcol]]>37 | is.na(df[[gapcol]]),station] <- NA
#   return(df)
# }


# #Develop ARIMAX model with Fourier series added as external regressors to model the seasonal pattern
# #and precipitation as a second external regressor.ARIMA (in R at least) cannot model seasonal periods > 350.
# #Inspired from https://robjhyndman.com/hyndsight/longseasonality/
# #and https://robjhyndman.com/hyndsight/forecasting-weekly-data/
# ts1KA41 <- msts(precip1KA41sub[, 'Flow'],seasonal.periods=365.25)
# bestfit <- list(aicc=Inf)
# for(i in 1:50) #Select the number of Fourier series by minimizing AIC 
# {
#   fit <- auto.arima(ts1KA41, seasonal=FALSE, xreg=cbind(fourier(ts1KA41, K=i), precip1KA41sub[, 'Precip']))
#   if(fit$aicc < bestfit$aicc){
#     print(i)
#     bestfit <- fit #Keep model with lowest AIC
#   }else {break};
# }
# kr <- KalmanSmooth(ts1KA41, bestfit$model) #Impute missing values with Kalman Smoother (from https://stats.stackexchange.com/questions/104565/how-to-use-auto-arima-to-impute-missing-values)
# id.na <- which(is.na(ts1KA41) & precip1KA41sub$gap<274) #Limit to gaps < 9 months
# pred <- ts1KA41
# for (i in id.na)
#   pred[i] <- bestfit$model$Z %*% kr$smooth[i,]
# precip1KA41sub[id.na,'Flow'] <- pred[id.na] #Replace NA values in the time series by predicted values
# precip1KA41sub[precip1KA41sub$Flow<0.05 & !is.na(precip1KA41sub$Flow),'Flow'] <- 0 #For values < 0 and improbably low values, replace with 0
# ggplot(precip1KA41sub[,], aes(x=Date, y=Flow)) + geom_point() + #Plot result
#   geom_point(data=precip1KA41sub[id.na,], color='red')
# 
# ggplot(precip1KA41sub[,], aes(x=Date, y=Flow, color=gap)) + geom_point() +
#   scale_colour_distiller(name='Gap (years)',palette='Spectral',breaks=seq(0,365.25,30),
#                          limits=c(min(precip1KA41sub$gap),max(precip1KA41sub$gap)))#Plot result
# 
# impute_preds2 <- merge(impute_preds, precip1KA41sub[,c('Date','Flow')], by='Date', all.x=T)
# impute_preds2[impute_preds2$Date>mindate & impute_preds2$Date<maxdate,'1KA41'] <- impute_preds2[impute_preds2$Date>mindate & impute_preds2$Date<maxdate,'Flow']
# impute_preds2 <- impute_preds2[,-which(colnames(impute_preds2)=='Flow')]
# impute_preds2[impute_preds2$gap_1KA41>=37 & impute_preds2$gap_1KA41<274 & !is.na(impute_preds2$gap_1KA41), 'source_1KA41'] <- 'ARIMAX_interpolated'


##################################Gap plot################
#Compute number of valid years on record depending on the percentage of missing data tolerated to consider a year valid and the corresponding
#maximum gap.
rufidat_gapyear <- data.frame(ID=unique(rufidat_gapsummary$ID))
for (i in seq(0.5,0,-0.1)) {
  df<-ddply(rufidat_gapsummary[rufidat_gapsummary$gap_per<=i & rufidat_gapsummary$max_gap<=(i*365.25),], .(ID), summarise, gcount=length(hyear))
  rufidat_gapyear <- merge(rufidat_gapyear,df, by='ID',all.x=T)
}
colnames(rufidat_gapyear) <- c('ID',paste('gagecount_gap', seq(0.5,0,-0.1),sep="_"))

#Compute for a range of durations, the number of stations that have data for at least this duration, for each percentage gap threshold
rufidat_gapplot <- ldply(seq(1,50,1), function(y) {
  adply(rufidat_gapyear[,2:ncol(rufidat_gapyear)], 2, function(x) length(which(x>y)))})
rufidat_gapplot$minyr <- sort(rep(seq(1,50,1), ncol(rufidat_gapyear[,2:ncol(rufidat_gapyear)])))
rufidat_gapplot$threshgap <- as.numeric(substr(rufidat_gapplot$X1,15,18))

gapplot <- ggplot(rufidat_gapplot, aes(x=minyr, y=V1, color=as.factor(100-100*threshgap))) + 
  scale_x_continuous(name='Record length (years)', breaks=c(1,seq(5,60,5)), expand=c(0,0)) +
  scale_y_continuous(name='Number of stream gauges', limits=c(0,40), breaks=seq(0,40,5),expand=c(0,0))+
  scale_color_discrete(name="Minimum yearly record completeness (% of days)") +
  guides(color=guide_legend(ncol=4)) +
  geom_line(size=1.5) +
  theme_bw() + 
  theme(text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position=c(0.65,0.85),
        legend.background = element_blank(),
        panel.grid.minor=element_blank(),
        legend.key.size = unit(3,"line"),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
png(file.path(outdir,'gapplot20180703.png'),width = 12, height=12,units='in',res=300)
print(gapplot)
dev.off()

plotseries <- function(gage){ #Make a graph of a time series highlighting deleted, interpolated, used and non-used data
  print(gage)
  genv <- gagesenv[gagesenv$RGS_No==gage,]
  #Generate FlowScreen time series
  gts<- create.ts(predsmelt[predsmelt$ID==gage,])  #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
  gts_sel <- merge(gts, rufidat_ycount, by='ID', all.x=T)
  gts_sel <- merge(gts_sel, rufidat_gapsummary[,c('ID','hyear','gap_per')], by=c('ID','hyear'))
  print(genv$RGC_Loc)
  #Make raw time series plot
  rawsgplot <-ggplot() +
    geom_point(data=gts_sel, aes(x=Date, y=Flow), color='#bf812d', size=0.5)+ 
    geom_point(data=rufidat_deleted[rufidat_deleted$ID==gage,], aes(x=Date, y=Flow),color='#e31a1c', size=0.5)+
    geom_point(data=rufidat_clean[rufidat_clean$ID==gage,], aes(x=Date, y=Flow),color='#9ecae1', size=0.5)+
    geom_point(data=rufidat_clean[rufidat_clean$ID==gage & 
                                    rufidat_clean$hyear %in%  as.data.frame(rufidat_gapsummary)[rufidat_gapsummary$ID==gage & 
                                                                                                  rufidat_gapsummary$gap_per<=0.1,'hyear'],],
               aes(x=Date, y=Flow),color='#045a8d', size=0.5)+
    scale_y_sqrt(expand=c(0,0),limits=c(0,max(gts$Flow)+1))+
    scale_x_date(date_breaks = "2 years", date_labels = "%Y", limits=as.Date(c('1954-01-01','2018-06-01')), expand=c(0,0)) + 
    labs(y='Discharge (m3/s)', 
         title=bquote(paste(.(gage)-.(as.character(genv$RGS_Loc)) ~'River at' ~ .(as.character(genv$RGS_Name)),'.')))+
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = 'none',
          text= element_text(size=8,lineheight = 0),
          plot.title= element_text(size=8,margin = margin(t = 0, r = 0, b = 0, l = 0)),
          plot.margin=unit(c(-0,0,-0,-0), "cm"),
          axis.title.x= element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.y=element_text(vjust=-0.5),
          panel.grid.minor= element_blank(),
          panel.border = element_blank(),
          axis.line = element_line()
    )
  p <- ggplot_gtable(ggplot_build(rawsgplot))
  lay= t(c(rep(1,4),2))
  png(file.path(outdir,paste(gage,'raw_sg.png',sep="_")),width = 6.5, height=4.5,units='in',res=300)
  print(grid.arrange(p, legendts, ncol = 11, layout_matrix = lay))
  dev.off()
}

plotflowscreen <- function(gage, div,thrs){ #make graphs of FlowScreen package outputs
  print(gage)
  genv <- gagesenv[gagesenv$RGS_No==gage,]
  #Generate FlowScreen time series
  gts<- create.ts(predsmelt[predsmelt$ID==gage,])  #Cannot run ts on multiple gages. Need to first subset by gage, then run ts.
  gname <- paste(genv$RGS_Loc," river at ",genv$RGS_Name,sep="")
  gts$Flow <- gts$Flow/div
  
  #Compute and output flowScreen metrics and plots
  try({
    res <- metrics.all(gts,NAthresh=thrs)
    ginfo <- data.frame(StationID=genv$RGS_No, StnName=gname, ProvState='Rufiji Basin',Country='Tanzania',
                        Lat=genv$POINT_Y, Long=genv$POINT_X, Area=genv$WsArea, RHN='RBWB')
    png(file.path(outdir,paste(gage,'screenb.png',sep="_")),width = 20, height=12,units='in',res=300)
    screen.summary(res, type="b", StnInfo=ginfo)
    dev.off()
    png(file.path(outdir,paste(gage,'screenl.png',sep="_")),width = 20, height=12,units='in',res=300)
    screen.summary(res, type="l", StnInfo=ginfo)
    dev.off()
    png(file.path(outdir,paste(gage,'screenh.png',sep="_")),width = 20, height=12,units='in',res=300)
    screen.summary(res, type="h", StnInfo=ginfo)
    dev.off()
  })
}
for (g in unique(predsmelt$ID)) {
  plotseries(g)
  plotflowscreen(g,div=1,thrs=0.5)
}

