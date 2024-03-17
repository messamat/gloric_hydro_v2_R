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

#------ digitform --------------------------------
#Get the number of significant digits based on the log10 of the median of each metric, if number >=10, no sig digit
digitform <- function(df, cols, extradigit=0, inplace=F) {
  if (!inplace) {
    df <- copy(df)
  }
  df[, (cols) := sapply(.SD, function(x) {
    as.character(round(x,
                       digits=ifelse(
                         median(as.numeric(x),na.rm=T)==0,
                         0,
                         ifelse(floor(log10(abs(median(as.numeric(x), na.rm=T))))>0,
                                0,abs(floor(log10(abs(median(as.numeric(x), na.rm=T)))))))
                       +extradigit)) 
    
  }, simplify=F),
  .SDcols=cols]
  
  return(df)
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


#------ fill_dt_dates ----------------------------------------------------------
#Make sure that the date column is complete for every year, either bounding
#by the actual date limits of the record or the 1st January of the year at the start
#and December 31st of the last year at the end
fill_dt_dates <- function(in_dt, full_yrs, date_col='date') {
  if (!('year' %in% names(in_dt))) {
    in_dt[, year := format(get(date_col), '%Y')]
  }
  
  yrs_to_keep <- in_dt[, unique(year)]
  
  if (full_yrs) {
    fullyrs_dt <- in_dt[, list(
      V1 = seq(
        as.IDate(paste0(min(year),
                        '-01-01')),
        as.IDate(paste0(max(year), 
                        '-12-31')),
        by='day')
    )] 
  } else {
    fullyrs_dt <- in_dt[, list(
      V1=seq(
        min(get(date_col)), 
        max(get(date_col)),
        by='day')
    )] 
  }
  
  fullyrs_dt[, year := as.integer(format(V1, '%Y'))] 
  setnames(fullyrs_dt, old='V1', new=date_col)
  
  in_dt <- merge(in_dt, fullyrs_dt, 
                 by=c(date_col,'year'), all.x=T, all.y=T) %>%
    .[year %in% yrs_to_keep,]
  
  return(in_dt)
}
#------ transform_scale_vars ---------------------------------------------------
transform_scale_vars <- function(in_dt, value_col='value', var_col=NULL, 
                                 inplace=FALSE, samp_frac=NULL,
                                 min_boxcox=-1, max_boxcox=2) {
  if (!inplace) {
    dt <- copy(in_dt)
  } else {
    dt <- in_dt
  }
  
  #Make sure 
  dt[, stat_min := min(get(value_col), na.rm=T), by=var_col] %>%
    .[, value_floored := fifelse(stat_min<0, get(value_col)-stat_min, get(value_col)),
      by=var_col] %>%
    .[, stat_non0min_floored := .SD[value_floored>0, min(value_floored, na.rm=T)],
      by=var_col] %>%
    .[, 
      value_pos := fifelse(stat_min<=0,
                           value_floored + stat_non0min_floored/2,
                           value_floored),
      by=var_col]
  
  #Transform all variables through monotonous trans to approach normal distribution
  #Get BoxCox lambda values rounded to the nearest 0.5 for simplicity/reproduceability 
  #(i.e., these values are less likely to change than finer ones)
  if (!is.null(samp_frac)) {
    samp <- function(x) {
      base::sample(x, samp_frac*length(x))
    }
  } else {
    samp <- function(x) {
      x
    }
  }
  
  if (!is.null(var_col)) {
    dt[
      , bc_lambda := 0.5*round(
        Rfast::bc(samp(value_pos), low=min_boxcox, up=max_boxcox)/0.5), 
      by=var_col]
    
    dt[, value_trans := fifelse(
      bc_lambda == 0,
      log(value_pos),
      ((value_pos^bc_lambda)-1)/bc_lambda
    ), by=var_col]
  } else {
    bc_lambda <- 0.5*round(
      Rfast::bc(samp(dt$value_pos), low=min_boxcox, up=max_boxcox)/0.5)
    
    if (bc_lambda[[1]] == 0) {
      dt[, value_trans := log(value_pos), 
               by=var_col]
    } else {
      dt[, value_trans := ((value_pos^bc_lambda[[1]])-1)/bc_lambda[[1]],
               by=var_col]
    }
  }
  
  if (is.null(var_col)) {
    trans_mean <- mean(dt$value_trans, na.rm=T)
    trans_sd <- sd(dt$value_trans, na.rm=T)
  }
  
  #Z-scale stats
  dt[
    , value_scaled := scale(value_trans, center=TRUE, scale=TRUE), 
    by=var_col]
  
  return_list <- list()
  
  if (!inplace) {
    return_list$dt <- dt
  }
  
  if (is.null(var_col)) {
    return_list$min_val <- dt$stat_min[[1]]
    return_list$min_floor <- dt$stat_non0min_floored[[1]]/2
    return_list$trans_mean <- trans_mean
    return_list$trans_sd <- trans_sd
    return_list$bc_lambda <- bc_lambda[[1]]
  }
  
  if (is.null(var_col)) {
    dt[, (value_col) := value_scaled] %>%
      .[, `:=`(stat_min=NULL,
               value_floored=NULL,
               stat_non0min_floored=NULL,
               value_pos=NULL,
               value_trans=NULL),]
  }
  
  return(return_list)
}


#------ prettydend_classes ----------------------------------------------------
#Make a nice looking dendogram based on a clustering output
prettydend_classes <- function(in_hclust, colorder=NULL, 
                               in_colors=NULL, in_labels=NULL,
                               in_kclass=7, classnames = NULL,
                               order_clusters=F) {
  
  classr <- dendextend::cutree(in_hclust, k=in_kclass, 
                               order_clusters_as_data = order_clusters)
  classr_df <- data.frame(ID=names(classr), gclass=classr) 
  
  if (!is.null(in_labels)) {
    in_hclust$labels <- in_labels
  }
  
  if (is.null(colorder)) colorder = 1:in_kclass
  
  # Choose the appropriate value for h based on the heights
  chosen_h <- in_hclust$height[length(in_hclust$height) - (in_kclass-1)]
  dendname_cut <- as.dendrogram(in_hclust) %>%
    cut(h=chosen_h)
  
  #Get a basic dendrogram with the right colors
  ggdendro_p <- dendname_cut$upper %>%
    dendextend::set("branches_lwd", 1) %>% #Set is also present in data.table, so important to use dendextend:: as package order varies in tar_make
    dendextend::color_branches(k=in_kclass, col=in_colors[colorder]) %>%
    dendextend::color_labels(k=in_kclass, col='white')  %>% #Create background labels in white to create a halo
    dendextend::as.ggdend(.)
  
  #Get cutoff length
  dend_segs <- ggdendro_p$segments  #[segment(dendr)$y>chosen_h,]
  hori_segs <- dend_segs[dend_segs$y==dend_segs$yend,]
  
  #Correct label positions that got mangled through "as.ggdend"
  new_leaves <- dend_segs[dend_segs$yend < min(hori_segs$yend),]
  names(new_leaves) <- paste0(names(new_leaves), 'leaves')
  ggdendro_p$labels$x <- new_leaves$xendleaves 
  ggdendro_p$labels$y <- min(hori_segs$yend)
  
  #Get labels to superimpose upon halo
  new_classlabels <- ggdendro_p$labels
  
  if (!is.null(classnames)) {
    classr_df <- merge(classr_df, classnames, by='gclass')  %>%
      as.data.table %>%
      .[, gclass := NULL] %>%
      setnames('classnames', 'gclass')
    new_classlabels$label <- paste("Class", classnames$classnames)
  } 
  
  new_classlabels$col <- in_colors[colorder]
  
  #Background rectangle for labels
  rect_df <- data.frame(
    xmin=min(dend_segs$x)-5,
    xmax=max(dend_segs$x)+5,
    ymin=0.5*min(hori_segs$yend), 
    ymax=0.95*min(hori_segs$yend)
  )
  
  #Format ggplot
  ggdendro_p_format <- ggplot(ggdendro_p) +
    scale_y_reverse(name="Euclidean distance")   + 
    geom_rect(data=rect_df, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill='white') +
    geom_text(data = new_classlabels, 
              aes(x = x, y = 0.94*min(hori_segs$yend), label = label, 
                  colour = col, size = cex), 
              hjust=0, angle = 0) +
    coord_flip(ylim=c(limits=c(max(dend_segs$yend), 0.75*min(hori_segs$yend))), #min(hori_segs$yend)
               clip='on') +
    theme_classic() +
    theme(plot.margin = margin(0,1,0,0, 'cm'),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  # ggdendro_p_format$layers[[3]]$aes_params$angle <- 0
  # ggdendro_p_format$layers[[3]]$aes_params$hjust<- 0.75
  ggdendro_p_format$layers[[3]] <- NULL
  
  return(list(classes=classr_df, 
              plot=ggdendro_p_format))
}
#------ rundiagnose_clustering ------------------------------------------------
rundiagnose_clustering <- function(in_mat, in_dist, in_method, min_nc, max_nc,
                                   dist_name="Euclidean distance"
) {
  print(in_method)
  hclust_res <- hclust(in_dist, method=in_method)
  
  #Compute cophenetic correlation --------------------------------------------
  cophcor <- cor(in_dist, cophenetic(hclust_res))
  
  dist_cophcor_dt <- merge(
    reshape2::melt(as.matrix(in_dist)),
    reshape2::melt(as.matrix(cophenetic(hclust_res))),
    by=c('Var1', 'Var2')) %>%
    setnames(c('value.x', 'value.y'), 
             c(dist_name, "Cophenetic dissimilarity"))
  
  #Plot cophenetic correlation------------------------------------------------
  p_cophcor <- ggplot(dist_cophcor_dt, 
                      aes(x=get(dist_name), y=`Cophenetic dissimilarity`)) +
    geom_point() +
    geom_smooth(method='lm') +
    annotate('text', x = 0.5, y=0.1,
             label=paste('Cophenetic correlation =', 
                         round(cophcor, 2))) +
    coord_cartesian(
      expand=F, 
      ylim=c(0, max(dist_cophcor_dt$`Cophenetic dissimilarity`)+0.05)) +
    theme_classic()
  
  #Mantel test ---------------------------------------------------------------
  clust_mantel <- vegan::mantel(
    xdis = in_dist,
    ydis = cophenetic(hclust_res),
    method = "pearson",
    permutations = 999
  )
  
  #Graph scree plot ----------------------------------------------------------
  scree_dt <- data.table(height=hclust_res$height,
                         groups=length(hclust_res$height):1)
  
  p_scree <- ggplot(scree_dt,
                    aes(x=groups, y=height)) +
    geom_point() +
    geom_line() + 
    scale_x_log10(breaks=c(1,5,10,20,50,100,max(scree_dt$groups))) +
    theme_bw()
  
  #Agglomerative coefficient -------------------------------------------------
  if (in_method !='median') {
    ag_coef <- coef.hclust(hclust_res)
  } else {
    ag_coef <- NA
  }
  
  #NbClust to determine number of clusters -----------------------------------
  # vc_methods <- c("kl","ch", "hartigan","ccc", "scott","marriot","trcovw", 
  #                "tracew","friedman", "rubin", "cindex", "db", "silhouette", 
  #                "duda", "beale", "ratkowsky", "ball", "ptbiserial", "pseudot2", 
  #                "gap", "frey", "mcclain", # "gamma", 
  #                "gplus", "tau", "dunn", 
  #                "hubert", "sdindex", "dindex", "sdbw", "alllong")
  # 
  # nbclust_tests <- lapply(vc_methods, function(vc_method){
  #   print(vc_method)
  #   tryCatch({
  #     en.nb <- NbClust(in_mat, distance = "euclidean", min.nc = min_nc,
  #                      max.nc = max_nc, method =in_method, 
  #                      index = vc_method)
  #     return(as.numeric(en.nb$Best.nc[1]))
  #   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # })
  
  nbclust_tests <- NbClust(data=in_mat, distance='euclidean',
                           #index=idx_totest,
                           min.nc=min_nc, max.nc=max_nc, method=in_method)
  
  return(list(
    hclust=hclust_res,
    cophcor = cophcor,
    p_cophcor = p_cophcor,
    ag_coef = ag_coef,
    scree_dt = scree_dt,
    p_scree = p_scree,
    clust_mantel = clust_mantel,
    nbclust_tests = nbclust_tests
  ))
}

#------ cuttree_and_visualize --------------------------------------------------
cuttree_and_visualize <- function(in_mat, in_hclust, in_kclass, in_colors, 
                                  in_meltdt, id_col, 
                                  value_col='value', variable_col='variable',
                                  aspect_col = 'aspect', 
                                  classnames= NULL, colorder=NULL) {
  #Get dendrogram and class membership
  dendo_format <- prettydend_classes(in_hclust = in_hclust, 
                                    in_kclass=in_kclass, 
                                    in_colors=in_colors,
                                    classnames=classnames,
                                    colorder=colorder)
  
  class_stats <- dendo_format$classes %>%
    as.data.table %>%
    .[, classn := .N, by='gclass'] %>%
    merge(in_meltdt, ., by.x=id_col, by.y='ID') %>%
    .[, gclass := factor(gclass, levels=sort(unique(gclass)))] %>%
    .[variable %in% colnames(in_mat),]
  
  #Test classification significance ------------------------------------------
  mat_dt <- in_mat %>%
    as.data.table %>%
    .[, (id_col) := rownames(in_mat)] %>%
    merge(dendo_format$classes, 
          by.x=id_col, by.y='ID') 
  
  #Check if classes are multivariate normal
  mshapirotest_classes <- lapply(
    mat_dt[, .N, by=gclass][N>=20, gclass],
    function(in_gclass) {
      #print(in_gclass)
      subdt <- mat_dt[gclass==in_gclass,]
      
      #Replace NA values
      subdt[, (names(subdt)) := sapply(
        .SD, function(x) suppressWarnings(zoo::na.aggregate(x)), simplify=F)]
      
      mat <- as.matrix(subdt[
        , -c(id_col, names(subdt)[!subdt[, sapply(.SD, uniqueN)] > 1]),
        with=F])
      rownames(mat) <- subdt[, get(id_col)]
      
      
      varcormat <- cor(mat)
      varcormat[lower.tri(varcormat)] <- NA
      
      vartoremove <- reshape2::melt(varcormat) %>%
        as.data.table %>%
        .[Var1 != Var2,] %>%
        .[abs(value) > 0.98, as.character(Var2)]
      
      if (length(vartoremove)>1) {
        sub_mat <- mat[, -which(colnames(mat) %in% vartoremove)]
      } else {
        sub_mat <- mat
      }
      
      shapiro_test <- t(sub_mat) %>%
        mvnormtest::mshapiro.test()
      
      return(list(
        gclass=in_gclass,
        shapiro_test=shapiro_test)
      )
    })
  
  #Test significant difference by stat
  # set up model
  diff_letters <- lapply(unique(class_stats$variable), function(var) {
    sub_dat <- class_stats[variable==var,]
    
    emmeans(lm(as.formula(paste(value_col, '~ gclass')), 
               data = sub_dat),
            specs = as.formula('~ gclass') 
    ) %>%
      cld(adjust = "sidak",
          Letters = letters,
          alpha = 0.05) %>%
      as.data.table %>%
      .[, variable := var]  %>%
      setnames('.group', 'grp_letter')
  }) %>% 
    rbindlist %>%
    .[, grp_letter := str_trim(grp_letter)] %>%
    merge(class_stats[!duplicated(paste0(variable, gclass)),
                      .(variable, aspect, gclass, classn)], 
          by=c('variable', 'gclass')) 
  
  class_stats[, facet_labels := paste0()]
  
  #Get box plot
  p_cluster_boxplot <- ggplot(
    class_stats, 
    aes(x=gclass, y=get(value_col), color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA, alpha=0.75) +
    geom_text(data=diff_letters, 
              aes(x=gclass, label=grp_letter, y=emmean), 
              color='black', size=3, alpha=0.5) +
    scale_y_continuous(name='Metric value') +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') + 
    facet_wrap(factor(get(variable_col),
                      levels=levels(class_stats[[variable_col]]))
               ~get(aspect_col), 
               scales='free', ncol=4) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.line = element_line(color='darkgrey'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(margin = margin(0,0,0.1,0, "cm")))
  
  return(list(
    class_dt = class_stats,
    p_dendo = dendo_format$plot,
    p_boxplot = p_cluster_boxplot,
    class_stats = diff_letters
  )
  )
}

#------ selectformat_predvars -----------------
#' Select + format predictor variables
#'
#' Select candidate predictor variables to use in subsequent modelling and create
#' a table of formatted names.
#'
#' @param inp_riveratlas_meta path to metadata table for river atlas variables
#' @param in_gaugestats data.table of formatted gauging station summary statistics and hydro-environmental attributes.
#' 
#' 
#' @return data.table of selected predictor variable codes, names, category, 
#' attribute, sources, references, etc. Most sources are from RiverATLAS 
#' technical documentation available at https://www.hydrosheds.org/page/hydroatlas.
#' 
#' @export
selectformat_predvars <- function(inp_riveratlas_meta, in_gaugestats) {
  #---- List predictor variables ----
  predcols<- c(
    #monthlydischarge_preds,
    'UPLAND_SKM',
    'slo_dg_cav',
    'slo_dg_uav',
    'dis_m3_pyr',
    'dis_m3_pmn',
    'dis_m3_pmx',
    'dis_m3_pvar',
    'dis_m3_pvaryr',
    'run_mm_cyr',
    'runc_ix_cyr', #runoff coefficient (runoff/precipitation)
    'sdis_ms_uyr', #specific discharge
    'sdis_ms_umn',
    'inu_pc_umn',
    'inu_pc_umx',
    'inu_pc_cmn',
    'lka_pc_cse',
    'lka_pc_use',
    #'dor_pc_pva', #anthropogenic - degree of regulation
    'ele_pc_rel',
    'slp_dg_cav',
    'slp_dg_uav',
    'clz_cl_cmj',
    'snw_pc_uyr',
    'snw_pc_cyr',
    'snw_pc_cmx',
    'glc_cl_cmj',
    
    'pnv_cl_cmj',
    'for_pc_use',
    'for_pc_cse',
    'gla_pc_use',
    'gla_pc_cse',
    'prm_pc_use',
    'prm_pc_cse',
    'swc_pc_uyr',
    'swc_pc_cyr',
    'swc_pc_cmn',
    'lit_cl_cmj',
    'kar_pc_use',
    'kar_pc_cse',
    
    'cly_pc_cav',
    'cly_pc_uav',
    'slt_pc_cav',
    'slt_pc_uav',
    'snd_pc_cav',
    'snd_pc_uav',
    
    'ari_ix_cav',
    'ari_ix_uav',
    'bio1_dc_cav',
    'bio2_dc_cav',
    'bio3_dc_cav',
    'bio4_dc_cav',
    'bio5_dc_cav',
    'bio6_dc_cav',
    'bio7_dc_cav',
    'bio10_dc_cav',
    'bio11_dc_cav',
    'bio12_mm_cav',
    'bio13_mm_cav',
    'bio14_mm_cav',
    'bio15_mm_cav',
    
    'bio1_dc_uav',
    'bio2_dc_uav',
    'bio3_dc_uav',
    'bio4_dc_uav',
    'bio5_dc_uav',
    'bio6_dc_uav',
    'bio7_dc_uav',
    'bio10_dc_uav',
    'bio11_dc_uav',
    'bio12_mm_uav',
    'bio13_mm_uav',
    'bio14_mm_uav',
    'bio15_mm_uav'
  )
  
  #---- Associate HydroATLAS column names with variables names ----
  
  #Get predictor variable names
  metaall <- readxl::read_xlsx(inp_riveratlas_meta,
                               sheet='Overall') %>%
    setDT
  
  metascale <- readxl::read_xlsx(inp_riveratlas_meta,
                                 sheet='scale') %>%
    setDT %>%
    setnames(c('Key','Spatial representation'),
             c('Keyscale', 'Spatial.representation'))
  
  metastat <- readxl::read_xlsx(inp_riveratlas_meta,
                                sheet='stat') %>%
    setDT %>%
    setnames(c('Key','Temporal or statistical aggregation or other association'),
             c('Keystat', 'Temporal.or.statistical.aggregation.or.other.association'))
  
  meta_format <- as.data.table(expand.grid(`Column(s)`=metaall$`Column(s)`,
                                           Keyscale=metascale$Keyscale,
                                           Keystat=metastat$Keystat)) %>%
    .[metaall, on='Column(s)'] %>%
    .[metascale, on = 'Keyscale'] %>%
    .[metastat, on = 'Keystat',
      allow.cartesian=TRUE]
  
  meta_format[, `:=`(
    unit = substr(`Column(s)`, 5, 6),
    varcode = paste0(gsub('[-]{3}', '', `Column(s)`),
                     Keyscale,
                     fifelse(grepl("[0-9]",Keystat),
                             str_pad(Keystat, 2, side='left', pad='0'),
                             Keystat)),
    varname = paste(Attribute,
                    Spatial.representation,
                    Temporal.or.statistical.aggregation.or.other.association))]
  
  #Add newly generated variables to meta_format (variable labels)
  addedvars <- data.table(varname=c('Precipitation catchment Annual min/max',
                                    'Discharge watershed Annual min/max',
                                    'Discharge watershed Annual min/average',
                                    'Elevation catchment average - watershed average',
                                    'Runoff coefficient catchment Annual average',
                                    'Specific discharge watershed Annual average',
                                    'Specific discharge watershed Annual min',
                                    'Drainage area',
                                    'Groundwater table depth catchment average'),
                          varcode=c('pre_mm_cvar',
                                    'dis_m3_pvar', 
                                    'dis_m3_pvaryr',
                                    'ele_pc_rel',
                                    'runc_ix_cyr',
                                    'sdis_ms_uyr', 
                                    'sdis_ms_umn',
                                    'UPLAND_SKM',
                                    'gwt_m_cav'
                          )
  )
  
  oldcolnames <- c('Spatial.representation',
                   'Temporal.or.statistical.aggregation.or.other.association',
                   'Source Data')
  newcolnames <- c('Spatial representation',
                   'Temporal/Statistical aggreg.',
                   'Source')
  
  predcols_dt <- merge(data.table(varcode=predcols),
                       rbind(meta_format, addedvars, fill=T),
                       by='varcode', all.x=T, all.y=F)   %>%
    setnames(oldcolnames, newcolnames) %>%
    setorder(Category, Attribute,
             `Spatial representation`, `Temporal/Statistical aggreg.`)
  
  #Format table
  predcols_dt[varcode=='UPLAND_SKM', `:=`(
    Category = 'Physiography',
    Attribute= 'Drainage Area',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='',
    Source = 'HydroSHEDS',
    Citation = 'Lehner & Grill 2013'
  )]
  
  predcols_dt[varcode=='dis_m3_pvar', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='dis_m3_pvaryr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='runc_ix_cyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Runoff coefficient',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2, WorldClim v2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='sdis_ms_uyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='sdis_ms_umn', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='mn',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='ele_pc_rel', `:=`(
    Category = 'Physiography',
    Attribute= 'Elevation',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='(cav-uav)/uav',
    Source = 'EarthEnv-DEM90',
    Citation = 'Robinson et al. 2014'
  )]
  
  predcols_dt[varcode=='cmi_ix_cvar', `:=`(
    Category = 'Climate',
    Attribute= 'Climate Moisture Index',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WorldClim v2 & Global-PET v2',
    Citation = 'Fick et al. 2017',
    varname = 'Climate moisture index catchment monthly mn/mx'
  )]
  
  predcols_dt[varcode=='cmi_ix_uvar', `:=`(
    Category = 'Climate',
    Attribute= 'Climate Moisture Index',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WorldClim v2 & Global-PET v2',
    Citation = 'Fick et al. 2017',
    varname = 'Climate moisture index watershed monthly mn/mx'
  )]
  
  predcols_dt[varcode=='gwt_m_cav', `:=`(
    Category = 'Hydrology',
    Attribute= 'Groundwater table depth',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='av',
    Source = 'Global Groundwater Map',
    Citation = 'Fan et al. 2013'
  )]
  
  #Remove duplicates (which were created with keystat meaning different things e.g; 09 meaning september, 2009, class 9)
  predcols_dtnodupli<- predcols_dt[!(
    (Category == "Climate" &
       grepl('Class.*', `Temporal/Statistical aggreg.`)) |
      (Category == "Landcover" &
         !grepl('(Class|Spatial).*', `Temporal/Statistical aggreg.`))
  ),]  %>%
    unique(by='varcode')
  
  return(predcols_dtnodupli)
}

#------ comp_derivedvar -----------------
#' Compute derived variables
#'
#' Format and compute derived a set of environmental variables from input
#' dataset that contains
#' \href{https://www.hydrosheds.org/page/hydroatlas}{HydroATLAS} variables for
#' the purpose of predicting intermittency.
#'
#' @param dt data.table which contains rows as records and HydroATLAS
#'   environmental variables as columns.
#' @param copy (logical) whether to make a copy of \code{dt} (if TRUE) or to
#'   modify it in place (if FALSE)
#'
#' @details This function only formats and compute derived variables deemed
#'   relevant for the intermittency analysis. \cr
#'   \cr
#'   Steps include: \cr
#'   1. Convert some -9999 that should be 0 after checking them
#'   2. Count the number of -9999 values per column
#'   3. Convert the rest of the -9999 values to NA
#'   4. Compute a set of derived variables based on existing variables in
#'   HydroATLAS
#'   (e.g. \code{pre_mm_cvar= fifelse(pre_mm_cmx==0, 0, pre_mm_cmn/pre_mm_cmx)})
#'   5. Correct scaling from RiverATLAS
#'   (e.g. Degree of regulation is out of 10,000 in HydroATLAS)
#'
#' ```
#' @source Linke, S., Lehner, B., Dallaire, C. O., Ariwi, J., Grill, G., Anand,
#'   M., ... & Tan, F. (2019). Global hydro-environmental sub-basin and river
#'   reach characteristics at high spatial resolution. Scientific Data, 6(1),
#'   1-15.
#'
#' @export
comp_derivedvar <- function(in_dt, copy=FALSE) {
  if (copy) {
    in_dt2 <- copy(in_dt)
  } else {
    in_dt2 <- in_dt
  }
  
  
  #---- Inspect and correct -9999 and NA values ----
  print('Inspect and correct -9999 values')
  #check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
  in_dt2[snw_pc_cyr == -9999, snw_pc_cyr:=0]
  in_dt2[snw_pc_cmx == -9999, snw_pc_cmx:=0]
  
  #Places with very high slope value (in Greenland) have -9999 - replace with high value
  #in_dt2[, max(slp_dg_uav)]
  in_dt2[slp_dg_cav == -9999, `:=`(slp_dg_cav = 750,
                                   slp_dg_uav = 650)]
  
  #bio5 is -9999 in a few places in Greenland. Set at lowest existing value
  # in_dt2[bio5_dc_cav>-9999, min(bio5_dc_cav)]
  # in_dt2[bio5_dc_uav != -9999, min(bio5_dc_uav)]
  in_dt2[bio5_dc_cav == -9999, bio5_dc_cav := -950]
  in_dt2[bio5_dc_uav == -9999, bio5_dc_uav:= -9500]
  
  print('Number of NA values per column')
  colNAs<- in_dt2[, lapply(.SD, function(x) sum(is.na(x)))]
  print(colNAs)
  
  print('Number of -9999 values per column')
  col9999<- in_dt2[, lapply(.SD, function(x) sum(x==-9999))]
  print(col9999)
  
  #-9999 in cly_pc_cav, slt, and snd are places with no soil mask (urban areas, lakes, glaciers, etc.)
  
  
  #Define column groups
  gladcols <- unlist(lapply(c('cav', 'uav'),function(s) {
    lapply(c('wloss', 'wdryp', 'wwetp', 'whfrq', 'wseas', 'wperm', 'wfresh',
             'wwpix', 'wdpix'),
           function(v) {
             paste0(v, '_pc_', s)
           })
  })
  )
  
  sgcols <- c('cly_pc_cav','slt_pc_cav', 'snd_pc_cav',
              'cly_pc_uav','slt_pc_uav', 'snd_pc_uav')
  
  #Bioclim columns in celsius degrees
  biocolsdc <- unlist(lapply(c('cav', 'uav'),
                             function(s) paste0('bio', c(1:7,10,11), '_dc_', s)
  ))
  biocolsnegative <- grep('bio[189][01]*_dc_uav', biocolsdc, value=T)
  
  # cmi_cmcols <- paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0))
  # cmi_umcols <- paste0('cmi_ix_u', str_pad(1:12, width=2, side='left', pad=0))
  
  #Convert -9999 to NAs
  for (j in which(sapply(in_dt2,is.numeric))) { #Iterate through numeric column indices
    if (!(j %in% which(names(in_dt2) %in% c(gladcols, sgcols, 'wet_cl_cmj')))) {
      set(in_dt2,which(in_dt2[[j]]==-9999),j, NA) #Set those to NA if -9999
    }
  }
  
  #Scale variables based on HydroATLAS v1.0 documentation and v1.0.9 processing
  in_dt2[, `:=`(
    ari_ix_cav = ari_ix_cav/1000,
    ari_ix_uav = ari_ix_uav/1000,
    lka_pc_cse = lka_pc_cse/10,
    lka_pc_use = lka_pc_use/10
  )]
  
  in_dt2[, (biocolsdc) := lapply(.SD, function(x) (x/100)), .SDcols = biocolsdc]
  in_dt2[, (biocolsnegative) := lapply(.SD, function(x) x-100),
         .SDcols=biocolsnegative]
  
  #---- Compute derived predictor variables ----
  print('Compute derived predictor variables')
  in_dt2[, `:=`(
    #min/max monthly watershed discharge
    dis_m3_pvar=fifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx),
    #min monthly/average yearly watershed discharge
    dis_m3_pvaryr=fifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),
    #runoff coefficient (runoff/precipitation)
    runc_ix_cyr = fifelse(bio12_mm_cav==0, 0, run_mm_cyr/bio12_mm_cav),
    #Specific discharge
    sdis_ms_uyr = dis_m3_pyr/UPLAND_SKM,
    sdis_ms_umn = dis_m3_pmn/UPLAND_SKM
  )]
  return(in_dt2)
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
                                 breaks=c(0, 100, 500, 1000, 5000, 10000),
                                 labels=c(0, 100, 500, 1000, 5000, 10000)),
    bio14_mm_uav  = scale_x_sqrt(expand=c(0,0),
                                 breaks = c(0, 10, 50, 100, 200, 400),
                                 labels=c(0, 10, 50, 100, 200, 400)),
    cly_pc_uav = scale_x_continuous(labels=scales::percent_format(scale=1),
                                    expand=c(0,0)),
    cmi_ix_uyr = scale_x_continuous(),
    dis_m3_pyr = scale_x_log10(breaks=c(0.1, 1, 10, 100, 1000),
                               labels=c(0.1, 1, 10, 100, 1000),
                               expand=c(0,0)),
    dor_pc_pva = scale_x_continuous(labels=scales::percent_format(scale=1),
                                    expand=c(0,0)),
    for_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=scales::percent_format(scale=1),
                              expand=c(0,0)),
    gla_pc_use = scale_x_continuous(labels=scales::percent_format(scale=1),
                                    expand=c(0,0)),
    kar_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=scales::percent_format(scale=1),
                              expand=c(0,0)),
    lka_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=scales::percent_format(scale=1),
                              expand=c(0,0)),
    pet_mm_uyr = scale_x_continuous(expand=c(0,0)),
    sdis_ms_uyr = scale_x_sqrt(breaks=c(0.01, 0.05, 0.1, 0.2),
                               expand=c(0,0)),
    snw_pc_uyr = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=scales::percent_format(scale=1),
                              expand=c(0,0)),
    run_mm_cyr = scale_x_continuous(expand=c(0,0)),
    swc_pc_uyr = scale_x_continuous(labels=scales::percent_format(scale=1),
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
    ire_pc_use = scale_x_continuous(labels=scales::percent_format(scale=1),
                                    expand=c(0,0))
  ) %>%
    .[(names(.) %in% names(in_df)) & names(.) %in% varstoplot]
  #Only keep those variables that are actually in df and that we want to plot
  
  scales_y <- unlist(rep(list(scale_y_continuous(expand=c(0,0))),
                         labels = scales::scientific_format(),
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
      scale_fill_manual(values=c('#35978f', '#bf812d'))
    
  } else if (vartoplot == "glc_pc_u16") {
    rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
                       by=glc_pc_u16]
    gclz <- in_gaugedt[,.N/in_gaugedt[,.N],by=glc_pc_u16]
    bindclz <- rbind(rivclz, gclz, idcol='source')%>%
      setnames(c( 'source', vartoplot, 'density'))
    
    penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
      geom_bar(aes(fill=as.factor(source)), stat='identity',
               position = 'identity', alpha=1/2, width=.6) +
      scale_fill_manual(values=c('#35978f', '#bf812d'))
    #
    #     penvhist <- ggplot(in_gaugedt, aes_string(x=vartoplot)) +
    #       geom_histogram(data=in_rivdt, aes(weight = LENGTH_KM),
    #                      fill='#2b8cbe', alpha=0.5, bins=101) +
    #       geom_histogram(fill='#dd3497', alpha=0.5, bins=101)
    
  } else {
    penvhist <- ggplot(in_gaugedt, aes_string(x=vartoplot)) +
      geom_density(data=in_rivdt, aes(weight = LENGTH_KM),
                   fill='#35978f', alpha=0.5) +
      geom_density(fill='#bf812d', alpha=0.5) +
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
    "bio1_dc_uav", "bio7_dc_uav", "bio11_dc_uav",
    "bio12_mm_uav", "bio14_mm_uav", "clz_cl_cmj",
    "ari_ix_uav", "dis_m3_pyr", "sdis_ms_uyr", "UPLAND_SKM",
    "lka_pc_use", "snw_pc_uyr", "kar_pc_use", "for_pc_use") #, "glc_pc_u16")
  
  if ("dis_m3_pyr" %in% varstoplot_hist) {
    setDT(in_rivernetwork)[, dis_m3_pyr := dis_m3_pyr + 1]
    setDT(in_gaugepred)[, dis_m3_pyr := dis_m3_pyr + 1]
  }
  
  in_gaugepred[, IRpredcat_full := fifelse(predprob1>=0.5, 1, 0)]
  
  #Get legend
  pleg <- ggplot(in_gaugepred, aes(x=dis_m3_pyr, fill=factor(IRpredcat_full))) +
    geom_density(alpha=1/2) +
    scale_fill_manual(values=c('#35978f', '#bf812d'),
                      name = 'Dataset',
                      labels=c('Global non-perennial river reaches (predicted)',
                               'Gauging stations')) +
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
  grid::grid.newpage()
  do.call("grid.arrange", list(grobs=penvhist_grobs, nrow=5))
}

#


#------ reckless_KLdiv  -------------------------------------------------------
#Adapted from rags2ridges to make it faster by not checking for symmetry
reckless_KLdiv <- function(Mtest, Mref, Stest, Sref, symmetric = FALSE){
  if (!inherits(Mtest, "numeric")){
    stop("Input (Mtest) is of wrong class")
  }
  else if (!inherits(Mref, "numeric")){
    stop("Input (Mref) is of wrong class")
  }
  else if (length(Mtest) != length(Mref)){
    stop("Mtest and Mref should be of same length")
  }
  else if (!is.matrix(Stest)){
    stop("Input (Stest) is of wrong class")
  }
  else if (!is.matrix(Sref)){
    stop("Input (Sref) is of wrong class")
  }
  # else if (!isSymmetric(Stest)){ #isSymmetric takes a lot of time
  #   stop("Stest should be symmetric")
  # }
  # else if (!isSymmetric(Sref)){
  #   stop("Sref should be symmetric")
  # }
  else if (dim(Stest)[1] != length(Mtest)){
    stop("Column and row dimension of Stest should correspond to length Mtest")
  }
  else if (dim(Sref)[1] != length(Mref)){
    stop("Column and row dimension of Sref should correspond to length Mref")
  }
  else if (!inherits(symmetric, "logical")){
    stop("Input (symmetric) is of wrong class")
  }
  else {
    # Evaluate KL divergence
    KLd <- (sum(diag(solve(Stest) %*% Sref)) +
              t(Mtest - Mref) %*% solve(Stest) %*% (Mtest - Mref) -
              nrow(Sref) - log(det(Sref)) + log(det(Stest)))/2
    
    # Evaluate (original) symmetric version KL divergence
    if (symmetric){
      KLd <- KLd + (sum(diag(solve(Sref) %*% Stest)) +
                      t(Mref - Mtest) %*% solve(Sref) %*% (Mref - Mtest) -
                      nrow(Sref) - log(det(Stest)) + log(det(Sref)))/2
    }
    
    # Return
    return(as.numeric(KLd))
  }
}

#------ reckless_ridgeP  -------------------------------------------------------
#Adapted from rags2ridges to make it faster by not checking for symmetry
reckless_ridgeP <- function(S, lambda, type = "Alt", 
                            target = rags2ridges::default.target(S)){
  ##############################################################################
  # - NOTES:
  # - When type = "Alt" and target is p.d., one obtains the
  #   van Wieringen-Peeters type I estimator
  # - When type = "Alt" and target is null-matrix, one obtains the
  #   van Wieringen-Peeters type II estimator
  # - When target is not the null-matrix it is expected to be p.d. for the
  #   vWP type I estimator
  # - The target is always expected to be p.d. in case of the archetypal I
  #   estimator
  # - When type = "Alt" and target is null matrix or of form c * diag(p), a
  #   rotation equivariant estimator ensues. In these cases the expensive
  #   matrix square root can be circumvented
  ##############################################################################
  
  # if (!isSymmetric(S)) {
  #   stop("S should be a symmetric matrix")
  # }
  if (lambda <= 0) {
    stop("lambda should be positive")
  }
  else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
    stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
  }
  
  # Calculate Ridge estimator
  # Alternative estimator
  if (type == "Alt"){
    # if (!isSymmetric(target)) {
    #   stop("Shrinkage target should be symmetric")
    # } else
    if (dim(target)[1] != dim(S)[1]) {
      stop("S and target should be of the same dimension")
    } else {
      P <- rags2ridges:::.armaRidgeP(S, target, lambda)
    }
    dimnames(P) <- dimnames(S)
  }
  
  # Archetypal I
  if (type == "ArchI"){
    if (lambda > 1){
      stop("lambda should be in (0,1] for this type of Ridge estimator")
    } 
    # else if (!isSymmetric(target)){
    #   stop("Shrinkage target should be symmetric")
    else if (dim(target)[1] != dim(S)[1]){
      stop("S and target should be of the same dimension")
    } else if (any(eigen(target, symmetric = TRUE,
                         only.values = TRUE)$values <= 0)){
      stop("Target should always be p.d. for this type of ridge estimator")
    } else {
      P <- solve((1 - lambda) * S + lambda * solve(target))
    }
  }
  
  # Archetypal II
  if (type == "ArchII"){
    P <- solve(S + lambda * diag(nrow(S)))
  }
  
  # Set class and return
  attr(P, "lambda") <- lambda
  class(P) <- c("ridgeP", class(P))
  return(rags2ridges::symm(P))
}


#------ get_KLdiv --------------------------------------------------------------
#Adapted from rags2ridges]
get_KLdiv <- function(comp_mat, ref_mat=NULL, cov_ref=NULL, mean_ref=NULL) {
  #Compute multivariate covariance and mean of ref data
  if (is.null(cov_ref)) {
    cov_ref  <-  rags2ridges::covML(as.matrix(ref_mat))
  }
  
  if (is.null(mean_ref)) {
    mean_ref <- colMeans(as.matrix(ref_mat))
  }
  
  #Compute multivariate covariance and mean of comparison data
  cov_comp  <- rags2ridges::covML(comp_mat)
  mean_comp <- colMeans(comp_mat)
  
  ## Regularize singular Cov1
  P <- reckless_ridgeP(cov_comp, nrow(comp_mat))
  cov_compreg<- solve(P)
  
  ## Obtain KL divergence
  KLdiv_out <-  reckless_KLdiv(mean_comp, mean_ref, cov_compreg, cov_ref)
  
  return(list(
    cov_ref = cov_ref,
    mean_ref = mean_ref,
    KLdiv = KLdiv_out
  ))
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

plot_anthropo_stats <- function(in_gmeta_formatted, 
                                vlines = list(dor=2, crop=25, pop=100, built=1),
                                export = T,
                                fig_outdir = NULL) {
  #Degree of regulation plot ---------------------------------------------------
  dt_format_dor <- melt_anthropo_stats(in_dt=in_gmeta_formatted,
                                       fieldroot = 'dor') 
  
  p_dor <- ggplot(dt_format_dor, 
                  aes(x=value, y=year, fill = n, 
                      group=year, height = after_stat(count))) +
    geom_density_ridges(stat='density', scale=15, alpha=1/2,
                        rel_min_height = 0.0005)  +
    geom_vline(xintercept = vlines$dor)+
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
  
  p_crop <- ggplot(dt_format_crop, 
                   aes(x=value, y=year, fill = n, 
                       group=year, height = after_stat(count))) +
    geom_density_ridges(stat='density', scale=8, alpha=1/2,
                        rel_min_height = 0.001)  +
    geom_vline(xintercept = vlines$crop)+
    scale_x_continuous(
      name = 'Crop extent upstream (%)',
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
  
  p_pop <- ggplot(dt_format_pop, 
                  aes(x=value, y=year, fill = n, 
                      group=year, height = after_stat(count))) +
    geom_density_ridges(stat='density', scale=5, alpha=1/2,
                        rel_min_height = 0.001)  +
    geom_vline(xintercept = vlines$pop)+
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
  
  p_built <- ggplot(dt_format_built, 
                    aes(x=value, y=year, fill = n, 
                        group=year, height = after_stat(count))) +
    geom_density_ridges(stat='density', scale=10, alpha=1/2,
                        rel_min_height = 0.001)  +
    geom_vline(xintercept = vlines$built)+
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
  
  #Assemble and export patchwork -------------------------
  p_patchwork <- (p_crop + p_dor+ plot_layout(axes='collect'))/
    (p_pop + p_built + plot_layout(axes='collect'))
  
  if (export & !is.null(fig_outdir)) {
    ggsave(file.path(fig_outdir, paste0('anthropo_plot',
                                        format(Sys.Date(), '%Y%m%d'), '.png')),
           p_patchwork,
           width = 20, height = 20, units='cm'
    )
  }
  
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
  
  ngauges <- gauges_sub[n_years >= 15, .N, by=year]  
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
# in_grdc_no=4208549
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
  # merge(in_rivice[order(-river_ice_fraction) & !duplicated(date),
  #                 -'mean_elv',with=F],
  #       by=c('grdc_no', 'date'),
  #       all.x=T, all.y=F
  # )
  
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
  potential_npr <- q_dt_attri[, .SD[Qobs<=0.001,.N]>0]
  
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
# in_row <- data_for_qc[potential_npr==T & integer_perc<0.95,][grdc_no=='4208549',]
# in_data <- qread(in_row[, qs_path][[1]])
# in_nearg_cols <- in_row[, near_gcols_sel][[1]]
# arima_split=F
# plot_fit = F
# 
# detect_outliers_ts(in_data, in_nearg_cols, arima_split=F, plot_fit=F)

detect_outliers_ts <- function(in_data, in_nearg_cols, run_arima = T,
                               arima_split=F, plot_fit=F) {
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
  
  in_data[, date := as.Date(date)] %>%
    .[, jday := as.numeric(format(date, '%j'))]
  
  #Run ARIMA model to detect potential outliers---------------------------------
  #start <- Sys.time()
  if (run_arima) {
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
    
    #Find differences from forecast outside of the julian day 95% interval
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
  } else {
    in_data[, arima_outlier_rolldiff := NA]
  }
  
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
    arima_model = ifelse(run_arima, qARIMA$mod, NA),
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
# in_outliers_output_dt <- tar_read(q_outliers_flags)
# in_no <- "4208230"
# in_outliers_path<- in_outliers_output_dt[grdc_no==in_no, out_qs][[1]]
# n_freq=365.25
# in_floor=0

compute_metastatistics_util <- function(in_outliers_path, in_no) {
  #print(in_no)
  in_dt_edit_path <- paste0(
    tools::file_path_sans_ext(in_outliers_path),
    '_edit.csv')
  in_dt_clean <- fread(in_dt_edit_path, select=c('date', 'Qobs', 'year'))
  
  if (in_dt_clean[, sum(!is.na(Qobs) & (Qobs > 0))] > 1) {
    in_dt_clean <- fill_dt_dates(in_dt=in_dt_clean, 
                                 full_yrs=T, 
                                 date_col='date') 
    
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
        edited_data_path = in_dt_edit_path,
        Q50 = in_dt_clean[
          year %in% .SD[missingdays_edit_interp10 < 15,year],
          quantile(Qobs_interp, 0.5, na.rm=T)]
      )] 
    
    return(metastats_dt)
  }
}

#------ compute metastatistics wrapper -----------------------------------------
# in_outliers_output_dt <- tar_read(q_outliers_flags)
# unique(in_outliers_output_dt$grdc_no)

compute_metastatistics_wrapper <- function(in_outliers_output_dt,
                                           in_gaugep_dt) {
  out_dt <- unique(in_outliers_output_dt, by=c('grdc_no'))[,
                                                           compute_metastatistics_util(
                                                             in_outliers_path=out_qs,
                                                             in_no=grdc_no),
                                                           by=grdc_no
  ] %>% 
    merge(in_gaugep_dt[, list(y_geo=y_geo, grdc_no=as.character(grdc_no))],
          by='grdc_no', all.y=F)
  
  return(out_dt)
}

#------ analyze_metastats -----------------------------------------------
#in_metastats_dt <- tar_read(metastats_dt)

analyze_metastats <- function(in_metastats_dt,
                              export = T,
                              fig_outdir = figdir) {
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
    scale_x_continuous(name="Minimum number of years of valid records") +
    scale_y_continuous(name="Maximum number of missing daily records per year after interpolation") +
    scale_color_gradientn(colors=c("#440154FF","#404788FF", '#2D708EFF',
                                   '#20A387FF',  "#55C667FF")) + 
    facet_wrap(~max_interp,
               labeller = as_labeller(c(
                 `0` = 'Missing data interpolation: none',
                 `5` = 'Missing data interpolation: max. 5-day gaps',
                 `7` = 'Missing data interpolation: max. 7-day gaps',
                 `10`= 'Missing data interpolation: max. 10 days gaps')
               )
    )+
    theme_bw() +
    theme(legend.position='none',
          text=element_text(size=14))
  
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
  
  if (export & !is.null(fig_outdir)) {
    ggsave(file.path(fig_outdir, paste0('cleaned_ngauges_plot',
                                        format(Sys.Date(), '%Y%m%d'), '.png')),
           plot_ngauges,
           width = 20, height = 20, units='cm'
    )
  }
  
  return(list(
    plot_ngauges = plot_ngauges,
    metastats_meanyrs = metastats_meanyrs,
    metastats_nyears = metastats_nyears
  ))
}

#------ compute_noflow_hydrostats_util -------------------------------------
#Hydrostatistics from 
#Sauquet, E., Shanafield, M., Hammond, J., Sefton, C., Leigh, C., & Datry, T. (2021). 
#Classification and trends in intermittent river flow regimes in Australia, 
#northwestern Europe and USA: a global perspective. Journal of Hydrology, 126170. 
#https://doi.org/10.1016/j.jhydrol.2021.126170

compute_noflow_hydrostats_util <- function(in_dt, 
                                           in_lat,
                                           q_thresh = 0.001) {
  
  dt <- copy(in_dt)
  
  #Prepare data for computing statistics -----------------------------------------
  #Uniquely identify no-flow period
  dt[, noflow_period := rleid(Qobs_interp <= q_thresh)] %>%
    .[, noflow_periodyr := rleid(Qobs_interp <= q_thresh), by=year] %>%
    .[Qobs_interp > q_thresh, noflow_period := NA] %>%
    .[!is.na(noflow_period), noflow_period_dur := .N, by = noflow_period] %>%
    .[!is.na(noflow_period), noflow_periodyr_dur := .N, by = .(noflow_periodyr, year)]
  
  #Identify continuous blocks of 5 years (for d80)
  whole_5yrblock <- dt[order(date) & !duplicated(year),] %>%
    .[, year_rleid :=  cumsum(c(1, diff(year) != 1))] %>%
    .[, year_5blockid := ceiling(seq_along(year)/5), by='year_rleid'] %>%
    .[, year_5blockidn := .N, by=c('year_rleid', 'year_5blockid')] %>%
    .[year_5blockidn==5, blockid := 10*year_rleid+year_5blockid]
  
  #Compute also partial 5-year blocks if there is no whole 5-year block (for d80)
  partial_5yrblock <- dt[order(date) & !duplicated(year),]%>%
    .[, year_5blockid_partial := ceiling(seq_along(year)/5)] %>%
    .[, year_5blockidn_partial := .N, by=c('year_5blockid_partial')] %>%
    .[year_5blockidn_partial==5, blockid := year_5blockid_partial]
  
  if (whole_5yrblock[!is.na(blockid), .N]>0) {
    dt <- merge(dt, whole_5yrblock[, .(year, blockid)], by='year', all.x=T)
  } else {
    dt <- merge(dt, partial_5yrblock[, .(year, blockid)], by='year', all.x=T)
  }
  
  #Convert each day i with no flow into an angular (ti) and represent it
  #by a unit vector with rectangular coordinates (cos(ti); sin(ti))
  dt[!is.na(noflow_period), `:=`(cos_t=cos(2*pi*(jday-1)/(diny(year)-1)),
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
    driest_6mocenter <- in_dt[get(jday_col) < 366 #can lead to artefacts if drying was infrequent but occurred during that time
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
    
    return(in_dt[order(date),])
  }
  dt <- identify_drywet6mo(in_dt=dt)
  
  #Compute seasonal recession
  compute_seasonal_recession <- function(in_dt, 
                                         q_col='Qobs_interp', date_col='date',
                                         jday_col = 'jday') {
    #Make sure that dates for missing data are included so that the rolling window
    #doesn't include non-contiguous periods
    in_dt <- fill_dt_dates(in_dt, full_yrs=TRUE)
    
    #Compute long-term 30-day moving average Q by julian day
    #and re=organize this statistic by water year such that the first day of the 
    #water year is the day is the maximum long-term daily discharge
    #the period at the beginning of the year is thus the falling limb
    falling_limb <- in_dt[, mean30dQ := frollmean(Qobs_interp, 31, align='center')] %>%
      .[, list(mean30jdQ =  mean(mean30dQ, na.rm=T)), by='jday'] %>%
      .[c(which(jday >= .[which.max(mean30jdQ), jday]), 
          which(jday < .[which.max(mean30jdQ), jday])),] %>%
      .[, ix := .I]
    
    #Identify the seasonal recession comprised between the 90th quantile and median
    #whose duration is Drec
    seasonal_recession <- falling_limb[
      !seq(1, min(which(mean30jdQ <= quantile(mean30jdQ, 0.9)))-1),] %>%
      .[!seq(min(which(mean30jdQ <= median(falling_limb$mean30jdQ))), .N),]
    
    ggplot(falling_limb, aes(x=ix, y=mean30jdQ)) +
      geom_point() +
      geom_point(data=seasonal_recession,
                 color='red')
    
    return(seasonal_recession)
  }
  sea_rec <- compute_seasonal_recession(in_dt=dt)
  
  
  #---- Median duration of runoff events -----------------------------------------
  #Compute baseflow
  compute_baseflow_gustard <- function(in_dt, q_col='Qobs_interp', 
                                       date_col = 'date', plot=F) {
    #Make sure that dates for missing data are included so that the rolling window
    #doesn't include non-contiguous periods
    in_dt <- fill_dt_dates(in_dt, full_yrs=TRUE)
    
    #Divide record in 5-day blocks
    in_dt[order(date), block5 := c(rep(seq(1,floor(.N/5)), each=5), 
                                   rep(NA, .N %% 5))]
    #Get minimum discharge in each block
    in_dt[, blockminQ := min(get(q_col)), by=block5]
    
    #Identify the day of minimum discharge in each block. When multiple days
    #have the same minimum discharge, take the median day of minimum value within the block
    median_blockmindate_ix <- in_dt[!is.na(block5),
                                    .SD[, which(get(q_col) == blockminQ)],
                                    by=block5]  %>%
      .[, list(blockminQ_medix = floor(median(V1))), by=block5]
    
    in_dt <- merge(in_dt, median_blockmindate_ix, by='block5') %>%
      .[, blockminQdate := date==.SD[blockminQ_medix, date], by=block5]
    
    #Identify turning points. Those days of minimum flow within each 5-day block
    #whose discharge*0.9 is equal or less than the minimum flow in the neighboring
    #5-day blocks
    turning_points <- in_dt[, blockminQ[[1]], by=block5] %>%
      .[, list(
        block5 = block5,
        blockminQdate = TRUE,
        turning_ptQ = fifelse(
          ((0.9*V1) <= data.table::shift(V1, n=1)) 
          & ((0.9*V1) <= data.table::shift(V1, n=-1)), V1, -999)
      )] %>%
      .[turning_ptQ == -999, turning_ptQ := NA]
    
    in_dt <- merge(in_dt, turning_points, 
                   by=c('block5', 'blockminQdate'), all.x=T)
    
    #Interpolate between turning points
    in_dt[, bfQ := as.vector(
      forecast::na.interp(turning_ptQ, lambda=NULL, linear=T))] %>%
      .[bfQ > get(q_col), bfQ := get(q_col)]
    
    bfp <- ggplot(in_dt, aes(x=date, y=get(q_col), group=format(date, "%Y"))) +
      geom_line() + 
      geom_line(aes(y=bfQ), color='red') + 
      scale_y_sqrt() +
      theme_bw() 
    
    if (plot) {
      print(bfp)
    }
    
    return(in_dt)
  }
  
  compute_medianDr <- function(in_dt, q_col='Qobs_interp', bf_col='bfQ',
                               date_col = 'date') {
    in_dt[, runoffQ := get(q_col) - get(bf_col)]
    in_dt[, maxrunoff := max(runoffQ), by=year]
    in_dt[, maxrunoff_date := .SD[min(which(runoffQ==maxrunoff)), get(date_col)] ,
          by=year]
    in_dt[, maxrunoff_date_recessionyr := maxrunoff_date + 365]
    
    runoff_enddates_dt <- lapply(unique(in_dt$maxrunoff_date), function(d) {
      #print(d)
      data.table(
        maxrunoff_date = d,
        runoffevent_enddate = in_dt[date %in% seq(d+1, d + 365, by='day'),
                                    .SD[min(which(runoffQ<=(maxrunoff/2))), date]
        ]
      )
    }) %>% 
      rbindlist %>%
      .[!is.na(runoffevent_enddate),] #If end of record or missing data the following year.
    
    maxrunoff_events_dt <- merge(in_dt,  runoff_enddates_dt, 
                                 by.x='date', by.y='maxrunoff_date', all.x=F)
    
    r_event_p <- ggplot(in_dt, aes(x=date, y=runoffQ, group=year)) + 
      geom_line(color='darkgrey') +
      geom_point(data=maxrunoff_events_dt, color='red') +
      geom_point(data=in_dt[date %in% maxrunoff_events_dt$runoffevent_enddate,], 
                 color='blue') +
      theme_bw() 
    #ggplotly(r_event_p)  
    
    medianDr <- median(maxrunoff_events_dt[, runoffevent_enddate-maxrunoff_date, 
                                           by=year]$V1)
    
    return(medianDr)
  }
  
  #Compute statistics ------------------------------------------------------------
  q_stats <- dt[, list(
    #Intermittence
    f0 = sum(Qobs_interp <= q_thresh)/.N,
    
    #Duration
    meanD = mean(
      .SD[!is.na(noflow_period) & !duplicated(noflow_period), noflow_period_dur]),
    medianD = median(
      .SD[!is.na(noflow_period) & !duplicated(noflow_period), noflow_period_dur]),
    sdD = sd(
      .SD[!is.na(noflow_period) & !duplicated(noflow_period), noflow_period_dur]),
    
    #d80 is computed, excluding all years not within continuous blocks of 5 years 
    #in chronological order (i.e., no re-ordering within continuous periods of 
    #record to minimize or maximize d80)
    d80 = .SD[!is.na(blockid), 
              fifelse(sum(!is.na(noflow_periodyr_dur)) > 0,
                      max(noflow_periodyr_dur, na.rm=T), 0)
              , by=blockid][, ceiling(quantile(V1, 0.8))],
    
    #Frequency
    meanN = mean(
      .SD[, uniqueN(noflow_period, na.rm=T), by=year]$V1),
    medianN = median(
      .SD[, uniqueN(noflow_period, na.rm=T), by=year]$V1),
    sdN = sd(
      .SD[, uniqueN(noflow_period, na.rm=T), by=year]$V1),
    
    #Timing
    theta = atan2(mean(sin_t, na.rm=T),
                  mean(cos_t, na.rm=T) 
    ),
    r = sqrt(mean(cos_t, na.rm=T)^2 
             + mean(sin_t, na.rm=T)^2),
    
    #Rate of change
    Ic = ((quantile(Qobs_interp, 0.9)- quantile(Qobs_interp, 0.01))/ #concavity index
            (quantile(Qobs_interp, 0.99)- quantile(Qobs_interp, 0.01))),
    
    #Proportion of zero-flow days with max monthly temperature under 0
    Fper = .SD[tmax<0, sum(Qobs_interp <= q_thresh)]/
      sum(Qobs_interp <= q_thresh),
    
    FperM10 = .SD[tmax<(-10), sum(Qobs_interp <= q_thresh)]/
      sum(Qobs_interp <= q_thresh),
    
    #50th and 90th quantile of PDSI during no-flow
    PDSIdiff = .SD[is.na(noflow_period), quantile(PDSI, 0.5, na.rm=T)]
    -.SD[!is.na(noflow_period), quantile(PDSI, 0.5, na.rm=T)],
    
    P90PDSI = .SD[!is.na(noflow_period), quantile(PDSI, 0.9, na.rm=T)]
  )]
  
  #Adjust theta for stations in southern hemisphere. To avoid discontinuities, 
  #shift the seasons progressively up to 23.5 degrees. 
  #After 23.5, shift by a full half year 
  #i.e. Meteorological summer starts December 1st instead of June 1st
  q_stats[, theta := fifelse(theta < 0,
                             theta + 2*pi,
                             theta)] %>%
    .[, theta :=  fifelse(in_lat >= 0,
                          theta,
                          (theta - pi*max(-1, (in_lat/23.5)))%%(2*pi)
    )]
  
  #Compute r x cos(theta) and r x sin(theta) for the actual classification (see Sauquet et al. 2021)
  q_stats[, `:=`(
    r_cos_theta = r*cos(theta),
    r_sin_theta = r*sin(theta)
  )]
  
  #Compute seasonal predictability of no-flow events (Sd6)
  q_stats$Sd6 <- dt[, ym := format(date, '%Y%m')] %>%
    .[, any(!is.na(noflow_period)), by=.(year, ym, dry_6mo)] %>%
    .[, sum(V1, na.rm=T), by=.(year, dry_6mo)] %>% #Compute number of months with zero-flows for each wet and dry period and year
    .[, 1-(.SD[!dry_6mo,mean(V1)]/.SD[dry_6mo,mean(V1)])]
  
  # Rate of change	Drec	Seasonal recession time scale (Catalogne, 2012) (day)
  q_stats[, Drec := nrow(sea_rec)]
  
  # BFI	Baseflow index computed with the smoothed minima method introduced by the 
  #Institute of Hydrology (1980) (dimensionless)
  #See "Gustard, A., & Demuth, S. (2008). Manual on low-flow estimation and 
  #prediction. Operational Hydrology Report No. 50 (World Meteorological Organization (WMO), Ed.). Opera."s
  q_stats$bfi <- compute_baseflow_gustard(dt)[, mean(bfQ/Qobs_interp, na.rm=T)]  
  
  # medianDr	Median duration of runoff event (*) (day)
  q_stats$medianDr <- compute_baseflow_gustard(dt) %>%
    compute_medianDr
  
  return(q_stats)
}


#------ compute_noflow_hydrostats_wrapper -------------------------
# in_metastats_dt <- tar_read(metastats_dt)
# in_metastats_analyzed <- tar_read(metastats_analyzed)
# in_gaugep_dt <- tar_read(gaugep_dt)
# max_interp_sel = 5
# max_miss_sel = 0
# min_nyears = 15
# q_thresh = 0.001

compute_noflow_hydrostats_wrapper <- function(in_metastats_analyzed,
                                              in_metastats_dt,
                                              in_gaugep_dt,
                                              max_interp_sel = 10,
                                              max_miss_sel = 0,
                                              min_nyears = 15,
                                              q_thresh = 0.001) {
  interp_fn <- paste0('missingdays_edit',
                      fifelse(max_interp_sel >0,
                              paste0('_interp', max_interp_sel),
                              '')
  )
  
  gauges_sel_no <- in_metastats_analyzed$metastats_nyears[
    max_interp==max_interp_sel & max_miss==max_miss_sel & nyears>=min_nyears,
    grdc_no]
  
  
  qstats_dt <- future_lapply(gauges_sel_no, function(in_no) {
    #print(in_no)
    #For given gauge, get years with less than the maximum number of missing days
    #given maximum interpolation period
    meta_sub_yrs <- in_metastats_dt[grdc_no==in_no
                                    & get(interp_fn)<=max_miss_sel,]
    
    g_lat <- in_metastats_dt[grdc_no==in_no, unique(y_geo)]
    
    #Read data, keeping only years with number of missing days under threshold
    q_dt <- fread(meta_sub_yrs$edited_data_path[[1]]) %>%
      .[year %in% meta_sub_yrs$year, 
        .(grdc_no, date, jday, Qobs, year, tmax, PDSI)] %>%#subset columns for speed
      .[min(which(!is.na(Qobs))):max(which(!is.na(Qobs))),] %>% #remove leading and trailing NAs
      .[!duplicated(date),] %>%
      fill_dt_dates(full_yrs=T, date_col='date')
    
    #
    q_dt[is.na(Qobs), `:=`(grdc_no = as.integer(in_no),
                           jday = as.integer(format(date, '%j'))
    )] %>%
      .[is.na(Qobs) & year>1960, `:=`(tmax = nafill(tmax, type='locf'),
                                      PDSI = nafill(PDSI, type='locf')
      )]
    
    #Interpolate discharge, keeping only interpolation periods under threshold
    na_interp_dt_custom(in_dt=q_dt,
                        in_var='Qobs')
    #Fill NA periods < 10 days on the edges of the contiguous periods of record
    q_dt[, Qobs_interp := nafill(Qobs_interp, type='locf')]
    q_dt[, Qobs_interp := nafill(Qobs_interp, type='nocb')]
    
    q_dt[NAperiod_n > max_interp_sel, Qobs_interp := NA]
    
    
    out_stats <- compute_noflow_hydrostats_util(in_dt = q_dt, 
                                                in_lat = g_lat,
                                                q_thresh = 0.001)
    out_stats$grdc_no <- in_no
    
    return(out_stats)
  }) %>% rbindlist(fill=T)
  
  return(qstats_dt)
}


#------ preformat_hydrostats ---------------------------------------------------
#in_hydrostats <- tar_read(noflow_hydrostats)

preformat_hydrostats <- function(in_hydrostats) {
  #Set statistics order and weight
  hydrostats_order <- data.table(
    variable = c("f0", 
                 "meanN", "medianN","sdN",   
                 "meanD", "medianD", "sdD", "d80",
                 "Drec", "Ic", "bfi", "medianDr",
                 "Fper","FperM10", "PDSIdiff", "P90PDSI",
                 "theta", "r", "Sd6", 'r_cos_theta', 'r_sin_theta'),
    aspect = c('Intermittence',
               rep('Frequency', 3),
               rep('Duration', 4),
               rep('Rate of change', 4),
               rep('Climate dependence', 4),
               rep('Timing', 5)),
    weight = c(1,
               rep(1/3, 3),
               rep(1/4, 4),
               rep(1/4, 4),
               rep(1/4, 4),
               rep(1/3, 3)
    )
  ) %>% 
    .[, var_order := .I] %>%
    .[, `:=`(variable = factor(variable, levels=variable),
             aspect = factor(aspect, levels=unique(aspect)))]
  
  #Pre-format statistics -------------------------------------------------------
  id_dt <- in_hydrostats[f0>0, .(grdc_no, .I)]
  hydrostats_raw <- in_hydrostats[f0>0, -'grdc_no', with=F] %>%
    .[, I := .I]
  
  #Assign 0 sD when only one drying event
  hydrostats_raw[is.na(sdD), sdD := 0]
  
  #Check distribution and melt
  hydrostats_distrib_p <- ggplot(melt(hydrostats_raw), aes(x=value)) +
    geom_density() +
    facet_wrap(~variable, scales='free')
  
  hydrostats_meltattri <- melt(hydrostats_raw, id.vars = 'I') %>%
    merge(hydrostats_order, ., by='variable', all.x=T) 
  
  
  hydrostats_trans <- transform_scale_vars(in_dt=hydrostats_meltattri, 
                                           value_col='value', var_col='variable',
                                           inplace=FALSE)$dt
  
  hydrostats_distrib_scaled_p <- ggplot(hydrostats_trans, aes(x=value_scaled)) +
    geom_density() +
    facet_wrap(~variable, scales='free')
  
  #Cast to matrix
  hydrostats_scaled_cast <- dcast(
    data = hydrostats_trans,
    formula = I~factor(variable, levels=hydrostats_order$variable),
    value.var = 'value_scaled')
  
  #Impute PDSIq50 and q90 with a median
  hydrostats_scaled_cast[
    is.na(PDSIdiff),
    `:=`(PDSIdiff = median(hydrostats_scaled_cast$PDSIdiff, na.rm=T),
         P90PDSI = median(hydrostats_scaled_cast$P90PDSI, na.rm=T)
    )]
  
  hydrostats_mat_sub <- as.matrix(hydrostats_scaled_cast[,-c('I'), with=F])
  row.names(hydrostats_mat_sub) <- id_dt[order(I), grdc_no]
  
  # The variables r × sin(θ) and r × cos(θ) were used instead of r and θ to avoid 
  # an artificial break in winter induced by angles.
  
  return(list(
    dt = merge(hydrostats_trans, id_dt, by='I'),
    mat = hydrostats_mat_sub)
  )
}

#------ cluster_gauges --------------------------------------------
# in_hydrostats_preformatted <- tar_read(noflow_hydrostats_preformatted)
# in_colors = class_colors
# stats_sel = c('f0', 'medianN', 'sdN',
#               'medianD', 'sdD',
#               'Ic', 'bfi', 'medianDr',
#               'Fper', 'FperM10',
#               'PDSIdiff', 'P90PDSI',
#               'r_cos_theta', 'r_sin_theta')

cluster_noflow_gauges_full <- function(in_hydrostats_preformatted,
                                       stats_sel, in_colors) {
  hydrostats_dt <- in_hydrostats_preformatted$dt
  hydrostats_order <- in_hydrostats_preformatted$dt[
    !duplicated(variable),
    .(variable, aspect, weight, var_order)]
  
  hydrostats_mat <- in_hydrostats_preformatted$mat
  hydrostats_mat_sub <- hydrostats_mat[
    , which(colnames(hydrostats_mat) %in% stats_sel)]
  
  #Correlation among variables  ------------------------------------------------
  var_cor <- cor(hydrostats_mat, method='spearman', use="pairwise.complete.obs")
  
  p_varscor <- ggcorrplot(var_cor, #method = "circle", 
                          hc.order = FALSE, hc.method = 'average',
                          type = "upper", lab=T, lab_size =3,
                          digits=2, insig='blank',
                          outline.color = "white") +
    scale_fill_distiller(
      name=str_wrap("Correlation coefficient Spearman's rho", 20),
      palette='RdBu', 
      limits=c(-1, 1), 
      breaks=c(-1, -0.5, 0, 0.5, 1)) +
    theme(legend.position = c(0.8, 0.3))
  
  #Compute  Euclidean distance based on correlation coefficients and variable weights
  hydro_dist <- cluster::daisy(
    hydrostats_mat_sub, 
    metric = "euclidean") %>%
    as.dist
  #weights = hydrostats_order$weight #not worth bothering with weighting given 
  #that weights only range from 0.25 to 0.33
  
  #
  #Cluster stations based on UPGMA or Ward's---------------------------------
  algo_list <- list('average', 'median', 'ward.D2')
  hclust_reslist <- lapply(
    algo_list, function(in_method) {
      rundiagnose_clustering(
        in_mat = hydrostats_mat_sub,
        in_dist = hydro_dist,
        in_method = in_method,
        min_nc = 5,
        max_nc = 15)
    })
  names(hclust_reslist) <- algo_list
  
  # hclust_reslist$average$cophcor
  # hclust_reslist$median$cophcor
  # hclust_reslist$ward.D2$cophcor
  # 
  # hclust_reslist$average$p_scree
  # hclust_reslist$average$nbclust_tests$Best.nc
  # hclust_reslist$ward.D2$p_scree
  # hclust_reslist$ward.D2$nbclust_tests$Best.nc
  
  
  #Define class colors------------------------------------------------------------
  #in_colors <- c("#176c93","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#7a5614","#6baed6","#00441b", '#e41a1c') #9 classes with darker color (base blue-green from Colorbrewer2 not distinguishable on printed report and ppt)
  # classcol <- c('#999900', '#728400', '#008F6B', '#005E7F', '#4A4A4A', '#A1475D',
  #               '#756200', '#C96234', '#4782B5', "#00441b",'#984ea3', '#e41a1c',
  #               '#4daf4a','#ff7f00')
  # classcol_temporal <- c(,'#377eb8',,'#666666','#a65628')
  
  
  #Make table of gauge classes and good looking dendogram-----------------------
  mat_forbp <- hydrostats_mat[
    , -which(colnames(hydrostats_mat) %in% c('r_cos_theta', 'r_sin_theta'))]
  
  nclass_list <- c(6, 9)
  cluster_analyses_avg <- lapply(nclass_list, function(kclass) {
    cuttree_and_visualize(
      in_mat =  mat_forbp,
      in_hclust = hclust_reslist$average$hclust,
      in_kclass = kclass,
      in_colors = in_colors, 
      in_meltdt = hydrostats_dt,
      id_col = 'grdc_no')
  })
  names(cluster_analyses_avg) <- paste0('ncl', nclass_list)
  
  nclass_list <- c(6, 7,8, 9)
  cluster_analyses_ward2 <- lapply(nclass_list, function(kclass) {
    cuttree_and_visualize(
      in_mat =  mat_forbp,
      in_hclust = hclust_reslist$ward.D2$hclust,
      in_kclass = kclass,
      in_colors = in_colors, 
      in_meltdt = hydrostats_dt,
      id_col = 'grdc_no')
  })
  names(cluster_analyses_ward2) <- paste0('ncl', nclass_list)
  
  # cluster_analyses_avg$ncl6$class_dt[, classn[[1]], by=gclass]
  # cluster_analyses_avg$ncl6$p_boxplot
  # cluster_analyses_avg$ncl9$class_dt[, classn[[1]], by=gclass]
  # cluster_analyses_avg$ncl9$p_boxplot
  # 
  # cluster_analyses_ward2$ncl6$class_dt[, classn[[1]], by=gclass]
  # cluster_analyses_ward2$ncl6$p_boxplot
  # cluster_analyses_ward2$ncl8$class_dt[, classn[[1]], by=gclass]
  # cluster_analyses_ward2$ncl8$p_boxplot
  
  #Return objects --------------------------------------------------------------
  return(list(
    hclust_reslist_all = hclust_reslist,
    cluster_analyses = cluster_analyses_ward2,
    hydro_dist = hydro_dist,
    chosen_hclust = 'ward.D2',
    kclass = 9,
    p_varscor = p_varscor
  ))
}

#------ Export gauge classes ---------------------------------------------------
# in_cluster_postanalysis <- tar_read(sel_cluster_postanalysis)
# in_path_gaugep = tar_read(path_gaugep)
# out_shp_root=file.path(resdir, 'gaugep_classstats_ward')
export_gauges_classes <- function(in_cluster_postanalysis, 
                                  in_path_gaugep,
                                  out_shp_root) {
  
  out_shp <- paste0(out_shp_root,
                    in_cluster_postanalysis$class_dt[, length(unique(gclass))],
                    '.shp')
  gaugep <- terra::vect(dirname(in_path_gaugep), layer=basename(in_path_gaugep))

  class_dt <- in_cluster_postanalysis$class_dt
  class_cast <- data.table::dcast(class_dt, grdc_no+gclass~variable, value.var = 'value')
  
  gaugep_classstats <- terra::merge(gaugep, class_cast, by='grdc_no', all.x=F)
  
  terra::writeVector(gaugep_classstats, out_shp, overwrite=T)
  return(out_shp)
}


#------ Plot hydrograph --------------------------------------------------
# in_dt <- q_dt_bind
# value_col <- 'Qobs'
# date_col <- 'date'
# lat_col <- 'y_geo'
# in_color <- '#4d004b'
# smoothing_window = 5

plot_hydrograph <- function(in_dt, value_col, date_col, lat_col, back_col,
                            q_thresh, in_color, 
                            noflow_qthresh=0.001, smoothing_window = 7) {
  
  in_dt[, yrmean := mean(get(value_col), na.rm=T),
        by=format(get(date_col), '%Y')]
  
  if (!('jday' %in% names(in_dt))) {
    in_dt[, jday := format(get(date_col), '%j')]
  }
  
  if (!is.null(lat_col)) {
    in_dt[, cal_doy := ifelse(
      y_geo > 0, 
      format(as.Date(as.numeric(jday)-1, 
                     origin=paste0(year-1, '-01-01')),
             "%m-%d"),
      format(as.Date(as.numeric(jday)-1,
                     origin=paste0(year-1, '-07-01')),
             "%m-%d"))]
  } else {
    in_dt[, cal_doy := format(as.Date(as.numeric(jday)-1, 
                                      origin=paste0(year-1, '-01-01')),
                              "%m-%d")]
  }
  
  #Compute statistics on long-term daily flow (average, min, max, Q10, Q25, Q75, Q90) across all stations and years for each class
  padding_ix <- c(rev(seq(365,1)[seq(smoothing_window%/%2)]), 
                  1:365, 
                  seq(smoothing_window%/%2))
  
  classflowstats <- in_dt[,list(
    classmeanfull = mean(get(value_col), na.rm=T), 
    classmean = mean(get(value_col)/yrmean, na.rm=T),
    classQ75 = quantile(get(value_col)/yrmean, .25, na.rm=T),
    classQ25 = quantile(get(value_col)/yrmean, .75, na.rm=T),
    classQ50 = quantile(get(value_col)/yrmean, .5, na.rm=T),
    classQ90 = quantile(get(value_col)/yrmean, .10, na.rm=T),
    classQ10 = quantile(get(value_col)/yrmean, .90, na.rm=T),
    classmax = max(get(value_col)/yrmean, na.rm=T),
    classmin = min(get(value_col)/yrmean, na.rm=T),
    classsd = sd(get(value_col)/yrmean, na.rm=T)
  ), by=cal_doy] %>%
    .[padding_ix,] %>%
    .[, (names(.)[-1]) := sapply(
      .SD, function(x) frollmean(x, n=smoothing_window, na.rm=T, align='center'), 
      simplify=F),
      .SDcols = names(.)[-1]] %>%
    .[!is.na(classmean),]
  
  #Compute daily no-flow probability by calendar day and scale it to plot in the
  #same panel as the flow stats
  classnoflowstats <- in_dt[!is.na(get(value_col)),list(
    prob_noflow=.SD[get(value_col)<noflow_qthresh,.N]/.N
  ), by=cal_doy] %>%
    .[padding_ix,] %>%
    .[, prob_noflow_scaled := frollmean(
      prob_noflow*max(classflowstats$classQ10, na.rm=T),
      n=smoothing_window, na.rm=T, align='center')]%>%
    .[!is.na(prob_noflow_scaled ),]
  
  
  unit_scale_axis <- function(y) {
    (y - min(y)) / (max(y)-min(y))
  }
  
  #Facetted standardized average yearly hydrograph for each class + Q90-Q10 ribbon
  classhydro_facet <- ggplot(as.data.frame(classflowstats), 
                             aes(x=as.Date(cal_doy, format="%m-%d")+181)) + 
    #geom_ribbon(aes(ymin=ifelse(classmean-2*classsd>=0,classmean-2*classsd,0), ymax=classmean+2*classsd,
    #                fill=factor(classnames)),alpha=0.3) +
    geom_ribbon(aes(ymin=classQ90, ymax=classQ10), fill=in_color, alpha=0.3) +
    geom_ribbon(aes(ymin=classQ75, ymax=classQ25), fill=in_color, alpha=0.5) +
    geom_line(aes(y=classQ50), color=in_color, linewidth=1) + 
    geom_line(data=classnoflowstats, aes(y=prob_noflow_scaled), 
              color='black', linewidth=1.2, alpha=0.5) + 
    # scale_color_manual(name='Hydrologic class',values=classcol[1:kclass]) +
    # scale_fill_manual(name='Hydrologic class',values=classcol[1:kclass]) +
    scale_y_sqrt(
      name=expression(frac('Daily discharge', 'Mean annual discharge')),
      expand=c(0,0), limits=c(0,NA),
      sec.axis = sec_axis(
        name = 'No-flow probability',
        trans = ~ unit_scale_axis(.),
        breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5 , 0.75, 1))
    ) + 
    scale_x_date(name=expression(Date~frac(North, South)), 
                 date_breaks = "3 month", date_labels='%b',
                 sec.axis = sec_axis(trans = ~ . - 181),
                 expand=c(0,0)) + 
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank(),
          text = element_text(size=12),
          axis.text.y.left = element_text(color=in_color))
  
  return(classhydro_facet)
}
#------ Plot class hydrograph --------------------------------------------------
# in_class_dt <- class_dt_sub
# #in_no <- '3649455'
# max_interp_sel = 5
# max_miss_sel = 0

plot_class_hydrograph <- function(in_class_dt, 
                                  in_metastats_dt,
                                  max_interp_sel,
                                  max_miss_sel,
                                  in_color,
                                  noflow_qthresh=0.001) {
  
  interp_fn <- paste0('missingdays_edit',
                      fifelse(max_interp_sel >0,
                              paste0('_interp', max_interp_sel),
                              '')
  )
  
  q_dt_bind <- lapply(in_class_dt$grdc_no, function(in_no) {
    print(in_no)
    #For given gauge, get years with less than the maximum number of missing days
    #given maximum interpolation period
    meta_sub_yrs <- in_metastats_dt[grdc_no==in_no
                                    & get(interp_fn)<=max_miss_sel,]
    
    #Read data, keeping only years with number of missing days under threshold
    q_dt <- fread(meta_sub_yrs$edited_data_path[[1]],
                  select = c('grdc_no', 'date', 'Qobs', 'tmax', 'year')) %>% #subset columns for speed
      .[year %in% meta_sub_yrs$year,] %>%
      .[min(which(!is.na(Qobs))):max(which(!is.na(Qobs))),] %>% #remove leading and trailing NAs
      .[!duplicated(date),] %>%
      .[, grdc_no := as.character(grdc_no)] 
    
    return(q_dt)
  }) %>% 
    rbindlist %>%
    merge(in_metastats_dt[!duplicated(grdc_no), .(grdc_no, y_geo)],
          by='grdc_no', all.y=F)
  
  out_facet <- plot_hydrograph(
    in_dt = q_dt_bind,
    value_col = 'Qobs',
    date_col = 'date',
    lat_col = 'y_geo',
    in_color = in_color,
    noflow_qthresh = noflow_qthresh, 
    smoothing_window = 7
  )
  
  return(out_facet)
}

#------ Plot class hydrograph wrapper --------------------------------------------------
# in_cluster_postanalysis <- tar_read(sel_cluster_postanalysis)
# in_metastats_dt <- tar_read(metastats_dt)
# max_interp_sel = 5
# max_miss_sel = 0
# noflow_qthresh=0.001

plot_class_hydrograph_wrapper <- function(in_cluster_postanalysis,
                                          in_metastats_dt,
                                          max_interp_sel,
                                          max_miss_sel,
                                          noflow_qthresh=0.001
) {
  
  class_dt <- in_cluster_postanalysis$class_dt
  kclass <- class_dt[, length(unique(gclass))]

  k_colors <- in_cluster_postanalysis$p_dendo$layers[[4]]$data %>%
    as.data.table %>%
    .[order(label), col]
  
  hydrograph_facets <- lapply(seq(kclass), function(in_class) {
    print(in_class)
    
    class_dt_sub <- class_dt[gclass == in_class,] %>%
      .[!duplicated(grdc_no), .(grdc_no, gclass, classn)]
    
    out_facet <- plot_class_hydrograph(in_class_dt = class_dt_sub,
                                       in_metastats_dt = in_metastats_dt,
                                       max_interp_sel = max_interp_sel,
                                       max_miss_sel = max_miss_sel,
                                       in_color = k_colors[in_class],
                                       noflow_qthresh=noflow_qthresh)
    return(out_facet)
  })
  
  out_plotgrid <- wrap_plots(hydrograph_facets, 
                             ncol=floor(sqrt(length(hydrograph_facets))),
                             axes='collect') +
    plot_layout(axis_titles = "collect") + 
    plot_annotation(tag_levels = '1')  & 
    theme(plot.tag = element_text(size = 10))
  
  return(out_plotgrid)
}

#------ Analyze cluster sensitivity --------------------------------------------
# in_noflow_clusters <- tar_read(noflow_clusters)
# in_cluster_postanalysis <- tar_read(sel_cluster_postanalysis)
# in_hydrostats_preformatted <- tar_read(noflow_hydrostats_preformatted)
# classnames <- manual_class_order
# stats_sel = hydrostats_sel
analyze_cluster_sensitivity <- function(in_noflow_clusters,  
                                        in_cluster_postanalysis,
                                        in_hydrostats_preformatted,
                                        stats_sel,
                                        classnames,
                                        export = T,
                                        fig_outdir = NULL) {
  #Permute each variable
  class_dt <- in_cluster_postanalysis$class_dt
  kclass <- class_dt[, length(unique(gclass))]
  
  base_clust_reslist <- in_noflow_clusters$hclust_reslist_all[[
    in_noflow_clusters$chosen_hclust]]
  
  base_hclust_cut <- base_clust_reslist$hclust %>%
    dendextend::cutree(k=kclass, order_clusters_as_data = FALSE)
  
  hydrostats_mat_sub <- in_hydrostats_preformatted$mat %>%
    .[, which(colnames(.) %in% stats_sel)]
  
  cols_ix <- seq(dim(hydrostats_mat_sub)[2])
  rows_ix <- seq(dim(hydrostats_mat_sub)[1])
  
  #---------------- Compute variable importance in clustering ------------------
  get_col_ARI <- function(in_mat, col_ix, rows_ix, 
                          in_method, in_kclass, comp_clust_cut) {
    mat_permut <- in_mat
    mat_permut[,col_ix] <- mat_permut[
      sample(rows_ix, size=length(rows_ix), replace=FALSE),col_ix] 
    
    hydro_dist <- cluster::daisy(
      mat_permut, 
      metric = "euclidean") %>%
      as.dist
    
    new_clust <- fastcluster::hclust(hydro_dist, method=in_method) %>%
      dendextend::cutree(k=in_kclass, order_clusters_as_data = FALSE)
    
    out_ari <- aricode::ARI(comp_clust_cut, new_clust)
  }
  
  hydrostats_ari <- future_lapply(cols_ix, function(col_to_permute) {
    ari_list <- replicate(500, get_col_ARI(in_mat=hydrostats_mat_sub,
                                           col_ix=col_to_permute,
                                           rows_ix=rows_ix,
                                           in_method=in_noflow_clusters$chosen_hclust,
                                           in_kclass=kclass,
                                           comp_clust_cut=base_hclust_cut))
    return(data.table(
      var = colnames(hydrostats_mat_sub)[[col_to_permute]],
      ari = ari_list
    )) 
  }) %>% rbindlist %>%
    .[, mean_ari := mean(ari), by=var] %>%
    .[, var := factor(var, levels=.SD[order(mean_ari), unique(var)])] %>%
    merge(unique(class_dt[, .(variable, aspect)]), 
          by.x='var', by.y='variable')
  
  
  #Plot variable importance
  p_varimp <- ggplot(hydrostats_ari, aes(x=var, y=ari, fill=aspect)) +
    geom_violin(alpha=0.7, draw_quantiles=0.5, color=NA) +
    geom_text(data=hydrostats_ari[!duplicated(var),], 
              aes(y=mean_ari, label=round(mean_ari,2)),
              #vjust=-.5, 
              alpha=0.75) +
    scale_y_continuous(name=
      'More important                                           Less important
      Adjusted Rand Index (ARI) after metric permutation'
    ,
    breaks=seq(0.2,1,0.2), limits=c(0.2,1), expand=c(0,0)) +
    scale_x_discrete(name='Hydrologic metric') +
    scale_fill_manual(name='Flow regime aspect',
                      values=c('#543005', '#80cdc1', '#bf812d', 
                               '#35978f', '#dfc27d')) +
    #scale_fill_viridis(name='Flow regime aspect', discrete=T, option='D') +
    coord_flip() +
    theme_classic() +
    theme(panel.grid.major = element_line(),
          text=element_text(size=12),
          axis.title.x= element_text(hjust = 0.5))
  
  #Save plot-----------------
  if (export & !is.null(fig_outdir)) {
    ggsave(file.path(fig_outdir, paste0('p_varimp_',
                                        format(Sys.Date(), '%Y%m%d'), '.png')),
           (p_varimp),
           width = 15, height = 15, units='cm'
    )
  }
  #---------------- Compute stability against single-gauge removal  ------------
  #See how cluster stability evoldes with increasing number of clusters
  boot_clust_klist <- lapply(seq(3,kclass), function(k) {
    fpc::clusterboot(
      data= hydrostats_mat_sub,
      B=100,
      bootmethod = 'boot',
      clustermethod = hclustCBI,
      k=k, 
      cut="number", 
      method=in_noflow_clusters$chosen_hclust,
      showplots=FALSE
    )
  })
  names(boot_clust_klist) <- paste0('cl', seq(3, kclass))
  
  #Compute cluster stability for chosen number of clusters and format it
  boot_clust <- fpc::clusterboot(
    data= hydrostats_mat_sub,
    B=100,
    bootmethod = 'boot',
    clustermethod = hclustCBI,
    k=kclass, 
    cut="number", 
    method=in_noflow_clusters$chosen_hclust,
    showplots=FALSE
  )
  
  base_hclust_cut_ordered <- base_clust_reslist$hclust %>%
    dendextend::cutree(k=kclass, order_clusters_as_data = T) %>%
    as.data.table(keep.rownames=T) %>%
    setnames(c('grdc_no', 'gclass_ordered')) %>%
    merge(
      in_noflow_clusters$cluster_analyses[[paste0('ncl', kclass)]]$class_dt[
        !duplicated(grdc_no), .(grdc_no, gclass)],
      by='grdc_no')
  
  
  out_bootstats <- base_hclust_cut_ordered[
    , .N, by=c('gclass', 'gclass_ordered')] %>%
    .[order(gclass_ordered), `:=`(
      bootmean = boot_clust$bootmean,
      bootbrd = boot_clust$bootbrd,
      bootrecover = boot_clust$bootrecover
    )] 
  
  if (!is.null(classnames)) {
    out_bootstats <- out_bootstats[order(gclass),] %>%
      .[, `:=`(gclass = classnames,
               gclass_ordered = NULL)] %>%
      .[order(gclass),]
  }
  
  #clusters with a stability value less than 0.6 should be considered unstable. 
  #Values between 0.6 and 0.75 indicate that the cluster is measuring a pattern 
  #in the data, but there isn’t high certainty about which points should be 
  #clustered together
  
  return(list(
    var_imp_ari = hydrostats_ari,
    gboot_clust = out_bootstats,
    boot_clust_klist = boot_clust_klist
  ))
}

#------ Analyze gauge representativeness ---------------------------------------
# in_gaugep_dt <- tar_read(gaugep_dt)
# in_noflow_clusters <- tar_read(noflow_clusters)
# in_predvars <- tar_read(predvars)
# in_gires_dt <- tar_read(gires_dt)

analyze_gauge_representativeness <- function(in_noflow_clusters,
                                             in_gaugep_dt,
                                             in_predvars,
                                             in_gires_dt,
                                             export = T,
                                             fig_outdir = NULL
) {
  #Get gauges class and attributes
  kclass <- in_noflow_clusters$kclass
  in_gaugep_dt[, grdc_no := as.character(grdc_no)]
  
  gclass_dt <- in_noflow_clusters$cluster_analyses[[
    paste0('ncl', kclass)]]$class_dt %>%
    .[!duplicated(grdc_no), .(grdc_no, gclass)] %>%
    merge(in_gaugep_dt[, -c('UPLAND_SKM', 'LENGTH_KM'), with=F], 
          by='grdc_no') %>%
    merge(in_gires_dt, by='HYRIV_ID', all.y=F)
  
  #Analyze relationship between predicted probability of intermittence
  #and classes
  giresprob_class_p <- ggplot(melt(gclass_dt, id.vars=c('gclass','grdc_no'), 
                                   measure.vars = c('predprob1', 'predprob30')),
                              aes(x=gclass, y=value, fill=variable)) +
    geom_boxplot() +
    scale_fill_discrete(
      name=str_wrap('Predicted mean annual probability of no-flow', 25),
      labels=c('At least 1 day/year', 'At least 30 day/year')
    ) +
    theme_classic()
  
  #names(gclass_dt)
  
  #Compute density distributions of gauges and entire network of non-perennial rivers
  print('Make distribution plot')
  envhist <- layout_ggenvhist(
    in_rivernetwork = in_gires_dt[predprob1>=0.5,],
    in_gaugepred = gclass_dt,
    in_predvars = in_predvars)
  
  #Scale data for computing distributional distance ----------------------------
  cols_to_analyze <- c(
    "bio1_dc_uav", "bio7_dc_uav", "bio11_dc_uav",
    "bio12_mm_uav", "bio14_mm_uav", "slp_dg_uav",
    "ari_ix_uav", "dis_m3_pyr", "sdis_ms_uyr", "UPLAND_SKM",
    "lka_pc_use", "snw_pc_uyr", "kar_pc_use", "for_pc_use") %>%
    .[. %in% names(in_gires_dt)]
  
  #Scale network data ----------------------------------------------------------
  print('Scale variables')
  net_scaling_parameters <- lapply(cols_to_analyze, function(in_col) {
    #print(in_col)
    transform_out <- transform_scale_vars(in_dt=in_gires_dt,
                                          value_col=in_col,
                                          var_col=NULL, 
                                          inplace=F,
                                          samp_frac = 0.1,
                                          min_boxcox=-2,
                                          max_boxcox=2)
    in_gires_dt[, (in_col) := transform_out$dt[[in_col]]]
    
    return(transform_out[-1])
  })
  scaling_parameters_dt <- data.table::rbindlist(net_scaling_parameters) %>%
    .[, col := cols_to_analyze]
  
  #Scale gauges data -----------------------------------------------------------
  gclass_dt_trans <- copy(gclass_dt)
  
  transform_column_wpars <- function(in_dt, in_col, pars) {
    floor_par <- ifelse(pars$min_val <= 0, pars$min_floor - pars$min_val, 0)
    in_dt <- in_dt[, 
                   (in_col) := (BoxCox(get(in_col) + floor_par, 
                                       lambda = pars$bc_lambda) 
                                - pars$trans_mean) / pars$trans_sd]
  }
  
  # Use lapply to iterate over cols_to_analyze
  purrr::walk(cols_to_analyze, function(in_col) {
    pars <- scaling_parameters_dt[col == in_col,]
    transform_column_wpars(in_dt=gclass_dt_trans, in_col, pars)
  })
  
  #Compute overall Wasserstein distance for each variable ----------------------
  print('Compute overall Wasserstein distance for each variable ')
  gires_varmeans <- in_gires_dt[predprob1>=0.5,
                                lapply(.SD, mean), 
                                .SDcols = cols_to_analyze]
  
  # calculating standardized bias and Wasserstein distance for each variable
  bias_dt <- gclass_dt_trans %>%
    data.table::melt(measure.vars = cols_to_analyze) %>%
    merge(as.data.table(t(gires_varmeans), keep.rownames = T), 
          by.x='variable', by.y='rn') %>%
    .[, list(bias = mean(value-V1, na.rm=T)), by=variable]
  
  wasser_dt <- lapply(cols_to_analyze, function(in_col) {
    print(in_col)
    if (gclass_dt_trans[!is.na(get(in_col)), .N]>0) {
      #Take a sample of full dataset
      n_segs <- in_gires_dt[predprob1>=0.5 & !is.na(get(in_col)), .N]
      net_samp <- sample(n_segs, n_segs/10)
      
      #Compute Wasserstein distance
      wassd <- waddR::wasserstein_metric(
        x=in_gires_dt[predprob1>=0.5 & !is.na(get(in_col)),][net_samp,get(in_col)],
        y=gclass_dt_trans[!is.na(get(in_col)), get(in_col)],
        p=2)
    } else {
      wassd <- NA
    }
    return(data.table(variable=in_col, wasser=wassd))
  }) %>% rbindlist
  
  univar_dist <- merge(bias_dt, wasser_dt, by='variable') %>%
    merge(in_predvars[, .(varcode, Attribute)],
          by.x='variable', by.y='varcode') %>%
    .[, Attribute := str_trim(
      gsub("(BIO[1-9]{1,2}[ ][-][ ])|([(]BIO5[-]BIO6[)])",
           "", Attribute))] %>%
    .[, Attribute := factor(Attribute, 
                            levels=.SD[order(wasser),Attribute])]
  
  #Plot Wasserstein distance -------
  p_wasser <- ggplot(univar_dist, aes(x=Attribute,
                                      y=wasser, color=bias<0)) +
    geom_segment(aes(xend=Attribute, yend=0.2), alpha=0.75) +
    geom_point(size=4, alpha=0.75) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    scale_y_continuous(breaks=c(0.2, 0.6, 1.0, 1.4)) +
    scale_color_manual(name= 'Mean bias',
                       labels = c('Negative', 'Positive'),
                       values=c('#01665e', '#8c510a')) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = c(0.8, 0.5),
          axis.title.y = element_blank(),
          text=element_text(size=12))
  
  p_wasser_envhist <- p_wasser + envhist + plot_layout(design= "ABBBBBB
  ABBBBBB
  ABBBBBB
  ABBBBBB"
  )
  
  #Optimize gains in multivariate gauge-network similarity -----------------------
  #For each ungauged river reach, compute the change in Kullback-Leibler 
  #divergence between the gauge multivariate normal distribution and that of the 
  #network if a gauge was placed on that reach
  #https://cfwp.github.io/rags2ridges/
  print('Compute Kullback-Leibler divergence ')
  in_comp_mat <- as.matrix(gclass_dt_trans[, cols_to_analyze, with=F])
  in_ref_mat <- as.matrix(in_gires_dt[predprob1>=0.5, cols_to_analyze, with=F])
  
  KLdiv_current <- get_KLdiv(
    comp_mat=rbind(in_comp_mat, in_ref_mat[1,]),
    ref_mat=as.matrix(in_gires_dt[predprob1>=0.5, cols_to_analyze, with=F])
  )
  
  start=Sys.time()
  KLdiv_marginal <- in_gires_dt[
    predprob1>=0.5,
    get_KLdiv(
      comp_mat = rbind(
        in_comp_mat,
        as.matrix(.SD)),
      cov_ref = KLdiv_current$cov_ref,
      mean_ref = KLdiv_current$mean_ref)$KLdiv,
    by='HYRIV_ID',
    .SDcols = cols_to_analyze] %>%
    .[, list(
      HYRIV_ID,
      KLdiv_diff = V1 - KLdiv_current$KLdiv)]
  print(Sys.time()-start)
  
  
  # Save plots and return results -----------------------------------------------
  if (export & !is.null(fig_outdir)) {
    ggsave(file.path(fig_outdir, paste0('p_wasser_envhist_',
                                        format(Sys.Date(), '%Y%m%d'), '.pdf')),
           (p_wasser_envhist),
           width = 35, height = 20, units='cm'
    )
  }
  
  return(list(
    plot_giresprob_class = giresprob_class_p,
    plot_wasser_envhist = p_wasser_envhist,
    univar_dist =  univar_dist,
    KLdiv_current = KLdiv_current,
    KLdiv_marginal = KLdiv_marginal
  ))    
}



#------ Analyze environmental correlates ---------------------------------------
# in_gaugep_dt <- tar_read(gaugep_dt)
# in_noflow_clusters <- tar_read(noflow_clusters)
# in_predvars <- tar_read(predvars)
#in_gires_dt <- tar_read(gires_dt)
#in_colors = class_colors

analyze_environmental_correlates <- function(in_cluster_postanalysis,
                                             in_gaugep_dt,
                                             in_predvars,
                                             in_gires_dt,
                                             in_colors,
                                             export = T,
                                             fig_outdir = NULL
) {
  #Read and compute derived variables for global river network
  print('Read network')
  
  #Get gauges class and attributes
  class_dt <- in_cluster_postanalysis$class_dt
  kclass <- class_dt[, length(unique(gclass))]
                     
  in_gaugep_dt[, grdc_no := as.character(grdc_no)]
  
  gclass_dt <- class_dt %>%
    .[!duplicated(grdc_no), .(grdc_no, gclass)] %>%
    merge(in_gaugep_dt[, -c('UPLAND_SKM', 'LENGTH_KM'), with=F], 
          by='grdc_no') %>%
    merge(in_gires_dt, by='HYRIV_ID', all.y=F)
  
  
  gclass_dt[, gclass_n := .N, by=gclass][
    , gclass_weight := max(gclass_n)/gclass_n ]
  
  setnames(gclass_dt, in_predvars$varcode, in_predvars$varname, skip_absent = T)
  
  #----------------- BUILD TREE ------------------------------------------------
  #Train a classification tree
  tree_formula <- as.formula(paste0("gclass~`",
                                    paste(names(gclass_dt)[names(gclass_dt)
                                                           %in% in_predvars$varname],
                                          collapse='`+`'),
                                    "`")
  )
  
  split.fun <- function(x, labs, digits, varlen, faclen)
  {
    # replace commas with spaces (needed for strwrap)
    labs <- gsub(",", " ", labs)
    for(i in 1:length(labs)) {
      # split labs[i] into multiple lines
      labs[i] <- paste(strwrap(labs[i], width=30), collapse="\n")
    }
    labs
  }
  
  rpart_mod1 <- rpart(formula=tree_formula, data=gclass_dt,
                      weights = gclass_dt$gclass_weight,
                      control = rpart.control(cp = 0.008)
  )
  
  tree_plot <- rpart.plot(rpart_mod1, type=3, extra=8,
                          uniform=T, tweak=1.8,
                          box.palette = as.list(in_colors),
                          split.fun=split.fun)
  
  #----------------- CREATE STATS TABLE ----------------------------------------
  gclass_dt_sub<-  gclass_dt[
    , c("grdc_no",
        "gclass",
        "predprob1",
        "Drainage area",
        "Natural Discharge pour point Annual average",
        "BIO11 - Mean Temperature of Coldest Quarter catchment Average",
        "Snow Cover Extent catchment Maximum or Annual maximum" ,
        "BIO14 - Precipitation of Driest Month catchment Average",
        "BIO15 - Precipitation Seasonality (Coefficient of Variation) catchment Average",
        "Global Aridity Index watershed Average",
        "Terrain Slope watershed Average",
        "Sand Fraction in Soil 0-100 cm watershed Average"),
    with=F
  ] 
  
  env_tab <- melt(rbind(gclass_dt_sub,
                        copy(gclass_dt_sub)[, gclass := 'All']),
                  id.vars=c('gclass', 'grdc_no')) %>%
    .[, list(mean = mean(value, na.rm=T), 
             q10 = quantile(value, 0.1, na.rm=T),
             q90 = quantile(value, 0.9, na.rm=T),
             classn = length(unique(grdc_no))),
      by=.(variable, gclass)] %>%
    digitform(cols=c('mean', 'q10', 'q90'),
              extradigit = 1, inplace=F) %>%
    .[, stat_format := paste0(mean, ' (', q10, '-', q90, ')')] %>%
    dcast(gclass+classn~variable, value.var = 'stat_format')
  
  #----------------- BUILD BOXPLOT ------------------------------------------------
  #Get box plot
  bp_theme <- function(){
    theme_bw() %+replace%
      theme(legend.position = 'none',
            axis.line = element_line(color='darkgrey'),
            axis.title.y = element_text(
              angle = 90,
              margin = margin(t = 0, r = 0, b = 0, l = 0)),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            text=element_text(size=12),
            strip.background = element_blank(),
            strip.text = element_text(margin = margin(0,0,0.1,0, "cm")))
  }
  
  
  bp_uparea <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, y=`Drainage area`, color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_log10(name=expression('Drainage area -'~km^2)) +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  bp_dis <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, y=`Natural Discharge pour point Annual average`, color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_log10(name=expression('Average discharge -'~m^3~s^-1),
                  breaks=c(0.01, 0.1, 1, 10, 100, 1000)) +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  bp_coldtemp <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, 
        y=`BIO11 - Mean Temperature of Coldest Quarter catchment Average`, 
        color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(name=str_wrap(
      'Mean T of coldest quarter in catchment - °C', 30)) +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  bp_bio14 <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, 
        y=`BIO14 - Precipitation of Driest Month catchment Average`, 
        color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_sqrt(name=str_wrap(
      'Precipitation of driest month in catchment - mm', 30)) +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  bp_bio15 <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, 
        y=`BIO15 - Precipitation Seasonality (Coefficient of Variation) catchment Average`, 
        color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(name=str_wrap(
      'Precipitation seasonality (CV) in catchment', 30)) +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  bp_slo <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, 
        y=`Terrain Slope watershed Average`/10, 
        color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(name='Terrain slope - degrees') +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  bp_sand <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, 
        y=`Sand Fraction in Soil 0-100 cm watershed Average`, 
        color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(name=str_wrap(
      'Sand fraction in watershed soil (0-100 cm) - %', 30)) +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  bp_predprob <- ggplot(
    gclass_dt_sub, 
    aes(x=gclass, 
        y=predprob1, 
        color=factor(gclass))) +
    geom_jitter(size=0.4, alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(name=str_wrap(
      'Predicted probability of intermittence', 30)) +
    scale_x_discrete(name='Class') +
    scale_color_manual(values=in_colors) +
    #geom_hline(yintercept=0) +
    coord_cartesian(clip='off') +
    bp_theme()
  
  
  bp_mosaic <- (bp_uparea + bp_dis + bp_bio14 + bp_bio15 +
                  bp_coldtemp + bp_predprob + bp_slo + bp_sand + 
                  plot_layout(ncol = 2) 
  )
  
  gclass_dt_sub[, range(`Drainage area`), by=gclass]
  
  #Save and return------------------------------------
  if (export & !is.null(fig_outdir)) {
    ggsave(file.path(fig_outdir, paste0('p_boxplot_env_ward_D2',
                                        format(Sys.Date(), '%Y%m%d'), '.png')),
           (bp_mosaic),
           width = 15, height = 25, units='cm'
    )
    
    pdf(file=file.path(fig_outdir, paste0('p_tree_env_ward_D2',
                                          format(Sys.Date(), '%Y%m%d'), '.pdf')),
        width = 20, height = 10)
    rpart.plot(rpart_mod1, type=3, extra=8,
               uniform=T, tweak=1.15,
               box.palette = as.list(in_colors),
               split.fun=split.fun)
    dev.off()
    
    fwrite(env_tab, 
           file.path(fig_outdir,paste0('envstats_',
                                       format(Sys.Date(), '%Y%m%d'), '.csv')
           ))
  }
  
  return(list(
    p_boxplot = bp_mosaic,
    p_tree = tree_plot,
    tab = env_tab
  ))
}


