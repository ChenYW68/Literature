###########################################################################
#                                   Create Y_ts/X_ts
###########################################################################
ParYtsXts <- function(data = NULL,
                      include = list(YearMonth = c(201506),
                                     YEAR = c(2015),
                                     MONTH = c(1:12),
                                     DAY = c(1:31)),
                      remove = list(YearMonth = c(0),
                                    YEAR = c(0),
                                    MONTH = c(0),
                                    DAY = c(0)),
                      Y = "REAL_PM25",
                      X = c("CMAQ_PM25_30", "REAL_TEMP",
                            "REAL_PRES",
                            "REAL_RAIN", "REAL_HUMI",
                            "REAL_IRAIN", "REAL_LON_WIND",
                            "REAL_LAT_WIND", "REAL_DEWP"),
                      date_time = "DATE_TIME",
                      siteid = "SITEID",
                      Hs = NULL,
                      Z_ts = NULL)
{
  if (is.null(data)) { stop("Must provide data.\n")}
  setDT(data)
  data_base <- data %>%
    dplyr::filter(
      YEAR_MONTH %in% include$YearMonth,
      # # YEAR %in% include$YEAR,
      # MONTH %in% include$MONTH,
      # DAY %in% include$DAY,
      # YEAR_MONTH %nin% remove$YearMonth,
      # YEAR %nin% remove$YEAR,
      # MONTH %nin% remove$MONTH,
      # DAY %nin% remove$DAY,
      # SITEID %nin% c(3, 57, 62, 63, 8, 11,18,72)
    ) %>% setorderv(cols = c(siteid, date_time))
  
  data_base <- setDF(data_base)
  #data_base <-subset(data_base, with(data_base, !(DAY %in% c(30, 31) & MONTH == 8)))
  setDT(data_base)
  # %>%
  #   dplyr::select(SITEID, CITY, CITY_NAME
  #                 , DATE_TIME
  #                 , YEAR, MONTH, DAY
  #                 , YEAR_MONTH
  #                 , REAL_PM25
  #                 , CMAQ_PM25
  #                 , CMAQ_PM25_30
  #                 , NA.Spline
  #                 , NA.Kriging
  #                 , NA.Kalman
  #                 , CMAQ_ID
  #                 , LON, LAT
  #                 , LON_X, LAT_Y)
  
  setDF(data_base)
  data_base_true <- data_base
  NA.check <- is.na(data_base[, Y])
  if(length(NA.check[NA.check == T])>0)
  {
    data_base[, Y] = if_else(is.na(data_base[, Y])
                             , data_base$NA.Kriging
                             , data_base[, Y])
  }
  Nt <- length(unique(data_base[, date_time]))
  SiteId <- unique(data_base[, siteid]) %>% sort()
  n <- length(SiteId)
  data_base$T_index <- data_base_true$T_index <- 1:Nt
  Y_ts <- Y_ts_true <- matrix(NA, nrow = Nt, ncol = n)
  # setDT(PM25_CMAQ)
  
  ######################################################################
  #                         create  data structure
  ######################################################################
  p1 <- length(X) + 1
  if(!is.null(Z_ts)){
    if(length(dim(Z_ts))>2){
      p2 <- dim(Z_ts)[1]
      Hz = sapply(seq_len(dim(Z_ts)[1]), function(p)
        Hs %*% Z_ts[p,, ]
        , simplify="array") %>% aperm(c(3, 1, 2))
      p <- p1 +p2
    }else{
      p2 <- 1
      Hz = Hs %*% Z_ts[p,, ]
      p <- p1 +p2 
    }
  }else{
    p <- p1
  }
  X_ts <- array(NA, dim = c(p, n, Nt)
                , dimnames = list(as.character(c("Intercepts", X)), c(SiteId),
                                  as.character(unique(data_base[, date_time]))
                ))
  X_ts[1,,] = 1
  for(t in 1:Nt)
  {
    Y_ts[t, ] <- data_base[data_base$T_index == t,
                           colnames(data_base) == Y]
    Y_ts_true[t, ] <- data_base_true[data_base_true$T_index == t,
                                     colnames(data_base_true) == Y]
    for(k in 2:p1)
    {
      
      X_ts[k, , t]  <- data_base[data_base$T_index == t,
                                 colnames(data_base) == X[k - 1]]
    }
  }
  if(!is.null(Z_ts)){
    if(length(dim(Z_ts))>2){ 
      for(j in (p1 + 1): p)
      {
        X_ts[j, , ]  <- Hz[j - p1,,]
      }}else{
        X_ts[p1 + 1, , ]  <- Hz
      }
  }
  if(dim(X_ts)[1]>2){
    for(k in 2:dim(X_ts)[1])
    {
      X_ts[k,,] <- matrix(scale(as.vector(X_ts[k,,]))[, 1],
                          nrow = nrow(X_ts[k,,]), ncol = Nt)
    }}
  # colnames(X_ts)
  colnames(Y_ts) <- colnames(Y_ts_true) <- SiteId
  rownames(Y_ts) <- as.character(unique(data_base[, date_time]))
  rownames(Y_ts_true) <- as.character(unique(data_base_true[, date_time]))
  
  Y_X_ts <- list(Y_ts = Y_ts, Y_ts_true = Y_ts_true, X_ts = X_ts)
  return(Y_X_ts)
}
