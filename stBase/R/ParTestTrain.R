######################################################################
######################################################################
Test_Train_Fun <- function(Site = Site, Data_Str, 
                           Yts_Xts, Z_ts = NULL,
                           Object = "Object", num)
{
  Hs = Data_Str$Hs;
  Y_ts_true = Yts_Xts$Y_ts_true
  Y_ts = Yts_Xts$Y_ts;
  X_ts = Yts_Xts$X_ts
  setDF(Site)
  Obj <- sort(unique(as.character(Site[, Object])))[num]
  cat("\n\nTest object(dataset): ", Obj, "...\n\n")
  ######################################################################
  ######################################################################
 setDF(Site)
  Test.ID <- unique(Site[as.character(Site[, Object]) == Obj,
                         "SITEID"]) %>% as.character()
  Test.Hs <- Hs[rownames(Hs) %in% Test.ID,]
  Train.Hs <- Hs[rownames(Hs) %nin% Test.ID,]

  Test.Y_ts_true <- Y_ts_true[, colnames(Y_ts_true) %in% Test.ID]
  Test.Y_ts <- Y_ts[, colnames(Y_ts) %in% Test.ID]
  Test.X_ts <- X_ts[, dimnames(X_ts)[[2]] %in% Test.ID, ]

  Train.Y_ts_true <- Y_ts_true[, colnames(Y_ts_true) %nin% Test.ID]
  Train.Y_ts <- Y_ts[, colnames(Y_ts) %nin% Test.ID]
  Train.X_ts <- X_ts[, dimnames(X_ts)[[2]] %nin% Test.ID, ]

  #####################################################################
  if(is.null(Z_ts))
  {
    Train <- list(
      n = ncol(Train.Y_ts)
      , Nt = nrow(Train.Y_ts)
      , N.BAUs = Data_Str$N.BAUs
      , Y_ts_true = Train.Y_ts_true 
      , Y_ts = Train.Y_ts
      , X_ts = Train.X_ts
      , BAUs.Dist = Data_Str$BAUs.Dist
      , Adj.Mat = Data_Str$G
      , Hs = Train.Hs
    )
    Test <- list(
        Y_ts_true = Test.Y_ts_true
      , Y_ts = Test.Y_ts
      , X_ts = Test.X_ts
      , H = Test.Hs
      , Object = Obj
    )
  }else{
    Train <- list(
      n = ncol(Train.Y_ts)
      , Nt = nrow(Train.Y_ts)
      , N.BAUs = Data_Str$N.BAUs
      , Y_ts_true = Train.Y_ts_true 
      , Y_ts = Train.Y_ts
      , X_ts = Train.X_ts
      , Z_ts = Z_ts
      , BAUs.Dist = Data_Str$BAUs.Dist
      , Adj.Mat = Data_Str$G
      , Hs = Train.Hs
    )
    Test <- list(
      Y_ts_true = Test.Y_ts_true
      , Y_ts = Test.Y_ts
      , X_ts = Test.X_ts
      , Z_ts = Z_ts
      , H = Test.Hs
      , Object = Obj
    ) 
  }
  
  Train.Test = list(Train = Train
             , Test = Test
             , Object = Obj)
  return(Train.Test)
}
######################################################################
######################################################################
Total_Data <- function(CMAQ_PM25, Data_Str, X = "CMAQ_PM25_30",
                       YearMonth = c(201511, 201512))
{
  CMAQ_PM25_Test <- CMAQ_PM25 %>% filter(YEAR_MONTH %in% YearMonth
                                         # ,day(DATE_TIME)  %in% c(1, 2)
  )

  ###########################################################################
  Nt <- length(unique(CMAQ_PM25_Test$DATE_TIME))
  CMAQ_PM25_Test$T_index <- 1:Nt
  CMAQ_ID <- unique(CMAQ_PM25_Test$CMAQ_ID)
  N.BAUs <- length(CMAQ_ID)
  CMAQ_N <- length(CMAQ_ID)
  CMAQ_X_ts <- array(NA, dim = c(2, N.BAUs, Nt))
  CMAQ_X_ts[1, , ] = 1
  setDF(CMAQ_PM25_Test)
  for(t in 1:Nt)
  {
    # CMAQ_X_ts[2, , t] <- dcast(CMAQ_PM25_Test[T_index == t, ]
    #                            , . ~ CMAQ_ID, value.var = "CMAQ_PM25_30"
    # )[1, 2:(CMAQ_N + 1)]  %>% as.numeric()

    CMAQ_X_ts[2, , t] <- CMAQ_PM25_Test[T_index == t,
                         colnames(CMAQ_PM25_Test) == X]
  }
  ###########################################################################
  Test <- list(
                H = Data_Str$H
                , Nt = Nt
                , X_ts = CMAQ_X_ts
                , CMAQ_N = CMAQ_N
                , CMAQ_PM25_Test = CMAQ_PM25_Test
              )
  Total = list(Test = Test)
  return(Total)
}
# GSD = Total_Data(CMAQ_PM25, Data_Str, YearMonth = c(201511, 201512))

