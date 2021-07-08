spMixCall <- function( Site, Yts_Xts,  
                       Z_ts = NULL, Data_Str,
                       Object = "Object", Obj.Seq = 1,
                       Tab = "spCali", Total = FALSE,
                       prior = NULL, para = NULL,
                       true.para = NULL, Database = FALSE,
                       parallel = FALSE, verbose = TRUE,
                       verbose.VB = FALSE, heavy.tail = FALSE,
                       method = c("ensemble"),
                       response.transf = c("sr"),# c("normal", "sr", "log")
                       covariate.transf = c("normal"),
                       Ensemble.size = 100, 
                       ds = 1e-2, cs = 0.5,
                       ct = 1, IS.size = 200,  
                       Thresh = c(1, 1, 1),
                       Remove.Cores = 5, N.Chunk = 1,
                       itMax = 1e2,  tol.vb = 1e-2,
                       tol.real = 1e-3, 
                       seed = 1234)
{
  call = match.call()
  if(is.null(Site)){stop("Must provide Site data.\n")}
  if(is.null(Yts_Xts)){stop("Must provide response and covariates data.\n")}
  if(is.null(Data_Str)){stop("Must provide Data_Str data.\n")}
  if(is.null(prior))
  {
    ####################################################################
    p <- dim(Yts_Xts$X_ts)[1]
    prior <- list(
      beta = list(mu = rep(0, p), Sigma2 = 1e5*diag(p))
      , tau2 = list(a = 1, b = 1)
      , theta1 = list(mu = 1e-2, Sigma2 =  1e5)
    )
  }
  
  if(is.null(para))
  {
    para <- list(
      beta = list(E_beta = rep(0, p), Sigma2 = diag(p))
      , theta1 = list(E_theta1 = 1e-3, Sigma2 = 1)
      , alpha = list(E_alpha = 1, Sigma2 = diag(1))
      , k = list(E_k = 5, a = 1e-3, b = 1e1)
      , k0 = list(E_k0 = 8, a = 1e-3, b = 1e1)
      , theta2 = list(E_theta2 = 1, a = 1e-3, b = cs*max(Data_Str$BAUs.Dist)) #
      , tau2 = list(E_tau2 = 0.05, a = 1, b = 1)
    )
  }
  ######################################################################
  ######################################################################
  if(is.null(true.para))
  {
    if(is.null(Z_ts)){
      true.para <- list(betaX = rep(NA, dim(Yts_Xts$X_ts)[1]),
                        Obs.tau2 = NA, alpha = NA,
                        Proc.tau2 = NA,  Proc0.tau2 = NA,
                        a = NA, b = NA,
                        theta = c(NA, NA), k0 = NA, k = NA)
      
    }else{
      true.para <- list(betaX = rep(NA, dim(Yts_Xts$X_ts)[1]),
                        betaZ = rep(NA, dim(Z_ts)[1]),
                        Obs.tau2 = NA, alpha = NA,
                        Proc.tau2 = NA,  Proc0.tau2 = NA,
                        a = NA, b = NA,
                        theta = c(NA, NA), k0 = NA, k = NA)
    }
  }
  ######################################################################
  ######################################################################
  if(Total)
  {
    # GSD = Total_Data(CMAQ_PM25, Data_Str, Yts_Xts$Y_ts,
    #                  Yts_Xts$X_ts, YearMonth = c(YearMonth))
    Train <- list(
      n = ncol(Yts_Xts$Y_ts)
      , Nt = nrow(Yts_Xts$Y_ts)
      , N.BAUs = Data_Str$N.BAUs
      , Y_ts_true = Yts_Xts$Y_ts_true
      , Y_ts = Yts_Xts$Y_ts
      , X_ts = Yts_Xts$X_ts
      , BAUs.Dist = Data_Str$BAUs.Dist
      , Adj.Mat = Data_Str$G
      , Hs = Data_Str$Hs
    )
    
    ######################################################################
    #                       Whether to store in the database
    ######################################################################
    Tab_Name <- paste0(Tab, '_'
                       , if_else(month(Sys.Date()) > 9
                                 , as.character(month(Sys.Date()))
                                 , paste0("0",  as.character(month(Sys.Date()))))
                       , "_", if_else(day(Sys.Date()) > 9
                                      , as.character(day(Sys.Date()))
                                      , paste0("0",  as.character(day(Sys.Date()))))
                       # , "_", if_else(hour(Sys.time())>9
                       #                , as.character(hour(Sys.time()))
                       #                , paste0("0",  as.character(hour(Sys.time()))))
                       
    )
    
    if(Database){
      DSN_01 <- odbcConnect("DSN_01", uid = "myname", pwd = "mypwd"
                            , believeNRows = FALSE, case = "toupper")
  
      
      sqlDrop(DSN_01, paste0(Tab_Name), errors = F)
      database = list(DSN = DSN_01, Table = Tab_Name)
    }else{ database = list(DSN = NULL, Table = Tab_Name)}
    ######################################################################
    ######################################################################
    cat("Oject(dataset): Full dataset ...\n")
    CV.Re <- spMixVBEnKs(data = Train, test = NULL, prior = prior, 
                         para = para, true.para = true.para,
                         database = database, Object = "Total", 
                         parallel = parallel, verbose = verbose, 
                         verbose.VB = verbose.VB, 
                         heavy.tail = heavy.tail, method = method,
                         response.transf = response.transf,
                         covariate.transf = covariate.transf,
                         Ensemble.size = Ensemble.size,
                         ds = ds, cs = cs, ct = ct, IS.size = IS.size,
                         Thresh = Thresh, Remove.Cores = Remove.Cores,
                         N.Chunk = N.Chunk, itMax = itMax, seed = seed,
                         tol.vb = tol.vb, tol.real = tol.real)
    CV.Re$Call = call
    class(CV.Re) <- "spMixCall"
  }else{
    CV.Re <- list()
    setDF(Site)
    # n0 <- length(unique(Site[, Object]))
    # if(Object == "Object"){n0 <- n0 -1}
    for(num in Obj.Seq)
    {
      GSD = Test_Train_Fun(Site, Data_Str, Yts_Xts, Z_ts = Z_ts, Object, num)
      Train = GSD$Train
      Test =  GSD$Test
      # Train$Y_ts = sqrt(Train$Y_ts)
      # Train$X_ts = sqrt(Train$X_ts)
      # Test$X_ts = sqrt(Test$X_ts)
      ######################################################################
      #                       Whether to store in the database
      ######################################################################
      Tab_Name <- paste0(Tab, '_'
                         , if_else(month(Sys.Date()) > 9
                                   , as.character(month(Sys.Date()))
                                   , paste0("0",  as.character(month(Sys.Date()))))
                         , "_", if_else(day(Sys.Date()) > 9
                                        , as.character(day(Sys.Date()))
                                        , paste0("0",  as.character(day(Sys.Date()))))
                         # , "_", if_else(hour(Sys.time())>9
                         #                , as.character(hour(Sys.time()))
                         #                , paste0("0",  as.character(hour(Sys.time()))))
                         , "_", GSD$Object
      )
      if(Database){
        DSN_01 <- odbcConnect("DSN_01", uid = "myname", pwd = "mypwd"
                              , believeNRows = FALSE, case = "toupper")
        
        
        sqlDrop(DSN_01, paste0(Tab_Name), errors = F)
        database = list(DSN = DSN_01, Table = Tab_Name)
      }else{ database = list(DSN = NULL, Table = Tab_Name)}
      ######################################################################
      ######################################################################
      R <- spMixVBEnKs(data = Train, test = Test,
                       prior = prior,  para = para,
                       true.para = true.para, database = database,
                       Object = GSD$Object,  parallel = parallel,
                       verbose = verbose, verbose.VB = verbose.VB,
                       heavy.tail = heavy.tail,
                       method = method, response.transf = response.transf,
                       covariate.transf = covariate.transf,
                       Ensemble.size = Ensemble.size,
                       ds = ds, cs = cs, ct = ct, IS.size = IS.size,
                       Thresh = Thresh,  N.Chunk = N.Chunk,
                       Remove.Cores = Remove.Cores, seed = seed,
                       itMax = itMax, tol.vb = tol.vb,
                       tol.real = tol.real)
      CV.Re[[num]] <- R
      
      # FILE <- paste0("./data/Generate_Data/Test/")
      # test.result <- R$test.result
      # ParaEst.list <- R$temp
      # FinalEst <- R$FinalEst
      # EnKs <- R$Ks
      # PIU <- R$PIU
      # IS <- R$IS
      # save(test.results, file = paste0(FILE, R$City_Name, "_Test.Rdata"))
      # save(ParaEst.list, file = paste0(FILE, R$City_Name, "_ParaEst.Rdata"))
      # save(FinalEst, file = paste0(FILE, R$City_Name, "_FinalEst.Rdata"))
      # save(EnKs, file = paste0(FILE, R$City_Name, "_EnKs.Rdata"))
      # save(PIU, file = paste0(FILE, R$City_Name, "_PIU.Rdata"))
      # save(IS, file = paste0(FILE, R$City_Name, "_IS.Rdata"))
      # rm(test.result, ParaEst.list, FinalEst, EnKs,PIU, IS, R)
    }
    # R = 1;
  }
  return(CV.Re)
}

######################################################################
print.spMixCall = function(x, ...) {
  if(length(x) == 0){stop("The dimension of R is 1.\n")}
  # print warning if maxit reached
  
  
  # print call
  cat("\nThere were", paste0(length(x) ,
                             " Objects that have been cross-verified by leaving one object out."),
      "\n\n")
  # print summry
  CV <- EST <- NULL
  for(i in 1:length(x))
  {
    da <-  x[[i]][["Final_Pred"]]
    index <- which(colnames(da) == "K0")
    index.obj <- which(colnames(da) == "Object")
    if(length(index) ==0)
    {
      index <- which(colnames(da) == "theta1")
    }
    EST <- rbind(EST, da[, c(index.obj, 1:index)])
    CV <- rbind(CV, da[, c(-c(1:index))])
  }
  
  EST[, -1] <- round(EST[, -1], 3)
  setorder(CV,  Object);setorder(EST,  Object);
  CV[, "Object"] <- as.character(CV[, "Object"])
  da0 <- setDF(CV[1, ])
  for(j in 1:length(CV))
  {
    da0[1, j] <- ".."
  }
  index.obj <- which(colnames(CV) == "Object")
  CV <- rbind(CV, da0, cbind(data.frame(Object = "Avg"),
                             t(round(colMeans(CV[, -index.obj]), 3))))%>%
    setcolorder("Object")
  # print(EST)
  #
  # print(CV)
  return(list(EST = as.data.frame(EST), CV = as.data.frame(CV)))
}
