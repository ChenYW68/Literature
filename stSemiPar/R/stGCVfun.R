stGCVfun <- function(X){
  source("E:/Literature/semiBase/R/util.R")
  Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
  source("./R/stSemiPar.R")
  # source("./R/stSemiPar_WLS.R")
  # source("./R/PSTVB_Packages.R")
  if(method %nin% c("WLS")){

 # call <- function(X){
 #   tryCatch(expr = {
 #     return(
    # assign("Q", as.matrix(data$Qst), envir = .GlobalEnv)
    fit <- stSemiPar(y_ts = data$Y_ts, 
                   x_ts = data$X_ts, 
                   z_ts = data$Z_ts,
                   loc = data$loc,
                   time = data$time,
                   method = method,
                   Kernel = Kernel,
                   h = H[X, ],
                   # prob = c(as.vector(H[X, 4]), prob[2]),
                   nuUnifb = nuUnifb,
                   nu = nu,
                   nThreads = nThreads,
                   nIter = nIter)
    # )
 #   }, error = function(e){
 #     return(NA)
 #     
 #   }
 #   )
 # }
 #   fit = NA
 #   while(is.na(fit)){
 #     fit <- call(X)
 #   }
  }else{
    # for (X in 1:50) {
    #  cat("X = ", X, "\n")
    fit <- stSemi_WLS(y_ts = data$Y_ts,
                      x_ts = data$X_ts,
                      z_ts = data$Z_ts,
                      loc = data$loc,
                      time = data$time,
                      # prob = c(as.vector(H[X, 4]), prob[2]),
                      Kernel = Kernel,
                      h =  as.vector(H[X]),
                      nThreads = nThreads,
                      nIter = nIter)
     # }
  }
  # plot(fit$fit.value, fit$y_ts)
  y.fitted.error <- spT.validation(fit$fit.value, fit$y_ts)[c(2, 3, 9)]
  y.fitted.error <- matrix(y.fitted.error, nrow = 1)
  colnames(y.fitted.error) <- paste0("y", c("_RMSE", "_MAE", "_CRPS"))
  
  
  da.Beta <- NULL
  if(!is.null(data$X_ts)){
    for (i in 1:(dim(data$X_ts)[1])) {
      da.Beta <- cbind(da.Beta, fit$beta[i, 1])
    }
    colnames(da.Beta) <- paste0("Beta_", 1:(dim(data$X_ts)[1]))
  }
  da.theta <- bias.theta <- NULL
  if(!is.null(data$theta)){
    for(i in 1:ncol(data$theta)){
      theta.error <- spT.validation(fit$theta$alpha[, i], data$theta[, i])[c(2, 3, 9)]
      da.theta <- c(da.theta, as.vector(theta.error))
      # bias.theta <- c(bias.theta, c(sum((fit$theta$alpha[, i] -  data$theta[, i])^2),
      #                               sd((fit$theta$alpha[, i] -  data$theta[, i])^2)
      #                               ))
      bias.theta <- c(bias.theta, sfsmisc::integrate.xy(x = data$time, 
                                                        fx = (fit$theta$alpha[, i] -  
                                                                data$theta[, i])^2)               
      )
    }
    da.theta <- matrix(da.theta, nrow = 1)
    # rownames(da.theta) <- paste0("theta.", 1:ncol(data$theta))
    colnames(da.theta) <- paste0(rep(paste0("theta_", 1:ncol(data$theta)),
                                     each = 3), c("_RMSE", "_MAE", "_CRPS"))
    
    bias.theta <- matrix(bias.theta, nrow = 1)
    colnames(bias.theta) <- paste0(rep(paste0("theta_", 1:ncol(data$theta)),
                                     each = 1), c("_MISE"))
    
  }
  
  # S <- fit$theta$S
  # if(is.null(data$Vc)){
  #   C1 <- fit$C
  #   C2 <- fit$C
  #   C1 = as(C1, "sparseMatrix") 
  #   C2 = as(C2, "sparseMatrix") 
  #   for (i in 1:(nrow(data$Y_ts) - 1)) {
  #     C2 <- Matrix::bdiag(C2, C1)
  #   }
  # }else{
    # C2 <- fit$C  %*% fit$C%*% fit$C
  # }
  
  
  GCV = y.fitted.error[1]^2/(mean(diag(1 - fit$theta$St)))^2 #%*% C2
  
  da <- as.data.frame(cbind(da.Beta, da.theta, bias.theta, y.fitted.error))
  rownames(da) <- NULL
  da <- list(summary = data.frame(H0 = H[X, 1], 
                                  H1 = H[X, 2], 
                   # covHs = H[X, 2],
                   # covHt = H[X, 3],
                   # taper_s = H[X, 4],
                    da, #Var.s = fit$Var.s,
                    sigma.sq.1 = fit$sigma.sq.s[1],
                    sigma.sq.2 = fit$sigma.sq.s[2],
                    range1 = fit$range_nugget[1],
                    range2 = fit$range_nugget[2],
                    tau_sq = fit$range_nugget[3],
                    GCV = GCV,
                    nIter = fit$nIter),
             GCV = GCV, 
             Cs = fit$Cs,
             alpha = fit$theta$alpha,
             Vt = fit$Vt,
             fit.value = fit$fit.value,
             y_ts = fit$y_ts
             )
  return(da)
}


GCVparaSemi <- function(data, H = expand.grid(meanH = seq(1E-2, 1,, 5),
                                              covHs = seq(1E-2, 1,, 5),
                                              covHt = seq(1E-2, 1,, 5)), 
                        method = "WI", Kernel = c(0, 0),
                        nuUnifb = 1, nu = 0.5, nThreads = 10,
                        nIter = 10)
{
  
  cl <- parallel::makeCluster(getOption("cl.cores", 
                                        ceil(nThreads/2)))#ceil(nThreads/2)
  parallel::clusterEvalQ(cl, library(plyr))
  parallel::clusterEvalQ(cl, library(dplyr))
  parallel::clusterEvalQ(cl, library(data.table))
  parallel::clusterEvalQ(cl, library(Matrix))
  parallel::clusterEvalQ(cl, library(fields))
  parallel::clusterEvalQ(cl, library(Hmisc))
  # parallel::clusterEvalQ(cl, library(semiBase))
  # assign("Q_temp", data$Q_temp, envir = .GlobalEnv)"Q_temp",
  parallel::clusterExport(cl = cl
                          , varlist = c("data",  "H"
                                        , "method", "Kernel"
                                        , "nuUnifb", "nu"
                                        , "nThreads", "nIter"
                                        , "spT.validation")
                          , envir = environment())
  temp <- parallel::parLapply(cl, X = 1:nrow(H), fun = stGCVfun)
  temp <- do.call("rbind", temp)
  # temp <- data.table::rbindlist(temp)
  stopCluster(cl)
  
  # library(snowfall)
  # # 并行初始化
  # sfInit(parallel = TRUE, cpus = n.omp.threads)
  # 
  # sfLibrary(plyr)
  # sfLibrary(dplyr)
  # sfLibrary(semiBase)
  # sfExport("data", "H", "k.fold"
  #           , "profile", "Kernel"
  #           , "GeomVariable_r"
  #           , "n.omp.threads"
  #           , "spSemiPara"
  #           , "Round"
  #           , "spT.validation")         # 载入依赖的对象
  # temp <- sfLapply(1:nrow(H), fun = CVparallel)
  # # 结束并行，返还内存等资源
  # sfStop()
  
  return(temp)
}