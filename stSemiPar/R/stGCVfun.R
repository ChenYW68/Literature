stGCVfun <- function(X){
  source("E:/Literature/semiBase/R/util.R")
  Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
  source("./R/stSemiPar.R")
  fit <- stSemiPar(y_ts = data$Y_ts, 
                   x_ts = data$X_ts, 
                   z_ts = data$Z_ts,
                   loc = data$loc,
                   Vc =  data$Vc,
                   time = data$time,
                   h = c(H[X, 1], H[X, 2], H[X, 3]),
                   Inde = Inde,
                   method = method,
                   nIter = nIter)
  # plot(fit$fit.value, fit$y_ts)
  y.fitted.error <- spT.validation(fit$fit.value, fit$y_ts)
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
      theta.error <- spT.validation(fit$theta$alpha[, i], data$theta[, i])
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
    # C2 <- fit$C  %*% fit$C
  # }
  
  
  GCV = y.fitted.error[1]^2/(mean(diag(1 - fit$theta$S)))^2 #%*% C2
  
  da <- as.data.frame(cbind(da.Beta, da.theta, bias.theta, y.fitted.error))
  rownames(da) <- NULL
  da <- data.frame(meanH = H[X, 1], 
                   covHs = H[X, 2], 
                   covHt = H[X, 3], 
                   da, GCV = GCV, 
                   nIter = fit$nIter)
  
  return(da)
}


GCVparaSemi <- function(data, H = expand.grid(covH = seq(1E-2, 1,, 5),
                             meanH = seq(1E-2, 1,, 5)), 
                        Inde = T, method = 1, nIter = 10, 
                        nThreads = 5)
{
  
  cl <- parallel::makeCluster(getOption("cl.cores", nThreads))
  parallel::clusterEvalQ(cl, library(plyr))
  parallel::clusterEvalQ(cl, library(dplyr))
  parallel::clusterEvalQ(cl, library(data.table))
  parallel::clusterEvalQ(cl, library(Matrix))
  # parallel::clusterEvalQ(cl, library(semiBase))
  
  parallel::clusterExport(cl = cl
                          , varlist = c("data", "H", "Inde"
                                        , "method" ,"nIter"
                                        , "spT.validation")
                          , envir = environment())
  temp <- parallel::parLapply(cl, X = 1:nrow(H), fun = stGCVfun)
  # temp <- do.call("rbind", temp)
  temp <- data.table::rbindlist(temp)
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