source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("./R/spT.validation.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
# Y_ts <- siMuIncF(n = 100, Nt = 10,
#                  x.0 = c(0),
#                  y.0 = c(0),
#                  delta = 0.1,
#                  para = list(tau.sq = 0.001, Phis = 0.3,
#                              nu = 1, sigma.sq.s = 1,
#                              sigma.sq.t = 1,
#                              Phit = 0.8,
#                              rho = 0.1,
#                              beta = c(1, 5)),
#                  nRatio = 0.8)
# range(Y_ts$D)
source("E:/Literature/semiBase/R/util.R")
# Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
DSN_01 <- odbcConnect(
                      "DSN_01",
                      uid = "myname",
                      pwd = "mypwd",
                      believeNRows = FALSE,
                      case = "toupper"
                    )
seed <- 1:50
Result <- NULL
Inde <- F
method <- 1
n <- 100
Nt <- 20
Phis <- 0.8
Phit <- 0.8
delta <- 0.1
tab <- paste0(Inde, "_", method, "_", n, "_", Nt, "_", 
              substr(Phis, 3, 3), "_",  
              substr(Phit, 3, 3)) 
alpha.est <- matrix(NA, nrow = 2*Nt, ncol = length(seed))
nIter = 20
Max <- vector()
Min <- vector()
for(iter in 1:length(seed))
{
  set.seed(seed[iter])
  simDa <- siMuIncF(n = n, Nt = Nt, 
                    x.0 = c(0),
                    y.0 = c(0), 
                    delta = delta,
                    para = list(tau.sq = 1, 
                                Phis = as.numeric(Phis), 
                                nu = 1, 
                                sigma.sq.s = 1,
                                sigma.sq.t = 1, 
                                Phit = as.numeric(Phit),
                                rho = 0.1,
                                beta = c(1, 5)))
  library(ape)
  # range(simDa$Y_ts)
  # M[iter] <- 0
  # w <- simDa$Vs
  # diag(w) <- 0
  # for(t in 1:Nt){
  #   M[iter] <-  M[iter] + max(eigen(w)$value)/min(eigen(w)$value)
  #     # Moran.I(simDa$W_ts[, t], 
  #     #                             (w),
  #     #       alternative = "g")$observed
  #   
  # }
  # M[iter] <- M[iter]/Nt
  # cat("iter = ", iter, "\n")
  # w = gpuR::eigen(simDa$Vc)$value
  # Max[iter] <- max(w)
  # Min[iter] <- min(w)
  
  
  
  Train <- list(Y_ts = simDa$Y_ts,
                X_ts = simDa$X_ts,#rep(1, nrow(TrainData)),simDa$train.X
                Z_ts = simDa$Z_ts,
                loc = simDa$loc,
                theta = simDa$theta,
                time = simDa$time,
                Vc = NULL#(simDa$Vc)
  )
  data <- Train
  H = expand.grid(meanH = seq(1E-2, 2E-1,, 5),
                  covHs = seq(1E-2, 2E-1,, 5)
                  # ,covHt = seq(1E-3, 2E-1,, 5)
                  )
  H$covHt <- H$covHs
  library(profvis)
  # profvis({
  temp <- GCVparaSemi(data, H = H,
                      Inde = Inde,
                      method = method,
                      nIter = nIter,
                      nThreads = 10)
  # })
  temp0 <- temp[which.min(temp$GCV), ]
  # start_time <- Sys.time()
  # profvis({
  fit <- stSemiPar(y_ts = data$Y_ts,
                   x_ts = data$X_ts,
                   z_ts = data$Z_ts,
                   loc = data$loc,
                   Vc = data$Vc,
                   time = data$time,#c(0.1, 0.1),#
                   h = c(temp0$meanH, temp0$covHs, temp0$covHt),#c(temp0$meanH, temp0$covHs, temp0$covHt), #c(temp0$covH, temp0$meanH),
                   Inde = Inde,
                   method = method,
                   nIter = nIter)
 # })
  # end_time <- Sys.time()
  # end_time - start_time
  
  pdf(file = paste0("./figure", "/semiTemp", ".pdf"),
      width = 10, height = 7)
  par(mfrow = c(2, 2))
  plot(simDa$time, simDa$theta[, 1],
       ylim = c(min(simDa$theta[, 1], fit$theta$alpha[, 1]),
                max(simDa$theta[, 1], fit$theta$alpha[, 1])))
  lines(simDa$time, fit$theta$alpha[, 1], col = "red")
  plot(simDa$time, simDa$theta[, 2],
       ylim = c(min(simDa$theta[, 2], fit$theta$alpha[, 2]),
                max(simDa$theta[, 2], fit$theta$alpha[, 2])))
  lines(simDa$time, fit$theta$alpha[, 2], col = "red")
  plot(fit$y_ts, fit$fix.effect.fit + fit$theta$y.fit,
       cex = 0.5, col = "black", pch = 20)
  dev.off()


  if(!is.null(data$theta)){
    for(i in 1:ncol(data$theta)){
      alpha.est[((i - 1)*Nt + 1):((i)*Nt), iter] <- fit$theta$alpha[, i]
    }
  }

  Result = rbind(Result, temp0)
  rownames(Result) <- NULL
  if (iter == 1) {
    sqlDrop(DSN_01, tab, errors = F)
  }
  sqlSave(
          DSN_01,
          as.data.frame(temp0),
          tab,
          append = TRUE,
          colnames = FALSE,
          rownames = FALSE,
          safer = TRUE,
          fast = TRUE
        )
  print(round(Result[, c(1:5, 8, 11:13)], 4))
  cat("..........................................\n")
  print(round(colMeans(Result[, c(1:5, 8, 11:13)]), 4))
  
}
save(alpha.est, file = paste0("./data/",
                              tab, "_alpha_est.RData"))


# system.time({
#   fit <- stSemiPar(y_ts = Y_ts$Y_ts,
#                    x_ts = Y_ts$X_ts,
#                    z_ts = Y_ts$Z_ts,
#                    time = Y_ts$time,
#                    Inde = F, nIter = 1e1,
#                    h = c(temp0$covH,
#                          temp0$meanH),
#                    method = 2)
# })
# par(mfrow = c(3, 3))
# plot(simDa$time, simDa$theta[, 1],
#      ylim = c(min(simDa$theta[, 1], fit$theta$alpha[, 1]),
#               max(simDa$theta[, 1], fit$theta$alpha[, 1])))
# lines(simDa$time, fit$theta$alpha[, 1], col = "red")
# plot(simDa$time, simDa$theta[, 2],
#      ylim = c(min(simDa$theta[, 2], fit$theta$alpha[, 2]),
#               max(simDa$theta[, 2], fit$theta$alpha[, 2])))
# lines(simDa$time, fit$theta$alpha[, 2], col = "red")
# plot(fit$y_ts, fit$fix.effect.fit + fit$theta$y.fit)
# 
# spT.validation(fit$y_ts, fit$fix.effect.fit + fit$theta$y.fit)
# spT.validation(simDa$theta[, 1], fit$theta$alpha[, 1])
# spT.validation(simDa$theta[, 2], fit$theta$alpha[, 2])
# fit$beta
mean(Max)
mean(Min)
