


source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("./R/spT.validation.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
source("./R/stSemiPar_OLS.R")

# load("./data/testData.RData")
# source("E:/Literature/semiBase/R/util.R")
# Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
# W1 <- theta_stWang( y_ts = fix.residuals,
#                     z_ts = z_ts,
#                     Time = time,
#                     Q = Q, S = S,
#                     Kernel = Kernel,
#                     h = h, nuUnifb = 1,
#                     nu = 0.5,
#                     nThreads = 10)
# source("E:/Literature/semiBase/R/util.R")
# W2 = stWang(y_ts = fix.residuals,
#             z_ts = z_ts,
#             Time = time,
#             Q = Q,
#             theta_St = S,
#             h = h)
# 
# all.equal(W1$alpha[, 1], W2$alpha[, 1])
# all.equal(W1$S, W2$S)
# all.equal(W1$St, W2$St)
# all.equal(W1$y.fit, W2$y.fit)


n <- 100
Nt <- 20
Phis <- 0.8
sigma.sq.s <- 0.5
prob = c(1, 1e0)
simDa <- siMuIncF(n = n, Nt = Nt,
                  x.0 = c(0),
                  y.0 = c(0),
                  delta = 0.1,
                  para = list(
                    Phis = as.numeric(Phis),
                    nu = c(0.5, 0.5),
                    sigma.sq.s = c(sigma.sq.s),
                    sigma.sq.t = c(1),
                    # Phit = as.numeric(Phit),
                    # rho = 0.1,
                    # tau.sq = 0,
                    beta = c(1, 5)))
# range(simDa$Y_ts)
# range(simDa$Vs)
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

source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")

Train <- list(Y_ts = simDa$Y_ts,
              X_ts = simDa$X_ts,#rep(1, nrow(TrainData)),simDa$train.X
              Z_ts = simDa$Z_ts,
              loc = simDa$loc,
              theta = simDa$theta,
              time = simDa$time
)
data <- Train
# H = expand.grid(meanH = seq(min(simDa$time[-1])/2,
#                             max(simDa$time)/10,, 5),
#                 covHs = seq(5E-2, 3E-1,, 5)
#                 # covHt = seq(1E-3, 2E-1,, 5)
# )
# H$covHt <- H$meanH
# H <- as.matrix(H)
library(profvis)
# profvis({

# start_time <- Sys.time()
# profvis({

# for(k in 1:25){
# k = 2
M = c("WI", "WDt", "WDst", "WDstR")
method <- M[4] 
source("E:/Literature/semiBase/R/util.R")
fit <- stSemiPar(y_ts = data$Y_ts,
                 x_ts = data$X_ts,
                 z_ts = data$Z_ts,
                 loc = data$loc,
                 time = data$time,#c(1e-1, 1e-1, 1e-1),#
                 Kernel = 0,
                 h = c(1e-2, 1e-1, 1e-1),
                 prob = c(1, 2),
                 method = method,
                 nIter = 20)
# fit$beta
# fit$beta
# fit$theta$alpha
# fit1$theta$alpha



# fit <- stSemi_OLS(y_ts = data$Y_ts,
#                  x_ts = data$X_ts,
#                  z_ts = data$Z_ts,
#                  loc = data$loc,
#                  time = data$time,#c(0.1, 0.1),#
#                  h = c(1e-2, 1e-1, 1e-1),#c(temp0$meanH, temp0$covHs, temp0$covHt),#c(temp0$meanH, temp0$covHs, temp0$covHt), #c(temp0$covH, temp0$meanH),
#                  prob = c(1, 1),
#                  nIter = nIter)
# })
# end_time <- Sys.time()
# end_time - start_time
# fit$theta$alpha
# fit$Theta$alpha
{
  pdf(file = paste0("./figure", "/semiTemp", ".pdf"),
      width = 10, height = 10)
  par(mfrow = c(3, 3))
  for(i in 1:dim(data$Z_ts)[1]){
    plot(simDa$time, simDa$theta[, i],
         ylim = c(min(simDa$theta[, i], fit$theta$alpha[, i]),
                  max(simDa$theta[, i], fit$theta$alpha[, i])))
    lines(simDa$time, fit$theta$alpha[, i], col = "red")
  }

  # plot(simDa$time, simDa$theta[, 2],
  #      ylim = c(min(simDa$theta[, 2], fit$theta$alpha[, 2]),
  #               max(simDa$theta[, 2], fit$theta$alpha[, 2])))
  # lines(simDa$time, fit$theta$alpha[, 2], col = "red")
  plot(as.vector(fit$y_ts), as.vector(fit$fit.value),
       cex = 0.5, col = "black", pch = 20)
  plot(as.vector(fit$y_ts)-as.vector(fit$fit.value),
       cex = 0.5, col = "black", pch = 20)
  # dev.off()
  # 
  # pdf(file = paste0("./figure", "/Cov_st", ".pdf"),
  #     width = 8, height = 8)
  # par(mfrow = c(2, 2))
  {
    brk <- 10
    size = 0.8
    library(plot.matrix)
    if(method %in% c("WDst", "WDstR")){
      d <- range((-simDa$D))
      plot((-simDa$D), border = NA, 
           main = "Distance.matrix", 
           breaks= (round(seq(d[1], d[2],, brk), 2)), 
           # breaks= seq(M1, M2,, bk), 
           col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
           key=list(side=4, cex.axis=size), fmt.key="%.2f")
      
      M1 <- min(simDa$Vs, fit$Cs$ModCov$Cov, fit$Cs$Cov)
      M2 <- max(simDa$Vs, fit$Cs$ModCov$Cov, fit$Cs$Cov)
      plot(simDa$Vs, border = NA, main = "Vs", breaks=NULL, 
           # breaks= seq(M1, M2,, bk), 
           col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
           key=list(side=4, cex.axis=size), fmt.key="%.2f")
      
      plot(fit$Cs$ModCov$Cov, border = NA, main = "Vs.est",
           # breaks= seq(M1, M2,, bk),  #
           breaks=NULL, 
           col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
           key=list(side=4, cex.axis=size), fmt.key="%.2f")  
      plot(fit$Cs$Cov, border = NA, main = "Vs.taper", 
           # breaks= seq(M1, M2,, bk),  #
           breaks = NULL, 
           col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
           key=list(side=4, cex.axis=size), fmt.key="%.2f") 
      
      # image.plot(simDa$Vs, main = "Vs")
      # image.plot(fit$Cs$ModCov$Cov, main = "Vs.est")
      # image.plot(fit$Cs$Cov, main = "Vs.taper")
    }
    brk <- Nt*Nt
    if(method %nin% c("WI")){
      M1 <- min(simDa$Vt, fit$Ct$ModCov$Cov, fit$Ct$Cov)
      M2 <- max(simDa$Vt, fit$Ct$ModCov$Cov, fit$Ct$Cov)
      plot(simDa$Vt, border = NA, main = "Vt", 
           # breaks= seq(M1, M2,, bk),  #
           breaks=NULL, 
           col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
           key=list(side=4, cex.axis=size), 
           fmt.key="%.2f")
      plot(fit$Ct$ModCov$Cov, border = NA, main = "Vt.est",
           # breaks= seq(M1, M2,, bk),  #
           breaks = NULL,  
           col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
           key = list(side = 4, cex.axis=size), fmt.key="%.2f")
      # plot(fit$Ct$Cov, border = NA, main = "Vt.taper",
      #      # breaks= seq(M1, M2,, bk),  #
      #      breaks=NULL, 
      #      col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
      #      key=list(side=4, cex.axis=size), fmt.key="%.2f")
      # image.plot(simDa$Vt, main = "Vt")
      # image.plot(fit$Ct$ModCov$Cov, main = "Vt.est")
      # image.plot(fit$Ct$Cov, main = "Vt.taper")
    }
  }
  dev.off()
}
