rm(list=ls())
source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel1.R")
source("./R/spT.validation.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
source("./R/stSemiPar_WLS.R")
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
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
# install.packages("E:/Literature/semiBase_1.0.zip", repos = NULL, type = "win.binary")
# library(semiBase)
DSN_01 <- odbcConnect(
  "DSN_01",
  uid = "myname",
  pwd = "mypwd",
  believeNRows = FALSE,
  case = "toupper"
)

##############################################################
seed <- 1:50
n <- 20
Nt <- 20
M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
##############################################################
Phis <- "0.2"
sigma.sq.s <- 0.5
##############################################################
prob = c(1.0e0, 1.5e0)
Kernel <- c(0, 0)
nIter = 20
set.seed(seed[39])
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
# library(ape)
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
              time = simDa$time
)
data <- Train
time <- rep(data$time, each = nrow(data$Y_ts))
fit <- mgcv::gam(as.vector(data$Y_ts)~as.vector(data$X_ts[1,,]) +
                   s(time, k =15), drop.intercept = T
                 # ,method="REML"
                 )
fit$coefficients[1]
dev.off()
plot(fit)
# dev.off()
# plot(as.vector(data$Y_ts),fit$fitted.values)
