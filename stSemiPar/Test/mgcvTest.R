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
