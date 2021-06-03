source("./R/SetSimuModel.R")
source("./R/PSTVB_Packages.R")
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
n =50
Nt = 20
prob = c(1, 2)
h <- c(5e-2, 1e-1, 1e-1)
simDa <- siMuIncF(n = n, Nt = Nt,
                  x.0 = c(0),
                  y.0 = c(0),
                  delta = 0.1,
                  para = list(
                    Phis = as.numeric(0.8),
                    nu = c(0.5, 0.5),
                    sigma.sq.s = c(1e-1, 1e-1),
                    sigma.sq.t = c(0, 1),
                    # Phit = as.numeric(Phit),
                    # rho = 0.1,
                    # tau.sq = 0,
                    beta = c(1, 5)))

# Vs <- Matern(d = simDa$D, range = 0.5,
#              smoothness = 0.5,
#              phi = 1)
# # Vt <- diag(Nt)
# 
# Vt <- rdist(simDa$time , simDa$time)
# phi.t <- 0.5
# Vt <- phi.t^abs(row(Vt) - col(Vt))/(1 - phi.t^2)
# 
# ## Cross covariance
# Vc <- kronecker(Vs, Vt)
# # Time <- seq(0, 1,, Nt)
# 
# W_ts <- t(matrix(Matrix::crossprod(Matrix::chol(Vc),
#                                    rep(rnorm(n * Nt))),
#                  # + rnorm(n * Nt, 0, sd = sqrt(0.1))
#                  nrow = Nt, ncol = n))
# 
# simDa$Y_ts <- simDa$thetaF[1,,] + W_ts
# Q = solve(Vc)


Cs <- semiCovs(simDa$Y_ts - simDa$thetaF[1,,], 
               Coord = simDa$loc,     
               Kernel = 0,
               h = h[2], 
               prob = prob[1],
               nuUnifb = 1,
               nu = 0.8,
               nThreads = 10)
# h[1] <- 1E-3
Ct <- semiCovt(simDa$Y_ts - simDa$thetaF[1,,],
               Time = simDa$time, 
               Kernel = 0,
               h = c(h[3], h[3]),
               prob = prob[2],
               nuUnifb = 1,
               nu = 0.8,
               nThreads = 10)
# Cov.t <- Ct$ModCov$Cov  
# Matrix::diag(Ct$ModCov$Cov) <- Matrix::diag(Ct$ModCov$Cov)/ Cs$Var.s
Q <- kronecker(Matrix::solve(diag(n)),  #diag(n)
               Matrix::solve(Ct$ModCov$Cov))

Train <- list(Y_ts = simDa$Y_ts,
              X_ts = simDa$X_ts,#rep(1, nrow(TrainData)),simDa$train.X
              Z_ts = simDa$Z_ts,
              loc = simDa$loc,
              theta = simDa$theta,
              time = simDa$time
)
data <- Train
n <- nrow(data$Y_ts)
Nt <- ncol(data$Y_ts)
theta <- theta_WI(y_ts = data$Y_ts,
                   z_ts = data$Z_ts,
                   Time = data$time,
                   Q = diag(Nt), #solve(Y_ts$Vt)；diag(ncol(Y_ts$Y_ts))
                   h = h[1])
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
W1 = stWang_R(y_ts = data$Y_ts,
               z_ts = data$Z_ts,
               Time = data$time,
               Q = Q,
               theta_St = theta$St,
               Kernel = 0,
               h = h[1])

start_time <- Sys.time()
W2<- theta_stWang( y_ts = data$Y_ts,
                   z_ts = data$Z_ts,
                   Time = data$time,
                   Q = Q,
                   S = theta$St,
                   Kernel = 0,
                   h = h[1], 
                   nThreads = 10)
end_time <- Sys.time()
end_time - start_time

all.equal(W1$alpha[, 1], W2$alpha[, 1])
all.equal(W1$S, W2$S)
all.equal(W1$St, W2$St)
all.equal(W1$y.fit, W2$y.fit)

source("E:/Literature/semiBase/R/util.R")
start_time <- Sys.time()
wang <- theta_Wang_Space_Inde(y_ts = data$Y_ts,
                       z_ts = data$Z_ts,
                       Time = data$time,
                       Q = solve(Ct$ModCov$Cov),
                       S0 = theta$St,
                       Kernel = 0,
                       h = h[1])
end_time <- Sys.time()
end_time - start_time 
all.equal(W2$y.fit, wang$y.fit)

y.fit <- matrix(W2$S %*% as.vector(t(data$Y_ts)), 
                 nrow = n, ncol = Nt)
all.equal(as.vector((y.fit)), as.vector(wang$y.fit))
plot(as.vector(y.fit), as.vector((simDa$Y_ts)))
points(theta$S%*% as.vector(t(data$Y_ts)), 
       as.vector(simDa$Y_ts), col = "red")



M1 <- min(simDa$theta, W2$alpha[, 1], theta$alpha[, 1], wang$alpha[, 1])
M2 <- max(simDa$theta, W2$alpha[, 1], theta$alpha[, 1], wang$alpha[, 1])



# W$S[1:5, 1:5]
# Sk[1:5, 1:5]
all.equal(W2$alpha[, 1], wang$alpha[, 1])
all.equal(W2$S, wang$S)
all.equal(W2$St, wang$St)
all.equal(simDa$theta[, 1], theta$alpha[, 1])
all.equal(simDa$theta[, 1], wang$alpha[, 1])
all.equal(simDa$theta[, 1], W2$alpha[, 1])


source("E:/Literature/semiBase/R/util.R")
theta1 <- semPar.space.time( y_ts = simDa$Y_ts,
                            z_ts = simDa$Z_ts,
                            Time = simDa$time,
                            Q = Q,
                            S = theta$St,
                            Kernel = 0,
                            h = h[1],
                            nuUnifb = 0,
                            nu = 0,
                            nThreads = 1,
                            method = "WDt")
theta2 <- semPar.space.time( y_ts = simDa$Y_ts,
                            z_ts = simDa$Z_ts,
                            Time = simDa$time,
                            Q = Q,
                            S = theta$St,
                            Kernel = 0,
                            h = h[1],
                            nuUnifb = 0,
                            nu = 0,
                            nThreads = 1,
                            method = "WDstR")
all.equal(theta1$alpha[, 1], theta2$alpha[, 1])
source("./R/stSemiPar_WLS.R")
source("E:/Literature/semiBase/R/util.R")
wls <- stSemi_WLS(y_ts = simDa$Y_ts, 
                       x_ts = simDa$X_ts, 
                       z_ts = simDa$Z_ts,
                       time = simDa$time,
                       loc = simDa$loc,
                       prob = c(1, 1),
                       Kernel = c(0, 0),
                       h = c(1e-1, 1e-1, 1e-1),
                       nuUnifb = 1,
                       nu = 0.8,
                       nThreads = 10,
                       nIter = 10)
plot(wls$fit.value, simDa$Y_ts)

plot(simDa$time, simDa$theta, type = "o", ylim = c(M1, M2))
lines(simDa$time, W2$alpha[, 1], col = "blue", lwd = 2)
lines(simDa$time, theta$alpha[, 1], col = "red", lwd = 1.5)
lines(simDa$time, wang$alpha[, 1], col = "green", lwd = 2.0)
lines(simDa$time, wls$theta$alpha[, 1], col = "black", lwd = 2)
