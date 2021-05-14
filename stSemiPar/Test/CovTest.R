source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
#source("./R/spT.validation.R")

Y_ts <- siMuIncF(n = 100, Nt = 10, 
                 x.0 = c(0),
                 y.0 = c(0), 
                 delta = 0.1,
                 para = list(tau.sq = 0.001, 
                             Phis = 0.5, 
                             nu = 1, 
                             sigma.sq.s = 1,
                             sigma.sq.t = 1, 
                             Phit = 0.8,
                             rho = 0.1,
                             beta = c(0, 0)), 
                 nRatio = 0.8)

source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
Cov1 = semiCov(Y_ts$W_ts,
               Time = Y_ts$time,
               Kernel = 0,
               h = c(1e-1, 0.2),
               nuUnifb = 0,
               nu = 0,
               nThreads = 10)

Cov1$Cmat
# Cov$mCov
#
plot(diag(Y_ts$Vt), diag(Cov1$Cmat))
plot(as.vector(Y_ts$Vt), as.vector(Cov1$Cmat))
points(as.vector(Y_ts$Vt), as.vector(Cov1$mCov),
       col = "red")

source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
cov <- semiCovs(Y_ts$W_ts, 
                Y_ts$loc,     
                Kernel = 0,
                h = 0.1, 
                d = NULL,
                prob = 0.3,
                nuUnifb = 0,
                nu = 0,
                nThreads = 10)
cov

n <- nrow(Y_ts$W_ts)
Cmat <- matrix(0, nrow = n, ncol= n)
cv_index <- as.data.frame(which(upper.tri(Cmat), arr.ind = T))
setorderv(cv_index, c("row", "col"))

cv_index.0 <- as.matrix(cv_index)
Cmat[cv_index.0] <- cov$C


d <- fields::rdist( Y_ts$loc,  Y_ts$loc)
d <- c(as.vector(d[upper.tri(d)]))
thresh <- quantile(d, probs = 0.3)
index.upp <- which(d >= thresh)
index.low <- which(d < thresh)






# Y_ts$Vt
# eigen(Y_ts$Vt)$values
# eigen(Cov$Cmat)$values

# Cmat <- matrix(NA, nrow = Nt, ncol= Nt)
# cv_index <- as.data.frame(which(upper.tri(Cmat), arr.ind = T))
# setorderv(cv_index, c("row", "col"))
# cv_index <- as.matrix(cv_index)
# Cmat[cv_index] <- a$Cov$Cov

source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
# b = stSemiVary( y = Y_ts$Y_ts, Z = Y_ts$Z_ts[1,,], 
#              Time = 0:(ncol(Y_ts$W_ts) - 1),
#              Q = solve(Cov$Cmat),
#              Kernel = 0, h = 0.1, 
#              nThreads = 10)
source("./R/stSemiPar.R")
system.time({
fit <- stSemiPar(y_ts = Y_ts$Y_ts, 
                  x_ts = Y_ts$X_ts, 
                  z_ts = Y_ts$Z_ts,
                  time = Y_ts$time,
                  Inde = F, nIter = 1e1,
                  h = c(1e-1, 8e-2),
                  method = 2)
})
par(mfrow = c(3, 3))
plot(Y_ts$time,Y_ts$theta1, 
     ylim = c(min(Y_ts$theta1, fit$theta$alpha[, 1]), 
              max(Y_ts$theta1, fit$theta$alpha[, 1])))
lines(Y_ts$time, fit$theta$alpha[, 1], col = "red")
plot(Y_ts$time, Y_ts$theta2, 
     ylim = c(min(Y_ts$theta2, fit$theta$alpha[, 2]), 
              max(Y_ts$theta2, fit$theta$alpha[, 2])))
lines(Y_ts$time, fit$theta$alpha[, 2], col = "red")
plot(fit$y_ts, fit$fix.effect.fit + fit$theta$y.fit)

spT.validation(fit$y_ts, fit$fix.effect.fit + fit$theta$y.fit)
spT.validation(Y_ts$theta1, fit$theta$alpha[, 1])
spT.validation(Y_ts$theta2, fit$theta$alpha[, 2])
fit$beta


n <- nrow(Y_ts$Y_ts)
Nt <- ncol(Y_ts$Y_ts)
Da <- data.frame(y_ts = as.vector(((Y_ts$Y_ts[, 1:Nt]))),
                 x1 = as.vector((Y_ts$X_ts[1, , 1:Nt])),
                 x2 = as.vector((Y_ts$X_ts[2, , 1:Nt])),
                 z1 = as.vector((Y_ts$Z_ts[1, , 1:Nt])),
                 z2 = as.vector((Y_ts$Z_ts[2, , 1:Nt])),
                 t = rep(1:Nt, each = n)/Nt)
Fit <- mgcv::gam(y_ts ~ -1 + x1 + x2 + s(t, by = z1, k = 5)+ 
                   s(t, by = z2), data = Da)
plot(Fit)
plot(Da$y_ts, Fit$fitted.values)
spT.validation(Da$y_ts, Fit$fitted.values)
Fit$coefficients[1:2]