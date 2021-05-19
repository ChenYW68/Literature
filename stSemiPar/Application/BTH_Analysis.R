load("./data/BTH_Yts_Xts.Rdata")
load("./data/Site.Rdata")
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
# b = stSemiVary( y = Y_ts$Y_ts, Z = Y_ts$Z_ts[1,,], 
#              Time = 0:(ncol(Y_ts$W_ts) - 1),
#              Q = solve(Cov$Cmat),
#              Kernel = 0, h = 0.1, 
#              nThreads = 10)
source("./R/stSemiPar.R")
col <- c(3:7)
Nt <- 90
range(dist(Site[, 6:7]))
fit <- stSemi_OLS(y_ts = t(sqrt(Yts_Xts$Y_ts[1:Nt, ])), 
                 x_ts = Yts_Xts$X_ts[c(1, 2), , 1:Nt], 
                 z_ts = Yts_Xts$X_ts[col, , 1:Nt],
                 loc = as.matrix(Site[, 4:5]),
                 Vc = 1,  prob = c(1e-2, 1e-1),
                 time = (0:(Nt - 1))/(Nt - 1),
                 Inde = F, nIter = 5e0,
                 h = c(1e-1, 1e-1, 1e-1),
                 method = 1)
fit$beta
par(mfrow = c(2, 2))
image.plot(fit$Cs$Cov)
image.plot(fit$Cs$ModCov$Cov)
image.plot(fit$Ct$Cov)
image.plot(fit$Ct$ModCov$Cov)

Nc <- dimnames(Yts_Xts$X_ts[col, , 1:Nt])[[1]]
par(mfrow = c(3, 3))
plot(fit$y_ts^2, fit$fit.value^2)
for (k in 1:length(Nc)) {
  plot(fit$time, fit$Theta$alpha[, k], 
       col = "red", xlab = "time",
       ylab = Nc[k], type = "l")
}
# Fit.y <- fit$fix.effect.fit + fit$theta$y.fit
spT.validation(t(Yts_Xts$Y_ts[1:Nt, ]), fit$fit.value^2)
spT.validation(t(Yts_Xts$Y_ts[1:Nt, ]), fit$fit.value.2^2)


n <- ncol(Yts_Xts$Y_ts)
Da <- data.frame(y_ts = as.vector(t(sqrt(Yts_Xts$Y_ts[1:Nt, ]))),
                 x1 = as.vector((Yts_Xts$X_ts[2, , 1:Nt])),
                 z2 = as.vector((Yts_Xts$X_ts[3, , 1:Nt])),
                 z3 = as.vector((Yts_Xts$X_ts[4, , 1:Nt])),
                 z4 = as.vector((Yts_Xts$X_ts[5, , 1:Nt])),
                 z5 = as.vector((Yts_Xts$X_ts[6, , 1:Nt])),
                 z6 = as.vector((Yts_Xts$X_ts[7, , 1:Nt])),
                 t = rep(1:Nt, each = n)/Nt)
colnames(Da) <- c("y_ts", dimnames(Yts_Xts$X_ts)[[1]][2:7], "t")

Fit <- mgcv::gam(y_ts ~ 1 + CMAQ_PM25_30 + s(t, by = REAL_LON_WIND)+ 
            s(t, by = REAL_TEMP)+ s(t, by = REAL_PRES)+ 
            s(t, by = REAL_DEWP)+ s(t, by = REAL_LAT_WIND), data = Da)
plot(Fit)
Fit$coefficients[1:2]
plot(Da$y_ts^2, Fit$fitted.values^2)
spT.validation(Da$y_ts^2, Fit$fitted.values^2)
