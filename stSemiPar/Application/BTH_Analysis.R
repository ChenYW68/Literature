load("./data/BTH_Yts_Xts.Rdata")
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
# b = stSemiVary( y = Y_ts$Y_ts, Z = Y_ts$Z_ts[1,,], 
#              Time = 0:(ncol(Y_ts$W_ts) - 1),
#              Q = solve(Cov$Cmat),
#              Kernel = 0, h = 0.1, 
#              nThreads = 10)
source("./R/stSemiPar.R")
col <- c(3:7)
Nt <- 92
fit <- stSemiPar(y_ts = t(sqrt(Yts_Xts$Y_ts[1:Nt, ])), 
                 x_ts = Yts_Xts$X_ts[1:2, , 1:Nt], 
                 z_ts = Yts_Xts$X_ts[col, , 1:Nt],
                 time = (0:(Nt - 1))/(Nt - 1),
                 Inde = F, nIter = 5e0,
                 h = c(1e-1, 1e-1),
                 method = 1)
fit$beta
Nc <- dimnames(Yts_Xts$X_ts[col, , 1:Nt])[[1]]
par(mfrow = c(3, 3))
plot(fit$y_ts, fit$fix.effect.fit + fit$theta$y.fit)
for (k in 1:length(Nc)) {
  plot(fit$time, fit$theta$alpha[, k], 
       col = "red", xlab = "time",
       ylab = Nc[k], type = "l")
}
Fit.y <- fit$fix.effect.fit + fit$theta$y.fit
spT.validation(t(Yts_Xts$Y_ts[1:Nt, ]), Fit.y^2)

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
plot(Da$y_ts, Fit$fitted.values)
spT.validation(Da$y_ts^2, Fit$fitted.values^2)
