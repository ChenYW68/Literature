source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
Nt <- 60
y = sqrt(t(Yts_Xts$Y_ts[1:Nt, ]))
time <- 0:(Nt - 1)#1:Nt(0:(Nt - 1))/(Nt - 1)
Nt <- ncol(y)
m <- rowMeans(y)
for(i in 1:nrow(y)){
  y[i, ] <- y[i, ] - m[i]
}

PCA = sp.FPCA(y_ts = y, time = time, nharm = 2, 
              lambda = 1e4)
cumsum(PCA$pca_fd$varprop)
PCA$pca_fd$values[1:3]
PCA$PCA[, 2] %*% PCA$PCA[, 2]
fda::matplot(time,PCA$PCA, xlab = 'time', ylab = 'PCs',
             lwd = 2, lty = 1, cex.lab = 1.5, cex.axis = 1.5,
             type = 'l')
# legend(0.5,-0.05,c('PC1','PC2'),col=1:3,lty=0.1,lwd=1)
title('Principle Component Functions')



y = (data$Y_ts)
Nt <- ncol(data$Y_ts)
m <- rowMeans(y)
for(i in 1:nrow(y)){
  y[i, ] <- y[i, ] - m[i]
}
time <- 0:(Nt - 1)
PCA = sp.FPCA(y_ts = y, time = time, #data$time
              nharm = 3, lambda = 1e2)
cumsum(PCA$pca_fd$varprop)
PCA$pca_fd$values[1:3]
PCA$PCA[, 1] %*% PCA$PCA[, 1]
fda::matplot(time,PCA$PCA, xlab = 'time', ylab = 'PCs',
             lwd = 2, lty = 1, cex.lab = 1.5, cex.axis = 1.5,
             type = 'l')
# legend(0.5,-0.05,c('PC1','PC2'),col=1:3,lty=0.1,lwd=1)
title('Principle Component Functions')