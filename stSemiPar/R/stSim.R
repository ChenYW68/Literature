rm(list=ls())
source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel1.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
# install.packages("E:/semiBase_1.0.tar.gz", repos = NULL, type = "source")

source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
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
n <- 100
Nt <- 20
M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
##############################################################
# G1: 0.05, 0.08; 0.1, 0.3; 0.5, 0.8; 1, 1.5; 0.5, 0.05
Phis <- c("0.5", "0.05")
sigma.sq.s <- c(2.5, 0.5)  #c(1, 0.5)
##############################################################
Kernel <- c(3, 0)
nIter <- 30
##############################################################
for(l in c(2, 4, 5)){
  # l=5
  Result <- NULL
  method <- M[l]
  ##############################################################
  # tab <- paste0(method, "I_", substr(sigma.sq.s[1], 3, 3),"_", 
  #               n, "_", Nt,S "_", substr(Phis[1], 3, 3))
  tab <- paste0(method, "_ED_", n, "_", Nt, "_", "CG_5")
  # , substr(Phis, 4, 4))
  alpha.est <- matrix(NA, nrow = 2*Nt, ncol = length(seed))
  # Max <- vector()
  # Min <- vector()
  ##############################################################
  ##############################################################
  for(iter in 1:length(seed))
  {
    # iter = 1
    set.seed(seed[iter])
    simDa <- siMuIncF(n = n, Nt = Nt, 
                      para = list(
                        Phis = as.numeric(Phis), 
                        nu = rep(0.5, length(sigma.sq.s)), 
                        sigma.sq.s = c(sigma.sq.s),
                        nugget = c(0.5), 
                        # Phit = as.numeric(Phit),
                        # rho = 0.1,
                        # tau.sq = 0, 
                        beta = c(1, 0),
                        J = 2))
    
    data <- list(Y_ts = simDa$Y_ts,
                 X_ts = simDa$X_ts,#rep(1, nrow(TrainData)),simDa$train.X
                 Z_ts = simDa$Z_ts,
                 loc = simDa$loc,
                 theta = simDa$theta,
                 time = simDa$time
                 , Qst = simDa$Qst
    )
    # assign("Q", as.matrix(data$Qst), envir = .GlobalEnv)
    H = unique(round(seq(5e-2, 2e-1,, 3), 3))
    
    H = as.matrix(expand.grid(H, H))
     # H = as.matrix(data.frame(H, H))
    library(profvis)
    # profvis({
    start_time <- proc.time()
    temp <- GCVparaSemi(data, H = H,
                        method = method,
                        Kernel = Kernel,
                        nuUnifb = 1, 
                        nu = 0.5, 
                        nIter = nIter,
                        nThreads = 10)
    # })
    # temp0 <- temp[which.min(temp$GCV), ]
    GCV <- NULL
    for(i in 1:nrow(H)){
      GCV <- c(GCV, temp[[(nrow(H) + i)]])
    }
    index <- which.min(GCV)
    temp0 <- round(temp[[index]], 3)
    # start_time <- Sys.time()
    # profvis({
    
    # for(k in 1:25){
    # k = 2
    # M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
    # method <- M[6]
    # if(method %nin% c("WLS")){
    #   # start_time <- Sys.time()
    #   # source("E:/Literature/semiBase/R/util.R")
    #   # source("./R/stSemiPar.R")
      # fit <- stSemiPar(y_ts = data$Y_ts,
      #                  x_ts = data$X_ts,
      #                  z_ts = data$Z_ts,
      #                  loc = data$loc,
      #                  time = data$time,#1e-1,#
      #                  Kernel = Kernel,
      #                  h = H[1, ],#as.matrix(temp0[, 1:2]),
      #                  method = method,
      #                  nThreads = 15,
      #                  nIter = nIter)
    #   # end_time <- Sys.time()
    #   # end_time - start_time
    #   # fit$beta
    # }else{
    #   # fit$theta$alpha
    #   fit <- stSemi_WLS(y_ts = data$Y_ts,
    #                     x_ts = data$X_ts,
    #                     z_ts = data$Z_ts,
    #                     loc = data$loc,
    #                     time = data$time,
    #                     prob = c(temp0$taper_s, prob[2]),
    #                     Kernel = Kernel,
    #                     h =  c(temp0$meanH, temp0$covHs, temp0$covHt),
    #                     nThreads = 10,
    #                     nIter = nIter)
    # }
    # 
    # })
    end_time <- proc.time()
    
    # fit$theta$alpha
    # fit$Theta$alpha
    {
      pdf(file = paste0("./figure", "/semiTemp", ".pdf"),
          width = 10, height = 10)
      par(mfrow = c(3, 3))
      for(i in 1:dim(data$Z_ts)[1]){
        # plot(simDa$time, simDa$theta[, i],
        #      ylim = c(min(simDa$theta[, i], fit$theta$alpha[, i]),
        #               max(simDa$theta[, i], fit$theta$alpha[, i])))
        # lines(simDa$time, fit$theta$alpha[, i], col = "red")
        
        plot(simDa$time, simDa$theta[, i],
             ylim = c(min(simDa$theta[, i], temp[[nrow(H)*3 + index]][, i]),
                      max(simDa$theta[, i], temp[[nrow(H)*3 + index]][, i])))
        lines(simDa$time, temp[[nrow(H)*3 + index]][, i], col = "red")
      }
      
      # plot(simDa$time, simDa$theta[, 2],
      #      ylim = c(min(simDa$theta[, 2], fit$theta$alpha[, 2]),
      #               max(simDa$theta[, 2], fit$theta$alpha[, 2])))
      # lines(simDa$time, fit$theta$alpha[, 2], col = "red")
      # plot(as.vector(fit$y_ts), as.vector(fit$fit.value),
      #      cex = 0.5, col = "black", pch = 20)
      # plot(as.vector(fit$y_ts)-as.vector(fit$fit.value),
      #      cex = 0.5, col = "black", pch = 20)
      
      
      plot(as.vector(temp[[nrow(H)*6 + index]]), as.vector(temp[[nrow(H)*5 + index]]),
           cex = 0.5, col = "black", pch = 20)
      plot(as.vector(temp[[nrow(H)*6 + index]])-as.vector(temp[[nrow(H)*5 + index]]),
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
        # d <- range((-simDa$D))
        # plot((-simDa$D), border = NA, 
        #      main = "Distance.matrix", 
        #      breaks= NULL,#(round(seq(d[1], d[2],, brk), 2)), 
        #      # breaks= seq(M1, M2,, bk), 
        #      col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
        #      key=list(side=4, cex.axis=size), fmt.key="%.2f")
        if((method %in% c("WEC_t", "WEC_st", "WEC_stw"))&(as.numeric(Phis[1])!=0)){
          
          plot(as.matrix(simDa$Cs), border = NA, main = "Cs.true",
               # breaks= seq(M1, M2,, bk),  #
               breaks=NULL,
               col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10),
               key=list(side=4, cex.axis=size), fmt.key="%.2f")
          plot(as.matrix(temp[[nrow(H)*2 + index]]), border = NA, main = "Cs.est",
               # breaks= seq(M1, M2,, bk),  #
               breaks=NULL,
               col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10),
               key=list(side=4, cex.axis=size), fmt.key="%.2f")
          
          plot(eigen(simDa$Cs)$value, eigen(temp[[nrow(H)*2 + index]])$value)
          
          # M1 <- min(min(simDa$Vt),  min(fit$Vt))
          # M2 <- max(max(simDa$Vt),  max(fit$Vt))
          # plot(as.matrix(simDa$Vt), border = NA, main = "Vs", breaks=NULL,
          #      # breaks= seq(M1, M2,, bk),
          #      col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10),
          #      key=list(side=4, cex.axis=size), fmt.key="%.2f")
          # 
          # plot(as.matrix(fit$Vt), border = NA, main = "Vs.est",
          #      # breaks= seq(M1, M2,, bk),  #
          #      breaks=NULL,
          #      col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10),
          #      key=list(side=4, cex.axis=size), fmt.key="%.2f")
          
          
          # plot(simDa$time, simDa$Vt[, 1], col = "black")
          # lines(simDa$time, simDa$Vt[, 2], col = "red")
          # lines(simDa$time, simDa$Vt[, 3], col = "green")
          matplot(simDa$time, as.matrix(simDa$Vt),xlab='time',ylab='PCs.True',
                  lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')
          
          matplot(simDa$time, as.matrix(temp[[nrow(H)*4 + index]]),xlab='time',ylab='PCs.Est',
                  lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')
          # legend(0.5,-0.05,c('PC1','PC2'),col=1:3,lty=0.1,lwd=1)
          title('Principle Component Functions')
          
          # image.plot(simDa$Vs, main = "Vs")
          # image.plot(fit$Cs$ModCov$Cov, main = "Vs.est")
          # image.plot(fit$Cs$Cov, main = "Vs.taper")
        }
        brk <- Nt*Nt
        if(method %in% c("WI1")){
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
      # cat("K = ", k, "H = ", as.vector(H[k, ]), "\n\n\n")
      # }
    }
    
    
    if(!is.null(data$theta)){
      for(i in 1:ncol(data$theta)){
        # alpha.est[((i - 1)*Nt + 1):((i)*Nt), iter] <- fit$theta$alpha[, i]
        alpha.est[((i - 1)*Nt + 1):((i)*Nt), iter] <- temp[[nrow(H)*3 + 1]][, i]
      }
    }
    temp0$Iter <- iter
    temp0$elapsed <- round((end_time - start_time)[[3]], 3)
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
    # print(round(Result[, c(1:3, 7, 10:14, 16:18)], 4))# c(1:5, 8, 10:14)])
    # cat("..................", tab, "..................\n")
    # print(round(colMeans(Result[, c(1:3, 7, 10:14, 16:18)]), 4))
    
    print(round(Result[, c(1:4, 8, 17:19)], 4))# c(1:5, 8, 10:14)])
    cat("..................", tab, "..................\n")
    print(round(colMeans(Result[, c(1:4, 8, 17:19)]), 4))
    
    cat("bias = ", mean(Result$Beta_1) - 1, "sd = ", sd(Result$Beta_1), "\n")
    if(iter == length(seed)){
      save(alpha.est, file = paste0("./data/", tab, "_alpha_est.RData"))
    }
  }
}


