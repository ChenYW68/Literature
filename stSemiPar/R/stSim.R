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
n <- 50
Nt <- 20
M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
##############################################################
Phis <- "0.8"
sigma.sq.s <- 0.5
##############################################################
prob = c(1.0e0, 1.5e0)
Kernel <- c(0, 0)
nIter = 20
##############################################################
for(l in c(5)){
  # l=5
  Result <- NULL
  method <- M[l]
  ##############################################################
  
  tab <- paste0(method, "DD_", substr(sigma.sq.s, 3, 3),"_", 
                n, "_", Nt, "_", substr(Phis, 3, 3)) 
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
    H = expand.grid(meanH = seq(5e-2, 2e-1,, 10),
                    covHs = seq(1E-1, 1E0,, 5)
                     # covHs = seq(1E-2, 1E-1,, 5)
                    # covHt = seq(1E-3, 2E-1,, 5)
                     ,taper_s = prob[1]#seq(0.2, 1, , 3)
    )
    # H$covHs <- H$meanH
    H$covHt <- H$meanH
    H <- as.data.table(H)
    setcolorder(H, c("meanH", "covHs", "covHt", "taper_s"))
    
    # H$taper_s <- rep(prob[1], 5)#H$covHs
    H <- as.matrix(H)
    library(profvis)
    # profvis({
    temp <- GCVparaSemi(data, H = H,
                        prob = prob,
                        method = method,
                        Kernel = Kernel,
                        nuUnifb = 1, 
                        nu = 0.5, 
                        nIter = nIter,
                        nThreads = 15)
     # })
    temp0 <- temp[which.min(temp$GCV), ]
    # start_time <- Sys.time()
    # profvis({
    
    # for(k in 1:25){
    # k = 2
    # M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
    # method <- M[6]
    if(method %nin% c("WLS")){
      fit <- stSemiPar(y_ts = data$Y_ts,
                       x_ts = data$X_ts,
                       z_ts = data$Z_ts,
                       loc = data$loc,
                       time = data$time,#c(1e-1, 1e-1, 1e-1),#
                       Kernel = Kernel,
                       h = c(temp0$meanH, temp0$covHs, temp0$covHt),
                       prob = c(temp0$taper_s, prob[2]),
                       method = method,
                       nIter = nIter)
      # fit$beta
    }else{
      # fit$theta$alpha
      fit <- stSemi_WLS(y_ts = data$Y_ts,
                        x_ts = data$X_ts,
                        z_ts = data$Z_ts,
                        loc = data$loc,
                        time = data$time,
                        prob = c(temp0$taper_s, prob[2]),
                        Kernel = Kernel,
                        h =  c(temp0$meanH, temp0$covHs, temp0$covHt),
                        nThreads = 10,
                        nIter = nIter)
    }
    
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
        d <- range((-simDa$D))
        plot((-simDa$D), border = NA, 
             main = "Distance.matrix", 
             breaks= NULL,#(round(seq(d[1], d[2],, brk), 2)), 
             # breaks= seq(M1, M2,, bk), 
             col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
             key=list(side=4, cex.axis=size), fmt.key="%.2f")
        if((method %in% c("WEC_st", "WEC_stw", "WLS"))&(as.numeric(Phis)!=0)){
          
          # plot(simDa$Vc, border = NA, main = "Vst.true",
          #      # breaks= seq(M1, M2,, bk),  #
          #      breaks=NULL, 
          #      col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
          #      key=list(side=4, cex.axis=size), fmt.key="%.2f")
          # plot(fit$C, border = NA, main = "Vst.est",
          #      # breaks= seq(M1, M2,, bk),  #
          #      breaks=NULL, 
          #      col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
          #      key=list(side=4, cex.axis=size), fmt.key="%.2f")
          # 
          
          
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
      # cat("K = ", k, "H = ", as.vector(H[k, ]), "\n\n\n")
      # }
    }
    
    
    if(!is.null(data$theta)){
      for(i in 1:ncol(data$theta)){
        alpha.est[((i - 1)*Nt + 1):((i)*Nt), iter] <- fit$theta$alpha[, i]
      }
    }
    temp0$Iter <- iter
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
    print(round(Result[, c(1:6, 10, 13, 16)], 4))# c(1:5, 8, 10:14)])
    cat("..................", tab, "..................\n")
    print(round(colMeans(Result[, c(1:6, 10, 13, 15)]), 4))
    
    if(iter == length(seed)){
      save(alpha.est, file = paste0("./data/", tab, "_alpha_est.RData"))
    }
  }
}


