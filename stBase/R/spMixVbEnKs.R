spMixVBEnKs <- function(data, test = NULL, prior = NULL,
                        para = NULL,  true.para = NULL,
                        database = NULL, Object = "Total",
                        parallel = FALSE, verbose = TRUE,
                        verbose.VB = FALSE, heavy.tail = FALSE,
                        method = c("ensemble"),
                        response.transf = c("normal", "sr", "log"),
                        covariate.transf = c("normal", "sr", "log"),
                        Ensemble.size = 100,
                        ds = 1e-3, cs = 0.5,  ct = 1, df = 2,
                        IS.size = 200, Thresh = rep(1, 3),
                        Remove.Cores = 3,  N.Chunk = 1,
                        itMax = 1e2, tol.vb = 1e-2, 
                        seed = NULL,
                        tol.real = 1e-3)
{
  if(!is.null(seed)){set.seed(seed)}
  if (is.null(data)) {stop("Must provide data.\n")}
  ######################################################################
  ######################################################################
  if(response.transf == "sr"){ data$Y_ts = sqrt(data$Y_ts)}
  if(response.transf == "log"){data$Y_ts = log(data$Y_ts)}
  if(covariate.transf == "sr"){
    data$X_ts[-1,,] = sqrt(data$X_ts[-1,,])
    if(!is.null(test)){
      test$X_ts[-1,,] = sqrt(test$X_ts[-1,,])
    }
    if(!is.null(data$Z_ts)){
      data$Z_ts[-1,,] = -sqrt(data$Z_ts[-1,,])}
  }
  if(covariate.transf == "log"){
    data$X_ts[-1,,] = log(data$X_ts[-1,,])
    if(!is.null(test)){
      test$X_ts[-1,,] = log(test$X_ts[-1,,])
    }
  }
  if(covariate.transf == "sc")
  {
    for(k in 2:dim(data$X_ts)[1])
    {
      data$X_ts[k,,] = matrix(scale(as.vector(data$X_ts[k,,]))[, 1],
                              nrow = nrow(data$X_ts[k,,]), ncol = data$Nt)
      if(!is.null(test)){
        test$X_ts[k,,] =  matrix(scale(as.vector(test$X_ts[k,,]))[, 1],
                                 nrow = nrow(test$X_ts[k,,]), ncol = data$Nt)}
      if(!is.null(data$Z_ts))
      {
        data$Z_ts[k,,] = matrix(scale(as.vector(data$Z_ts[k,,]))[, 1],
                                nrow = nrow(data$Z_ts[k,,]), ncol = data$Nt)
      }
    }
    
  }
  ######################################################################
  IS <- list(theta2 = list(sample = NULL, weight = NULL),
             K = list(sample = NULL, weight = NULL),
             K0 = list(sample = NULL, weight = NULL),
             Thresh = Thresh)
  ######################################################################
  ######################################################################
  criterion = TRUE;iter = 0;initial.para <- para; RowName <- vector()
  ###########################################################
  if(heavy.tail&(length(para$Obs.tau2$E_tau2)==1))
  {
    para$Obs.tau2$E_tau2 = rep(para$Obs.tau2$E_tau2, data$n)
  }
  if(heavy.tail & is.null(para$obs_tau2)){
    para$obs_tau2 = para$Obs.tau2$E_tau2
  }
  if(is.null(para$Q$E_Q)){
    para$Q$E_Q <-  data$Adj.Mat + para$k$E_k^2 *  base::diag(data$N.BAUs)
    para$Q0$E_Q0 <- data$Adj.Mat + para$k0$E_k0^2 * base::diag(data$N.BAUs)
  }
  
  # detect missing data  ----------------------------------------------------
  ###########################################################
  para$v0$E_v0 <- rep(0, data$N.BAUs)
  if(!is.null(test)){
    test.data <- test$Y_ts_true
    test.index <- which(is.na(test$Y_ts_true), arr.ind = T)}
  # test$Y_ts[test.index] <- 0
  
  miss.index <- which(is.na(data$Y_ts), arr.ind = T)
  if(nrow(miss.index) >0)
  {
    cat("There are some missing data...\n")
    X_ts = sapply(seq_len(nrow(miss.index)), function(s) 
      (data$X_ts[, miss.index[s, 2] , miss.index[s, 1]]) %*% 
        para$beta$E_betaX
      , simplify = "matrix")
    if(!is.null(data$Z_ts)){
      X_ts = sapply(seq_len(nrow(miss.index)), function(s)
        data$Hs[miss.index[s, 2], ] %*% t(data$Z_ts[, , miss.index[s, 1]]) %*% 
          para$beta$E_betaZ, simplify="array") + X_ts
    }
    Fill_Miss_Value <- data$Y_ts[miss.index] <- X_ts
  }
  spTaper <- fields::Wendland(data$BAUs.Dist
                              , theta = max(data$BAUs.Dist)*cs
                              , dimension = 1, k = 1)
  bandKernel <- spTaper
  # bandKernel[bandKernel > 0] <- 1
  para.ks = IniEnKs(data = data, para = para,
                    Ensemble.size = Ensemble.size,
                    heavy.tail = heavy.tail,
                    spTaper = spTaper, 
                    bandKernel = bandKernel,
                    ds = ds, ct = ct)
  sp = para.ks$sp
  Ks <- TaEnKs(para.ks)
  
  rm(para.ks)
  
  log.lik0 <- loglik(data = data, para = para, Ks = Ks,
                     ds = ds, sp = sp, heavy.tail = heavy.tail)
  
  ######################################################################
  while (criterion & iter <= itMax)
  {
    S <- covSfun.Gaussian(data$Nt, Ks = Ks)
    T0 <- proc.time()
    Re.spVB <- spVB.Gaussian(data, Ks = Ks, S = S,
                             prior = prior,  para = para,
                             IS = IS, parallel = parallel,
                             sp = sp, ds = ds, VB.err = tol.vb,
                             IS.size = IS.size, N.Chunk = N.Chunk,
                             Remove.Cores = Remove.Cores,
                             verbose.VB = verbose.VB,
                             iter = iter)
    PIU <- Re.spVB$para
    PIU$v0$E_v0 <- Ks$Xf[1, ]
    IS <- Re.spVB$IS
    t1 <- proc.time()
    Ks <- TaEnKs(para.ks = IniEnKs(data = data, para = PIU,
                                   Ensemble.size = Ensemble.size,
                                   heavy.tail = heavy.tail,
                                   spTaper = spTaper, 
                                   bandKernel = bandKernel,
                                   ds = ds, ct = ct))
    t2 <- proc.time()
    rm(S)
    if(verbose){
      cat("TaEnKs algorithm takes time: \n")
      print(t2 - t1)
      cat("................................................................. \n")
    }
    
    #TRAIN DATA
    if(nrow(miss.index)>0)
    {
      W_ts = sapply(seq_len(nrow(miss.index)), function(s)
        data$Hs[miss.index[s, 2], ] %*% 
          as.vector(t(Ks$Xf[miss.index[s, 1], ]))
        , simplify="array")
      
      X_ts = sapply(seq_len(nrow(miss.index)), function(s) 
        (data$X_ts[, miss.index[s, 2] , miss.index[s, 1]]) %*% 
          PIU$beta$E_betaX
        , simplify = "matrix")
      if(!is.null(data$Z_ts)){
        X_ts = sapply(seq_len(nrow(miss.index)), function(s)
          data$Hs[miss.index[s, 2], ] %*% t(data$Z_ts[, , miss.index[s, 1]]) %*% 
            para$beta$E_betaZ
          , simplify="array") + X_ts
      }
      Fill_Miss_Value <- data$Y_ts[miss.index] <- X_ts + W_ts
    }else{
      Fill_Miss_Value = NULL;
    }
    if(method == c("ensemble"))
    {
      n.Enseble = dim(Ks$EnXf)[3]
      Y.fit <- array(0, dim = c(data$Nt, data$n, n.Enseble),
                     dimnames = list(base::rownames(data$Y_ts),
                                     base::colnames(data$Y_ts),
                                     paste0("En.", 1:n.Enseble)
                     ))
      for(t0 in 1:data$Nt)
      {
        W_ts = apply(PIU$alpha$E_alpha*data$Hs, 1, FUN = '%*%', Ks$EnXf[t0 + 1, ,])
        X_ts <- X_ts_Transf(data$Nt, data$X_ts, PIU$beta$E_betaX)[t0, ] 
        if(!is.null(data$Z_ts))
        {
          X_ts <-  X_ts +
            PIU$alpha$E_alpha*(X_ts_Transf(data$Nt, data$Z_ts, 
                                           PIU$beta$E_betaZ) %*% t(data$Hs))[t0, ]
        }
        X <- matrix(X_ts, nrow = data$n, ncol = n.Enseble)
        Y.fit[t0, ,]  <- X + t(W_ts)
        # if(heavy.tail)
        # {
        #   # tau2 <- rgamma(data$n*n.Enseble, shape = PIU$Obs.tau2$a, rate = PIU$Obs.tau2$b)
        #   Y.fit[t0, ,]  <- X + t(W_ts)# +
        #     # matrix(as.vector(spCalibration::rmvn(1,  mu = rep(0, data$n*n.Enseble),
        #     #        L = sp$linalg$cholesky(diag(1/tau2)))),
        #     #        nrow = data$n, ncol = n.Enseble)
        # }else{
        #   Y.fit[t0, ,]  <- X + t(W_ts) #+ matrix(rnorm(data$n*n.Enseble,
        #                                          #       0, sqrt(PIU$Obs.tau2$E_tau2)),
        #                                          # nrow = data$n,
        #                                          # ncol = n.Enseble)
        # }
      }
      Y.fit.mean <- apply(Y.fit, c(1, 2), mean)
      if(response.transf == "sr") { 
        spT <- spT.validation(data$Y_ts_true, Y.fit.mean^2, Y.fit^2, F)
      }else if(response.transf == "log"){
        spT <- spT.validation(data$Y_ts_true,
                              exp(Y.fit.mean), exp(Y.fit), F)
      }else{spT <- spT.validation(data$Y_ts_true,
                                  (Y.fit.mean), (Y.fit), F)}
      
      Y.fit <- Y.fit.mean
    }else{
      Y.fit <- matrix(NA, nrow =  data$Nt, ncol = data$n)
      for(t0 in 1:data$Nt)
      {
        W_ts <- PIU$alpha$E_alpha*data$Hs %*% Ks$Xf[t0 + 1, ]
        X_ts <- X_ts_Transf(data$Nt, data$X_ts, PIU$beta$E_betaX)[t0, ]
        if(!is.null(data$Z_ts))
        {
          X_ts <-  X_ts +
            PIU$alpha$E_alpha*(X_ts_Transf(data$Nt, data$Z_ts, 
                                           PIU$beta$E_betaZ) %*% t(data$Hs))[t0, ]
        }
        Y.fit[t0, ] <- X_ts + W_ts
      }
      spT <- spT.validation(data$Y_ts_true, Y.fit, NULL, F)
    }
    if(response.transf == "sr") { Y.fit <- Y.fit^2
    }else if(response.transf == "log"){ Y.fit <- exp(Y.fit)
    }else{Y.fit <- Y.fit}
    
    print(spT)
    CRPS1 <- spT[16] %>% as.vector() %>% as.numeric()
    ES1 <- spT[17] %>% as.vector() %>% as.numeric()
    
    y_ts <- data$Y_ts_true
    
    Er_da <- as.vector(Y.fit - y_ts)
    RMSE1 <- round(sqrt(mean(Er_da^2, na.rm = T)), 3)
    
    MB1   = round(mean(Er_da, na.rm = T), 3)
    NMB1   = round(sum(Er_da, na.rm = T)*100/sum(y_ts, na.rm = T), 3)
    NME1  = round(sum(abs(Er_da), na.rm = T)*100/sum(y_ts, na.rm = T), 3)
    if(verbose){
      # cat("\n ........................\n")
      cat(" Fitting RMSE = ", RMSE1, "; \n"
          , "Mean Bias (MB) = ", MB1, "; \n"
          , "Normalized Mean Bias(%) = ", NMB1,"; \n"
          , "Normalized Mean Error(%) = ", NME1,"; \n")
      cat("................................................................. \n")
    }
    if(!is.null(test))
    {
      N <- ncol(test$Y_ts_true)
      # Nt <- nrow(test$Y_ts)
      if(method == c("ensemble"))
      {
        n.Enseble = dim(Ks$EnXf)[3]
        Y.test <- array(0, dim = c(data$Nt, N, n.Enseble),
                        dimnames = list(base::rownames(test$Y_ts),
                                        base::colnames(test$Y_ts),
                                        paste0("En.", 1:n.Enseble)
                        ))
        for(t0 in 1:data$Nt)
        {
          W_ts = apply(PIU$alpha$E_alpha*test$H, 1, FUN = '%*%', Ks$EnXf[t0 + 1, ,])
          # X_ts <- X_ts_Transf(data$Nt, test$X_ts, PIU$beta$E_betaX)[t0, ]
          if(!is.null(data$Z_ts))
          {
            X_ts <-  X_ts +
              PIU$alpha$E_alpha*(X_ts_Transf(data$Nt, data$Z_ts, 
                                             PIU$beta$E_betaZ) %*% t(test$H))[t0, ]
          }
          # X <- matrix(X_ts, nrow = N, ncol = n.Enseble)
          # Y.test[t0, ,]  <- X + t(W_ts)
          
          beta.sample <- mvnfast::rmvn(n.Enseble,
                                       mu = PIU$beta$E_betaX, 
                                       sigma = PIU$beta$betaX.Sigma2) %>% t()
          
          Xbeta_ts = apply(t(test$X_ts[,, t0]), 1, FUN = '%*%', beta.sample)
          
          Y.test[t0, ,]  <-  t(Xbeta_ts + W_ts)+ rnorm(N*n.Enseble, 0, 
                                                       sqrt(PIU$Obs.tau2$E_tau2))
        }
        # Y.test.mean <- apply(Y.test, c(1, 2), mean)
        Y.test.mean <- apply(Y.test, c(1, 2), quantile, prob = 0.5)
        
        
        if(response.transf == "sr") { 
          spT <- spT.validation(test$Y_ts_true, Y.test.mean^2, Y.test^2, F)
        }else if(response.transf == "log"){
          spT <- spT.validation(test$Y_ts_true,
                                exp(Y.test.mean), exp(Y.test), F)
        }else{spT <- spT.validation(test$Y_ts_true,
                                    (Y.test.mean), (Y.test), F)}
        
        Y.test <- Y.test.mean
      }else{
        Y.test <- matrix(NA, ncol = N , nrow =  data$Nt)
        for(t0 in 1:data$Nt)
        {
          W_ts <- PIU$alpha$E_alpha*test$H %*% Ks$Xf[t0 + 1, ]
          X_ts <- X_ts_Transf(data$Nt, test$X_ts, PIU$beta$E_betaX)[t0, ]
          if(!is.null(data$Z_ts))
          {
            X_ts <-  X_ts +
              PIU$alpha$E_alpha*(X_ts_Transf(data$Nt,  data$Z_ts, 
                                             PIU$beta$E_betaZ) %*% t(test$H))[t0, ]
          }
          Y.test[t0, ] <- X_ts + W_ts + rnorm(N, 0, sqrt(PIU$Obs.tau2$E_tau2))
        }
        spT <- spT.validation(test$Y_ts_true, Y.test, NULL, F)
      }
      if(response.transf == "sr"){ Y.test <- Y.test^2}
      if(response.transf == "log"){ Y.test <- exp(Y.test)}
      if(covariate.transf == "sr"){x_ts <-(t(test$X_ts[2,,]))^2}
      if(covariate.transf == "log"){x_ts <- exp(t(test$X_ts[2,,]))}
      if(covariate.transf == "normal"){x_ts <-(t(test$X_ts[2,,]))}
      # if(covariate.transf == "sc"){x_ts <-(t(test$X_ts[2,,]))^2}
      
      print(spT)
      CRPS2 <- spT[15] %>% as.vector() %>% as.numeric()
      ES2 <- spT[16] %>% as.vector() %>% as.numeric()
      
      
      Er_da2 <- Y.test - test$Y_ts_true
      
      test.yts <- test$Y_ts_true %>% as.vector()
      index <- !is.na(test.yts)
      x_ts <- as.vector(x_ts)
      rho0 <-  round(mean(cor(test.yts[index], x_ts[index])), 3)
      rho <- round(mean(cor(Y.test[index], test.yts[index])), 3)
      
      Er_da <- as.vector(Er_da2)
      RMSE2 <- round(sqrt(mean(Er_da^2, na.rm = T)), 3)
      MB2   <- round(mean(Er_da, na.rm = T), 3)
      # NMB2   = round(mean(Er_da/test$Y_ts), 3)
      NMB2 <-  round(sum(Er_da, na.rm = T)*100/sum(test$Y_ts_true, na.rm = T), 3)  #100%
      NME2 <-  round(sum(abs(Er_da), na.rm = T)*100/sum(test$Y_ts_true, na.rm = T), 3)
      if(verbose){cat("\n Object:", Object, "\n\n")
        cat(" Testing RMSE = ", RMSE2, "; \n"
            , "Mean Bias (MB) = ", MB2, "; \n"
            , "Normalized Mean Bias(%) = ", NMB2, "; \n"
            , "Normalized Mean Error(%) = ", NME2,"; \n")
        cat("\n Pearson-Before:", rho0, "; Pearson-After:", rho, "\n")
      }
      # print(round(sqrt(colMeans(Er_da2^2)),3))
    }else{rho0 = rho = RMSE2 = MB2 = NMB2 = NME2 = CRPS2 = ES2 = NA}
    # if(verbose){cat("\n ........................\n")}
    # if(verbose){cat("covSfun function runing time: \n")
    # print(t2 - t1)}
    
    t3 <- proc.time()
    if(verbose){
      cat("\n................................................................. \n")
      cat("\nEach iteration takes time: \n")
      print(t3 - T0)}
    log.lik <- loglik(data = data, para = PIU, Ks = Ks,
                      ds = ds, sp = sp, heavy.tail = heavy.tail)
    
    # if(verbose){ cat("................................................................. \n")}
    
    log.lik.error <- abs((log.lik - log.lik0)/log.lik0)
    
    cat("\n Iteration:", iter + 1)
    cat("\n loglik:", log.lik.error,"\n\n")
    log.lik0 <- log.lik
    {
      if(is.null(data$Z_ts)){
        if(heavy.tail){
          result <- rbind(c(true.para$betaX, true.para$alpha,
                            true.para$a, true.para$b, 
                            true.para$theta[1], true.para$theta[2], 
                            true.para$k, true.para$Proc.tau2,
                            true.para$k0, true.para$Proc0.tau2, 0),
                          c(initial.para$beta$E_betaX, 
                            initial.para$alpha$E_alpha,
                            initial.para$Obs.tau2$a,
                            initial.para$Obs.tau2$b,
                            initial.para$theta1$E_theta1,
                            initial.para$theta2$E_theta2,
                            initial.para$k$E_k, 
                            initial.para$Proc.tau2$E_tau2, 
                            initial.para$k0$E_k0, 
                            initial.para$Proc0.tau2$E_tau2, 
                            0),
                          c(PIU$beta$E_betaX, PIU$alpha$E_alpha,
                            prior$Obs.tau2$a, prior$Obs.tau2$b, 
                            PIU$theta1$E_theta1,
                            PIU$theta2$E_theta2, PIU$k$E_k,
                            PIU$Proc.tau2$E_tau2,
                            PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                            iter + 1))
          colnames(result) <- c(dimnames(data$X_ts)[[1]],
                                "alpha", "a_tau2.prior", "b_tau2.prior",
                                "theta1", "theta2","K", "Proc_tau2" , "K0", 
                                "Proc0_tau2",  "iter")
        }else{
          result <- rbind(c(true.para$betaX, true.para$alpha,
                            true.para$Obs.tau2, 
                            true.para$theta[1], true.para$theta[2], 
                            true.para$k, true.para$Proc.tau2,
                            true.para$k0, true.para$Proc0.tau2, 0),
                          c(initial.para$beta$E_betaX, 
                            initial.para$alpha$E_alpha,
                            initial.para$Obs.tau2$E_tau2,
                            initial.para$theta1$E_theta1,
                            initial.para$theta2$E_theta2,
                            initial.para$k$E_k, 
                            initial.para$Proc.tau2$E_tau2, 
                            initial.para$k0$E_k0, 
                            initial.para$Proc0.tau2$E_tau2, 
                            0),
                          c(PIU$beta$E_betaX, PIU$alpha$E_alpha,
                            PIU$Obs.tau2$E_tau2, 
                            PIU$theta1$E_theta1,
                            PIU$theta2$E_theta2, PIU$k$E_k,
                            PIU$Proc.tau2$E_tau2,
                            PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                            iter + 1))
          colnames(result) <- c(dimnames(data$X_ts)[[1]],
                                "alpha", "E_tau2", "theta1", "theta2","K",
                                "Proc_tau2" , "K0", "Proc0_tau2",  "iter")
        }
      }else{
        if(heavy.tail){
          result <- rbind(c(true.para$betaX,  true.para$betaZ, 
                            true.para$alpha, true.para$a, 
                            true.para$b, true.para$theta[1],
                            true.para$theta[2], 
                            true.para$k, true.para$Proc.tau2,
                            true.para$k0, true.para$Proc0.tau2, 0),
                          c(initial.para$beta$E_betaX, 
                            initial.para$beta$E_betaZ, 
                            initial.para$alpha$E_alpha,
                            initial.para$Obs.tau2$a,
                            initial.para$Obs.tau2$b,
                            initial.para$theta1$E_theta1,
                            initial.para$theta2$E_theta2,
                            initial.para$k$E_k, 
                            initial.para$Proc.tau2$E_tau2, 
                            initial.para$k0$E_k0, 
                            initial.para$Proc0.tau2$E_tau2,  
                            0),
                          c(PIU$beta$E_betaX, PIU$beta$E_betaZ,
                            PIU$alpha$E_alpha,
                            prior$Obs.tau2$a, prior$Obs.tau2$b,
                            PIU$theta1$E_theta1,
                            PIU$theta2$E_theta2, 
                            PIU$k$E_k, PIU$Proc.tau2$E_tau2,
                            PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2, 
                            iter + 1))
          colnames(result) <- c(dimnames(data$X_ts)[[1]],
                                paste0("betaZ", 0:(length(PIU$beta$E_betaZ) - 1)),
                                "alpha", "a_tau2.prior", "b_tau2.prior", 
                                "theta1", "theta2", "K",
                                "Proc_tau2" , "K0", "Proc0_tau2", "iter")
        }else{
          result <- rbind(c(true.para$betaX,  true.para$betaZ, 
                            true.para$alpha, true.para$Obs.tau2,    
                            true.para$theta[1],
                            true.para$theta[2], 
                            true.para$k, true.para$Proc.tau2,
                            true.para$k0, true.para$Proc0.tau2, 0),
                          c(initial.para$beta$E_betaX, 
                            initial.para$beta$E_betaZ, 
                            initial.para$alpha$E_alpha,
                            initial.para$Obs.tau2$E_tau2,
                            initial.para$theta1$E_theta1,
                            initial.para$theta2$E_theta2,
                            initial.para$k$E_k, 
                            initial.para$Proc.tau2$E_tau2, 
                            initial.para$k0$E_k0, 
                            initial.para$Proc0.tau2$E_tau2,  
                            0),
                          c(PIU$beta$E_betaX, PIU$beta$E_betaZ,
                            PIU$alpha$E_alpha,
                            PIU$Obs.tau2$E_tau2, 
                            PIU$theta1$E_theta1,
                            PIU$theta2$E_theta2, 
                            PIU$k$E_k,
                            PIU$Proc.tau2$E_tau2,
                            PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2, 
                            iter + 1))
          colnames(result) <- c(dimnames(data$X_ts)[[1]],
                                paste0("betaZ", 0:(length(PIU$beta$E_betaZ) - 1)),
                                "alpha", "E_tau2", "theta1", "theta2", "K",
                                "Proc_tau2" , "K0", "Proc0_tau2", "iter")
        }
      }
      rownames(result) <- c("True:", "Init:", "Esti:")
    }
    if(verbose){
      print(round(result, 3))}
    para <- PIU
    ###########################################################
    # criterion <- (error > tol.real)
    criterion <- (log.lik.error > tol.real)
    
    ###########################################################
    {
      if(iter == 0)
      {
        criterion <- T
        if(is.null(data$Z_ts)){
          if(heavy.tail){
            temp <- cbind(rbind(c(true.para$betaX, true.para$alpha,
                                  true.para$a, true.para$b,
                                  true.para$theta, true.para$k,
                                  true.para$Proc.tau2,
                                  true.para$k0, true.para$Proc0.tau2, 
                                  rep(1, 16), iter = 0),
                                c(initial.para$beta$E_betaX,
                                  initial.para$alpha$E_alpha,
                                  initial.para$Obs.tau2$a,
                                  initial.para$Obs.tau2$b,
                                  initial.para$theta1$E_theta1,
                                  initial.para$theta2$E_theta2,
                                  initial.para$k$E_k, 
                                  initial.para$Proc.tau2$E_tau2, 
                                  initial.para$k0$E_k0, 
                                  initial.para$Proc0.tau2$E_tau2,
                                  rep(1, 16), iter = 0),
                                c(PIU$beta$E_betaX, PIU$alpha$E_alpha,
                                  prior$Obs.tau2$a, prior$Obs.tau2$b,
                                  PIU$theta1$E_theta1,
                                  PIU$theta2$E_theta2, 
                                  PIU$k$E_k, PIU$Proc.tau2$E_tau2,
                                  PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                  log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                  NMB1, NMB2, NME1, NME2, rho0, rho,
                                  CRPS1, CRPS2, ES1, ES2,
                                  round((t3 - T0)[[3]]), iter + 1)) %>% as.data.frame(),
                          Object)
            colnames(temp) <- c(dimnames(data$X_ts)[[1]],
                                "alpha", "a_tau2_prior", "b_tau2_prior", paste0("theta", 1:2),
                                "K", "Proc_tau2", "K0", "Proc0_tau2", "Iter_loglik_Error", "Fitting_RMSE",
                                "Testing_RMSE", "Fitting_MB", "Testing_MB",
                                "Fitting_NMB", "Testing_NMB", "Fitting_NME",
                                "Testing_NME", "Corr_Bef", "Corr_Aft",
                                "Fitting_CRPS", "Testing_CRPS", 
                                "Fitting_ES", "Testing_ES",
                                "elapsed", "iter", "Object")
          }else{
            temp <- cbind(rbind(c(true.para$betaX, true.para$alpha,
                                  true.para$Obs.tau2, 
                                  true.para$theta,
                                  true.para$k,
                                  true.para$Proc.tau2,
                                  true.para$k0,
                                  true.para$Proc0.tau2, 
                                  rep(1, 16), iter = 0),
                                c(initial.para$beta$E_betaX,
                                  initial.para$alpha$E_alpha,
                                  initial.para$Obs.tau2$E_tau2,
                                  initial.para$theta1$E_theta1,
                                  initial.para$theta2$E_theta2,
                                  initial.para$k$E_k, 
                                  initial.para$Proc.tau2$E_tau2, 
                                  initial.para$k0$E_k0, 
                                  initial.para$Proc0.tau2$E_tau2,
                                  rep(1, 16), iter = 0),
                                c(PIU$beta$E_betaX, PIU$alpha$E_alpha,
                                  PIU$Obs.tau2$E_tau2, PIU$theta1$E_theta1,
                                  PIU$theta2$E_theta2, 
                                  PIU$k$E_k, PIU$Proc.tau2$E_tau2,
                                  PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                  log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                  NMB1, NMB2, NME1, NME2, rho0, rho,
                                  CRPS1, CRPS2, ES1, ES2,
                                  round((t3 - T0)[[3]]), iter + 1)) %>% as.data.frame(),
                          Object)
            colnames(temp) <- c(dimnames(data$X_ts)[[1]],
                                "alpha", "E_tau2", paste0("theta", 1:2),
                                "K", "Proc_tau2" , "K0", "Proc0_tau2", "Iter_loglik_Error", "Fitting_RMSE",
                                "Testing_RMSE", "Fitting_MB", "Testing_MB",
                                "Fitting_NMB", "Testing_NMB", "Fitting_NME",
                                "Testing_NME", "Corr_Bef", "Corr_Aft",
                                "Fitting_CRPS", "Testing_CRPS", 
                                "Fitting_ES", "Testing_ES",
                                "elapsed", "iter", "Object")
          }
        }else{
          if(heavy.tail){    
            temp <- cbind(rbind(c(true.para$betaX, true.para$betaZ, 
                                  true.para$alpha, true.para$a,true.para$b,
                                  true.para$theta,
                                  true.para$k,
                                  true.para$Proc.tau2,
                                  true.para$k0,
                                  true.para$Proc0.tau2, 
                                  rep(1, 16), iter = 0),
                                c(initial.para$beta$E_betaX,
                                  initial.para$beta$E_betaZ,
                                  initial.para$alpha$E_alpha,
                                  initial.para$Obs.tau2$a,
                                  initial.para$Obs.tau2$b,
                                  initial.para$theta1$E_theta1,
                                  initial.para$theta2$E_theta2,
                                  initial.para$k$E_k, 
                                  initial.para$Proc.tau2$E_tau2, 
                                  initial.para$k0$E_k0, 
                                  initial.para$Proc0.tau2$E_tau2,
                                  rep(1, 16), iter = 0),
                                c(PIU$beta$E_betaX, 
                                  PIU$beta$E_betaZ, PIU$alpha$E_alpha,
                                  prior$Obs.tau2$a, prior$Obs.tau2$b,
                                  PIU$theta1$E_theta1,
                                  PIU$theta2$E_theta2, PIU$k$E_k,
                                  PIU$Proc.tau2$E_tau2,
                                  PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                  log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                  NMB1, NMB2, NME1, NME2, rho0, rho,
                                  CRPS1, CRPS2, ES1, ES2,
                                  round((t3 - T0)[[3]]), iter + 1)) %>% as.data.frame(),
                          Object)
            colnames(temp) <- c(paste0("betaX", 0:(length(PIU$beta$E_betaX) - 1)),
                                paste0("betaZ", 0:(length(PIU$beta$E_betaZ) - 1)),
                                "alpha", "a_tau2_prior", "b_tau2_prior", paste0("theta", 1:2),
                                "K", "Proc_tau2" , "K0", "Proc0_tau2", "Iter_loglik_Error", "Fitting_RMSE",
                                "Testing_RMSE", "Fitting_MB", "Testing_MB",
                                "Fitting_NMB", "Testing_NMB", "Fitting_NME",
                                "Testing_NME", "Corr_Bef", "Corr_Aft",
                                "Fitting_CRPS", "Testing_CRPS", 
                                "Fitting_ES", "Testing_ES",
                                "elapsed", "iter", "Object")
          }else{   
            temp <- cbind(rbind(c(true.para$betaX, true.para$betaZ, 
                                  true.para$alpha,
                                  true.para$Obs.tau2, true.para$theta,
                                  true.para$k,
                                  true.para$Proc.tau2,
                                  true.para$k0,
                                  true.para$Proc0.tau2, 
                                  rep(1, 16), iter = 0),
                                c(initial.para$beta$E_betaX,
                                  initial.para$beta$E_betaZ,
                                  initial.para$alpha$E_alpha,
                                  initial.para$Obs.tau2$E_tau2,
                                  initial.para$theta1$E_theta1,
                                  initial.para$theta2$E_theta2,
                                  initial.para$k$E_k, 
                                  initial.para$Proc.tau2$E_tau2, 
                                  initial.para$k0$E_k0, 
                                  initial.para$Proc0.tau2$E_tau2,
                                  rep(1, 16), iter = 0),
                                c(PIU$beta$E_betaX, 
                                  PIU$beta$E_betaZ, PIU$alpha$E_alpha,
                                  PIU$Obs.tau2$E_tau2, PIU$theta1$E_theta1,
                                  PIU$theta2$E_theta2, PIU$k$E_k,
                                  PIU$Proc.tau2$E_tau2,
                                  PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                  log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                  NMB1, NMB2, NME1, NME2, rho0, rho,
                                  CRPS1, CRPS2, ES1, ES2,
                                  round((t3 - T0)[[3]]), iter + 1)) %>% as.data.frame(),
                          Object)
            colnames(temp) <- c(dimnames(data$X_ts)[[1]],
                                paste0("betaZ", 0:(length(PIU$beta$E_betaZ) - 1)),
                                "alpha", "E_tau2", paste0("theta", 1:2),
                                "K", "Proc_tau2" , "K0", "Proc0_tau2", "Iter_loglik_Error", "Fitting_RMSE",
                                "Testing_RMSE", "Fitting_MB", "Testing_MB",
                                "Fitting_NMB", "Testing_NMB", "Fitting_NME",
                                "Testing_NME", "Corr_Bef", "Corr_Aft",
                                "Fitting_CRPS", "Testing_CRPS", 
                                "Fitting_ES", "Testing_ES",
                                "elapsed", "iter", "Object")
          }
        }
        RowName[1] <- "True:"
        RowName[2] <-  "Initi:"
        RowName[3] <-  "Iter1:"
        temp[, 1:(ncol(temp) - 3)] <- round(temp[, 1:(ncol(temp) - 3)], 5)
        if(!is.null(database$DSN))
        {
          sqlSave(database$DSN, temp, database$Table,
                  append = TRUE, colnames = FALSE,
                  rownames = FALSE, safer = TRUE, 
                  fast = TRUE)
        }
      }else{
        if(is.null(data$Z_ts)){
          if(heavy.tail){
            temp0 <- cbind(rbind(NULL, c(PIU$beta$E_betaX, PIU$alpha$E_alpha,
                                         prior$Obs.tau2$a, prior$Obs.tau2$b,
                                         PIU$theta1$E_theta1,
                                         PIU$theta2$E_theta2, 
                                         PIU$k$E_k, PIU$Proc.tau2$E_tau2,
                                         PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                         log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                         NMB1, NMB2, NME1, NME2, rho0, rho,
                                         CRPS1, CRPS2, ES1, ES2,
                                         round((t3 - T0)[[3]]), iter + 1))%>% as.data.frame(),
                           Object)
            colnames(temp0) <- c(dimnames(data$X_ts)[[1]],
                                 "alpha", "a_tau2_prior", "b_tau2_prior", paste0("theta", 1:2),
                                 "K", "Proc_tau2" , "K0", "Proc0_tau2", "Iter_loglik_Error",
                                 "Fitting_RMSE", "Testing_RMSE", "Fitting_MB",
                                 "Testing_MB", "Fitting_NMB", "Testing_NMB",
                                 "Fitting_NME", "Testing_NME",
                                 "Corr_Bef", "Corr_Aft",
                                 "Fitting_CRPS", "Testing_CRPS", 
                                 "Fitting_ES", "Testing_ES",
                                 "elapsed", "iter", "Object")
          }else{
            temp0 <- cbind(rbind(NULL, c(PIU$beta$E_betaX, PIU$alpha$E_alpha,
                                         PIU$Obs.tau2$E_tau2, PIU$theta1$E_theta1,
                                         PIU$theta2$E_theta2, 
                                         PIU$k$E_k, PIU$Proc.tau2$E_tau2,
                                         PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                         log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                         NMB1, NMB2, NME1, NME2, rho0, rho,
                                         CRPS1, CRPS2, ES1, ES2,
                                         round((t3 - T0)[[3]]), iter + 1))%>% as.data.frame(),
                           Object)
            colnames(temp0) <- c(dimnames(data$X_ts)[[1]],
                                 "alpha", "E_tau2", paste0("theta", 1:2),
                                 "K", "Proc_tau2" , "K0", "Proc0_tau2", "Iter_loglik_Error",
                                 "Fitting_RMSE", "Testing_RMSE", "Fitting_MB",
                                 "Testing_MB", "Fitting_NMB", "Testing_NMB",
                                 "Fitting_NME", "Testing_NME",
                                 "Corr_Bef", "Corr_Aft",
                                 "Fitting_CRPS", "Testing_CRPS", 
                                 "Fitting_ES", "Testing_ES",
                                 "elapsed", "iter", "Object")
          }
        }else{
          if(heavy.tail){
            temp0 <- cbind(rbind(NULL, c(PIU$beta$E_betaX, PIU$beta$E_betaZ, 
                                         prior$Obs.tau2$a, prior$Obs.tau2$b,
                                         PIU$Obs.tau2$b, PIU$theta1$E_theta1,
                                         PIU$theta2$E_theta2, PIU$k$E_k,
                                         PIU$Proc.tau2$E_tau2,
                                         PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                         log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                         NMB1, NMB2, NME1, NME2, rho0, rho,
                                         CRPS1, CRPS2, ES1, ES2,
                                         round((t3 - T0)[[3]]), iter + 1))%>% as.data.frame(),
                           Object)
            colnames(temp0) <- c(paste0("betaX", 0:(length(PIU$beta$E_betaX) - 1)),
                                 paste0("betaZ", 0:(length(PIU$beta$E_betaZ) - 1)),
                                 "alpha", "a_tau2_prior", "b_tau2_prior",paste0("theta", 1:2),
                                 "K", "Proc_tau2" , "K0", "Proc0_tau2", "Iter_loglik_Error",
                                 "Fitting_RMSE", "Testing_RMSE", "Fitting_MB",
                                 "Testing_MB", "Fitting_NMB", "Testing_NMB",
                                 "Fitting_NME", "Testing_NME",
                                 "Corr_Bef", "Corr_Aft",
                                 "Fitting_CRPS", "Testing_CRPS", 
                                 "Fitting_ES", "Testing_ES",
                                 "elapsed", "iter", "Object")
          }else{
            temp0 <- cbind(rbind(NULL, c(PIU$beta$E_betaX, PIU$beta$E_betaZ, 
                                         PIU$alpha$E_alpha,
                                         PIU$Obs.tau2$E_tau2, PIU$theta1$E_theta1,
                                         PIU$theta2$E_theta2, PIU$k$E_k,
                                         PIU$Proc.tau2$E_tau2,
                                         PIU$k0$E_k0, PIU$Proc0.tau2$E_tau2,
                                         log.lik.error, RMSE1, RMSE2, MB1, MB2,
                                         NMB1, NMB2, NME1, NME2, rho0, rho,
                                         CRPS1, CRPS2, ES1, ES2,
                                         round((t3 - T0)[[3]]), iter + 1))%>% as.data.frame(),
                           Object)
            colnames(temp0) <- c(dimnames(data$X_ts)[[1]],
                                 paste0("betaZ", 0:(length(PIU$beta$E_betaZ) - 1)),
                                 "alpha", "E_tau2", paste0("theta", 1:2),
                                 "K", "Proc_tau2" , "K0", "Proc0_tau2", 
                                 "Iter_loglik_Error",
                                 "Fitting_RMSE", "Testing_RMSE", "Fitting_MB",
                                 "Testing_MB", "Fitting_NMB", "Testing_NMB",
                                 "Fitting_NME", "Testing_NME",
                                 "Corr_Bef", "Corr_Aft",
                                 "Fitting_CRPS", "Testing_CRPS", 
                                 "Fitting_ES", "Testing_ES",
                                 "elapsed", "iter", "Object")
          }
          
        }
        temp0[, 1:(ncol(temp0) - 3)] <- round(temp0[, 1:(ncol(temp0) - 3)], 5)
        if(!is.null(database$DSN))
        {
          sqlSave(database$DSN, temp0, database$Table,
                  append = TRUE, colnames = FALSE,
                  rownames  = FALSE,
                  safer = TRUE, fast = TRUE)
        }
        temp <- rbind(temp, temp0)
        RowName[iter + 3] <- paste0("Iter", iter + 1, ":")  
      }
      iter <- iter + 1
    }
  }
  rownames(temp) <- RowName
  if(!is.null(test))
  {
    test.result <- as.data.frame(test$Y_ts_true)
    base::rownames(test.result) <- NULL
    # base::colnames(test.result) <- paste0("S", colnames(test.result))
    Y.test <- as.data.frame(Y.test)
    # test.result$DATE_TIME <- Y.test$DATE_TIME <-
    #   as.Date(base::rownames(test$Y_ts))
    test.result$Flag <- TRUE
    Y.test$Flag <- FALSE
    base::colnames(Y.test) <-  base::colnames(test.result)
    test.result <- rbind(test.result, Y.test)
    
    if(!file.exists("./1_CrossValidation")){
      dir.create("./1_CrossValidation")
      dir.create("./1_CrossValidation/data")
    }
    Ks = list(EnXs = Ks$EnXf, Xs = Ks$Xf)
    save(PIU, Ks, test, #test.result, spT,
         file = paste0("./1_CrossValidation/data/", 
                       database$Table, "_Summary.RData"))
    # cat("--------KS-------------", max(Ks$Xf),"\n")
    # Re = list(Test.Object = Object
    #           , FinalEst = round(result, 4)
    #           , PIU = PIU
    #           # , IS = IS
    #           , Ks = list(EnXs = Ks$EnXf, Xs = Ks$Xf)
    #           , Pred.list = test.result
    #           , Final_Pred = temp0
    #           , Prior = prior
    #           , ParaEst = temp
    #           , Miss_Infor = list(Index = miss.index,
    #                               Fill_Value = Fill_Miss_Value)
    #           , Fit_data = Y.fit
    #           , test.data = test.data
    # )
    Re <- 0
    # cat("\n .....KS.............. \n")
  }else{
    Re = list(FinalEst = round(result, 4)
              , PIU = PIU, IS = IS
              , Ks = list(EnXs = Ks$EnXf, Xs = Ks$Xf)
              , ParaEst = temp
              , Prior = prior
              , response.transf = response.transf
              , covariate.transf = covariate.transf
              , Fitting_RMSE = RMSE1
              , Fitting_MB = MB1
              , Fitting_NMB = NMB1
              , Fitting_NME = NME1
              , Corr_Bef = rho0
              , Fitting_CRPS = CRPS1
              , Fitting_ES = ES1
              , Object = Object
              , loglik = log.lik
              , Iterations = iter
              , itMax = itMax
    )
    if(!file.exists("./2_Full_Data_Model")){
      dir.create("./2_Full_Data_Model")
      dir.create("./2_Full_Data_Model/data")
    }
    save(Re, file = paste0(paste0("./2_Full_Data_Model/data/"), database$Table, ".RData"))
  }
  return(Re)
}
