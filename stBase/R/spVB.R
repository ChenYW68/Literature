spVB.Gaussian <- function(data = NULL, Ks = NULL,
                          S = NULL, prior = NULL,
                          para = NULL, IS = NULL,
                          parallel = FALSE,
                          ds = 1e-2, sp = NULL,
                          VB.err = 1e-1, 
                          IS.size = 100,
                          Remove.Cores = 1, 
                          N.Chunk = 1,
                          verbose.VB = FALSE,
                          method = c("IS"),
                          alpha.C = -1e2,
                          iter = 0)
{
  options(warn = -1)
  if(is.null(sp)) { sp <- import("scipy") }
  if(is.null(Remove.Cores)) {Remove.Cores = 1 }
  if(is.null(Remove.Cores)) { Remove.Cores = 1}
  if(is.null(N.Chunk)) {N.Chunk = 1}
  # CPU.Number = 5; N.Chunk = 1
  # if(verbose.VB){
  if(dplyr::between(iter , alpha.C, alpha.C + 3)){lower <- Upp <- 10}else{Upp = data$n*data$Nt/2;
  lower = 1e-1}
  cat("\n***************************************************************** \n")
  cat("                   Start to execute VB algorithm ! \n\n")
  # }
  # beta --------------------------------------------------------------------
  if(is.null(data$Z_ts))
  {
    X_CuSum <- XYXi_CuSum <- XYW_CuSum <- 0
    for(t in 1:data$Nt)
    {
      X_CuSum <- X_CuSum + data$X_ts[, ,t] %*% t(data$X_ts[, , t])
      XYXi_CuSum <- XYXi_CuSum + data$X_ts[, , t] %*%
        ( data$Y_ts[t, ] - para$alpha$E_alpha * data$Hs %*% Ks$Xf[t + 1, ])
    }
    # E_inveser_tau2 --------------------------------------------------------------------
    E_inveser_tau2 <- (para$Obs.tau2$a/para$Obs.tau2$b)
    post_betaX_sigma2 <- solve(E_inveser_tau2 * X_CuSum +
                                 solve(prior$beta$betaX.Sigma2))
    eta <- E_inveser_tau2 * XYXi_CuSum +
      solve(prior$beta$betaX.Sigma2) %*% prior$beta$betaX.mu
    post_betaX_mu <- post_betaX_sigma2 %*% eta %>% as.vector()
    if(verbose.VB){
      cat(".................................")
      cat("\npost expectation of beta: \n")
      print(round(post_betaX_mu, 4))
      cat("post_betaX_sigma2: \n")
      print(round(post_betaX_sigma2, 4))
    }
    rm(X_CuSum, XYXi_CuSum)
    Y_ts <- data$Y_ts - X_ts_Transf(data$Nt, data$X_ts, post_betaX_mu)
    # alpha -------------------------------------------------------------------
    if(iter <= alpha.C){
      alpha_eta <- 0
      for(t in 1:data$Nt)
      {
        alpha_eta <- alpha_eta + t(Ks$Xf[t + 1, ])%*% t(data$Hs) %*%
          (Y_ts[t, ])
      }
      alpha_eta <- as.vector(alpha_eta*E_inveser_tau2)
      post_alpha_sigma2 <- as.vector(solve(sp$trace(t(data$Hs) %*%
                                                      data$Hs %*% S$S11)*E_inveser_tau2 +
                                             1/prior$alpha$Sigma2))
      post_alpha_mu <- as.vector( post_alpha_sigma2 * (alpha_eta +
                                                         prior$alpha$mu/prior$alpha$Sigma2))
      if(verbose.VB){
        cat("\npost expectation and sigma2 of alpha: \n")
        print(c(post_alpha_mu, post_alpha_sigma2))
        # cat(".................................\n\n")
      }}else{
        post_alpha_mu = para$alpha$E_alpha;post_alpha_sigma2 =0}
    # tau2 --------------------------------------------------------------------
    a_tau2 = prior$Obs.tau2$a + data$n * data$Nt / 2
    t1 <- proc.time()
    {
      # system.time({
      f_tau2 = sapply(seq_len(data$Nt), function(t) {
        (t(data$Y_ts[t, ]) %*% data$Y_ts[t, ] -
           2 * gpuR::t(data$Y_ts[t, ]) %*% gpuR::t(data$X_ts[, , t]) %*% post_betaX_mu)[1] +
          sp$trace((gpuR::vclMatrix(data$X_ts[, , t]) %*%
                      gpuR::t(gpuR::vclMatrix(data$X_ts[, , t])) %*%
                      (post_betaX_sigma2 + post_betaX_mu %*% t(post_betaX_mu))
          )%>% as.matrix()) -
          (2 * gpuR::t(post_alpha_mu*data$Hs %*% Ks$Xf[t + 1, ]) %*%
             (Y_ts[t, ]))[1]
      }, simplify = "matrix") %>% sum()
      # })
    }
    f_tau2 <- f_tau2 + post_alpha_mu^2 * sp$trace((gpuR::t(gpuR::vclMatrix(data$Hs)) %*%
                                                     gpuR::vclMatrix(data$Hs) %*% S$S11)%>% as.matrix())
    if(verbose.VB){
      t2 <- proc.time()
      cat("\n\nEstimation of tau_2 takes time: \n")
      print(t2 - t1)}
    
    b_tau2 <-  prior$Obs.tau2$b + f_tau2/2
    E_tau2 <-  b_tau2 / (a_tau2 - 1)
    # E_tau2 <-  a_tau2 / b_tau2
    
  }else{
    # BetaZ -------------------------------------------------------------------
    Hv_ts <- sapply(seq(1, data$n), function(s)
      (apply(Ks$Xf[2:(data$Nt + 1), ], 1, FUN = '%*%',
             data$Hs[s, ])), simplify = "array")
    Z_ts <-  X_ts_Transf(Nt = data$Nt,
                         X_ts = data$Z_ts,
                         beta = para$beta$E_betaZ) %*% gpuR::t(data$Hs) 
    
    X_CuSum <- XYXi_CuSum <- 0
    for(t in 1:data$Nt)
    {
      X_CuSum <- X_CuSum +  data$X_ts[, ,t] %*% t(data$X_ts[, , t])
      XYXi_CuSum <- XYXi_CuSum +  data$X_ts[, , t] %*%
        ( data$Y_ts[t, ]  - para$alpha$E_alpha*Hv_ts[t, ] -
            para$alpha$E_alpha*Z_ts[t, ])
    }
    E_inveser_tau2 <- (para$Obs.tau2$a/para$Obs.tau2$b)
    post_betaX_sigma2 <- solve(E_inveser_tau2 * X_CuSum + 
                                 solve(prior$beta$betaX.Sigma2))
    
    eta <- E_inveser_tau2 * XYXi_CuSum + 
      solve(prior$beta$betaX.Sigma2) %*% prior$beta$betaX.mu
    
    post_betaX_mu <- post_betaX_sigma2 %*% eta %>% as.vector()
    
    if(verbose.VB){
      cat("..........................")
      cat("\npost expectation of betaX: \n")
      print(round(post_betaX_mu, 4))
      cat("post_betaX_sigma2: \n")
      print(round(post_betaX_sigma2, 4))
    }
    rm(X_CuSum, XYXi_CuSum)
    
    X_ts <-  X_ts_Transf(Nt = data$Nt, 
                         X_ts = data$X_ts,
                         beta = post_betaX_mu)
    # t1 <- proc.time()
    Hz = sapply(seq_len(data$Nt), function(t) 
      data$Hs %*% t(data$Z_ts[,, t])
      , simplify="array")  # dim(n, p, Nt)
    
    Z_CuSum <- ZYZi_CuSum <- 0
    for(t in 1:data$Nt)
    {
      Z_CuSum <- Z_CuSum +  t(Hz[, ,t]) %*% Hz[, ,t]
      ZYZi_CuSum <- ZYZi_CuSum +  para$alpha$E_alpha* t(Hz[, ,t]) %*%
        ( data$Y_ts[t, ]  - para$alpha$E_alpha*Hv_ts[t, ] -
            X_ts[t, ])
    }
    post_betaZ_sigma2 <- solve( (para$alpha$E_alpha^2) *
                                  E_inveser_tau2 * Z_CuSum + 
                                  solve(prior$beta$betaZ.Sigma2))
    
    eta <- E_inveser_tau2 * ZYZi_CuSum + solve(prior$beta$betaZ.Sigma2) %*% prior$beta$betaZ.mu
    
    post_betaZ_mu <- post_betaZ_sigma2 %*% eta %>% as.vector()
    
    if(verbose.VB){
      cat("\n\npost expectation of betaZ: \n")
      print(post_betaZ_mu)
      cat("post_betaZ_sigma2: \n")
      print(post_betaZ_sigma2)
      cat("..........................")
    }
    rm(Z_CuSum, ZYZi_CuSum)
    
    
    # step2.2: upate alpha ----------------------------------------------------
    
    Y_ts <- data$Y_ts - X_ts_Transf(data$Nt, data$X_ts, post_betaX_mu)
    
    HZ_tsBeta <- X_ts_Transf(data$Nt, data$Z_ts, post_betaZ_mu) %*% t(data$Hs)
    
    alpha1 = sapply(seq_len(data$Nt), function(t)
      t(Hv_ts[t, ] + HZ_tsBeta[t, ]) %*% Y_ts[t, ]
      , simplify="array") %>% sum()
    
    HvD <- sum(as.vector(gpuR::diag(gpuR::t(gpuR::vclMatrix(data$Hs)) %*% 
                                      gpuR::vclMatrix(data$Hs)%*% 
                                      gpuR::vclMatrix(S$S11) )))
    S.betaZ <- (gpuR::tcrossprod(post_betaZ_mu) + 
                  post_betaZ_sigma2)
    alpha2 = sapply(seq_len(data$Nt), function(t)
      2 * HZ_tsBeta[t, ] %*% Hv_ts[t, ]  + 
        sum(as.vector( gpuR::diag(gpuR::tcrossprod(t(Hz[,, t])) %*% S.betaZ)))
      , simplify="array")%>% sum() + HvD
    
    
    if(iter <= alpha.C){
      alpha_eta <- E_inveser_tau2 * alpha1
      D_alpha <- E_inveser_tau2* alpha2
      post_alpha_sigma2 <- solve(D_alpha + 1/prior$alpha$Sigma2)[1]
      post_alpha_mu <- as.vector( post_alpha_sigma2 * (alpha_eta +
                                                         prior$alpha$mu/prior$alpha$Sigma2))
      if(verbose.VB){
        cat("\npost expectation and sigma2 of alpha: \n")
        print(c(post_alpha_mu, post_alpha_sigma2))
        # cat(".................................\n\n")
      }}else{
        post_alpha_mu = para$alpha$E_alpha;post_alpha_sigma2 =0}
    # tau2 --------------------------------------------------------------------
    a_tau2 = prior$Obs.tau2$a + data$n * data$Nt / 2
    t1 <- proc.time()
    {
      # system.time({
      f_tau2 = sapply(seq_len(data$Nt), function(t) {
        (t(data$Y_ts[t, ]) %*% data$Y_ts[t, ] -
           2 * gpuR::t(data$Y_ts[t, ]) %*% gpuR::t(data$X_ts[, , t]) %*% post_betaX_mu)[1] +
          sp$trace((gpuR::vclMatrix(data$X_ts[, , t]) %*%
                      gpuR::t(gpuR::vclMatrix(data$X_ts[, , t])) %*%
                      (post_betaX_sigma2 + post_betaX_mu %*% 
                         t(post_betaX_mu))
          )%>% as.matrix())
      }, simplify = "matrix") %>% sum()
      # })
      
      f_tau2 = f_tau2  - 2 * post_alpha_mu * alpha1 +
        (post_alpha_mu^2 + post_alpha_sigma2) * alpha2
    }
    if(verbose.VB){
      t2 <- proc.time()
      cat("\nEstimation of tau_2 takes time: \n")
      print(t2 - t1)}
    b_tau2 <-  (prior$Obs.tau2$b + f_tau2/2)
    E_tau2 = b_tau2/(a_tau2 - 1)
  }
  if(verbose.VB){
    # cat(".................................\n\n")
    cat("post expectation of tau2: ", E_tau2, "\n a: ", a_tau2, "b: ", b_tau2)
    # cat(".................................\n")
  }
  
  # Q and M ------------------------------------------------------------------
  # rho.space <- fields::Wendland(data$BAUs.Dist
  #                       , theta = max(data$BAUs.Dist)*cs
  #                       , dimension = 1, k =1)

  ###########################################################################

  ###########################################################################
  M <- gpuR::vclMatrix(ds*exp(- (data$BAUs.Dist^2 / para$theta2$E_theta2)))*Ks$bandKernel
  # Q <- para$Q$E_Q
  Q_M <- gpuR::vclMatrix(para$Q$E_Q) %*% M
  
  # step4: upate Theta1's distribution
  f1_theta1 <- sum(as.vector(gpuR::diag(gpuR::t(M) %*% Q_M %*% gpuR::vclMatrix(S$S00))))
  f2_theta1 <- sum(as.vector(gpuR::diag(Q_M %*% gpuR::vclMatrix(S$S01))))
  ###########################################################################
  # M <- as(ds*exp(- (data$BAUs.Dist^2 / para$theta2$E_theta2))*Ks$bandKernel, "sparseMatrix")
  # Q <- as(para$Q$E_Q, "sparseMatrix")
  # Q_M <- Q %*% M
  # f1_theta1 <- Trace_Muti(as.matrix(Matrix::t(M) %*% Q_M), S$S00)
  # f2_theta1 <- Trace_Muti(as.matrix(Q_M), S$S01)
  
  #   sum(as.vector(gpuR::diag(gpuR::t(T_Mt_Q) %*% gpuR::vclMatrix(S$S01))))
  
  post_theta1_sigma2 <- solve(f1_theta1 + 1/prior$theta1$Sigma2) %>% as.vector()
  theta1_eta <- f2_theta1 + prior$theta1$mu/prior$theta1$Sigma2
  
  post_theta1_mu <- post_theta1_sigma2 %*% theta1_eta %>% as.vector()
  
  ###########################################################################
  # bandKernel <- Ks$bandKernel
  # D <- data$BAUs.Dist 
  # S00 <- S$S00
  # S01 <- S$S01
  # theta1 <- post_theta1_mu
  # sigTheta1 <- post_theta1_sigma2
  # save(bandKernel, Q, D, ds, S00, S01, theta1, sigTheta1, file = "./data/temp.RData")
  ###########################################################################
  
  if(verbose.VB){
    cat("\n\npost expectation and sigma2 of theta1: \n")
    print(c(post_theta1_mu, post_theta1_sigma2))
    # cat(".................................\n\n")
  }
  ###########################################################################
  ###########################################################################
  # lambda0 ------------------------------------------------------------------
  Q <- as(para$Q$E_Q, "sparseMatrix")
  # theta2 ----------------------------------------------------------------------
  if(method %in% c("IS")){
    BAUs.Dist <- data$BAUs.Dist
    bandKernel <- Ks$bandKernel
    {
      #cat("........................................................\n")
      G_theta2 <- vector()
      if(iter == 0){
        # rtheta2 <- runif(IS.size, para$theta2$a, para$theta2$b)
        # Da.RS <- data.frame(Rs = rtheta2,
        #                     Gw = dunif(rtheta2, para$theta2$a,
        #                                para$theta2$b, log = T))
        Da.RS <- data.frame(Rs = seq(para$theta2$a, para$theta2$b,
                                     (para$theta2$b - para$theta2$a)/(IS.size - 1)),
                            Gw = (1/IS.size))
      }else{
        Da.RS <- data.frame(Rs = IS$theta2$sample,
                            Gw = IS$theta2$weight)
      }
      #print(Da.RS[1:5,])
      
      theta2 <- Da.RS$Rs
      SIR <- nrow(Da.RS)
      if(IS$Thresh[1] == 1)
      {
        t1 <- proc.time()
        if(parallel == T){
          Chunk <- chunk(1:SIR, N.Chunk)
          G_theta2 <- NULL
          for(B in 1:N.Chunk){
            Seq <- Chunk[[B]]
            no_cores <- detectCores() - Remove.Cores
            cl <- makeCluster(getOption("cl.cores", no_cores))
            clusterEvalQ(cl,library(gpuR))
            clusterEvalQ(cl,library(stBase))
            clusterEvalQ(cl,library(Matrix))
            clusterExport(cl = cl
                          , varlist = c("BAUs.Dist", "ds", "bandKernel"
                                        , "Q", "theta2", "S"
                                        , "post_theta1_mu"
                                        , "post_theta1_sigma2"
                                        , "f_theta", "Trace_Muti")
                          , envir = environment())
            temp <- parallel::parLapply(cl, X = Seq, fun = f_theta)
            temp <- as.vector(do.call("rbind", temp))# * Da.RS$Gw[Seq]
            stopCluster(cl)
            G_theta2 <- c(G_theta2, temp)
          
          # G_theta2 = c(G_theta2, f_theta2(theta = theta2, ds = ds, 
          #                                  m = dim(para$Q$E_Q), 
          #                                  n = IS.size, 
          #                                  D = data$BAUs.Dist, 
          #                                  bandKernel = Ks$bandKernel, 
          #                                  Q = para$Q$E_Q, 
          #                                  S00 = S$S00, S01 = S$S01, 
          #                                  theta1 = post_theta1_mu, 
          #                                  sigTheta1 = post_theta1_sigma2, 
          #                                  nThreads = no_cores)$logLik)
          }
          
        }else{
          assign("S", S, envir = .GlobalEnv)
          assign("ds", ds, envir = .GlobalEnv)
          assign("bandKernel", bandKernel, envir = .GlobalEnv)
          assign("post_theta1_mu", post_theta1_mu, envir = .GlobalEnv)
          assign("post_theta1_sigma2", post_theta1_sigma2, envir = .GlobalEnv)
          assign("BAUs.Dist", BAUs.Dist, envir = .GlobalEnv)
          assign("Q", Q, envir = .GlobalEnv)
          assign("theta2", theta2, envir = .GlobalEnv)
          G_theta2 = sapply(X = seq_len(SIR), f_theta, simplify = "array")
        }
        t2 <- proc.time()
        if(verbose.VB){
          cat("\nEstimation of theta2 takes time: \n")
          print((t2 - t1))[[3]]}
        
        if_else(is.infinite(G_theta2), min(G_theta2[!is.finite(G_theta2)]), G_theta2)
        # importance --------------------------------------------------------------
        da <- importance(Da.RS$Rs, G_theta2, Da.RS$Gw, IS.err, iter = iter)
        
        E_theta2 <- da$EIS
        Da <- data.frame(x = da$x, Gw = da$Pro)# %>% dplyr::distinct(Rs, Gw)
        
        p <- ggplot(Da) + geom_point(aes(x = x, y = Gw)) +
          geom_vline(xintercept = E_theta2,
                     col = "red", linetype="dotted") +
          labs(x = TeX("$\\theta_2$"), y = "density") + theme_bw()
        print(p)
        IS$theta2$sample = da$Rs
        IS$theta2$weight = da$Gw
        IS$theta2$Final.sample = da$x
        IS$theta2$Final.weight = da$Pro
        
        
        Err <- abs(E_theta2 - para$theta2$E_theta2)/para$theta2$E_theta2
        if(verbose.VB){cat("theta2_w_err: ", Err)}
        if(Err < VB.err)
        {
          IS$Thresh[1] = 0
        }
        
      }else{
        E_theta2 = para$theta2$E_theta2
      }
    }
  }else{
    assign("S", S, envir = .GlobalEnv)
    assign("ds", ds, envir = .GlobalEnv)
    assign("bandKernel", bandKernel, envir = .GlobalEnv)
    assign("post_theta1_mu", post_theta1_mu, envir = .GlobalEnv)
    assign("post_theta1_sigma2", post_theta1_sigma2, envir = .GlobalEnv)
    assign("BAUs.Dist", data$BAUs.Dist, envir = .GlobalEnv)
    assign("Q", Q, envir = .GlobalEnv)
    # assign("theta2", theta2, envir = .GlobalEnv)
    op.theta2 <- optim(c(para$theta2$E_theta2),
                       f_theta_op, 
                       method = "L-BFGS-B",
                       lower = lower, #upper = 6,
                       control = list(REPORT = 1,
                                      trace = 0)
                       , hessian = F
    )
    E_theta2 <- (op.theta2$par) 
  }
  if(verbose.VB){
    cat("\npost expectation of theta2: ", E_theta2, "\n\n")
    # cat("........................................................\n")
  }

  
  #K --------------------------------------------------------------------
  # M <- ds*exp(- BAUs.Dist^2 / E_theta2)*Ks$bandKernel
  Nt <- data$Nt
  Adj.Mat = as(data$Adj.Mat, "sparseMatrix")
  
  Proc.tau2 <- para$Proc.tau2$E_tau2
  M <- post_theta1_mu * ds*exp(- data$BAUs.Dist^2 / E_theta2)*Ks$bandKernel
  if(method %in% c("IS")){
    {
      # cat("........................................................\n")
      G_k <- vector()
      if(iter == 0){
        # rK <- runif(IS.size, para$k$a, para$k$b)
        # Da.RS <- data.frame(Rs = rK,
        #                     Gw = dunif(rK, para$k$a, para$k$b, log = T))
        
        Da.RS <- data.frame(Rs = seq(para$k$a, para$k$b,
                                     (para$k$b - para$k$a)/(IS.size - 1)),
                            Gw = (1/IS.size))
      }else{
        Da.RS <- data.frame(Rs = IS$K$sample,
                            Gw = IS$K$weight)
      }
      Rs = Da.RS$Rs
      SIR <- nrow(Da.RS)
      if(IS$Thresh[2] == 1)
      {
        
        TeK1 <- sum(as.vector(gpuR::diag((S$S11))))
        TeK2 <- post_theta1_mu *Trace_Muti(M, S$S01)
        TeK3 <- (post_theta1_sigma2 + post_theta1_mu^2) *
          Trace_Muti(as.matrix(gpuR::t(M) %*%  gpuR::vclMatrix(M)), S$S00)
        TeK <- Proc.tau2 *(-TeK1/2 + TeK2 - TeK3/2)
        t1 <- proc.time()
        if(parallel == T){
          Chunk <- chunk(1:SIR, N.Chunk)
          G_k <- NULL
          for(B in 1:N.Chunk){
            Seq = Chunk[[B]]
            ############################################################
            no_cores <- detectCores() - Remove.Cores
            cl <- parallel::makeCluster(getOption("cl.cores", no_cores))
            clusterEvalQ(cl,library(gpuR))
            clusterEvalQ(cl,library(stBase))
            parallel::clusterExport(cl = cl
                                    , varlist = c("Rs", "TeK", "Adj.Mat"
                                                  , "Nt", "f_k")
                                    , envir = environment()
            )
            temp <- parallel::parLapply(cl, X = Seq, fun = f_k)
            temp <- as.vector(do.call("rbind", temp))# * Da.RS$Gw[Seq]
            stopCluster(cl)
            G_k <- c(G_k, temp)
            gc()
          }
        } else{
          assign("Rs", Rs, envir = .GlobalEnv)
          assign("TeK", TeK, envir = .GlobalEnv)
          # assign("Proc.tau2", Proc.tau2, envir = .GlobalEnv)
          assign("Adj.Mat", Adj.Mat, envir = .GlobalEnv)
          assign("Nt", Nt, envir = .GlobalEnv)
          G_k = sapply(X = seq_len(SIR), f_k, simplify = "array")
          #
        }
        t2 <- proc.time()
        if(verbose.VB){
          cat("Estimation of K takes time: \n")
          print((t2 - t1))[[3]]}
        # if_else(is.infinite(G_k), min(G_k[!is.finite(G_k)]), G_k)
        ############################################################
        da <- importance(Da.RS$Rs, as.vector(G_k), Da.RS$Gw, IS.err, iter = iter)
        E_k <- da$EIS
        Da <- data.frame(x = da$x, Gw = da$Pro)# %>% dplyr::distinct(Rs, Gw)
        
        p <- ggplot(Da) + geom_point(aes(x = x, y = Gw)) +
          geom_vline(xintercept = E_k, col = "red", linetype="dotted") +
          labs(x = TeX("$\\K$"), y = "density") + theme_bw()
        print(p)
        
        IS$K$sample = da$Rs
        IS$K$weight = da$Gw
        IS$K$Final.sample = da$x
        IS$K$Final.weight = da$Pro
        
        Err <- abs(E_k - para$k$E_k)/para$k$E_k
        if(verbose.VB){cat("K_w_err: ", Err, "\n")}
        if(Err < VB.err)
        {
          IS$Thresh[2] = 0
        }
      }else{
        E_k = para$k$E_k
      }
    }
  }else{
    TeK1 <- sum(as.vector(gpuR::diag((S$S11))))
    TeK2 <- post_theta1_mu *Trace_Muti(M, S$S01)
    TeK3 <- (post_theta1_sigma2 + post_theta1_mu^2) *
      Trace_Muti(as.matrix(gpuR::t(M) %*%  gpuR::vclMatrix(M)), S$S00)
    TeK <- Proc.tau2 *(-TeK1/2 + TeK2 - TeK3/2)
    assign("TeK", TeK, envir = .GlobalEnv)
    # assign("Proc.tau2", Proc.tau2, envir = .GlobalEnv)
    assign("Adj.Mat", Adj.Mat, envir = .GlobalEnv)
    assign("Nt", Nt, envir = .GlobalEnv)
    op.k <- optim(c(log(para$k$E_k)), f_k_op,
                  method = "Nelder-Mead",
                  # lower = -5, #upper = 6,
                  control = list(REPORT = 1,
                                 trace = 0)
                  , hessian = F
    )
    E_k <- exp(op.k$par) 
  }
  if(verbose.VB){
    cat("post expectation of K: ", E_k, "\n\n")
  }

  # tau2_proce --------------------------------------------------------------
  G = data$Adj.Mat + E_k^2 * diag(data$N.BAUs)
  # M = as(M, "sparseMatrix")
  # G_M <- as(G, "sparseMatrix") %*% M #gpuR::vclMatrix(M)
  G_M <- G %*% gpuR::vclMatrix(M)
  a.proc <- prior$Proc.tau2$a + data$N.BAUs*data$Nt/2
  b.proc <- prior$Proc.tau2$b + (Trace_Muti(G, S$S11) - 2*post_theta1_mu*
                                   Trace_Muti(as.matrix(G_M), S$S01) +
                                   (post_theta1_sigma2 + post_theta1_mu^2) *
                                   Trace_Muti(as.matrix(gpuR::t(M) %*% G_M), S$S00))/2
  
  # b.proc <- prior$Proc.tau2$b + (Trace_Muti(as.matrix(G), S$S11) - 2*post_theta1_mu*
  #                                  Trace_Muti(as.matrix(G_M), S$S01) +
  #                                  (post_theta1_sigma2 + post_theta1_mu^2) * 
  #                                  Trace_Muti(as.matrix(Matrix::t(M) %*% G_M), S$S00))/2
  Proc.tau2 <-  a.proc/(b.proc)
  if(verbose.VB){
    cat("post expectation of Proc.tau2: ", Proc.tau2, "\n\n")
    # cat("............................s............................\n")
  }
  
  # Proc.tau2 = 1
  # ###########################################################################
  # K0--------------------------------------------------------------------
  Mu0 <- gpuR::tcrossprod(Ks$Xf[1, ])
  Pt0 <- Ks$Pt[, , 1]
  Proc0.tau2 <-  para$Proc0.tau2$E_tau2
  trace_Pt0 <- Proc0.tau2*sum(as.vector(gpuR::diag(Pt0 + Mu0)))
  if(method %in% c("IS")){
    {
      # cat("........................................................\n")
      G_k0 <- vector()
      if(iter == 0){
        # rK0 <- runif(IS.size, para$k0$a, para$k0$b)
        # Da.RS <- data.frame(Rs = rK0,
        #                     Gw = dunif(rK0, para$k0$a, para$k0$b, log = T))
        Da.RS <- data.frame(Rs = seq(para$k0$a, para$k0$b,
                                     (para$k0$b - para$k0$a)/(IS.size - 1)),
                            Gw = (1/IS.size))
      }else{
        Da.RS <- data.frame(Rs = IS$K0$sample,
                            Gw = IS$K0$weight)
      }
      SIR <- nrow(Da.RS)
      Rs = Da.RS$Rs
      if(IS$Thresh[3] == 1)
      {
        
        t1 <- proc.time()
        if(parallel == T){
          Chunk <- chunk(1:SIR, N.Chunk)
          G_k0 <- NULL
          for(B in 1:N.Chunk){
            Seq = Chunk[[B]]
            ############################################################
            no_cores <- detectCores() - Remove.Cores
            cl <- parallel::makeCluster(getOption("cl.cores", no_cores))
            clusterEvalQ(cl,library(gpuR))
            #clusterEvalQ(cl,library(reticulate))
            parallel::clusterExport(cl=cl
                                    , varlist=c("Rs", "Adj.Mat"
                                                , "trace_Pt0"
                                                , "f_k0")
                                    , envir = environment()
            )
            temp <- parallel::parLapply(cl, X = Seq, fun = f_k0)
            temp <- as.vector(do.call("rbind", temp))# * Da.RS$Gw[Seq]
            stopCluster(cl)
            G_k0 <- c(G_k0, temp)
          }
          # ############################################################
          # # })
        }else{
          assign("Rs", Rs, envir = .GlobalEnv)
          assign("Adj.Mat", Adj.Mat, envir = .GlobalEnv)
          # assign("Proc0.tau2", Proc0.tau2, envir = .GlobalEnv)
          assign("trace_Pt0", trace_Pt0, envir = .GlobalEnv)
          G_k0 = sapply(X = seq_len(SIR), f_k0, simplify = "array")
          
        }
        t2 <- proc.time()
        if(verbose.VB){
          # cat("........................................................\n")
          cat("Estimation of K0 takes time: \n")
          print((t2 - t1))[[3]]}
        if_else(is.infinite(G_k0), min(G_k0[!is.finite(G_k0)]), G_k0)
        ############################################################
        da <- importance(Da.RS$Rs, G_k0, Da.RS$Gw, IS.err, iter = iter)
        E_k0 <- da$EIS
        Da <- data.frame(x = da$x, Gw = da$Pro)# %>% dplyr::distinct(Rs, Gw)
        
        p <- ggplot(Da) + geom_point(aes(x = x, y = Gw)) +
          geom_vline(xintercept = E_k0, col = "red", linetype = "dotted") +
          labs(x = TeX("$\\K_0$"), y = "density") + theme_bw()
        print(p)
        
        IS$K0$sample = da$Rs
        IS$K0$weight = da$Gw
        IS$K0$Final.sample = da$x
        IS$K0$Final.weight = da$Pro
        
        Err <- abs(E_k0 - para$k0$E_k0)/para$k0$E_k0
        if(verbose.VB){ cat("K0_w_err: ", Err, "\n")}
        if(Err < VB.err)
        {
          IS$Thresh[3] = 0
        }
      }else{
        E_k0 = para$k0$E_k0
      }
    }
  }else{
    assign("Adj.Mat", Adj.Mat, envir = .GlobalEnv)
    # assign("Proc0.tau2", Proc0.tau2, envir = .GlobalEnv)
    assign("trace_Pt0", trace_Pt0, envir = .GlobalEnv)
    op.k0 <- optim(c(log(para$k0$E_k0)),
                   f_k0_op,
                   method = "Nelder-Mead",
                   # lower = -5, #upper = 6,
                   control = list(REPORT = 1, 
                                  trace = 0)
                   , hessian = F
    )
    E_k0 <- exp(op.k0$par) 
  }
  if(verbose.VB){
    cat("post expectation of K0: ", E_k0, "\n\n")
    # cat("............................s............................\n")
  }
  
  G0 = data$Adj.Mat + E_k0^2 * diag(data$N.BAUs)
  
  a.proc0 <- prior$Proc0.tau2$a + data$N.BAUs/2
  b.proc0 <- prior$Proc0.tau2$b + Trace_Muti(G0, (Mu0 + Pt0))/2
  Proc0.tau2 <-  a.proc0/(b.proc0)
  if(verbose.VB){
    cat("post expectation of Proc0.tau2: ", Proc0.tau2, "\n\n")
    # cat("............................s............................\n")
  }
  # Proc0.tau2 = 1
  # result ------------------------------------------------------------------
  if(is.null(data$Z_ts))
  {
    para <- list(
      para = list(
        beta = list(E_betaX = post_betaX_mu
                    , betaX.Sigma2 = post_betaX_sigma2)
        , theta1 = list(E_theta1 = post_theta1_mu
                        , Sigma2 = post_theta1_sigma2)
        , alpha = list(E_alpha = post_alpha_mu
                       , Sigma2 = post_alpha_sigma2)
        , k = list(E_k = E_k, a = para$k$a, b = para$k$b)
        , Q = list(E_Q = Proc.tau2*G)
        , Proc.tau2 = list(E_tau2 = Proc.tau2,
                           a = a.proc, b = b.proc)
        , k0 = list(E_k0 = E_k0, a = para$k0$a, b = para$k0$b)
        , Q0 = list(E_Q0 = Proc0.tau2*G0)
        , Proc0.tau2 = list(E_tau2 = Proc0.tau2,
                            a = a.proc0, b = b.proc0)
        , theta2 = list(E_theta2 = E_theta2, a = para$theta2$a,
                        b = para$theta2$b)
        , Obs.tau2 = list(E_tau2 = E_tau2, a = a_tau2, b = b_tau2)
      )
      , IS = IS
    )
  }else{
    para <- list(
      para = list(
        beta = list(E_betaX = post_betaX_mu
                    , betaX.Sigma2 = post_betaX_sigma2
                    , E_betaZ = post_betaZ_mu
                    , betaZ.Sigma2 = post_betaZ_sigma2)
        , theta1 = list(E_theta1 = post_theta1_mu
                        , Sigma2 = post_theta1_sigma2)
        , alpha = list(E_alpha = post_alpha_mu
                       , Sigma2 = post_alpha_sigma2)
        , k = list(E_k = E_k, a = para$k$a, b = para$k$b)
        , Q = list(E_Q = Proc.tau2*G)
        , Proc.tau2 = list(E_tau2 = Proc.tau2,
                           a = a.proc, b = b.proc)
        , k0 = list(E_k0 = E_k0, a = para$k0$a, b = para$k0$b)
        , Q0 = list(E_Q0 = Proc0.tau2*G0)
        , Proc0.tau2 = list(E_tau2 = Proc0.tau2,
                            a = a.proc0, b = b.proc0)
        , theta2 = list(E_theta2 = E_theta2, a = para$theta2$a,
                        b = para$theta2$b)
        , Obs.tau2 = list(E_tau2 = E_tau2, a = a_tau2, b = b_tau2)
        
      )
      , IS = IS
    )}
  return(para)
}