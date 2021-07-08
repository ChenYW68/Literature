TaEnKs = function(para.ks)
{
  cat("***************************************************************** \n")
  cat("                Start to execute TaEnKs algorithm ! \n\n")
  {
    xf = array(NA, c(para.ks$Nt + 1, para.ks$N.BAUs
                     , para.ks$Ensemble.size))
    xp = array(NA, c(para.ks$Nt + 1, para.ks$N.BAUs
                     , para.ks$Ensemble.size))
    forecast = array(NA, c(para.ks$Nt + 1, para.ks$N.BAUs
                           , para.ks$Ensemble.size))
    VIF = 0.0
    d = 0:(para.ks$Nt)
    rho.time = taper(d, para.ks$ct + 1)
    names(rho.time) = paste(d)
   
    
    L_Q0_upper <- Matrix::solve(Matrix::chol(as(para.ks$Q0, "sparseMatrix")))  #para.ks$sp$linalg$cholesky
    L_Q_upper <- Matrix::solve(Matrix::chol(as(para.ks$Q, "sparseMatrix")))
    
    # L_Q0_upper <- Matrix::solve(Matrix::t(as(Matrix::Cholesky(as(para.ks$Q0, "sparseMatrix")), "sparseMatrix")))
    # L_Q_upper1 <- Matrix::solve(Matrix::t(as(Matrix::Cholesky(as(para.ks$Q, "sparseMatrix")), "sparseMatrix")))
                      
    
    # library(mvnfast)
    if(para.ks$heavy.tail){
      xf[1,,] = rMvn(para.ks$Ensemble.size,
                     mu = para.ks$mu0,
                     L = L_Q0_upper,
                     Cov = F)#para.ks$lambda0
    }else{
      xf[1,,] = rMvn(para.ks$Ensemble.size,
                     as.vector(para.ks$mu0), 
                     L_Q0_upper,
                     Cov = F)# 
    }
    
    #
  }
  ###############################################################
  Hs = gpuR::vclMatrix(para.ks$Hs)
  if(para.ks$Nt >= 1)
  {
    pb <- progress_bar$new(format = "|:bar| :current/:total (:percent in :elapsed)"
                           , total = para.ks$Nt
                           , clear = FALSE
                           # , width= para.ks$Nt
    )
    for (i in 2:(para.ks$Nt + 1))
    {
      pb$tick()
      Sys.sleep(1 / (2*para.ks$Nt))
      if(para.ks$heavy.tail){
        R <- para.ks$R * diag(para.ks$n)
        xi =  rMvn(para.ks$Ensemble.size,
                   mu = as.vector(para.ks$mu),
                   L = L_Q_upper,
                   Cov = F) #%>% t()data$lambda[1]para.ks$lambda[i - 1]
        # }
        e =  rMvn(para.ks$Ensemble.size,
                  mu = rep(0, para.ks$n),
                  L = chol(R))# %>% t()para.ks
        # v = rt(para.ks$Ensemble.size, df = para.ks$T.df - 2)
      }else{
        xi = rMvn(para.ks$Ensemble.size,
                  as.vector(para.ks$mu), L_Q_upper,
                  Cov = F)
        e = matrix(rnorm(para.ks$n*para.ks$Ensemble.size,
                         0, sqrt(para.ks$tau2)),
                   ncol = para.ks$Ensemble.size)
      }
      xp[1:(i -1), , ] = xf[1:(i -1), , ]
      # 预测
      forecast <-  mat_vector_multi(para.ks$Mt, xf[i - 1, , ])
      
      xp[i,,] = (forecast + xi) #先验
      yp =  para.ks$Hs %*% xp[i, , ] + e
      
      Phat = (gpuR::cov(gpuR::t(gpuR::vclMatrix(xp[i, , ]))) * para.ks$rho.space)
      
      HP <- Hs %*% Phat %*% gpuR::t(Hs)
      if(para.ks$heavy.tail){
        HP = HP + gpuR::vclMatrix(R)#/para.ks$tau2[i - 1]
      }else{
        HP = HP + gpuR::vclMatrix(para.ks$R)
      }
      Sk <- gpuR::t(Hs) %*% gpuR::solve(HP)
      
      for (j in max(1, i - para.ks$ct):i)
      {
        rho = rho.time[paste(i - j)] * para.ks$rho.space
        
        Phat.ji = rho * gpuR::cov(gpuR::t(gpuR::vclMatrix(xp[j, , ]))
                                  , gpuR::t(gpuR::vclMatrix(xp[i, , ])))
        
        # if(para.ks$heavy.tail){Phat.ji = Phat.ji*(para.ks$T.df - 2)/para.ks$T.df}
        K = Phat.ji %*% Sk
        
        xf[j, , ] = (xp[j, , ] + K %*%
                       gpuR::vclMatrix(para.ks$y[i - 1,] - yp)) %>% as.matrix()
        rm(K, Phat.ji, rho)
      }
      rm(Sk, Phat)
    }
    # cat("........................................\n")
  }
  rm(Hs);gc()
  # cat("........................................\n")
  # cat("   Smoother has being over;\nBeing going to calculate the Covariance matrix of ensemble! \n")
  ###############################################################
  ###############################################################
  
  # system.time({
  Sigma2.Smoother = sapply(seq_len(para.ks$Nt + 1)
                           , FUN = function(t){
                             gpuR::cov(gpuR::t(gpuR::vclMatrix(xf[t, , ])))%>%
                               as.matrix() * rho.time[paste(0)] * para.ks$rho.space
                           }, simplify = "array")
  # diag(Sigma2.Smoother[,,1]) <- diag(Sigma2.Smoother[,,1])
  
  Pt_t_1 = sapply(X = seq_len(para.ks$Nt)
                  , FUN = function(t) {
                    gpuR::cov(gpuR::t(gpuR::vclMatrix(xf[t, , ]))
                              , gpuR::t(gpuR::vclMatrix(xf[t + 1, , ])))%>%
                      as.matrix() * rho.time[paste(1)] * para.ks$rho.space
                  }, simplify = "array")
  
  # })
  # all.equal(Sigma2.Smoother, P1)
  Mu.Smoother = apply(xf, 1:2, mean)
  # Mu.pri = apply(xp, 1:2, mean)
  rm(xp, forecast)
  gc()
  # cat("           TaEnks'th calculation is over!!! \n")
  # cat("................................................................. \n")
  return(list(
    # EnXp = xp
    EnXf = xf,
    # Xp = Mu.pri
    Xf = Mu.Smoother,
    Pt = Sigma2.Smoother,
    Pt_t_1 = Pt_t_1,
    rho.time = rho.time,
    rho.space = para.ks$rho.space,
    bandKernel = para.ks$bandKernel
  )
  )
}


