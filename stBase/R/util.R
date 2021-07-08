# solve variance by row
colVar <- function(data){
  # colNames <- colnames(data)
  temp = NULL
  for(i in 1:ncol(data)){
    temp <- rbind(temp, var(data[, i]))
  }
  # colnames(temp) <- colNames
  return(temp)
}
# transform spatial coordinates by WGS84
spCoords.transform <- function(loc, col = c("LON", "LAT"), 
                               colname = c("LON_X", "LAT_Y"), 
                               method = 2){
  nc <- ncol(loc)
  if(method == 1){
    d <- loc[, col]
    coordinates(d) <- col
    proj4string(d) <- CRS("+proj=longlat +datum=WGS84")# WGS 84
    CRS.new <- CRS("+proj=utm +zone=51 ellps=WGS84")
    sp.coord.trans <- spTransform(d, CRS.new)
    
    
    loc <- cbind(loc, sp.coord.trans@coords[,1], sp.coord.trans@coords[, 2])
  }else{
    LCC <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
    # transforming latitude and longitude using the Lambert conformal projection
    sp.coord <- SpatialPoints(coords = loc[, col],
                              proj4string = CRS("+proj=longlat +ellps=WGS84"))
    
    sp.coord.trans <- spTransform(sp.coord,  CRS(LCC))
    #   
    #   #  %>% setnames(c("LON", "LAT"), c("LON_X", "LAT_Y"))
    sp.coord.trans <- as.data.frame(sp.coord.trans)
    #   colnames(sp.coord.trans) <- colname
    # loc = as.data.frame(loc)
    loc <- cbind(loc, sp.coord.trans[, 1], sp.coord.trans[, 2])
    
  }
  
  # grid.coords$LON_X = sp.coord.trans@coords[,1]
  # grid.coords$LAT_Y = d.ch1903@coords[,2]
  colnames(loc)[(nc + 1):(nc + 2)] <- colname
  return(loc)
}

# solve quantile
quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}

# log likelihood for gamma distribution
log_gamma <- function(x, E, a, b, n)
{
  x = exp(x)
  loglike <- n * x[1]*log(x[2]) - n*log(gamma(x[1])) +
    (x[1] - 1) *sum(psigamma(a) - log(b)) - x[2] * sum(E)
  # sum(psigamma(a) - log(b))
  return(-loglike)
}

# sampling program for theta_2, zeta-squared and zeta_0 - squared
importance <- function(x, p, q, IS.err, iter = 0)
{
  p <-  exp(p/mean(abs(p))) #1.0*max(abs(p)) + p  #
  if(iter == 0){ Q <- p}else{
    B <- sum(p)
    Q <- p*(q)/B
  }
  A <- sum( (t(x)%*%(Q))) # numerator
  B <- sum(Q)
  k <- which.max(Q)[1]
  Pro = Q/B
  return(list(
    EIS = ifelse(iter == 0, A/B, x[k]),
    Rs = x , Gw = Pro,
    x = x , Pro = Pro))
}

# log likelihood of the HDCM in the first-level
loglik <- function(data, para, Ks, ds = 1e-2, sp,
                   heavy.tail = FALSE)
{
  # Mt = para$theta1$E_theta1 * ds* exp(- data$BAUs.Dist^2 / para$theta2$E_theta2)
  
  
  if(heavy.tail){
    Hv <- sapply(seq(1, data$n), function(s)
      para$alpha$E_alpha*(apply(Ks$Xf[2:(data$Nt + 1), ], 1,
                                FUN = '%*%', data$Hs[s, ])), simplify = "array")
    
    # epsilon = data$Y_ts -  X_ts_Transf(data$Nt, data$n,
    #                     data$X_ts, para$beta$E_beta) - Hv
    
    epsilon = data$Y_ts -  X_ts_Transf(Nt = data$Nt,
                                       X_ts = data$X_ts,
                                       beta = para$beta$E_betaX)
    if(!is.null(data$Z_ts)){
      epsilon = epsilon -
        X_ts_Transf(Nt = data$Nt, 
                    X_ts = data$Z_ts,
                    beta = para$beta$E_betaZ) %*% gpuR::t(data$Hs) 
    }
    
    Loglik.y.1 <-  - data$Nt * sum(log(para$Obs.tau2$E_tau2))
    Loglik.y.2 <- (para$Obs.tau2$E_tau2%*% colSums(epsilon^2))[1]
    # -2*sum(dgamma(para$Obs.tau2$E_tau2, para$Obs.tau2$a, para$Obs.tau2$b, log = T))
    # }
  }else{
    t = 2:(data$Nt + 1)
    Hv <- sapply(t, function(t) para$alpha$E_alpha * data$Hs %*% Ks$Xf[t, ],
                 simplify = "matrix") %>%  gpuR::t()
    
    Y_ts = data$Y_ts -  X_ts_Transf(data$Nt, 
                                    data$X_ts,
                                    para$beta$E_betaX) - Hv
    if(!is.null(data$Z_ts)){
      Y_ts <- Y_ts - X_ts_Transf(Nt = data$Nt,
                                 X_ts = data$Z_ts,
                                 beta = para$beta$E_betaZ) %*% gpuR::t(data$Hs) 
    }
    # ylikelihood
    Loglik.y.1 = data$n * data$Nt * log(para$Obs.tau2$E_tau2)
    
    t = 1:data$Nt
    Loglik.y.2 <- sapply(t, function(t) Y_ts[t, ] %*% Y_ts[t, ]
                         , simplify = "matrix") %>% sum()
    Loglik.y.2 <- Loglik.y.2/(para$Obs.tau2$E_tau2)
  }
  Loglik <- (- Loglik.y.1 - Loglik.y.2)/2
  return(Loglik)
}

# sample from a Gaussian process
rMvn <- function(n, mu = 0, L, Cov = T){
  # L is an upper triangular matrix
  p <- length(mu)
  if(any(is.na(match(dim(L), p)))){
    stop("Dimension problem!")
  }
  if(Cov){
    # sample by covariance matrix 
    R <- gpuR::t(gpuR::vclMatrix(rnorm(n * p), nrow = n, ncol = p) %*%
                   gpuR::vclMatrix(L) +
                   gpuR::vclMatrix(rep(mu, rep(n, p)),
                                   nrow = n, ncol = p)) %>%
      as.matrix()  # pxn
  }else{
    # sample by precision matrix ： L is an inverse from upper triangular matrix
    # R = (gpuR::vclMatrix((L)) %*% gpuR::vclMatrix(rnorm(p * n),
    #           nrow = p , ncol = n)+
    #        t(matrix(rep(mu, rep(n, p)), ncol = p))) %>%
    #   as.matrix()  # pxn
    # 
    R = (L %*% matrix(rnorm(n * p), nrow = p, ncol = n) +
          t(matrix(rep(mu, rep(n, p)), ncol = p))) %>%
    as.matrix()
  }
  return(R)
}

# tapering function 
taper <- function(distance, radius) Wendland(distance, radius, 1, 1)

#
chunk <- function(x,n)
{
  split(x, factor(sort(rank(x)%%n)))
}

#
X_ts_Transf <- function(Nt, X_ts, beta)
{
  t = seq_len(Nt)
  x_ts = sapply(t, function(t) gpuR::t(X_ts[, , t]) %*% beta
                , simplify = "matrix") %>% gpuR::t()
  return(x_ts)
}
mat_vector_multi <- function(mat, vec)
{
  Multi = apply(t(vec), 1, FUN = '%*%', as.matrix(gpuR::t(mat)))
  return(Multi)
}
################################################################################
IniEnKs <- function(data, para, heavy.tail = FALSE,
                    Ensemble.size = 100, 
                    spTaper = NULL, 
                    bandKernel = NULL,
                    ds = 1e-2, ct = 1)
{
  # library(reticulate)
  # conda_install("r-reticulate", "scipy")
  # sp <- import("scipy")
  # py_install("scipy")
  if(is.null(spTaper)){
    spTaper <- fields::Wendland(data$BAUs.Dist
                                  , theta = max(data$BAUs.Dist)*cs
                                  , dimension = 1, k = 1)
  }
  if(is.null(bandKernel)){
    bandKernel <- spTaper
    bandKernel[bandKernel > 0] <- 1
  }
  Mt = gpuR::vclMatrix(para$theta1$E_theta1 * ds *
                 exp(- data$BAUs.Dist^2 / para$theta2$E_theta2))*bandKernel
  
  # lamda <- eigen(C_Mt)$value
  if(para$theta1$E_theta1 >= 1e-3)
  {
    t1 <- proc.time()
    # lamda <- gpuR::eigen(Mt, 1, only.values = T)$values
    
    lamda <- as.vector(eigen(as.matrix(Mt), 1, 
                             only.values = T,
                             EISPACK = T)$values)
    t2 <- proc.time()
    
    if(is.complex(lamda[1]))
    {
      Lamda <- max(Mod(as.vector(lamda)))
      cat("\n ........................ \n")
      cat("Eigen of the Mt:", Lamda)
      # print(t2 - t1)
      cat("\n ........................ \n")
      if(Lamda >1) Mt <- Mt/(1.1 *  abs(Lamda))
    }else{
      Lamda <- max(lamda)
      if(Lamda >1) Mt <- Mt/(1.1 * abs(Lamda))
      cat("\n ........................ \n")
      cat("Eigen of the Mt:", Lamda)
      # print(t2 - t1)
      cat("\n ........................ \n")
    }
  }
  ###################################################################
  
  if(heavy.tail)
  {
    R = para$Obs.tau2$b/(para$Obs.tau2$a - 1)#para$obs_tau2#
  }else{
    R = diag(data$n)*(para$Obs.tau2$E_tau2)
  }
  Y_ts = data$Y_ts -  X_ts_Transf(Nt = data$Nt, 
                                  X_ts = data$X_ts,
                                  beta = para$beta$E_betaX)
  if(!is.null(data$Z_ts))
  {
    Y_ts =  Y_ts - X_ts_Transf(Nt = data$Nt,
                               X_ts = data$Z_ts,
                               beta = para$beta$E_betaZ) %*% gpuR::t(data$Hs) 
  }
  ###################################################################
  para.ks <- list(y = Y_ts
                  , Nt = data$Nt
                  , n = data$n
                  , N.BAUs = data$N.BAUs
                  , Mt = Mt                
                  , Q =  para$Q$E_Q            
                  , mu = rep(0, data$N.BAUs) 
                  , mu0 = rep(0, data$N.BAUs)
                  , Q0 = para$Q0$E_Q0            # 
                  , Hs = para$alpha$E_alpha*data$Hs   
                  , R =  R
                  , tau2 = para$Obs.tau2$E_tau2  # nugget effect 
                  , Ensemble.size = Ensemble.size
                  , BAUs.Dist = data$BAUs.Dist
                  , ct = ct
                  , heavy.tail = heavy.tail
                  , rho.space = spTaper
                  , bandKernel = bandKernel
  )
  rm(Mt, R, Y_ts, bandKernel);gc()
  return(para.ks)
}
################################################################################
covSfun.Gaussian <- function(Nt, Ks)
{
  temp = 0
  t = seq_len(Nt - 1) + 1
  # t1 <- proc.time()
  P1_t_1 = sapply(t, function(t)
    gpuR::tcrossprod(Ks$Xf[t, ]) +  Ks$Pt[, , t]
    , simplify="array")
  # t2 <- proc.time()
  # t2 - t1
  for(t in seq_along(t))
  {
    temp <- temp + P1_t_1[,, t]
  }
  
  S11 = temp + gpuR::tcrossprod(Ks$Xf[Nt + 1, ]) + Ks$Pt[, , Nt + 1]
  S00 = temp + gpuR::tcrossprod(Ks$Xf[1, ]) + Ks$Pt[, , 1]
  
  t = seq_len(Nt)
  
  temp10 = sapply(t, function(t)
    gpuR::tcrossprod(Ks$Xf[t + 1, ], Ks$Xf[t, ]) + gpuR::t(Ks$Pt_t_1[, , t])
    , simplify="array")
  temp01 = sapply(t, function(t)
    gpuR::tcrossprod(Ks$Xf[t, ], Ks$Xf[t + 1, ]) +  Ks$Pt_t_1[, , t]
    , simplify="array")
  S10 = S01 <- 0
  for(t in seq_along(t))
  {
    S10 <- S10 + temp10[,, t]
    S01 <- S01 + temp01[,, t]
  }
  
  return(S = list(S00 = S00, S11 = S11, S10 = S10, S01 = S01))
}

# no closed form parameters
# f_theta <- function(X)
# {
#   M <- gpuR::vclMatrix(ds*exp(- (BAUs.Dist^2 / theta2[X])))*bandKernel
#   Q_M <- (gpuR::vclMatrix(Q) %*% M)
#   f = post_theta1_mu * sum(as.vector(gpuR::diag(Q_M %*%
#                                                   gpuR::vclMatrix(S$S01)))) -
#     (post_theta1_sigma2 + post_theta1_mu^2) *
#     sum(as.vector(gpuR::diag( gpuR::t(M) %*% Q_M %*%
#                                 gpuR::vclMatrix(S$S00))))/2
#   rm(M);gc()
#   return(f)
# }

f_theta <- function(X)
{
  M <- as(ds*exp(- (BAUs.Dist^2 / theta2[X]))*bandKernel, 
          "sparseMatrix")
  Q_M <- Q %*% M
  # f = post_theta1_mu * sum(Matrix::diag(Q_M %*% 
  #                       (S$S01))) -
  #   (post_theta1_sigma2 + post_theta1_mu^2) *
  #   sum((Matrix::diag(crossprod(M, Q_M) %*%
  #                               (S$S00))))/2
  f = post_theta1_mu * Trace_Muti(as.matrix(Q_M), S$S01) -
    (post_theta1_sigma2 + post_theta1_mu^2) *
    Trace_Muti(as.matrix(Matrix::t(M) %*% Q_M),
                        S$S00)/2
  rm(M, Q_M);gc()
  return(f)
}





f_theta2 <- function(theta, ds, m, n, D, bandKernel, Q, 
                     S00, S01, theta1, sigTheta1, 
                     nThreads){
  
  storage.mode(theta) <- "double"
  storage.mode(ds) <- "double"
  storage.mode(m) <- "integer"
  storage.mode(n) <- "integer"

  storage.mode(D) <- "double"
  storage.mode(bandKernel) <- "double"
  storage.mode(Q) <- "double"
  storage.mode(S00) <- "double"
  storage.mode(S01) <- "double"
  storage.mode(theta1) <- "double"
  storage.mode(sigTheta1) <- "double"
  storage.mode(nThreads) <- "integer"
  
  return(theta2_fun_C(m, n, theta, ds,
                       D, bandKernel,
                       Q, S00, S01, 
                       theta1, sigTheta1, 
                       nThreads))
}




###########################################################################
f_k <- function(X)
{
  Matrix::diag(Adj.Mat) <- Matrix::diag(Adj.Mat) + Rs[X]^2
  # L_Q <- Matrix::Cholesky(Adj.Mat)
  # D_L_Q <- Matrix::determinant(L_Q, 
  #                              logarithm = TRUE)$modulus
  # rm(L_Q)
  
  D_L_Q <- Matrix::determinant(Adj.Mat, logarithm = TRUE)$modulus
  f <- Nt * D_L_Q/2 + Rs[X]^2*TeK
  rm(D_L_Q, TeK);gc()
  return(f)
}
###########################################################################
f_k0 <- function(X)
{
  Matrix::diag(Adj.Mat) <-  Matrix::diag(Adj.Mat) + Rs[X]^2
  # L_Q0 <- Matrix::Cholesky(Adj.Mat)
  # D_L_Q0 <- Matrix::determinant(L_Q0, logarithm = TRUE)$modulus
  # 
  # rm(L_Q0)#Mu0 +
  # f <-  D_L_Q0 - Rs[X]^2 * trace_Pt0/2
  # 
  
  D_L_Q0 <- Matrix::determinant(Adj.Mat, logarithm = TRUE)$modulus
  
  f <- D_L_Q0/2 - Rs[X]^2 * trace_Pt0/2
  rm(D_L_Q0);gc()
  return(f)
}

###################################################################
# 污染物质关联CMAQ
Species_Link_CMAQ <- function( Site , CMAQ_Table = CMAQ_PM25  # CMAQ 数据表
                               , CMAQ_Site
                               , Distance_Threshold  #选择加权的距离阈值
)
{ 
  
  library(spBayes)
  library(data.table)
  library(Hmisc)
  library(plyr)
  setDT(CMAQ_Table)
  D <- iDist(setDF(Site[, c("LON_X","LAT_Y")]), 
             setDF(CMAQ_Site[, c("LON_X","LAT_Y")]))
  colnames(D) <- unique(CMAQ_Site$CMAQ_ID)
  for(s in 1:nrow(D))
  {
    if(s == 1)
    {
      Index <- which.min( D[s,] <= Distance_Threshold )
      # Index <- which.min( D[s,])
      Index1 <- as.numeric(rownames(as.data.frame(Index)))
      temp <- CMAQ_Table[CMAQ_ID %in% Index1, ]
      
      # 计算权重
      Inv_Dist = 1/D[s, Index]
      inf_index <- which(is.infinite(Inv_Dist))
      
      if(length(inf_index) > 0 )
      {
        Inv_Dist[inf_index] <- 
          max(Inv_Dist[-inf_index])+0.05
      }
      W <- data.frame(CMAQ_ID = Index1, Inv_Dist = Inv_Dist 
                      , Weight = Inv_Dist/sum(Inv_Dist))
      
      # 关联权重 
      temp <- sqldf("select a.*, b.Weight
            from temp a
            left join W b
            on a.CMAQ_ID = b.CMAQ_ID")
      # 汇总
      temp = plyr::ddply(temp,
                         .(DATE_TIME),
                         plyr::summarize,
                         PM25_AVG_CMAQ = sum(CMAQ_PM25 * Weight),    
                         .progress = "text"
      )
      temp$LON_X = Site$LON_X[s]
      temp$LAT_Y = Site$LAT_Y[s]
    }else{
      Index <- which( D[s,] <= Distance_Threshold )
      # Index <- which.min( D[s,])
      Index1 <- as.numeric(rownames(as.data.frame(Index)))
      da <- CMAQ_Table[CMAQ_ID %in% Index1, ]
      
      # 计算权重
      Inv_Dist = 1/D[s, Index]
      inf_index <- which(is.infinite(Inv_Dist))
      
      if(length(inf_index) > 0 )
      {
        Inv_Dist[inf_index] <- 
          max(Inv_Dist[-inf_index])+0.05
      }
      W <- data.frame(CMAQ_ID = Index1, Inv_Dist = Inv_Dist 
                      , Weight = Inv_Dist/sum(Inv_Dist))
      
      # 关联权重 
      da <- sqldf("select a.*, b.Weight
            from da a
            left join W b
            on a.CMAQ_ID = b.CMAQ_ID")
      # 汇总
      da = plyr::ddply(da,
                       .(DATE_TIME),
                       plyr::summarize,
                       PM25_AVG_CMAQ = sum(CMAQ_PM25 * Weight),    
                       .progress = "text"
      )
      da$LON_X = Site$LON_X[s]
      da$LAT_Y = Site$LAT_Y[s]
      
      
      temp = rbind( temp, da)
    }
    # cat("已关联第", s, "个站点(关联点数:", length(Index),").\n")
  }
  return(temp)
}    
# 站点
# load("G:/Mirror/github/HDCM/data/Site.Rdata")
# load("./data/CMAQ_PM25.RData")
# CMAQ_Table = CMAQ_PM25
# load("G:/Mirror/github/HDCM/data/CMAQ_Site.RData")
# CMAQ_Site <- cmaq_site
# 
# 
# temp <- Species_Link_CMAQ( Site         # 站点
#                            , CMAQ_Table = CMAQ_PM25  # CMAQ 数据表
#                            , CMAQ_Site = CMAQ_Site
#                            , Distance_Threshold = 30  #选择加权的距离阈值
# )


# da <- temp %>%
#   filter(year(temp$DATE_TIME) %in% c(2015),
#          month(temp$DATE_TIME) %in% c(Month),
#          day(temp$DATE_TIME) %in% c(Day)) %>%
#   setorder(DATE_TIME, LON_X, LAT_Y)

# Predict.Location.trans[1:5,  ]
# load("./data/Model_Base_Table_Update.RData")
# da1 <- Model_Base_Table_Update %>%
#   filter(YEAR %in% c(2015),
#          MONTH %in% c(Month),
#          DAY %in% c(Day))%>%  
#   setorder(LON_X, LAT_Y,DATE_TIME)

# da2 <- Predict.Location.trans%>%
#   filter(year(Predict.Location.trans$DATE_TIME) %in% c(2015),
#          month(Predict.Location.trans$DATE_TIME) %in% c(Month),
#          day(Predict.Location.trans$DATE_TIME) %in% c(Day)) %>%
#   setorder(DATE_TIME, LON_X, LAT_Y)
# da[1:5,]
# da2[1:5, ]
# 
# 
# 
# plot(da[, 2], da1$CMAQ_PM25_30)
# 
# 
# 
# Model_Base_Table_2021 <- da1[1:5, -c("CMAQ_PM25_30")] %>%
#   left_join(da[1:5,], by = c("DATE_TIME", "LON_X","LAT_Y")) %>%
#   setnames("PM25_AVG_CMAQ", "CMAQ_PM25_30")
# 
# 
# 
# Model_Base_Table_2021 <- Model_Base_Table_Update[, -c("CMAQ_PM25_30")] %>%
#   left_join(temp, by = c("DATE_TIME", "LON_X","LAT_Y")) %>%
#   setnames("PM25_AVG_CMAQ", "CMAQ_PM25_30")
# 
# save(Model_Base_Table_2021, file = "./data/Model_Base_Table_2021.RData")


spT.validation <- function (z = NULL, zhat = NULL, zhat.Ens = NULL, names = FALSE) 
{
  VMSE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sum(u^2)/length(u), 4)
  }
  ## root mean square error
  RMSE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sqrt(sum(u^2)/length(u)), 4)
  }
  MAE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  MAPE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))/abs(z)
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u) * 100, 4)
  }
  ## normalised mean gross error
  NMGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    NMGE <- abs(c(x[[2]] - x[[1]]))/sum(x[[1]])
    round(mean(NMGE), 4)
  }  
  
  BIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  rBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u) * mean(z, na.rm = TRUE)), 4)
  }
  ## normalised mean bias
  nBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    round(sum(x, na.rm = T)*100/sum(z, na.rm = T), 4)
  }
  
  rMSEP <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat, na.rm = TRUE) - z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), 4)
  }
  ## correlation coefficient
  Coef <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    Coef <- suppressWarnings(cor(x[[2]], x[[1]])) ## when SD=0; will return NA
    round(Coef, 4)
  }
  
  ##  Coefficient of Efficiency
  COE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    COE <- 1 - sum(abs(x[[2]] - x[[1]])) / sum(abs(x[[1]] - mean(x[[1]])))
    round(COE, 4) 
  }
  # FAC2 <- function(z, zhat) {
  #   z <- z %>% as.data.frame() %>% as.vector()
  #   zhat <- zhat %>% as.data.frame() %>% as.vector()
  #   index <- !is.na(z)
  #   round(mean(zhat[index]/z[index]), 4)
  # }
  ## fraction within a factor of two
  FAC2 <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    ratio <- x[[2]]/x[[1]]
    len <- length(ratio)
    if (len > 0) {
      FAC2 <- length(which(ratio >= 0.5 & ratio <= 2)) / len
    } else {
      FAC2 <- NA
    }
    round(FAC2, 4)
  }
  ## mean gross error
  MGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    MGE <- round(mean(abs(x[[2]] - x[[1]])), 4)
  }
  
  ##  Index of Agreement
  IOA <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    
    LHS <- sum(abs(x[[2]] - x[[1]]))
    RHS <- 2 * sum(abs(x[[1]] - mean(x[[1]])))
    
    if (LHS <= RHS) IOA <- 1 - LHS / RHS else IOA <- RHS / LHS - 1
    round(IOA, 4)
  }
  
  
  CRPS <<- function(z, zhat){
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    
    Nt <- nrow(z)
    n <- ncol(z) 
    
    CRPS <- vector()
    for(p in 1:n)
    {
      test.index <- which(is.na(z[, p]))%>% as.vector()
      if(length(test.index) == 0 ) {
        z0 <- z[, p] %>% as.vector()
        zt <- zhat[, p] %>% as.vector()
        s <- sd(zt)
      }else{
        z0 <- z[-test.index, p] %>% as.vector()
        zt <- zhat[-test.index, p] %>% as.vector()
        s <- sd(zt)
      }
      CRPS[p] <- verification::crps(z0, cbind(zt, s))$CRPS
    }
    CRPS <- round(mean(CRPS, na.rm = T), 4) 
  } 
  
  CRPS.ES <- function(z, zhat, zhat.Ens) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    
    Nt <- nrow(z)
    n <- ncol(z) 
    # cat("Nt = ", Nt)
    if(!is.null(zhat.Ens))
    {
      # z <- as.matrix(z)
      test.index <- which(is.na(z), arr.ind = T)
      if(nrow(test.index) >0){ z[test.index] <- 0 }
      CRPS <- ES <- matrix(NA, nrow = nrow(z), ncol = n)
      for(t in 1:Nt)
      {
        CRPS[t, ] <- crps_sample(z[t, ], zhat.Ens[t,,])
        ES[t, ] <- es_sample(z[t, ], zhat.Ens[t,,])
      }
      CRPS[test.index] <-  ES[test.index] <- NA
      #crps
      CRPS <- round(mean(CRPS, na.rm = T), 4)
      ES <- round(mean(ES, na.rm = T), 4)
    }else{
      # z <- as.matrix(z)
      # zhat <- as.matrix(zhat)
      
      CRPS <- ES <- vector()
      for(t in 1:Nt)
      {
        test.index <- which(is.na(z[t, ]))%>% as.vector()
        if(length(test.index) == 0 ) {
          z0 <- z[t, ] %>% as.vector()
          zt <- zhat[t, ] %>% as.matrix()
        }else{
          z0 <- z[t, -test.index] %>% as.vector()
          zt <- zhat[t, -test.index] %>% as.matrix()
        }
        CRPS[t] <- mean(crps_sample(z0, zt))
        ES[t] <-  es_sample(z0, zt)
      }
      CRPS <- round(mean(CRPS, na.rm = T), 4)
      ES <- round(mean(ES, na.rm = T), 4)
    }
    return(list(CRPS = CRPS, ES = ES))
  }
  if (names == TRUE) {
    cat("##\n Mean Squared Error (MSE) \n Root Mean Squared Error (RMSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Normalized Mean Error(NME: %) \n Bias (BIAS) \n Relative Bias (rBIAS) \n Normalized Mean Bias(NMB: %) \n Relative Mean Separation (rMSEP) \n Correlation coefficient (Coef.) \n The fraction of predictions within a factor of two of observations (FAC2) \n Continuous ranked probability score (CRPS) based on sample mean and standard deviation \n CRPS based on empirical distribution function \n Energy score (ES) based on empirical distribution function \n##\n")
  }
  if((!is.null(z))&(!is.null(zhat))){
    out <- NULL
    out$MSE <- VMSE((z), (zhat))
    out$RMSE <- RMSE((z), (zhat))
    out$MAE <- MAE((z), (zhat))
    out$MAPE <- MAPE((z), (zhat))
    out$BIAS <- BIAS((z), (zhat))
    out$rBIAS <- rBIAS((z), (zhat))
    out$Coef <- Coef((z), (zhat))
    out$FAC2 <- FAC2(z, zhat)
    
    if(!is.null(zhat.Ens))
    {
      R <- CRPS.ES(z, zhat, zhat.Ens)
    }else{
      R <- CRPS.ES(z, zhat, zhat.Ens = NULL)
    }
    out$CRPS <- CRPS(z, zhat)
    out$CRPS.sample <- R$CRPS
    out$ES.sample <- R$ES
    out$NMGE <- NMGE((z), (zhat))
    out$MGE <- MGE(z, zhat)
    out$IOA <- IOA(z, zhat)
    out$NMB <- nBIAS((z), (zhat))
    out$rMSEP <- rMSEP((z), (zhat))
    out$COE <- COE(z, zhat)
    
    unlist(out)
  }else{
    return(0)
  }
}

Validation.Group.City <- function(table, col = c(2, 3), by = c(1)) 
{
  library(data.table)
  library(tidyr)
  library(plyr)
  spT.validation()
  if(length(col) != 2){col = col[1:2]}
  Da <- table[, c(by, col)]
  setDT(Da)
  by <- 1:length(by)
  col <- c((length(by) + 1):(length(by) + 2))
  # setnames(Da, by, paste0("G", length(by)))
  setnames(Da, col, c("Y", "Y.Pred"))
  # setDF(table)#table[table$Miss_Flag == "FALSE", ]
  result = ddply(Da,
                 .(CITY), #as.quoted(paste0("G", length(by)))
                 plyr::summarize,
                 RMSE = RMSE(Y, Y.Pred),
                 CRPS = CRPS(Y, Y.Pred),
                 Corr = Coef(Y, Y.Pred),
                 FAC2 = FAC2(Y, Y.Pred),
                 # RMSE = round(sqrt(mean((PM25.Pred - PM25)^2)),3),
                 # MB   = round(mean(PM25.Pred - PM25), 3),
                 # NMB   = round(sum(PM25.Pred - PM25)*100/sum(PM25), 3), 
                 # NME  = round(sum(abs(PM25.Pred - PM25))*100/sum(PM25),3),
                 # CORR_AFT  = round(cor(PM25.Pred, PM25), 3),
                 .progress = "text"
  ) %>% setorderv("CITY")
  # setDT(result)
  # setorderv(result, c("RMSE"), 1)
  result <- rbind(result, data.frame(CITY = "AVG"
                                     , RMSE = round(mean(result$RMSE), 4)
                                     , CRPS = round(mean(result$CRPS), 4)
                                     , Corr  = round(mean(result$Corr), 4)
                                     , FAC2 = round(mean(result$FAC2), 4)
                                     # , MB = round(mean(result$MB),3)
                                     # , NMB = round(mean(result$NMB),3)
                                     # , NME = round(mean(result$NME),3)
  )
  ) 
  return(result)
}
# results
test.HDCM <- function(test, Ks, PIU, seed = 1234){
 
  if(is.null(Ks$EnXs)){
    Ks$EnXs <- Ks$EnXf
  }
  N <- nrow(test$H)
  Nt <- dim(Ks$EnXs)[1] - 1
  n.Enseble <- dim(Ks$EnXs)[3]
  Y.test <- array(0, dim = c(Nt, N, n.Enseble),
                  dimnames = list(rownames(test$Y_ts),
                                  colnames(test$Y_ts),
                                  paste0("En.", 1:n.Enseble)
                  ))
  # Xbeta <- X_ts_Transf(Nt, test$X_ts, PIU$beta$E_betaX)
  # Fix.Var <- matrix(NA, nrow = Nt, ncol = N)
  
  H <- as(test$H, "sparseMatrix")
  s.e <- seq_len(n.Enseble)
  s.n <- seq_len(N)
  W_ts.sample <- matrix(0, nrow = N, ncol = n.Enseble)
  beta <- array(0, dim = c(length(PIU$beta$E_betaX), N, n.Enseble),
                dimnames = list(c(paste0("X", 1:length(PIU$beta$E_betaX))),
                                1:N,
                                paste0("En.", 1:n.Enseble)
                ))
  set.seed(seed)
  for(t in 1:Nt)
  {
    W_ts <- sapply(s.e,
                function(s.e) H %*% Ks$EnXs[t + 1, , s.e],
                simplify = "array")

    for(i in s.e){
      W_ts.sample[, i] <- as.vector(W_ts[[i]])
    }
    beta.sample <- mvnfast::rmvn(N*n.Enseble,
                                 mu = PIU$beta$E_betaX,
                                 sigma = PIU$beta$betaX.Sigma2)
    for(i in 1:length(PIU$beta$E_betaX)){
      if(nrow(beta.sample) == n.Enseble){
        beta[i,,] <- t(matrix(beta.sample[, i], nrow = n.Enseble, ncol = N))  
      }else{
        beta[i,,] <- matrix(beta.sample[, i], nrow = N, ncol = n.Enseble)
      }
    }

    betaX.sample <- sapply(s.n, function(s.n)
      t(test$X_ts[, s.n, t]) %*% beta[, s.n, ],
      simplify = "array")
    betaX.sample <- Matrix::t(betaX.sample[1, , ])

    y.pred.sample  <-  (betaX.sample + W_ts.sample +
                          rnorm(N*n.Enseble, 0,
                                sqrt(PIU$Obs.tau2$E_tau2)))
    y.pred.sample <- ifelse(y.pred.sample < 0, 0, y.pred.sample)
    Y.test[t, ,] <- ifelse(y.pred.sample^2 > 800, 800, y.pred.sample^2)
    
    
    # W_ts <- apply(test$H, 1, FUN = '%*%', Ks$EnXs[t + 1, ,])
    # beta.sample <- mvnfast::rmvn(n.Enseble,
    #                              mu = PIU$beta$E_betaX,
    #                              sigma = PIU$beta$betaX.Sigma2) %>% t()
    # Xbeta_ts <- apply(t(test$X_ts[,, t]), 1, FUN = '%*%', beta.sample)
    # temp.value <- t(Xbeta_ts + W_ts) + rnorm(N*n.Enseble, 0,
    #                                          sqrt(PIU$Obs.tau2$E_tau2))
    # temp.value <- ifelse(temp.value < 0, 0, temp.value)
    # Y.test[t, ,] <- ifelse(temp.value^2 > 800, 800, temp.value^2)
  }
  return(Y.test)
}
######################################################################
######################################################################
load_HDCM_fun <- function(file, tab, seed = 1234)
{
  # if(length(tab) == 4){
  #   load(paste0(file, "/data/", tab[1], "_", 
  #               tab[2], "_", tab[3], "_",
  #               tab[4], "_Summary.RData"))
  # }else{
  #   load(paste0(file, "/data/", tab[1], 
  #               tab[2], "_", tab[3], "_Summary.RData"))
  # }
  Tab <- list.files(paste0(file, "/data/"))
  HDCM <- Tab[grepl(tab[1], Tab) & grepl(paste0(tab[2], "_"), Tab) & grepl(tab[3], Tab)]
  load(paste0(file, "/data/", HDCM[1]))
  
  ######################################################################
  ######################################################################
  method <- "ensemble"
  
  Y.test <- test.HDCM(test = test, Ks = Ks, PIU = PIU, seed = seed)
  HDCM.test <- (apply(Y.test, c(1, 2), quant))
  ######################################################################
  ######################################################################
  Pred.L25 <- HDCM.test[1,,] %>% as.data.frame()
  Pred.Median <- HDCM.test[2,,] %>% as.data.frame()
  Pred.U95 <- HDCM.test[3,,] %>% as.data.frame()
  Pred.Mean <- (apply(Y.test, c(1, 2), mean)) %>% as.data.frame()
  
  # Err <- spT.validation(test$Y_ts_true, as.matrix(Pred.Median), NULL, T)
  # print(Err)
  
  Pred.L25$DATE_TIME <- rownames(Pred.L25) %>% as.Date()
  Pred.Median$DATE_TIME <- rownames(Pred.Median) %>% as.Date()
  Pred.U95$DATE_TIME <- rownames(Pred.U95) %>% as.Date()
  Pred.Mean$DATE_TIME <- rownames(Pred.Mean) %>% as.Date()
  
  ######################################################################
  ######################################################################
  REAL_PM25 <- test$Y_ts_true %>% as.data.frame()
  REAL_PM25$DATE_TIME <- rownames(REAL_PM25) %>% as.Date()
  
  REAL_PM25 <- gather(
    data = REAL_PM25,      #待转换的数据集名称
    key = "SITEID",       #转换后的分类字段名称（维度）
    value = "REAL_PM25" ,    #转换后的度量值名称
    -DATE_TIME 
  ) 
  
  Pred.L25 <- gather(
    data = Pred.L25,      #待转换的数据集名称
    key = "SITEID",       #转换后的分类字段名称（维度）
    value = "Pred.L25" ,    #转换后的度量值名称
    -DATE_TIME 
  ) 
  
  Pred.Median <- gather(
    data = Pred.Median,      #待转换的数据集名称
    key = "SITEID",       #转换后的分类字段名称（维度）
    value = "Pred.Median" ,    #转换后的度量值名称
    -DATE_TIME 
  ) 
  
  Pred.U95 <- gather(
    data = Pred.U95,      #待转换的数据集名称
    key = "SITEID",       #转换后的分类字段名称（维度）
    value = "Pred.U95" ,    #转换后的度量值名称
    -DATE_TIME 
  ) 
  
  Pred.Mean <- gather(
    data = Pred.Mean,      #待转换的数据集名称
    key = "SITEID",       #转换后的分类字段名称（维度）
    value = "Pred.Mean" ,    #转换后的度量值名称
    -DATE_TIME 
  ) 
  ######################################################################
  ######################################################################
  HDCM2.Ens <- Pred.L25 %>% left_join(Pred.Median, by = c("SITEID", "DATE_TIME")) %>% 
    left_join(Pred.U95, by = c("SITEID", "DATE_TIME")) %>% 
    left_join(REAL_PM25, by = c("SITEID", "DATE_TIME"))%>% 
    left_join(Pred.Mean, by = c("SITEID", "DATE_TIME"))
  HDCM2.Ens$CITY <- CITY.Name
  HDCM2.Ens$Model <- "HDCM.Ens"
  # save(HDCM2.Ens, file = paste0(file, "/data/HDCM2.Ens_", CITY.Name, 
  #                              ifelse(min(month) %in% c(201511, 201512, 201601),
  #                                     "_W", "_S"), ".RData"))
  return(HDCM2.Ens)
}
######################################################################
load_SVC_fun <- function(file, tab){
  load(paste0(file, "/data/", tab[1], "_", tab[2], ".RData"))
  SVC <- as.data.frame(SVC)
  SVC <- SVC[SVC$CITY %in% tab[3], ]
  
  # Err1 <- spT.validation(setDF(SVC1[, "True_REAL_PM25"]), 
  #                        setDF(SVC1[, "PM25.Pred"]),
  #                        NULL, F)
  # print(Err1)
  # Err <- spT.validation(SVC$True_REAL_PM25, SVC$PM25.Pred, NULL, F)
  # print(Err)
  
  # SVC1$Model <- "SVC1";
  # SVC1$Pred.Mean <- SVC1$PM25.Pred
  
  SVC$Model <- "SVC"
  SVC$Pred.Mean <- SVC$PM25.Pred
  
  
  # SVC1 <- SVC1[, c("DATE_TIME" , "SITEID", "PM25.L25", "PM25.Pred",
  #                  "PM25.U95", "PM25", "Pred.Mean", "CITY", "Model")]
  SVC <- SVC[, c("DATE_TIME" , "SITEID", "PM25.L25", "PM25.Pred",
                 "PM25.U95", "PM25", "Pred.Mean", "CITY", "Model")]
  
  colnames(SVC) <- c("DATE_TIME", "SITEID", "Pred.L25",
                     "Pred.Median", "Pred.U95", "REAL_PM25",
                     "Pred.Mean", "CITY","Model")
  # save(SVC, file = paste0(file, "./data/", to, ".RData"))
  return(SVC)
}
######################################################################
######################################################################
load_HDCM_SVC_data <- function(file, hdcm_tab, svc_tab, month, day, seed = 1234){
  data("SiteData", package = "stBase")
  ######################################################################
  ######################################################################
  CMAQ_PM25 <- Model_Base_Table_2021[, c("DATE_TIME", "SITEID",
                                         "CMAQ_PM25", "REAL_PM25","CITY", 
                                         "YEAR_MONTH",
                                         "DAY")] %>%
    filter(YEAR_MONTH %in% month, DAY %in% day,
           CITY %in% hdcm_tab[3]
    ) %>% dplyr::select(DATE_TIME, SITEID, CMAQ_PM25, REAL_PM25, CITY) %>%
    setnames("CMAQ_PM25", "Pred.Median")
  CMAQ_PM25$Pred.Mean <- CMAQ_PM25$Pred.Median
  CMAQ_PM25$Pred.L25 = CMAQ_PM25$Pred.L25 = CMAQ_PM25$Pred.U95 = NA
  CMAQ_PM25$Model <- "CMAQ"   
  
  HDCM.Ens <- load_HDCM_fun(file = file, tab = hdcm_tab, seed = seed)
  CMAQ_PM25 <- CMAQ_PM25 %>% setcolorder(base::colnames(HDCM.Ens))
  # colnames(CMAQ_PM2)
  ######################################################################
  ######################################################################
  HDCM.Ens <- HDCM.Ens %>%
    dplyr::filter(
      as.numeric(paste0(year(DATE_TIME), ifelse(month(DATE_TIME) < 10, 
                                                paste0("0", month(DATE_TIME)),
                                                month(DATE_TIME)))) %in% month,
      day(DATE_TIME) %in% day)
  ######################################################################
  ######################################################################
  SVC <- load_SVC_fun(file = file, tab = svc_tab) %>% dplyr::filter(
    as.numeric(paste0(year(DATE_TIME), ifelse(month(DATE_TIME) < 10, 
                                              paste0("0", month(DATE_TIME)),
                                              month(DATE_TIME)))) %in% month,
    day(DATE_TIME) %in% day)
  
  #-----------------------------------------------------
  CMAQ_PM25.Err1 <- spT.validation(CMAQ_PM25$REAL_PM25,  
                                   CMAQ_PM25$Pred.Median, 
                                   NULL, F)
  cat("\n CMAQ_PM25.Err.Median = \n")
  print(CMAQ_PM25.Err1)  
  
  HDCM.Ens.Err1 <- spT.validation(HDCM.Ens$REAL_PM25,  
                                  HDCM.Ens$Pred.Median, 
                                  NULL, F)
  cat("\n HDCM.Ens.Err.Median = \n")
  print(HDCM.Ens.Err1)
  HDCM.Ens.Err2 <- spT.validation(HDCM.Ens$REAL_PM25,  
                                  HDCM.Ens$Pred.Mean, 
                                  NULL, F)
  cat("\n HDCM.Ens.Err.Mean = \n")
  print(HDCM.Ens.Err2)
  #-----------------------------------------------------
  
  # SVC1.Err1 <- spT.validation(SVC1[, "REAL_PM25"],  
  #                             SVC1[, "Pred.Median"], 
  #                             NULL, F)
  # cat("SVC1.Err1 = \n")
  # print(SVC1.Err1)
  SVC.Err1 <- spT.validation(SVC[, "REAL_PM25"],  
                             SVC[, "Pred.Median"], 
                             NULL, F)
  cat("\n SVC.Err.Median = \n")
  print(SVC.Err1)
  
  #------------------------------------------------------------
  CMAQ_PM25 <- ddply(CMAQ_PM25
                     , .(CITY, DATE_TIME, Model)
                     , .fun = plyr::summarize
                     , REAL_PM25 = mean(REAL_PM25, na.rm = TRUE)
                     , Pred.L25 = mean(Pred.L25, na.rm = TRUE)
                     , Pred.U95 = mean(Pred.U95, na.rm = TRUE)
                     , Pred.Median = mean(Pred.Median, na.rm = TRUE)
                     , Pred.Mean = mean(Pred.Mean, na.rm = TRUE)
                     , .progress = "text"
  )
  CMAQ_PM25$RMSE <- as.vector(CMAQ_PM25.Err1[2]) %>% round(3)
  CMAQ_PM25$MAPE <- as.vector(CMAQ_PM25.Err1[4]) %>% round(3)
  CMAQ_PM25$Corr <- as.vector(CMAQ_PM25.Err1[7]) %>% round(3)
  CMAQ_PM25$FAC2 <- as.vector(CMAQ_PM25.Err1[8]) %>% round(3)
  CMAQ_PM25$CRPS <- as.vector(CMAQ_PM25.Err1[9]) %>% round(3)
  
  SVC <- ddply(SVC
               , .(CITY, DATE_TIME, Model)
               , .fun = plyr::summarize
               , REAL_PM25 = mean(REAL_PM25, na.rm = TRUE)
               , Pred.L25 = mean(Pred.L25, na.rm = TRUE)
               , Pred.U95 = mean(Pred.U95, na.rm = TRUE)
               , Pred.Median = mean(Pred.Median, na.rm = TRUE)
               , Pred.Mean = mean(Pred.Mean, na.rm = TRUE)
               , .progress = "text"
  )
  SVC$RMSE <- as.vector(SVC.Err1[2]) %>% round(3)
  SVC$MAPE <- as.vector(SVC.Err1[4]) %>% round(3)
  SVC$Corr <- as.vector(SVC.Err1[7]) %>% round(3)
  SVC$FAC2 <- as.vector(SVC.Err1[8]) %>% round(3)
  SVC$CRPS <- as.vector(SVC.Err1[9]) %>% round(3)
  
  
  HDCM.Ens <- ddply(HDCM.Ens
                    , .(CITY, DATE_TIME, Model)
                    , .fun = plyr::summarize
                    # , LON = mean(LON, na.rm = TRUE)
                    # , LAT = mean(LAT, na.rm = TRUE)
                    # , FAC2 = mean(Pred.Median/REAL_PM25, na.rm = TRUE)
                    , REAL_PM25 = mean(REAL_PM25, na.rm = TRUE)
                    , Pred.L25 = mean(Pred.L25, na.rm = TRUE)
                    , Pred.U95 = mean(Pred.U95, na.rm = TRUE)
                    , Pred.Median = mean(Pred.Median, na.rm = TRUE)
                    , Pred.Mean = mean(Pred.Mean, na.rm = TRUE)
                    , .progress = "text"
  )
  HDCM.Ens$RMSE <- as.vector(HDCM.Ens.Err1[2]) %>% round(3)
  HDCM.Ens$MAPE <- as.vector(HDCM.Ens.Err1[4]) %>% round(3)
  HDCM.Ens$Corr <- as.vector(HDCM.Ens.Err1[7]) %>% round(3)
  HDCM.Ens$FAC2 <- as.vector(HDCM.Ens.Err1[8]) %>% round(3)
  HDCM.Ens$CRPS <- as.vector(HDCM.Ens.Err1[9]) %>% round(3)
  return(list(CMAQ_PM25 = CMAQ_PM25, HDCM.Ens = HDCM.Ens, SVC = SVC))
}

######################################################################
Credible_Beijing <- function(file, hdcm_tab, svc_tab, month, day, seed)
{
  # if(length(hdcm_tab) == 4){CITY.Name <- hdcm_tab[4]}else{
    CITY.Name <- hdcm_tab[3]
  # }
  
  ######################################################################
  HDCM_SVC <- load_HDCM_SVC_data(file = file, hdcm_tab = hdcm_tab, day = day,
                                 month = month, seed = seed,
                                 svc_tab = svc_tab)
  
  # save(HDCM_SVC, file = paste0(file, "/data/HDCM_SVC_", CITY.Name,
  #                                   ifelse(min(month) %in% c(201511, 201512, 201601),
  #                                          "_W", "_S"), ".RData"))
  CMAQ_PM25 <- HDCM_SVC$CMAQ_PM25
  HDCM.Ens <- HDCM_SVC$HDCM.Ens
  SVC <- HDCM_SVC$SVC
  ######################################################################
  ######################################################################
  Real_PM25 <- CMAQ_PM25
  Real_PM25$Pred.Median = CMAQ_PM25$REAL_PM25
  Real_PM25$Model = "Observation"
  Real_PM25[, c(5, 6, 8, 9, 10, 11, 12, 13)] = NA
  CMAQ_PM25[, c(5, 6)] = NA
  ######################################################################
  Da <- rbind(Real_PM25, CMAQ_PM25, SVC, HDCM.Ens)
  ######################################################################
  #             plot  credible  interval
  ######################################################################
  label <- c("Observation", "CMAQ", "SVC with 95% interval", "HDCM with 95% interval")
  Bcol <- c("black", "grey80", "grey50", "red")
  ######################################################################
  Sys.setlocale("LC_TIME", "English")
  ######################################################################
  size = c(0.4, 0.8, 4)
  UP = max(Da$Pred.U95, na.rm = T) + 10
  scal = 60
  time <- unique(Da$DATE_TIME)
  # Da$fill <- ifelse(Da$Model %in% "HDCM1.Ens", "red",
  #                    ifelse(Da$Model %in% "HDCM2.Ens", "red",
  #                           ifelse(Da$Model %in% "SVC1", "gray",
  #                                  ifelse(Da$Model %in% "SVC2", "black",
  #                                         ifelse(Da$Model %in% "CMAQ", "transparent", "yellow")))))
  alpha <- c("1", "2", "3", "4")
  Da$alpha <- ifelse(Da$Model %in% "HDCM.Ens", alpha[4],
                     ifelse(Da$Model %in% "SVC", alpha[3],
                            ifelse(Da$Model %in% "CMAQ",
                                   alpha[2], alpha[1])))
  S = c("1", "2", "3", "4")
  Da$size <- ifelse(Da$Model %in% "HDCM.Ens", S[4],
                    ifelse(Da$Model %in% "SVC", S[3],
                           ifelse(Da$Model %in% "CMAQ",
                                  S[2], S[1])))
  
  Da$Model <- ifelse(Da$Model %in% "HDCM.Ens", "HDCM",
                     ifelse(Da$Model %in% "SVC", "SVC",
                            ifelse(Da$Model %in% "CMAQ", "CMAQ", "Observation")))
  
  
  Da$Model <- ordered(Da$Model, levels = c("Observation", "CMAQ", "SVC", "HDCM"))
  ls = c(1.3, 1.3, 1.3, 1.3)
  ######################################################################
  ######################################################################
  p <- ggplot(Da, aes(DATE_TIME, group  = Model)) +
    geom_ribbon(aes(ymin = Pred.L25, ymax = Pred.U95,
                    linetype = Model, fill = Model), 
                alpha = 0.3, size = size[1]) +
    geom_line(aes(y = Pred.Median, linetype = Model, col = Model, 
                  alpha = alpha, size = size)) +
    theme_bw() + 
    scale_colour_manual(name = '', values = c("gray40","gray50", 
                                              "#fcbba1", "#43a2ca"),
                        labels = label) +
    scale_fill_manual(name = '', values = c("transparent",
                                            "transparent",
                                            "#fcbba1",
                                            "#43a2ca"),
                      labels = label) +
    scale_alpha_manual(values = c(1, 1, 1, 1)) + 
    scale_size_manual(values = ls) + 
    scale_linetype_manual(name = '',values=c("dashed", "solid",
                                             "solid","solid"),
                          labels = label) +
    scale_x_continuous(expand = c(1e-2, 0)
                       , breaks  = unique(Da$DATE_TIME)[c(seq(1, 31, 7))]
                       , labels = c("Jun 01", "Jun 08","Jun 15","Jun 22","Jun 29")
                       # , labels = c("Nov 01", "Nov 08","Nov 15","Nov 22","Nov 29")
                       # , labels = c("Dec 01", "Dec 08","Dec 15","Dec 22","Dec 29")
    ) +
    scale_y_continuous(limits = c(0, UP), 
                       breaks  = seq(0, UP, scal), 
                       labels = seq(0, UP, scal)) +
    labs( x = "Date", fill = "",
          y = TeX("Observed and predicted PM$_{2.5}$ $( mg/m^3 )$")) +
    theme(axis.text = element_text(size = 18, colour = "black")
          ,axis.text.x = element_text(hjust = 0.25, size = 16, colour = "black") 
          , axis.title = element_text(size = 18, colour = "black")
          , legend.title = element_text(size = 18, colour = "black")
          , legend.text = element_text(size = 16, colour = "black")
          # , legend.title = element_blank()
          , legend.background = element_rect(colour = 'transparent'
                                             , fill = 'transparent')
          , legend.key.width = unit(5,"line")
          , panel.grid.major = element_blank()
          , panel.grid.minor = element_blank()
          , legend.position = c(0.4, 0.9)
          # , legend.margin = margin(t = -0.1, unit='cm')
          , strip.text =  element_text(size = 16, colour = "black")) +
    guides(col = guide_legend(override.aes = list(size = ls),
                              nrow = 2, byrow = TRUE),
           alpha = F, size = F)
  return(p)
}

Credible_Tangshan <- function(file, hdcm_tab, svc_tab, month, day, seed)
{
  # if(length(hdcm_tab) == 4){CITY.Name <- hdcm_tab[4]}else{
    CITY.Name <- hdcm_tab[3]
  # }
  
  ######################################################################
  HDCM_SVC <- load_HDCM_SVC_data(file = file, hdcm_tab = hdcm_tab, day = day,
                                 month = month, seed = seed,
                                 svc_tab = svc_tab)
  
  # save(HDCM_SVC, file = paste0(file, "/data/HDCM_SVC_", CITY.Name,
  #                                   ifelse(min(month) %in% c(201511, 201512, 201601),
  #                                          "_W", "_S"), ".RData"))
  CMAQ_PM25 <- HDCM_SVC$CMAQ_PM25
  HDCM.Ens <- HDCM_SVC$HDCM.Ens
  SVC <- HDCM_SVC$SVC
  ######################################################################
  #   2
  ######################################################################
  Real_PM25 <- CMAQ_PM25
  Real_PM25$Pred.Median <- CMAQ_PM25$REAL_PM25
  Real_PM25$Model <- "Observation"
  Real_PM25[, c(5, 6, 8, 9, 10, 11, 12, 13)] <- NA
  CMAQ_PM25[, c(5, 6)] <- NA
  ######################################################################
  Da <- rbind(Real_PM25, CMAQ_PM25, SVC, HDCM.Ens)
  ######################################################################
  #             plot  credible  interval
  ######################################################################
  label <- c("Observation", "CMAQ", "SVC with 95% interval", "HDCM with 95% interval")
  Bcol <- c("black", "grey80", "grey50", "red")
  ######################################################################
  Sys.setlocale("LC_TIME", "English")
  ######################################################################
  size = c(0.4, 0.8, 4)
  UP = max(Da$Pred.U95, na.rm = T) + 10
  scal = 200
  time <- unique(Da$DATE_TIME)
  # Da$fill <- ifelse(Da$Model %in% "HDCM1.Ens", "red",
  #                    ifelse(Da$Model %in% "HDCM2.Ens", "red",
  #                           ifelse(Da$Model %in% "SVC1", "gray",
  #                                  ifelse(Da$Model %in% "SVC2", "black",
  #                                         ifelse(Da$Model %in% "CMAQ", "transparent", "yellow")))))
  alpha <- c("1", "2", "3", "4")
  Da$alpha <- ifelse(Da$Model %in% "HDCM.Ens", alpha[4],
                     ifelse(Da$Model %in% "SVC", alpha[3],
                            ifelse(Da$Model %in% "CMAQ",
                                   alpha[2], alpha[1])))
  S = c("1", "2", "3", "4")
  Da$size <- ifelse(Da$Model %in% "HDCM.Ens", S[4],
                    ifelse(Da$Model %in% "SVC", S[3],
                           ifelse(Da$Model %in% "CMAQ",
                                  S[2], S[1])))
  
  Da$Model <- ifelse(Da$Model %in% "HDCM.Ens", "HDCM",
                     ifelse(Da$Model %in% "SVC", "SVC",
                            ifelse(Da$Model %in% "CMAQ", "CMAQ", "Observation")))
  
  
  Da$Model <- ordered(Da$Model, levels = c("Observation", "CMAQ", "SVC", "HDCM"))
  ls = c(1.3, 1.3, 1.3, 1.3)
  ######################################################################
  ######################################################################
  ggplot(Da, aes(DATE_TIME, group  = Model)) +
    geom_ribbon(aes(ymin = Pred.L25, ymax = Pred.U95,
                    linetype = Model, fill = Model), 
                alpha = 0.3, size = size[1]) +
    geom_line(aes(y = Pred.Median, linetype = Model, col = Model, 
                  alpha = alpha, size = size)) +
    theme_bw() + 
    scale_colour_manual(name = '', values = c("gray40","gray50", 
                                              "#fcbba1", "#43a2ca"),
                        labels = label) +
    scale_fill_manual(name = '', values = c("transparent",
                                            "transparent",
                                            "#fcbba1",
                                            "#43a2ca"),
                      labels = label) +
    scale_alpha_manual(values = c(1, 1, 1, 1)) + 
    scale_size_manual(values = ls) + 
    scale_linetype_manual(name = '',values=c("dashed", "solid",
                                             "solid","solid"),
                          labels = label) +
    scale_x_continuous(expand = c(1e-2, 0)
                       , breaks  = unique(Da$DATE_TIME)[c(seq(1, 31, 7))]
                       # , labels = c("Jun 01", "Jun 08","Jun 15","Jun 22","Jun 29")
                       # , labels = c("Nov 01", "Nov 08","Nov 15","Nov 22","Nov 29")
                       , labels = c("Dec 01", "Dec 08","Dec 15","Dec 22","Dec 29")
    ) +
    scale_y_continuous(limits = c(0, UP), 
                       breaks  = seq(0, UP, scal), 
                       labels = seq(0, UP, scal)) +
    labs( x = "Date", fill = "",
          y = TeX("Observed and predicted PM$_{2.5}$ $( μg/m^3 )$")) +
    theme(axis.text = element_text(size = 18, colour = "black")
          ,axis.text.x = element_text(hjust = 0.25, size = 16, colour = "black") 
          , axis.title = element_text(size = 18, colour = "black")
          , legend.title = element_text(size = 18, colour = "black")
          , legend.text = element_text(size = 16, colour = "black")
          # , legend.title = element_blank()
          , legend.background = element_rect(colour = 'transparent'
                                             , fill = 'transparent')
          , legend.key.width = unit(5,"line")
          , panel.grid.major = element_blank()
          , panel.grid.minor = element_blank()
          , legend.position = c(0.4, 0.9)
          # , legend.margin = margin(t = -0.1, unit='cm')
          , strip.text =  element_text(size = 16, colour = "black")) +
    guides(col = guide_legend(override.aes = list(size = ls),
                              nrow = 2, byrow = TRUE),
           alpha = F, size = F)
  return(p)
}



######################################################################
#                      plot  Cond.quantile 
######################################################################
Cond.quantile <- function (pred, obs, bins = NULL, city = NULL, 
                           thrs = c(0, 0), model = NULL, len = 6) 
{
  # old.par <- par(no.readonly = TRUE)
  # on.exit(par(old.par))
  if (!is.null(bins)) {
    if (min(bins) > min(obs) | max(bins) < max(obs)) {
      warning("Observations outside of bin range. \n")
    }
    if (min(bins) > min(pred) | max(bins) < max(pred)) {
      warning("Forecasts outside of bin range. \n")
    }
  }else {
    dat <- c(obs, pred)
    min.d <- min(dat)
    max.d <- max(dat)
    bins <- seq(floor(min.d), ceiling(max.d), length = len)
  }
  lo <- min(bins)
  hi <- max(bins)
  b <- bins[-length(bins)]
  labs <- b + 0.5 * diff(bins)
  # obs.cut <- cut(obs, breaks = bins, include.lowest = TRUE, 
  #                labels = labs)
  # obs.cut[is.na(obs.cut)] <- labs[1]
  # obs.cut <- as.numeric(as.character(obs.cut))
  frcst.cut <- cut(pred, breaks = bins, include.lowest = TRUE, 
                   labels = labs)
  frcst.cut[is.na(frcst.cut)] <- labs[1]
  frcst.cut <- as.numeric(as.character(frcst.cut))
  n <- length(labs)
  lng <- aggregate(obs, by = list(frcst.cut), length)
  med <- aggregate(obs, by = list(frcst.cut), median)
  q1 <- aggregate(obs, by = list(frcst.cut), quantile, 0.25)
  q2 <- aggregate(obs, by = list(frcst.cut), quantile, 0.75)
  q1$x[lng$x <= thrs[1]] <- NA
  q2$x[lng$x <= thrs[1]] <- NA
  q3 <- aggregate(obs, by = list(frcst.cut), quantile, 0.1)
  q4 <- aggregate(obs, by = list(frcst.cut), quantile, 0.9)
  q3$x[lng$x <= thrs[2]] <- NA
  q4$x[lng$x <= thrs[2]] <- NA
  X <- as.numeric(as.character(med$Group.1))
  quan.data <- rbind(data.frame(x = X, y = med$x, flag = "1", group = "Median")
                     ,data.frame(x = X, y = q1$x, flag = "2", group = "25th/75th Quantiles")
                     ,data.frame(x = X, y = q2$x, flag = "3", group = "25th/75th Quantiles")
                     ,data.frame(x = X, y = q3$x, flag = "4", group = "10th/90th Quantiles")
                     ,data.frame(x = X, y = q4$x, flag = "5", group = "10th/90th Quantiles")) %>%
    as.data.frame()
  quan.data$Model <- model
  quan.data$City <- city
  hist.data <- data.frame(z = frcst.cut, Model = model, City = city)
  return(list(quan.data = quan.data, hist.data = hist.data))
}

CondQuantile_Baoding <- function(file, hdcm_tab, svc_tab, month, day, seed)
{
  # if(length(hdcm_tab) == 4){CITY.Name <- hdcm_tab[4]}else{
    CITY.Name <- hdcm_tab[3]
  # }
  
  ######################################################################
  HDCM_SVC <- load_HDCM_SVC_data(file = file, hdcm_tab = hdcm_tab, day = day,
                                 month = month, seed = seed,
                                 svc_tab = svc_tab)
  
  
  # save(HDCM_SVC, file = paste0(file, "/data/HDCM_SVC_", CITY.Name,
  #                                   ifelse(min(month) %in% c(201511, 201512, 201601),
  #                                          "_W", "_S"), ".RData"))
  CMAQ_PM25 <- HDCM_SVC$CMAQ_PM25
  HDCM.Ens <- HDCM_SVC$HDCM.Ens
  SVC <- HDCM_SVC$SVC
  ######################################################################
  pred <- HDCM.Ens[, "Pred.Median"]
  obs <- HDCM.Ens[, "REAL_PM25"]
  # conditional.quantile(frcst, obs, main = "Sample Conditional Quantile Plot")
  Con.HDCM.Ens <- Cond.quantile(pred, obs, 
                                city =  unique(HDCM.Ens$CITY),
                                model = unique(HDCM.Ens$Model))
  
  
  pred <- SVC[, "Pred.Median"]
  obs <- SVC[, "REAL_PM25"] 
  # conditional.quantile(pred[, 1], obs[, 1], main = "Sample Conditional Quantile Plot")
  Con.SVC <- Cond.quantile(pred, obs,
                           city = unique(SVC$CITY),
                           model = unique(SVC$Model))
  
  
  pred <- CMAQ_PM25[, "Pred.Median"]
  obs <- CMAQ_PM25[, "REAL_PM25"]
  # conditional.quantile(pred, obs, main = "Sample Conditional Quantile Plot")
  Con.CMAQ_PM25 <- Cond.quantile(pred, obs,
                                 city = unique(CMAQ_PM25$CITY),
                                 model = unique(CMAQ_PM25$Model))
  
  
  # da1 <- rbind(Con.CMAQ_PM25$quan.data,
  #              # Con.HDCM1.Ens$quan.data,
  #              Con.HDCM.Ens$quan.data,
  #              # Con.SVC1$quan.data,
  #              Con.SVC$quan.data)
  # da2 <- rbind(Con.CMAQ_PM25$hist.data,
  #              Con.HDCM.Ens$hist.data,
  #              Con.SVC$hist.data)
  
  
  ######################################################################
  #                               plot
  ######################################################################
  library(ggsci)
  # load(paste0(file, CITY.Name, "_", month[1], "_cq.RData"))
  UP <- 180#max(da1$x, da1$y, na.rm = T) + 80
  da1 <- rbind(Con.CMAQ_PM25$quan.data,
               # Con.HDCM1.Ens$quan.data,
               Con.HDCM.Ens$quan.data,
               # Con.SVC1$quan.data,
               Con.SVC$quan.data)
  label <- as.character(unique(da1$group))
  da1$group <- ordered(da1$group,
                       levels=label)
  # da1$method <- ifelse(da2$model%in% "SVC1", "1",
  #                     ifelse(da2$model%in% "SVC2", "3",
  #                            ifelse(da2$model%in% "HDCM1.Ens", "2",
  #                                   ifelse(da2$model%in% "HDCM2.Ens", "4", "0"))))
  # da1$method <- ifelse(da1$model%in% "SVC1", 1,
  #                      ifelse(da1$model%in% "SVC2", 3,
  #                             ifelse(da1$model%in% "HDCM1.Ens", 2,
  #                                    ifelse(da1$model%in% "HDCM2.Ens", 4, 0))))
  da1$Model <- ifelse(da1$Model %in% "HDCM.Ens", "HDCM",
                      ifelse(da1$Model %in% "SVC", "SVC",
                             ifelse(da1$Model %in% "CMAQ", "CMAQ", "Observation")))
  
  da1$Model <- ordered(da1$Model,
                       levels=c("Observation", "CMAQ", "SVC", "HDCM"))
  
  
  Label <- as_labeller(c(`0` = "CMAQ" , `1` = "SVC", `2` = "HDCM"))
  # Label <- as_labeller(c(`0` = "CMAQ" , `1` = "SVC1" 
  #                        , `2` = "HDCM1"
  #                        , `3` = "SVC2", `4` = "HDCM2"))
  size <- c(16, 14)
  ls = c(0.8, 0.5, 0.5)
  {
    p <- ggplot() +  
      geom_line(da1, mapping =aes(x = x, y = y, 
                                  col = group, 
                                  group = flag,
                                  linetype = group,
                                  size = group))+
      geom_abline(intercept = 0, slope = 1, col = "gray", size = 0.5) +
      # geom_histogram(da2, binwidth = 70
      #                , position = "identity"
      #                ,mapping = aes(z, y=(..count..)/5
      #                )) +
      facet_wrap(~ Model, ncol = 3) +
      # facet_grid(city ~ method, scales = "free"
      #            , labeller = labeller(method = Label)
      # )+
      scale_x_continuous(limits = c(0, UP), 
                         expand = c(0, 0),
                         breaks  = seq(0, UP, 50),
                         labels = seq(0, UP, 50)) +
      scale_y_continuous(name = "",
                         limits = c(0, UP), 
                         expand = c(0, 0),
                         breaks  = c(8),
                         position = "right",
                         labels = ""
                         #   function(y) 
                         # { paste0(round(y*50, 0), "")
                         #   }                  
                         , sec.axis = sec_axis(~., 
                                               name = TeX("Observed PM$_{2.5}$ $( μg/m^3 )$"),
                                               breaks  = seq(0, UP, 50),
                                               labels = seq(0, UP, 50)
                                               # function(b) {
                                               #paste0(round(b,0), "")
                                               # }
                         )
      )  + scale_linetype_manual(name='', 
                                 values=c("solid", "longdash", "dotted"),
                                 labels = label) + 
      scale_size_manual(name='', 
                        values = ls, labels = label) +
      scale_colour_manual(name='', 
                          values = c("black", "black" , "black"), #, "#556B2F" , "blue"
                          labels = label)+
      theme_bw() + 
      labs(x = TeX("Predicted PM$_{2.5}$ $( μg/m^3 )$"),
           y.left = TeX("Observed PM$_{2.5}$ $( μg/m^3 )$")) +
      theme(axis.text = element_text(size = size[2], colour = "black")
            , axis.text.x  = element_text(angle = 0)
            , axis.title.x = element_text(size = size[1], colour = "black")
            , axis.title.y = element_text(size = size[1], colour = "black")
            , legend.title = element_text(size = size[1], colour = "black")
            , legend.text = element_text(size = size[2], colour = "black")
            # , legend.title = element_blank()
            , legend.background = element_rect(fill="transparent")
            , legend.key.width = unit(5,"line")
            , panel.grid.major = element_blank()
            , panel.grid.minor = element_blank()
            , legend.position = "top"#c(0.6, 0.1)
            , legend.margin = margin(t = 1, unit='cm')
            , axis.ticks.length.y.right = unit(-0, "cm")
            , strip.text =  element_text(size = size[2], colour = "black")
            # , axis.text.y.right  = element_text(vjust = -2,
            #                                     hjust = -300,
            #                              margin = margin(l = 15, r = 2))
      )  +  guides(linetype = guide_legend(override.aes = 
                                             list(size = ls), nrow=1, byrow=TRUE))
  }
  return(p)
}
CondQuantile_Cangzhou <- function(file, hdcm_tab, svc_tab, month, day, seed)
{
  # if(length(hdcm_tab) == 4){CITY.Name <- hdcm_tab[4]}else{
    CITY.Name <- hdcm_tab[3]
  # }
  
  ######################################################################
  HDCM_SVC <- load_HDCM_SVC_data(file = file, hdcm_tab = hdcm_tab, day = day,
                                 month = month, seed = seed,
                                 svc_tab = svc_tab)

  CMAQ_PM25 <- HDCM_SVC$CMAQ_PM25
  HDCM.Ens <- HDCM_SVC$HDCM.Ens
  SVC <- HDCM_SVC$SVC
  ######################################################################
  pred <- HDCM.Ens[, "Pred.Median"]
  obs <- HDCM.Ens[, "REAL_PM25"]
  # conditional.quantile(frcst, obs, main = "Sample Conditional Quantile Plot")
  Con.HDCM.Ens <- Cond.quantile(pred, obs, 
                                city =  unique(HDCM.Ens$CITY),
                                model = unique(HDCM.Ens$Model))
  
  
  pred <- SVC[, "Pred.Median"]
  obs <- SVC[, "REAL_PM25"] 
  # conditional.quantile(pred[, 1], obs[, 1], main = "Sample Conditional Quantile Plot")
  Con.SVC <- Cond.quantile(pred, obs,
                           city = unique(SVC$CITY),
                           model = unique(SVC$Model))
  
  
  pred <- CMAQ_PM25[, "Pred.Median"]
  obs <- CMAQ_PM25[, "REAL_PM25"]
  # conditional.quantile(pred, obs, main = "Sample Conditional Quantile Plot")
  Con.CMAQ_PM25 <- Cond.quantile(pred, obs,
                                 city = unique(CMAQ_PM25$CITY),
                                 model = unique(CMAQ_PM25$Model))
  ######################################################################
  #                               plot
  ######################################################################
  # CITY.Name <- "Handan"
  month = c(201511, 201512, 201601)
  # load(paste0(file, CITY.Name, "_", month[1], "_cq.RData"))
  da1 <- rbind(Con.CMAQ_PM25$quan.data,
               # Con.HDCM1.Ens$quan.data,
               Con.HDCM.Ens$quan.data,
               # Con.SVC1$quan.data,
               Con.SVC$quan.data)
  UP <- max(da1$x, da1$y, na.rm = T) + 20
  label <- as.character(unique(da1$group))
  da1$group <- ordered(da1$group,
                       levels=label)
  
  da1$Model <- ifelse(da1$Model %in% "HDCM.Ens", "HDCM",
                      ifelse(da1$Model %in% "SVC", "SVC",
                             ifelse(da1$Model %in% "CMAQ", "CMAQ", "Observation")))
  
  da1$Model <- ordered(da1$Model,
                       levels=c("Observation", "CMAQ", "SVC", "HDCM"))
  
  
  Label <- as_labeller(c(`0` = "CMAQ" , `1` = "SVC", `2` = "HDCM"))
  # Label <- as_labeller(c(`0` = "CMAQ" , `1` = "SVC1" 
  #                        , `2` = "HDCM1"
  #                        , `3` = "SVC2", `4` = "HDCM2"))
  size <- c(16, 14)
  ls = c(0.8, 0.5, 0.5)
 p <- ggplot() +  
    geom_line(da1, mapping =aes(x = x, y = y, 
                                col = group, 
                                group = flag,
                                linetype = group,
                                size = group))+
    geom_abline(intercept = 0, slope = 1, col = "gray", size = 0.5) +
    # geom_histogram(da2, binwidth = 70
    #                , position = "identity"
    #                ,mapping = aes(z, y=(..count..)/5
    #                )) +
    facet_wrap(~ Model, ncol = 3) +
    # facet_grid(city ~ method, scales = "free"
    #            , labeller = labeller(method = Label)
    # )+
    scale_x_continuous(limits = c(0, UP), 
                       expand = c(0, 0),
                       breaks  = seq(0, UP, 150),
                       labels = seq(0, UP, 150)) +
    scale_y_continuous(name = "",
                       limits = c(0, UP), 
                       expand = c(0, 0),
                       breaks  = c(8),
                       position = "right",
                       labels = ""
                       #   function(y) 
                       # { paste0(round(y*50, 0), "")
                       #   }                  
                       , sec.axis = sec_axis(~., 
                                             name = TeX("Observed PM$_{2.5}$ $( μg/m^3 )$"),
                                             breaks  = seq(0, UP, 150),
                                             labels = seq(0, UP, 150)
                                             # function(b) {
                                             #paste0(round(b,0), "")
                                             # }
                       )
    )  + scale_linetype_manual(name='', 
                               values=c("solid", "longdash", "dotted"),
                               labels = label) + 
    scale_size_manual(name='', 
                      values = ls, labels = label) +
    scale_colour_manual(name='', 
                        values = c("black", "black" , "black"), #, "#556B2F" , "blue"
                        labels = label)+
    theme_bw() + 
    labs(x = TeX("Predicted PM$_{2.5}$ $( μg/m^3 )$"),
         y.left = TeX("Observed PM$_{2.5}$ $( μg/m^3 )$")) +
    theme(axis.text = element_text(size = size[2], colour = "black")
          , axis.text.x  = element_text(angle = 0)
          , axis.title.x = element_text(size = size[1], colour = "black")
          , axis.title.y = element_text(size = size[1], colour = "black")
          , legend.title = element_text(size = size[1], colour = "black")
          , legend.text = element_text(size = size[2], colour = "black")
          # , legend.title = element_blank()
          , legend.background = element_rect(fill="transparent")
          , legend.key.width = unit(5,"line")
          , panel.grid.major = element_blank()
          , panel.grid.minor = element_blank()
          , legend.position = "top"#c(0.6, 0.1)
          , legend.margin = margin(t = 1, unit='cm')
          , axis.ticks.length.y.right = unit(-0, "cm")
          , strip.text =  element_text(size = size[2], colour = "black")
          # , axis.text.y.right  = element_text(vjust = -2,
          #                                     hjust = -300,
          #                              margin = margin(l = 15, r = 2))
    )  +  guides(linetype = guide_legend(override.aes = 
                                           list(size = ls), nrow=1, byrow=TRUE))
  return(p)
}

######################################################################
#                      plot  scatter 
######################################################################
Round <- function(x, n = 2)
{
  return(format(round(x, n), nsmall = n))
}

Scatter_Xingtai <- function(file, hdcm_tab, svc_tab, month, day, seed)
{
  # if(length(hdcm_tab) == 4){CITY.Name <- hdcm_tab[4]}else{
    CITY.Name <- hdcm_tab[3]
  # }
  
  ######################################################################
  HDCM_SVC <- load_HDCM_SVC_data(file = file, hdcm_tab = hdcm_tab, day = day,
                                 month = month, seed = seed,
                                 svc_tab = svc_tab)
  
  CMAQ_PM25 <- HDCM_SVC$CMAQ_PM25
  HDCM.Ens <- HDCM_SVC$HDCM.Ens
  SVC <- HDCM_SVC$SVC
  ######################################################################
  ######################################################################
  # 2
  ######################################################################
  ######################################################################
  Da <- rbind(CMAQ_PM25, SVC, HDCM.Ens)
  Da$Model <- ifelse(Da$Model %in% "HDCM.Ens", "HDCM",
                     ifelse(Da$Model %in% "SVC", "SVC",
                            ifelse(Da$Model %in% "CMAQ", "CMAQ", NA)))
  
  Da$method <- ifelse(Da$Model%in% "SVC", 1, ifelse(Da$Model%in% "HDCM", 2, 0))
  
  Da$RMSE = paste0("RMSE = ", Round((Da$RMSE)))
  Da$MAPE = paste0("MAPE = ", Round((Da$MAPE)))
  Da$Corr = paste0("Corr = ", Round((Da$Corr)))
  Da$FAC2 = paste0("FAC2 = ", Round((Da$FAC2)))
  Da$CRPS = paste0("CRPS = ", Round((Da$CRPS)))
  
  Da$crps_x <- ifelse(Da$Model %in% "HDCM", 0.85,
                      ifelse(Da$Model %in% "SVC", 0.70,
                             ifelse(Da$Model %in% "CMAQ", 0.75, NA)))
  
  Da$corr_x <- ifelse(Da$Model %in% "HDCM", 9,
                      ifelse(Da$Model %in% "SVC", 9,
                             ifelse(Da$Model %in% "CMAQ", 9, NA)))
  
  
  Da$fac2_x <- ifelse(Da$Model %in% "HDCM", 6,
                      ifelse(Da$Model %in% "SVC", 6,
                             ifelse(Da$Model %in% "CMAQ", 6, NA)))
  Up <- floor(max(Da$REAL_PM25, Da$Pred.Median, na.rm = T)) + 25
  Low <- 2#ceil(min(Da$REAL_PM25, Da$Pred.Median))
  x0 <- -12
  size = c(16, 14, 5)
  Label <- as_labeller(c(`0` = "CMAQ" ,`1` = "SVC", `2` = "HDCM"))
  library(viridis)
  Da$x = 45
  ######################################################################
  ######################################################################
  # 3   plot
  ######################################################################
  ######################################################################
  # library(grid)
  # grid.newpage();grid.draw(roundrectGrob(gp = gpar(lwd = NA)))
  {
    p <- ggplot(Da, aes(x = REAL_PM25, y = Pred.Median)) + 
      geom_point(size = 0.8) +
      # geom_point(aes(col = as.factor(month(DATE_TIME))), size = 0.8) +
      # geom_bin2d(binwidth = c(20, 20), shape = 10) +
      # stat_density_2d(alpha = 0.5, geom = "polygon", 
      #                 aes(fill = (after_stat(level)))) +
      # scale_fill_distiller(palette =  "Spectral", direction = -1)+
      
      # stat_density_2d(geom = "raster",
      #                 aes(fill = stat(density)),
      #                 contour = FALSE) +
      
      geom_abline(slope = 1, color = "black", size = 1) +
      geom_abline(slope = 2, color = "gray", size = 0.5) +
      geom_abline(slope = 0.5, color = "gray", size = .5) +
      # geom_text(y = Up*0.95,
      #           aes(x = x0 + x, label = RMSE,  group = Model),
      #           size = size[3]) +
      # geom_text(y = Up*0.87,
      #           aes(x = x0 + x + crps_x, label = CRPS, group = Model),
      #           size = size[3]) +
      # geom_text(y = Up*0.83,
      #           aes(x = x0 +  x + crps_x, label = MAPE, group = Model),
      #           size = size[3]) +
      
      geom_label(y = Up*0.9,
                 aes(x = x0 + x + fac2_x, label = FAC2, group = Model),
                 size = size[3], label.size = 0) +
      geom_label(y = Up*0.8,
                 aes(x = x0 + x + corr_x, label = Corr, group = Model),
                 size = size[3], label.size = 0) +
      
      annotate(geom="text",x = (Up - Low)*0.97/2,
               y = Up*0.90,
               angle = 60,
               label = "k = 2",
               color = "gray",
               size = 5) +
      annotate(geom="text", x = (Up)*0.93,
               y = Up*0.88,
               angle = 40,
               label = "k = 1",
               size = 5,
               fontface = 1) +
      annotate(geom="text",x = (Up)*0.90,
               y = Up*0.49,
               angle = 25,
               label = "k = 0.5",
               color = "gray",
               size = 5) +
      # facet_wrap(~ Model, ncol = 4) +
      facet_grid(~ method
                 , labeller = labeller(method = Label)
      ) +
      scale_x_continuous(limits = c(0, Up)
                         , expand = c(0, 0)
      ) +
      scale_y_continuous(limits = c(0, Up)
                         , expand = c(0, 0)
      ) +
      labs(color = "", x = TeX("Observed PM$_{2.5}$ ($μg/m^3$)")
           , y = TeX("Predicted PM$_{2.5}$ ($μg/m^3$)")) + theme_bw() + 
      theme( axis.title= element_text(size = size[1], colour = "black")
             , axis.text = element_text(size = size[2], colour = "black")
             # , legend.title = element_text(size = size[1], colour = "black")
             , legend.text= element_text(size = size[2], colour = "black")
             , legend.title = element_blank()
             # , legend.position="top"
             , legend.margin=margin(t = -0.1, unit='cm')
             , legend.background = element_rect(fill="transparent")
             , panel.grid.major = element_blank()
             , panel.grid.minor = element_blank()
             , legend.position = "top"
             # , legend.key.width = unit(1,"line")
             # , legend.key.height = unit(2,"line")
             , strip.text =  element_text(size = size[2], colour = "black")
      ) #+  guides(col = F)
  }
 return(p)
}
Scatter_Handan <- function(file, hdcm_tab, svc_tab, month, day, seed)
{
  # if(length(hdcm_tab) == 4){CITY.Name <- hdcm_tab[4]}else{
    CITY.Name <- hdcm_tab[3]
  # }
  
  ######################################################################
  HDCM_SVC <- load_HDCM_SVC_data(file = file, hdcm_tab = hdcm_tab, day = day,
                                 month = month, seed = seed,
                                 svc_tab = svc_tab)
  CMAQ_PM25 <- HDCM_SVC$CMAQ_PM25
  HDCM.Ens <- HDCM_SVC$HDCM.Ens
  SVC <- HDCM_SVC$SVC
  ######################################################################
  ######################################################################
  # 2
  ######################################################################
  ######################################################################
  Da <- rbind(CMAQ_PM25, SVC, HDCM.Ens)
  Da$Model <- ifelse(Da$Model %in% "HDCM.Ens", "HDCM",
                     ifelse(Da$Model %in% "SVC", "SVC",
                            ifelse(Da$Model %in% "CMAQ", "CMAQ", NA)))
  
  Da$method <- ifelse(Da$Model%in% "SVC", 1, ifelse(Da$Model%in% "HDCM", 2, 0))
  
  Da$RMSE = paste0("RMSE = ", Round((Da$RMSE)))
  Da$MAPE = paste0("MAPE = ", Round((Da$MAPE)))
  Da$Corr = paste0("Corr = ", Round((Da$Corr)))
  Da$FAC2 = paste0("FAC2 = ", Round((Da$FAC2)))
  Da$CRPS = paste0("CRPS = ", Round((Da$CRPS)))
  
  Da$crps_x <- ifelse(Da$Model %in% "HDCM", 0.85,
                      ifelse(Da$Model %in% "SVC", 0.70,
                             ifelse(Da$Model %in% "CMAQ", 0.75, NA)))
  Da$mape_x <- ifelse(Da$Model %in% "HDCM", 0.80,
                      ifelse(Da$Model %in% "SVC", 0.65,
                             ifelse(Da$Model %in% "CMAQ", 0.70, NA)))
  
  corr_x = 30
  Da$corr_x <- ifelse(Da$Model %in% "HDCM", corr_x,
                      ifelse(Da$Model %in% "SVC", corr_x,
                             ifelse(Da$Model %in% "CMAQ", corr_x, NA)))
  
  fac2_x = 20
  Da$fac2_x <- ifelse(Da$Model %in% "HDCM", fac2_x,
                      ifelse(Da$Model %in% "SVC", fac2_x,
                             ifelse(Da$Model %in% "CMAQ", fac2_x, NA)))
  
  
  
  Up <- floor(max(Da$REAL_PM25, Da$Pred.Median, na.rm = T)) + 60
  Low <- 10#ceil(min(Da$REAL_PM25, Da$Pred.Median))
  x0 <- 1
  size = c(16, 14, 5)
  Label <- as_labeller(c(`0` = "CMAQ",`1` = "SVC", `2` = "HDCM"))
  library(viridis)
  Da$x = 130
  ######################################################################
  ######################################################################
  # 3   plot
  ######################################################################
  ######################################################################
  # library(grid)
  # grid.newpage();grid.draw(roundrectGrob(gp = gpar(lwd = NA)))
  {
    p <- ggplot(Da, aes(x = REAL_PM25, y = Pred.Median)) + 
      geom_point(size = 0.8) +
      # geom_point(aes(col = as.factor(month(DATE_TIME))), size = 0.8) +
      # geom_bin2d(binwidth = c(20, 20), shape = 10) +
      # stat_density_2d(alpha = 0.5, geom = "polygon", 
      #                 aes(fill = (after_stat(level)))) +
      # scale_fill_distiller(palette =  "Spectral", direction = -1)+
      
      # stat_density_2d(geom = "raster",
      #                 aes(fill = stat(density)),
      #                 contour = FALSE) +
      
      geom_abline(slope = 1, color = "black", size = 1) +
      geom_abline(slope = 2, color = "gray", size = 0.5) +
      geom_abline(slope = 0.5, color = "gray", size = .5) +
      # geom_text(y = Up*0.95,
      #           aes(x = x0 + x, label = RMSE,  group = Model),
      #           size = size[3]) +
      # geom_text(y = Up*0.87,
      #           aes(x = x0 + x + crps_x, label = CRPS, group = Model),
      #           size = size[3]) +
      # geom_text(y = Up*0.83,
      #           aes(x = x0 +  x + crps_x, label = MAPE, group = Model),
      #           size = size[3]) +
      
      geom_label(y = Up*0.9,
                 aes(x = x0 + x + fac2_x, label = FAC2, group = Model),
                 size = size[3], label.size = 0) +
      geom_label(y = Up*0.8,
                 aes(x = x0 + x + corr_x, label = Corr, group = Model),
                 size = size[3], label.size = 0) +
      
      annotate(geom="text",x = (Up - Low)*0.97/2,
               y = Up*0.90,
               angle = 60,
               label = "k = 2",
               color = "gray",
               size = 5) +
      annotate(geom="text", x = (Up)*0.93,
               y = Up*0.88,
               angle = 40,
               label = "k = 1",
               size = 5,
               fontface = 1) +
      annotate(geom="text",x = (Up)*0.90,
               y = Up*0.49,
               angle = 25,
               label = "k = 0.5",
               color = "gray",
               size = 5) +
      # facet_wrap(~ Model, ncol = 4) +
      facet_grid(~ method
                 , labeller = labeller(method = Label)
      ) +
      scale_x_continuous(limits = c(0, Up)
                         , expand = c(0, 0)
      ) +
      scale_y_continuous(limits = c(0, Up)
                         , expand = c(0, 0)
      ) +
      labs(color = "", x = TeX("Observed PM$_{2.5}$ ($μg/m^3$)")
           , y = TeX("Predicted PM$_{2.5}$ ($μg/m^3$)")) + theme_bw() + 
      theme( axis.title= element_text(size = size[1], colour = "black")
             , axis.text = element_text(size = size[2], colour = "black")
             # , legend.title = element_text(size = size[1], colour = "black")
             , legend.text= element_text(size = size[2], colour = "black")
             , legend.title = element_blank()
             # , legend.position="top"
             # , legend.margin=margin(t = -0.1, unit='cm')
             , legend.background = element_rect(fill="transparent")
             , panel.grid.major = element_blank()
             , panel.grid.minor = element_blank()
             , legend.position = "top"
             # , legend.key.width = unit(1,"line")
             # , legend.key.height = unit(2,"line")
             , strip.text =  element_text(size = size[2], colour = "black")
      ) #+  guides(col = F)
  }
return(p)
}