spT.validation <- function(z, zhat, zhat.Ens = NULL, names = FALSE) 
{
  VMSE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sum(u^2)/length(u), 4)
  }
  ## root mean square error
  RMSE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sqrt(sum(u^2)/length(u)), 10)
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
  Coef <- function(z, zhat) {
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
  FAC2 <- function(z, zhat) {
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
  
  
  CRPS <- function(z, zhat){
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
  out <- NULL
  # out$MSE <- VMSE((z), (zhat))
  out$RMSE <- RMSE((z), (zhat))
  out$MAE <- MAE((z), (zhat))
  #out$MAPE <- MAPE((z), (zhat))
  # out$NMGE <- NMGE((z), (zhat))
  # out$BIAS <- BIAS((z), (zhat))
  # out$rBIAS <- rBIAS((z), (zhat))
  # out$NMB <- nBIAS((z), (zhat))
  # out$rMSEP <- rMSEP((z), (zhat))
  # out$Coef <- Coef((z), (zhat))
  # out$COE <- COE(z, zhat)
  # out$FAC2 <- FAC2(z, zhat)
  # out$MGE <- MGE(z, zhat)
  # out$IOA <- IOA(z, zhat)
  
  # if(!is.null(zhat.Ens))
  # {
  #   R <- CRPS.ES(z, zhat, zhat.Ens)
  # }else{
  #   R <- CRPS.ES(z, zhat, zhat.Ens = NULL)
  # }
  out$CRPS <- CRPS(z, zhat)
  # out$CRPS.sample <- R$CRPS
  # out$ES.sample <- R$ES
  unlist(out)
}

colVar <- function(data){
  colNames <- colnames(data)
  temp = NULL
  for(i in 1:ncol(data)){
    temp <- cbind(temp, var(data[, i]))
  }
  colnames(temp) <- colNames
  return(temp)
}

Rdist <- function(loc1, loc2, covModel = 0, phi = 1, nu = 0,
                  nuUnifb = 0, threads = 10){

  loc1 <- matrix(as.matrix(loc1), ncol = 2)
  loc2 <- matrix(as.matrix(loc2), ncol = 2)
  n1 <- nrow(loc1)
  n2 <- nrow(loc2)
  storage.mode(loc1) <- "double"
  storage.mode(loc2) <- "double" 
  storage.mode(covModel) <- "integer"
  storage.mode(phi) <- "double"
  storage.mode(nu) <- "double"
  storage.mode(nuUnifb) <- "integer"
  
  storage.mode(n1) <- "integer" 
  storage.mode(n2) <- "integer" 
  storage.mode(threads) <- "integer"
  Dvec <- RdistC(loc1, loc2, n1, n2, covModel, phi, nu, nuUnifb, threads)
  
  D <- Dvec$Dist %>% matrix(nrow = n1, ncol = n2, byrow = T)
  Corr <- Dvec$Corr %>% matrix(nrow = n1, ncol = n2, byrow = T)
  return(list(Dist = D, Corr = Corr, range = range(D)))
}

colVar <- function(data){
  colNames <- colnames(data)
  temp = NULL
  for(i in 1:ncol(data)){
    temp <- cbind(temp, var(data[, i]))
  }
  colnames(temp) <- colNames
  return(temp)
}
# mk.n.indx.list <- function(n.indx, n, m){
#     n.indx.list <- vector("list", n)
#     n.indx.list[1] <- NA
#     for(i in 2:n){
#         n.indx.list[[i]] <- n.indx[get.n.indx(i, m)]+1
#     }
#     n.indx.list
# }
Round <- function(x, n =2)
{
  return(format(round(x, n), nsmall = n))
}

RmkUIndx <- function(n, m, nn.indx, nn.indx.lu, search.type){
  
  n.indx <- (1+m)/2*m+(n-m-1)*m
  u.indx <- rep(0, n.indx)
  u.indx.lu <- rep(0, 2*n)
  ui.indx <- rep(0, n.indx)
  
  storage.mode(n) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(nn.indx) <- "integer"
  storage.mode(u.indx) <- "integer"
  storage.mode(u.indx.lu) <- "integer"
  storage.mode(ui.indx) <- "integer"
  storage.mode(nn.indx.lu) <- "integer"
  storage.mode(search.type) <- "integer"
  
  ptm <- proc.time()
  
  out <- mkUIndx(n, m, nn.indx, u.indx, u.indx.lu, ui.indx, nn.indx.lu, search.type)
  
  run.time <- proc.time() - ptm
  
  list("run.time"=run.time, "u.indx"=as.integer(u.indx), "u.indx.lu"=as.integer(u.indx.lu), "ui.indx"=as.integer(ui.indx))
}


RmkNNIndx <- function(coords, m, n.omp.threads=1){
  
  n <- nrow(coords)
  nIndx <- (1+m)/2*m+(n-m-1)*m
  nnIndx <- rep(0, nIndx)
  nnDist <- rep(0, nIndx)
  nnIndxLU <- matrix(0, n, 2)
  
  n <- as.integer(n)
  m <- as.integer(m)
  coords <- as.double(coords)
  nnIndx <- as.integer(nnIndx)
  nnDist <- as.double(nnDist)
  nnIndxLU <- as.integer(nnIndxLU)
  n.omp.threads <- as.integer(n.omp.threads)
  
  ptm <- proc.time()
  
  out <- mkNNIndx(n, m, coords, nnIndx, nnDist, nnIndxLU, n.omp.threads)
  
  run.time <- proc.time() - ptm
  
  list("run.time"=run.time, "nnIndx"=as.integer(nnIndx), "nnDist"=as.double(nnDist), "nnIndxLU"=nnIndxLU)
  
}

RmkNNIndxCB <- function(coords, m, n.omp.threads=1){
  
  n <- nrow(coords)
  nIndx <- (1+m)/2*m+(n-m-1)*m
  nnIndx <- rep(0, nIndx)
  nnDist <- rep(0, nIndx)
  nnIndxLU <- matrix(0, n, 2)
  
  n <- as.integer(n)
  m <- as.integer(m)
  coords <- as.double(coords)
  nnIndx <- as.integer(nnIndx)
  nnDist <- as.double(nnDist)
  nnIndxLU <- as.integer(nnIndxLU)
  n.omp.threads <- as.integer(n.omp.threads)
  
  ptm <- proc.time()
  
  out <- mkNNIndxCB(n, m, coords, nnIndx, nnDist, nnIndxLU, n.omp.threads)
  
  run.time <- proc.time() - ptm
  
  list("run.time"=run.time, "nnIndx"=as.integer(nnIndx), "nnDist"=as.double(nnDist), "nnIndxLU"=nnIndxLU)
}

get.n.indx <- function(i, m){
  i <- i-1
  if(i == 0){
    return(NA)
  }else if(i < m){
    n.indx.i <- i/2*(i-1)
    m.i <- i
    return((n.indx.i+1):((n.indx.i+1)+i-1))
  }else{
    n.indx.i <- m/2*(m-1)+(i-m)*m
    m.i <- m
    return((n.indx.i+1):((n.indx.i+1)+m-1))
  }
}

mk.n.indx.list <- function(n.indx, n, m){
  n.indx.list <- vector("list", n)
  n.indx.list[1] <- NA
  for(i in 2:n){
    n.indx.list[[i]] <- n.indx[get.n.indx(i, m)]+1
  }
  n.indx.list
}

parseFormula <-  function(formula, data, intercept=TRUE, justX=FALSE){
  
  # extract Y, X, and variable names for model formula and frame
  mt <- terms(formula, data=data)
  if(missing(data)) data <- sys.frame(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- mf$justX <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  if (!intercept){
    attributes(mt)$intercept <- 0
  }
  
  # null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)
  X <- as.matrix(X)         # X matrix
  xvars <- dimnames(X)[[2]] # X variable names
  xobs  <- dimnames(X)[[1]] # X observation names
  if (justX){
    Y <- NULL
  }
  else {
    Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
  }
  return(list(Y, X, xvars, xobs))
}

OutputBF <- function(m = 3, coords, 
                     cov.model = 1, 
                     sigmaSq = 1, phi = 1, 
                     nu = 0.5, nuUnifb = 0,
                     nnIndx, nnIndxLU,
                     nThreads = 10){
  coords <- as.matrix(coords)
  n <- nrow(coords)
 

  search.type.names <- c("brute", "cb")
  
  
  n <- as.integer(n)
  m <- as.integer(m)
  coords <- as.double(coords)
  
  cov.model <- as.integer(cov.model)
  nuUnifb <- as.integer(nuUnifb)
  nnIndx <- as.integer(nnIndx)
  nnIndxLU <- as.integer(nnIndxLU)
  
  sigmaSq <- as.double(sigmaSq)
  phi <- as.double(phi)
  nu <- as.double(nu)
  nThreads <- as.integer(nThreads)
  
  BF <- OputBF(n, m, coords, cov.model, nnIndx, nnIndxLU,
               sigmaSq, phi, nu, nuUnifb, nThreads)
  
  # Q <- Matrix::Matrix(BF$Q)
  # lowB <- matrix(BF$lowB, ncol = n)
  # invF <- matrix(BF$invF, ncol = n)
  
  return(list(neiB = BF$neiB, varF = BF$varF,
              invF = BF$invF, lowB = BF$lowB, 
              Q = BF$Q,
              nnIndx = nnIndx,
              nnIndxLU = nnIndxLU))
  
}
#section 2
ConExpeX <- function(X, Z, coords, 
                     Q, Kernel, h = 0.1, 
                     GeomVariable, nThreads = 10){
  n <- nrow(X)
  p <- ncol(X)
  storage.mode(Z) <- "double" 
  storage.mode(coords) <- "double" 
  # storage.mode(Weight) <- "double" 
  storage.mode(Q) <- "double"
  #storage.mode(invF) <- "double"
  
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(Kernel) <- "integer"
  storage.mode(h) <- "double" 
  
  storage.mode(GeomVariable) <- "integer"
  storage.mode(nThreads) <- "integer"
  
  return(ConExp(X, Z, coords,  
                Q, n, p, Kernel,
                h, GeomVariable, nThreads))
}

sFun <- function(Z, S, i, j, index1, ind){
  # m1 <-range(ind[[i]])
  return((Z[i, ]) %*% S[index1 == j, ind[[i]]])
}
sFun.Znull <- function(S, i, j, index1, ind){
  # m1 <-range(ind[[i]])
  return(S[index1 == j, ind[[i]]])
}

devFun <- function(y, S, i, j, index1, ind){
  # m1 <-range(ind[[i]])
  return(t(S[index1 == j, ind[[i]]] %*% y))
}

spLocLKest <- function(y, X, coord, 
                       Q,
                       fs.Density,
                       covModel,
                       h, nu,
                       nuUnifb = 0,
                       ad.width,
                       nThreads){
  # x.sd <- attr(X, "scaled:scale")
  coord <- as.matrix(cbind(coord))
  X <- as.matrix(cbind(X))
  ad.width <- as.matrix(cbind(ad.width))
  mm <- ncol(ad.width)
  n = nrow(X)
  p = ncol(X)
  storage.mode(nuUnifb) <- "integer"
  storage.mode(y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(coord) <- "double"
  storage.mode(Q) <- "double"
  storage.mode(fs.Density) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(covModel) <- "integer"
  storage.mode(h) <- "double" 
  storage.mode(nu) <- "double" 
  storage.mode(ad.width) <- "integer" 
  storage.mode(mm) <- "integer"
  storage.mode(nThreads) <- "integer"
  
  sLLE <- spatial_LLE(y, X, coord, Q, fs.Density, 
                      n, p, covModel, h, nu, nuUnifb,
                      ad.width, mm, nThreads)
  
  index1 <- rep(1:3, p)
  index2 <- rep(1:n, each = n)
  if(p==1){
    alpha <- cbind(sLLE$alpha[index1 == 1, ])
    lon.diff.coef <- cbind(sLLE$alpha[index1 == 2, ])
    lat.diff.coef <- cbind(sLLE$alpha[index1 == 3, ])
  }else{
    alpha <- t(sLLE$alpha[index1 == 1, ])
    lon.diff.coef <- t(sLLE$alpha[index1 == 2, ])
    lat.diff.coef <- t(sLLE$alpha[index1 == 3, ])
  }
 
  
  # coord.0 <- round(coord, 2)

  # rowName <- paste0("s(", coord.0[, 1],
  #                          ",", coord.0[, 2], ")")

  ind = split(1:length(index2), index2)
  
  S0 <- lapply(X = 1:n, sFun, Z = X, S = sLLE$S, j = 1, index1 = index1, ind)
  # lon.S <- lapply(X = 1:n, FUN = devFun, y = y, S = sLLE$S, j = 2, index1 = index1, ind)
  # lat.S <- lapply(X = 1:n, FUN = devFun, y = y, S = sLLE$S, j = 3, index1 = index1, ind)
  
  S0 <- do.call("rbind", S0)
  # lon.S <- do.call("rbind", lon.S)
  # lat.S <- do.call("rbind", lat.S)
  
  #   S0 <- rbind(S0, sLLE$S[index1 == 1, index2 == i])
  #   lon.S <- rbind(S.lon, sLLE$S[index1 == 1, index2 == i])
  #   lat.S <- rbind(S.lat, sLLE$S[index1 == 1, index2 == i])
  # } 
  if(is.null(colnames(X))){
    colName <- c("lon", "lat", paste0("alpha", 1:p))
  }else{
    colName <- c("lon", "lat", colnames(X)) 
  }
  
  alpha <- data.table(coord, alpha)
  lon.diff.coef <- data.table(coord, lon.diff.coef)
  lat.diff.coef <- data.table(coord, lat.diff.coef)
  # rownames(alpha) <- rownames(lon.diff.coef) <- rownames(lat.diff.coef) <- rowName
  colnames(alpha) <- colnames(lon.diff.coef) <- colnames(lat.diff.coef) <- colName
  # setDT();#setDT(lon.diff.coef);setDT(lat.diff.coef);
  
   
  fitted.values = S0 %*% y
  return(list(Alpha = sLLE$alpha, alpha = alpha,
              lon.diff.coef = lon.diff.coef,
              lat.diff.coef = lat.diff.coef,
              S = sLLE$S, S0 = S0, 
              fitted.values = fitted.values,
              residuals = y - fitted.values,
              Q = sLLE$Q
              # ,lon.S = lon.S,
              # lat.S = lat.S
              ))
}

spLocLinKernel <- function(y, X,
                        coord, Q,
                        fs.Density,
                        FitModel = TRUE,
                        pred.X = NULL,
                        pred.coord = NULL,
                        pred.nGrid = 0,
                        covModel = 2, 
                        h = 0.1, 
                        nu = 0,
                        nuUnifb = 0,
                        adapt.n = 0,
                        nThreads = 10){
  coord <- as.matrix(cbind(coord))
  if(FitModel){

     ad.width <- RANN::nn2(coord, coord, k = adapt.n + 1)$nn.idx - 1
    
    Fit = spLocLKest(y = y, X = X, coord = coord, 
                     Q = Q,
                     fs.Density = fs.Density,
                     covModel = covModel, 
                     h = h, nu = nu,
                     nuUnifb = nuUnifb,
                     ad.width = ad.width,
                     nThreads = nThreads)
  }else{
    Fit = NULL
  }
  # if((pred.nGrid > 0)| (!is.null(pred.coord))){
  if((pred.nGrid > 0)|(!is.null(pred.X)&(!is.null(pred.coord)))){
    if(is.null(pred.coord)){
      Min = c(min(coord[, 1]), min(coord[, 2]))
      Max = c(max(coord[, 1]), max(coord[, 2]))
      # }
      # cat(Min)
      pred.coord <- cbind(seq(Min[1], Max[1],, pred.nGrid), 
                          seq(Min[2], Max[2],, pred.nGrid))
      
      pred.coord <- expand.grid(x = pred.coord[, 1],
                                y = pred.coord[, 2])
    }
    coord <- as.matrix(coord)
    pred.coord <- as.matrix(pred.coord)
    # if(adapt.n>0){
    ad.width <- RANN::nn2(coord, pred.coord, 
                          k = adapt.n + 1)$nn.idx - 1
    mm <- ncol(ad.width)
    # }else{
    #   ad.width <- NA
    #   mm <- 0
    # }
    
    X <- as.matrix(cbind(X))
    if(!is.null(pred.X)){
      pred.X <- as.matrix(cbind(pred.X))
    }
    
    n = length(y)
    p = ncol(X)
    m = nrow(pred.coord)
    
    storage.mode(nuUnifb) <- "integer"
    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(coord) <- "double"
    storage.mode(Q) <- "double"
    storage.mode(fs.Density) <- "double"
    storage.mode(pred.coord) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(covModel) <- "integer"
    storage.mode(h) <- "double" 
    storage.mode(nu) <- "double" 
    storage.mode(ad.width) <- "integer"
    storage.mode(mm) <- "integer"
    storage.mode(nThreads) <- "integer"
    
    sLLE <- spatial_LLE_Pred(y, X, coord, Q,
                             fs.Density,
                             n, pred.coord, 
                             m, p, covModel, h, 
                             nu, nuUnifb, ad.width,
                             mm, nThreads)
    
    index1 <- rep(1:3, p)
    index2 <- rep(1:m, each = n)
    if(p==1){
      alpha <- cbind(sLLE$alpha[index1 == 1, ])
      lon.diff.coef <- cbind(sLLE$alpha[index1 == 2, ])
      lat.diff.coef <- cbind(sLLE$alpha[index1 == 3, ])
    }else{
      alpha <- t(sLLE$alpha[index1 == 1, ])
      lon.diff.coef <- t(sLLE$alpha[index1 == 2, ])
      lat.diff.coef <- t(sLLE$alpha[index1 == 3, ])
    }
    
    # pred.coord.0 <- pred.coord
    # pred.coord <- round(pred.coord, 2)
    rowName <- NULL
    
    # rowName <- paste0("s(", pred.coord[, 1],
    #                    ",", pred.coord[, 2], ")")
    
    
    # S = data.table(sLLE$S)
    ind = split(1:length(index2), index2)
    
    # cat(pred.X)
    # pred.X = NULL
    if(!is.null(pred.X)){
      S0 <- lapply(X = 1:m, FUN = sFun, Z = pred.X, S = sLLE$S, j = 1, index1 = index1, ind)
      S0 <- do.call("rbind", S0)
    }else{
      S0 <- lapply(X = 1:m, FUN = sFun.Znull, S = sLLE$S, j = 1, index1 = index1, ind)
    }
    
    # lon.S <- lapply(X = 1:m, FUN = devFun, y = y, S = sLLE$S, j = 2, index1 = index1, ind)
    # lat.S <- lapply(X = 1:m, FUN = devFun, y = y, S = sLLE$S, j = 3, index1 = index1, ind)
    
    
    
    # lon.S <- do.call("rbind", lon.S)
    # lat.S <- do.call("rbind", lat.S)
    
    if(is.null(colnames(X))){
      colName <- c("lon", "lat", paste0("alpha", 1:p))
    }else{
      colName <- c("lon", "lat", colnames(X)) 
    }
    
    alpha <- data.table(pred.coord, alpha)
    lon.diff.coef <- data.table(pred.coord, lon.diff.coef)
    lat.diff.coef <- data.table(pred.coord, lat.diff.coef)
    # rownames(alpha) <- rownames(lon.diff.coef) <- rownames(lat.diff.coef) <- rowName
    colnames(alpha) <- colnames(lon.diff.coef) <- colnames(lat.diff.coef) <- colName
    
    predict <- list(Alpha = sLLE$alpha, alpha = alpha
                    , lon.diff.coef = lon.diff.coef
                    , lat.diff.coef = lat.diff.coef
                    , S = sLLE$S, S0 = S0
                    # , lon.S = lon.S
                    # , lat.S = lat.S
                    , pred.coord = pred.coord)
  }else{
    predict = NULL
  }
  
  return(list(Fit.Model = Fit, predict = predict))
}

spLocPlot <- function(predModel, 
                      nCov = 1, Num = 10,
                      screen = c(60, -60),
                      cols = rainbow(10, alpha=0.5),
                      size= 20, simu = F,
                      path = 'H:/semiSpatailM/figure/',
                      FUN = function(x, y){
  sin(((x^2 + y^2))*pi)}){
  library(plot3D)
  alpha <- setDF(predModel$alpha)
  
  alpha <- data.table(alpha[, c(1, 2, nCov + 2)])
  if(!is.null(colnames(alpha))){
    Name <- colnames(alpha)[3]
  }else{
    Name = paste0("alpha", nCov + 2)
  }
 
  
  
  z.pre <- data.table::dcast(alpha, lon ~ lat, 
                             value.var = Name, 
                             mean) %>% as.matrix()
  lon = as.numeric(z.pre[, 1])
  z.pre <- z.pre[, -1]
  lat = as.numeric(colnames(z.pre))
  nr <- nrow(z.pre)
  nc <- ncol(z.pre)
  setDF(alpha)
  alpha$group = paste0(Name)
  setnames(alpha, Name, "alpha")
  # library(plot3D)
  # z <- outer(pred.coord[,1], pred.coord[,2], f1)
  
  
 
  if(simu){
    
    z <- data.table(lon = predModel$pred.coord[,1], 
                    lat = predModel$pred.coord[,2], 
                    alpha = FUN(predModel$pred.coord[,1],
                                predModel$pred.coord[,2]))
    z.true <- data.table::dcast(z, lon ~ lat, 
                                value.var = "alpha", fun = sum) %>% 
      as.matrix()
    z.true <- z.true[, -1]
    
    
    setDF(z)
    z$group = paste0("alpha", nCov, ".true")
    
    alpha <- rbind(alpha, z)
    alpha$group <- as.factor(alpha$group)
  }

  pdf(file = paste0(path, "/semiTemp", nCov, ".pdf"), width = 10, height = 5)
  if(simu){
    par(mfrow = c(1, 2))
    # contour2D(x = lon, y = lat, z = z.true, lwd = 1 ,
    #           levels = round(seq(min(z.true),
    #            max(z.true), , Num), 3))
    # contour2D(x = lon, y = lat, z = z.pre, lwd = 1
    #           ,levels = round(seq(min(z.pre),
    #                              max(z.pre), , Num), 3)
    #           )
    zcl = range(z.true)
    persp3D(x = lon, y = lat, z = z.true, 
            zlim = c(zcl[1], zcl[2]),
            phi = 20, contour =list(fill = cols, col= cols),
            colkey = list(length = 1, line.clab = 1.0,
                          # adj.clab = 1e-3,
                          # dist = -0.1,
                          # addlines = T,
                          # breaks = seq(zcl[1], zcl[2],, 10),
                          width = 1,
                          shift =0,
                          cex.axis = 0.8,
                          cex.clab = 0.85,
                          side = 4),
            lighting = TRUE, lphi = 90,
            clab = c("","true.fx"),
            bty = "f", plot = T)
    layout = c(2, 1)
   
  }else{
    par(mfrow = c(1, 2))
    z.true = z = NULL
    layout = c(1, 1)
    z.pre.index <- which(is.na(z.pre), arr.ind = T)
    z.pre[z.pre.index] <- 0
    # contour2D(x = lon, y = lat, z = z.pre, lwd = 2,
    #           xlab = "longitude", ylab = "latitude",
    #           levels = round(seq(min(z.pre), 
    #                              max(z.pre), , Num), 3))
    image2D(x = lon, y = lat, z = z.pre,  contour = list(col = "black", 
                                    labcex = 1, lwd = 3, alpha = .6),
            xlab = "longitude", ylab = "latitude",
            levels = round(seq(min(z.pre), 
                   max(z.pre), , Num), 3)
            )
    # zcl = range(z.pre)
    # persp3D(z = z.pre, zlim = c(zcl[1], zcl[2]),
    #         phi = 20,
    #         colkey = list(length = 1, width = 1,
    #                       shift = 0.15,
    #                       cex.axis = 0.8,
    #                       cex.clab = 0.85,
    #                       side = 4),
    #         lighting = TRUE, lphi = 90,
    #         clab = c("","predict.fx"), bty = "f", plot = T)
  }
  zcl = range(z.pre)
 
  persp3D(x = lon, y = lat, z = z.pre, zlim = c(zcl[1], zcl[2]),
          phi = 20, contour =list(fill = cols, col= cols),
          colkey = list(length = 1, line.clab = 1.0,
                        # adj.clab = 1e-3,
                        # dist = -0.1,
                        # addlines = T,
                        # breaks = seq(zcl[1], zcl[2],, 10),
                        width = 1,
                        shift =0,
                        cex.axis = 0.8,
                        cex.clab = 1,
                        side = 4),
          lighting = TRUE, lphi = 90,
          xlab = "longitude", ylab = "latitude",
          zlab = Name,
          clab = c("", paste0(Name, ": predict.f(s)")), bty = "f", plot = T)

  dev.off()
  
  # library(rgl)
  # zlim <- range(A$z.pre)
  # zlen <- zlim[2]-zlim[1] + 1
  # #Assign color values
  # colorlut <- terrain.colors(zlen, alpha = 0)
  # col <- colorlut[A$z.pre-zlim[1]+1]
  # open3d()
  # rgl.surface(x = A$lon, y = A$lat, z = A$z.pre, color = col, back = "lines") 
  # 
  
  
  # contour2D(x = lon, y = lat, z = z.pre, lwd = 2)
  # persp(z = z.true)
  # persp(z = z.pre)
  
  
  # lon.diff <- data.table(pred$lon.diff.coef)
  # lon.diff <- data.table::dcast(lon.diff[, c(1, 2, 3)], lon ~ lat, 
  #                               value.var = "alpha1", mean) %>% 
  #   as.matrix()
  # lon.diff <- lon.diff[, -1]
  
  
 
  # create gradient in x-direction
  
  # add as image with own color key, at bottom
  # image3D(z = zcl[1], colvar = lon.diff, add = T,
  #         colkey = list(length = 1, width = 1,
  #                       shift = -0.15,
  #                       cex.axis = 0.8,
  #                       cex.clab = 0.85,
  #                       side = 1),
  #         clab = c("","gradient"), plot = F)
  # # add contour
  # contour3D(z = zcl[1], colvar = lon.diff, add = TRUE,
  #           col = "black", plot = TRUE)

  # scatterplot3js(predModel$pred.coord[, 1],
  #                predModel$pred.coord[, 2], alpha[, c(nCov + 2)],
  #                phi = 40, #theta = 20,
  #                # color=rainbow(length(A$alpha)),
  #                # colkey = FALSE,
  #                cex = .3, size= .1,
  #                main = "Bivariate Normal")
  
   p <- wireframe(alpha~ lon*lat|group, data = alpha, 
            # groups = group,
            layout = layout,
            # scales = list(arrows = FALSE),
            drape = F,
            colorkey = F,
            # screen = list(z = 30, x = -60) , 
            ylab = list(label = "lat",
                               fontsize = size
                               , rot = 340)
            # # , xlab='blah1'
            , xlab = list(label = "lon",
                          fontsize = size
                          , rot = 45)
            , zlab = list(label = "f(s)",
                          fontsize = size
                          , rot = 95)
            , scales = list(arrows = F, cex = 1.3)
            , shade = T
            , screen = list(z = screen[1], x = screen[2])
            # # , colorkey = F
            ,  aspect = c(1, 0.8)
            # # , drape = F
             , col.regions = c("red", "blue")
            # , plot = F
            , par.strip.text = list(cex = 2, lines = 1.5))
  
  
  return(list(pred.coord = predModel$pred.coord,
              z.pre = z.pre,
              lon = lon, lat = lat,
              z.true = z.true,
              alpha = alpha,
              alpha.true = z,
              p = p))
}

semiQLM <- function(y, m, coords, covModel,
                    nnIndx, nnIndxLU, sigmaSq, phi, 
                    nu, nThreads){
  cov.model.names <- c("exponential","spherical","matern","gaussian")##order much match util.cpp spCor
  cov.model.indx <- which(covModel == cov.model.names)-1
  n <- length(y)
  N <- length(phi)
  storage.mode(y) <- "double" 
  storage.mode(n) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(N) <- "integer" 
  
  storage.mode(coords) <- "double"
  storage.mode(cov.model.indx) <- "integer"
  storage.mode(nnIndx) <- "integer" 
  storage.mode(nnIndxLU) <- "integer"
  
  storage.mode(sigmaSq) <- "double"
  storage.mode(phi) <- "double"
  storage.mode(nu) <- "double"
  
  storage.mode(nThreads) <- "integer"
  
  return(semiQLME(y, n, m, N, coords, cov.model.indx,
                  nnIndx, nnIndxLU, sigmaSq, phi, 
                  nu, nThreads) )
}

local_linear_kernel <- function(y, covZ,
                                covModel, h = 0.1,
                                nu = 0, nuUnifb = 0,
                                nThreads = 10){
  n <- length(y)
  storage.mode(y) <- "double" 
  storage.mode(covZ) <- "double" 
  
  storage.mode(n) <- "integer"
  storage.mode(covModel) <- "integer"
  storage.mode(h) <- "double" 
  storage.mode(nu) <- "double" 
  storage.mode(nThreads) <- "integer"
  storage.mode(nuUnifb) <- "integer"
  locLin = local_kernel_est(y, covZ, n, 
                            covModel, h, nu, 
                            nuUnifb, nThreads)
  
  
  return(list(S0 = t(matrix(locLin$S[1, ], n, n)),
              S1 = t(matrix(locLin$S[2, ], n, n)),
              alpha = locLin$alpha))
}
LocalPredict <- function(y, Z, TestZ,
                         covModel, h, nu, nuUnifb = 0, 
                         nThreads){
  n = length(y)
  nTest = length(TestZ)

  storage.mode(y) <- "double"
  storage.mode(Z) <- "double"
  storage.mode(TestZ) <- "double"

  storage.mode(nuUnifb) <- "integer"
  
  storage.mode(n) <- "integer"
  storage.mode(nTest) <- "integer"
  storage.mode(covModel) <- "integer"
  storage.mode(h) <- "double"
  storage.mode(nu) <- "double"
  storage.mode(nThreads) <- "integer"
  
  return(local_kernel_pred(y, Z, TestZ, n, nTest, 
                           covModel, h, nu, nuUnifb,
                           nThreads))
  
}
semiPredict <- function(model, newZ, newCoords){
  newZ <- as.matrix(newZ)
  newCoords <- as.matrix(newCoords)
  
  localPre <- spLocLinKernel(y = model$data$order$data$semiY,
                             X = model$data$order$data$Z,
                             coord = model$data$order$data$coords,
                             Q = model$data$order$neighbors$Q,
                             fs.Density = model$data$order$data$fs.Density,
                             pred.X = newZ,
                             pred.coord = newCoords,
                             FitModel = F,
                             pred.nGrid = 0,
                             covModel = model$model$mean.covModel,
                             h = model$model$meanH,
                             nu = model$model$var.nu, 
                             nuUnifb = model$model$nuUnifb, 
                             adapt.n = model$model$adapt.n,
                             nThreads = model$model$nThreads)
  return(list(S0 = localPre$predict$S0,
              S1 = localPre$predict$S,
              alpha = localPre$predict$alpha))
}


Predict1 <- function(model, newX = NULL, newZ = NULL, newCoords, newCovZ = NULL){
  if(!is.null(newX)){
    newX <- as.matrix(cbind(1, newX))}
  
  newCoords <- as.matrix(newCoords)
  n.0 <- nrow(newCoords)
  n <- nrow(model$data$order$data$coords)
  
  # if(model$model$center){
  #   for(i in 1:model$data$input$Pz){
  #     Rz <- range(newZ[, i])
  #     newZ[, i] <- (newZ[, i] - Rz[1])/(Rz[2]- Rz[1])
  #   }
  #   Rz <- range(newCovX)
  #   newCovX <- (newCovX - Rz[1])/(Rz[2]- Rz[1])
  # }
  if(!is.null(newCovZ)){
    newCovZ <- as.matrix(newCovZ)
    if(ncol(newCovZ) == 1){
      CovPre <- LocalPredict(y = model$data$order$residual^2, 
                             Z = model$data$order$data$covZ,
                             TestZ = newCovZ,  
                             covModel = model$model$var.covModel,
                             h = model$model$varH,
                             nu = model$model$var.nu,
                             nuUnifb = model$model$nuUnifb, 
                             nThreads = model$model$nThreads)
      
      var.pre <- CovPre$alpha[1,]
    }else{
      CovPre <- spLocLinKernel(y = model$data$order$residual^2,
                               X = model$data$order$data$covZ, 
                               coord = model$data$order$data$coords,
                               fs.Density = model$data$order$data$fs.Density,
                               pred.X = rep(1, length(newCovZ)),
                               pred.coord = newCoords,
                               FitModel = F,
                               pred.nGrid = 0, 
                               covModel = model$model$var.covModel, 
                               h = model$model$varH,
                               nu = model$model$var.nu,
                               nuUnifb = model$model$nuUnifb, 
                               adapt.n = model$model$adapt.n,
                               nThreads = model$model$nThreads)
      
      var.pre <- CovPre$predict$alpha$alpha1
    }
    
    # locpoly(newCovX, Fit$est$variance$sigmaSq, bandwidth = 0.25)
    # var.pre = local_linear_kernel(pre.residual^2, newCovX, model$coords, 
    #                               model$Kernel[2], 
    #                               model$varH, model$GeomVariable, 
    #                               model$n.omp.threads)$alpha[1, ]
    
    if(length(which(var.pre < 0))>= 1){
      var.pre[which(var.pre < 0)] <- min(var.pre[- which(var.pre < 0)])
    }
  }else{
    CovPre <- NULL
    var.pre <- rep(model$est$variance$sigmaSq[1], n.0)
  }
  nn.indx.0 <- RANN::nn2(model$data$order$data$coords, newCoords,
                         k = model$data$order$neighbors$n.neighbors)$nn.idx
  
  
  wPre <- vector()
  for(i in 1:n.0)
  {
    
    Ci <- sqrt(var.pre[i])*Rdist(newCoords[i,], 
                                 model$data$order$data$coords[nn.indx.0[i,],], 
                                 covModel = model$model$mean.covModel, 
                                 phi = model$model$phi, 
                                 nu = model$model$mean.nu,
                                 nuUnifb = model$model$nuUnifb, 
                                 threads = model$model$nThreads)$Corr
    
    siGma <- Rdist(model$data$order$data$coords[nn.indx.0[i,], ], 
                   model$data$order$data$coords[nn.indx.0[i,], ], 
                   covModel = model$model$mean.covModel, 
                   phi = model$model$phi, 
                   nu = model$model$mean.nu,
                   nuUnifb = model$model$nuUnifb, 
                   threads = model$model$nThreads)$Corr
    
    
    # Ci = sqrt(var.pre[i])*exp(-Rdist(newCoords[i,], 
    #                                  model$data$order$data$coords[nn.indx.0[i,],])$Dist/model$model$phi)
    
    # siGma <- exp(-Rdist(model$data$order$data$coords[nn.indx.0[i,], ],
    #                     model$data$order$data$coords[nn.indx.0[i,],])$Dist/model$model$phi)
    
    siGma <- diag(sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]])) %*% 
      siGma %*% diag(sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]]))
    # for(k1 in 1:model$data$order$neighbors$n.neighbors){
    #   siGma[k1, ] <- siGma[k1, ] * sqrt(model$data$order$sigmaSq[nn.indx.0[i, k1]])
    # }
    # for(k2 in 1:model$data$order$neighbors$n.neighbors){
    #   siGma[, k2] <- siGma[, k2] * sqrt(model$data$order$sigmaSq[nn.indx.0[i, k2]]) 
    # }
    
    wPre[i] <- as.vector((Ci * sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]])) %*% 
                           solve(siGma) %*% model$data$order$residual[nn.indx.0[i, ]]) 
  }
  
  if(!is.null(newZ))
  {
    newZ <- as.matrix(newZ)
    semiP <- semiPredict(model, newZ, newCoords)
    semiPred <- semiP$S0 %*% model$data$order$data$semiY
  }else{
    semiPred <- 0;
    semiP <- NULL
  }
  
  FixEffect <- 0
  if(!is.null(newX)){
    FixEffect <- newX %*% rbind(model$est$mean$intercept, model$est$mean$beta)
  } #+ wPre
  pred <-  FixEffect + semiPred + wPre
  # return(list(pred = pred))
  return(list(FinalPred = pred, FixEffect = FixEffect,
              semiPred = semiPred, RandmomPred = wPre, 
              varPred = var.pre,
              semiPreDa = semiP, CovPreDa = CovPre))
}

KrigPredict <- function(model, newCoords, newCovZ, verbose = TRUE){
  coords <- model$data$order$data$coords
  n <- nrow(coords)
  m <- model$data$order$neighbors$n.neighbors
  newCoords <- as.matrix(cbind(newCoords))
  N <- nrow(newCoords)
  nnIndx <- RANN::nn2(coords, newCoords, k = m)$nn.idx - 1
  wSamples <- model$data$order$residual
  covModel <- model$model$mean.covModel
  sigmaSq <- sqrt(model$data$order$sigmaSq)
  newCovZ <- as.matrix(cbind(newCovZ))
  if(ncol(newCovZ) == 1){
    CovPre <- LocalPredict(y = model$data$order$residual^2, 
                           Z = model$data$order$data$covZ,
                           TestZ = newCovZ, 
                           covModel = model$model$var.covModel,
                           h = model$model$varH,
                           nu = model$model$var.nu,
                           nuUnifb = model$model$nuUnifb, 
                           nThreads = model$model$nThreads)
    
    NewSigmaSq <- (CovPre$alpha[1, ])
    
  }else{
    CovPre <- spLocLinKernel(y = model$data$order$residual^2,
                             X = model$data$order$data$covZ, 
                             coord = model$data$order$data$coords,
                             fs.Density = model$data$order$data$fs.Density,
                             pred.X = rep(1, length(newCovZ)),
                             pred.coord = newCoords,
                             FitModel = F,
                             pred.nGrid = 0, 
                             covModel = model$model$var.covModel, 
                             h = model$model$varH,
                             nu = model$model$var.nu,
                             nuUnifb = model$model$nuUnifb, 
                             adapt.n = model$model$adapt.n,
                             nThreads = model$model$nThreads)
    
    NewSigmaSq <- (CovPre$predict$alpha$alpha1)
    
  }
  if(length(which(NewSigmaSq < 0))>= 1){
    NewSigmaSq[which(NewSigmaSq < 0)] <- min(NewSigmaSq[- which(NewSigmaSq < 0)])
  }
  NewSigmaSq <- sqrt(NewSigmaSq)
  
  phi <- model$model$phi
  nu <- model$model$var.nu
  nuUnifb <- model$model$nuUnifb
  nThreads <- model$model$nThreads
  
  storage.mode(coords) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(newCoords) <- "double"
  storage.mode(N) <- "integer"
  storage.mode(nnIndx) <- "integer"
  storage.mode(wSamples) <- "double"
  storage.mode(covModel) <- "integer"
  storage.mode(sigmaSq) <- "double"
  storage.mode(NewSigmaSq) <- "double"
  storage.mode(phi) <- "double"
  storage.mode(nu) <- "double"
  storage.mode(nuUnifb) <- "integer"
  storage.mode(nThreads) <- "integer"
  storage.mode(verbose) <- "integer"
  
  return(list(wPre = krigPred(coords, n, m, newCoords, N, nnIndx,
                     wSamples, covModel, sigmaSq,
                     NewSigmaSq, phi, nu, nuUnifb,
                     nThreads, verbose),
              NewSigmaSq = NewSigmaSq^2))
}

residual <- function(model, newCovZ = NULL, pred.coord = NULL,
                     pred.nGrid = 1e2, Num = 50, colvarName = "VarX",
                     a = 0, b = 1, FUN = identity, width = 10, height = 6,
                     path = getwd(), verbose = T){
  library(plot3D)
  if(is.null(pred.coord)){
    Min = c(min(model$data$order$data$coords[, 1]), min(model$data$order$data$coords[, 2]))
    Max = c(max(model$data$order$data$coords[, 1]), max(model$data$order$data$coords[, 2]))
    # }
    # cat(Min)
    pred.coord <- cbind(seq(Min[1], Max[1],, pred.nGrid), 
                        seq(Min[2], Max[2],, pred.nGrid))
    pred.coord <- expand.grid(x = pred.coord[, 1],
                          y = pred.coord[, 2])
  }
  
  # variance
  if(!is.null(newCovZ)){
    min.varX <-  quantile(newCovZ, 0.01)#max(min(Train$varX), 0.1)
    max.varX <-  quantile(newCovZ, 0.99)#max(Train$varX)
    newCovZ = seq(min.varX, max.varX, , pred.nGrid^2)
    
    # CovPre <- LocalPredict(y = Fit$data$order$residual^2, 
    #                        Z = Fit$data$order$data$covZ,
    #                        TestZ = newCovZ,  
    #                        covModel = Fit$model$var.covModel,
    #                        h = temp$varH,
    #                        nu = Fit$model$var.nu,
    #                        nuUnifb  = Fit$model$nuUnifb,
    #                        nThreads = Fit$model$nThreads)
    
    # var.pre <- CovPre$alpha[1,] 
  }else{
    newCovZ <- FUN(n^2, a, b)
  }
  Res <- KrigPredict(model = model, 
                       newCoords = pred.coord,
                       newCovZ = newCovZ, 
                       verbose = verbose)
  residual = data.table::data.table(lon = pred.coord[, 1],
                                    lat = pred.coord[, 2],
                                    residual = colMeans(Res$wPre$wMean))
  
  z.pre <- data.table::dcast(residual, lon ~ lat, 
                             value.var = "residual", 
                             mean) %>% as.matrix()
  lon = as.numeric(z.pre[, 1])
  z.pre <- z.pre[, -1]
  lat = as.numeric(colnames(z.pre))
  
  
  pdf(file = paste0(path, "/residual.pdf"), width = width, height = height)
  par(mfrow = c(1, 2))
  # variance
  ord <- model$data$order$ord
  residual.fit <- model$data$order$residual[order(ord),, drop = F]
  # pdf(file = paste0(path, "/variance.pdf"), width = 10, height = 6)
  # 
  plot(model$data$order$data$covZ, residual.fit, cex = 0.35, pch = 20,
       ylim = c(range(residual.fit, Fit$est$variance$sigmaSq)),
       ylab = "Residual and Scale (Sigma.sq) function", cex.axis = 0.8,
       xlab = colvarName, cex.lab = 1)
  # points(simDa$train.VarX, simDa$train.Var, cex = 1, pch = 20) 
  points(model$data$order$data$covZ, model$data$order$sigmaSq,
         cex = 0.5, col = "red", pch = 20)
  points(newCovZ, Res$NewSigmaSq, cex = .5, pch = 20) 
  abline(h = 0)
  # Residuals
  image2D(x = lon, y = lat, cex.axis = 0.8,
          z = z.pre,  contour = list(col = "black", 
               labcex = 1, lwd = 3, alpha = .6),
          xlab = "longitude", ylab = "latitude",
          main = "Residuals: w(s)",
          colkey = list(length = 1, line.clab = 1.0,
                        # adj.clab = 1e-3,
                        # dist = -0.1,
                        # addlines = T,
                        at = round(seq(min(z.pre), 
                                       max(z.pre), , Num), 2),
                        # at = round(seq(min(z.pre), max(z.pre), , Num), 1),
                        width = 1,
                        shift =0,
                        cex.axis = 1,
                        cex.clab = 1,
                        side = 4),
          levels = round(seq(min(z.pre), max(z.pre), , Num), 2)
  )
  # dev.off()
  
  dev.off()
  return(list(residual = residual,
              residual.mat = z.pre,
              newCovZ = newCovZ,
              NewSigmaSq = Res$NewSigmaSq,
              residual.err = Res$wPre$wSigma))
}


# residual.predict <- function(model, newCoords, newCovZ){
#     newCovZ <- as.matrix(cbind(newCovZ))
#     if(ncol(newCovZ) == 1){
#       CovPre <- LocalPredict(y = model$data$order$residual^2,
#                              Z = model$data$order$data$covZ,
#                              TestZ = newCovZ,
#                              covModel = model$model$var.covModel,
#                              h = model$model$varH,
#                              nu = model$model$var.nu,
#                              nuUnifb = model$model$nuUnifb,
#                              nThreads = model$model$nThreads)
# 
#       var.pre <- CovPre$alpha[1,]
#     }else{
#       CovPre <- spLocLinKernel(y = model$data$order$residual^2,
#                                X = model$data$order$data$covZ,
#                                coord = model$data$order$data$coords,
#                                fs.Density = model$data$order$data$fs.Density,
#                                pred.X = rep(1, length(newCovZ)),
#                                pred.coord = newCoords,
#                                FitModel = F,
#                                pred.nGrid = 0,
#                                covModel = model$model$var.covModel,
#                                h = model$model$varH,
#                                nu = model$model$var.nu,
#                                nuUnifb = model$model$nuUnifb,
#                                adapt.n = model$model$adapt.n,
#                                nThreads = model$model$nThreads)
# 
#       var.pre <- CovPre$predict$alpha$alpha1
#     }
# 
#     # locpoly(newCovX, Fit$est$variance$sigmaSq, bandwidth = 0.25)
#     # var.pre = local_linear_kernel(pre.residual^2, newCovX, model$coords,
#     #                               model$Kernel[2],
#     #                               model$varH, model$GeomVariable,
#     #                               model$n.omp.threads)$alpha[1, ]
# 
#     if(length(which(var.pre < 0))>= 1){
#       var.pre[which(var.pre < 0)] <- min(var.pre[- which(var.pre < 0)])
#     }
#   nn.indx.0 <- RANN::nn2(model$data$order$data$coords, newCoords,
#                          k = model$data$order$neighbors$n.neighbors)$nn.idx
# 
# 
#   wPre <- vector()
#   for(i in 1:nrow(newCovZ))
#   {
# 
#     Ci <- sqrt(var.pre[i])*Rdist(newCoords[i,],
#                                  model$data$order$data$coords[nn.indx.0[i,],],
#                                  covModel = model$model$mean.covModel,
#                                  phi = model$model$phi,
#                                  nu = model$model$mean.nu,
#                                  nuUnifb = model$model$nuUnifb,
#                                  threads = model$model$nThreads)$Corr
# 
#     siGma <- Rdist(model$data$order$data$coords[nn.indx.0[i,], ],
#                    model$data$order$data$coords[nn.indx.0[i,], ],
#                    covModel = model$model$mean.covModel,
#                    phi = model$model$phi,
#                    nu = model$model$mean.nu,
#                    nuUnifb = model$model$nuUnifb,
#                    threads = model$model$nThreads)$Corr
# 
# 
#     # Ci = sqrt(var.pre[i])*exp(-Rdist(newCoords[i,],
#     #                                  model$data$order$data$coords[nn.indx.0[i,],])$Dist/model$model$phi)
# 
#     # siGma <- exp(-Rdist(model$data$order$data$coords[nn.indx.0[i,], ],
#     #                     model$data$order$data$coords[nn.indx.0[i,],])$Dist/model$model$phi)
# 
#     siGma <- diag(sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]])) %*%
#       siGma %*% diag(sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]]))
#     # for(k1 in 1:model$data$order$neighbors$n.neighbors){
#     #   siGma[k1, ] <- siGma[k1, ] * sqrt(model$data$order$sigmaSq[nn.indx.0[i, k1]])
#     # }
#     # for(k2 in 1:model$data$order$neighbors$n.neighbors){
#     #   siGma[, k2] <- siGma[, k2] * sqrt(model$data$order$sigmaSq[nn.indx.0[i, k2]])
#     # }
# 
#     wPre[i] <- as.vector((Ci * sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]])) %*%
#                            solve(siGma) %*% model$data$order$residual[nn.indx.0[i, ]])
# 
# }
#  return(list(wPre = wPre, var.pre = var.pre))
# }




Predict <- function(model, newX = NULL, newZ = NULL, newCoords, newCovZ = NULL){
  if(!is.null(newX)){newX <- as.matrix(cbind(newX))}
  newCoords <- as.matrix(newCoords)
  n <- nrow(model$data$order$data$coords)
  # if(model$model$center){
  #   for(i in 1:model$data$input$Pz){
  #     Rz <- range(newZ[, i])
  #     newZ[, i] <- (newZ[, i] - Rz[1])/(Rz[2]- Rz[1])
  #   }
  #   Rz <- range(newCovX)
  #   newCovX <- (newCovX - Rz[1])/(Rz[2]- Rz[1])
  # }
  if(!is.null(newCovZ)){
    newCovZ <- as.matrix(cbind(newCovZ))
    CovPre <- KrigPredict(model = model, newCoords = newCoords, 
                newCovZ = newCovZ, verbose = TRUE)
    wPre <- colSums(CovPre$wPre$wMean)
    var.pre <- CovPre$NewSigmaSq
  }else{
    CovPre <- NULL
    var.pre <- wPre <- rep(model$est$variance$sigmaSq[1], n.0)
  }
  
  if(!is.null(newZ))
  {
    newZ <- as.matrix(newZ)
    semiP <- semiPredict(model, newZ, newCoords)
    semiPred <- semiP$S0 %*% model$data$order$data$semiY
  }else{
    semiPred <- 0;
    semiP <- NULL
  }

  FixEffect <- 0
  if(!is.null(newX)){
    if(ncol(newX)>1){
    FixEffect <- newX %*% rbind(model$est$mean$beta)
  }else{
    FixEffect <- newX * rbind(model$est$mean$beta)
  }
  }
  pred <-  FixEffect + semiPred + wPre + model$est$mean$intercept
  # return(list(pred = pred))
  return(list(FinalPred = pred, FixEffect = FixEffect,
              semiPred = semiPred, RandmomPred = wPre,
              varPred = var.pre,
              semiPreDa = semiP, CovPreDa = CovPre))
}




bivariate_local_kernel <- function(y, covZ, lon, lat,
                                   Kernel, h, 
                                   nThreads){
  n <- length(y)
  storage.mode(y) <- "double" 
  storage.mode(covZ) <- "double"
  storage.mode(lon) <- "double"
  storage.mode(lat) <- "double"
  storage.mode(Kernel) <- "integer"
  storage.mode(h) <- "double"
  storage.mode(nThreads) <- "integer"
  storage.mode(n) <- "integer"
  
  return(bivariate_local_kernel_est(y, covZ,lon, lat, n,
                                    Kernel, h, 
                                    nThreads))
}



SemiAlphaProfile <- function(y, Z, n, coords, B, varF, Q, nnIndx, nnIndxLU,
                  Kernel, meanH, GeomVariable, nThreads){
  n = length(y)
  storage.mode(y) <- "double"
  storage.mode(Z) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(B) <- "double"
  storage.mode(varF) <- "double"
  storage.mode(Q) <- "double"
  
  storage.mode(nnIndx) <- "integer" 
  storage.mode(nnIndxLU) <- "integer"
  
  storage.mode(Kernel) <- "integer"
  storage.mode(meanH) <- "double"
  
  storage.mode(GeomVariable) <- "integer"
  storage.mode(nThreads) <- "integer"
  locLin <- SemiAlpha(y, Z, n, coords, B, varF, Q, nnIndx, nnIndxLU,
                    Kernel, meanH,
                    GeomVariable,
                    nThreads)
 return(list(S0 = t(matrix(locLin$S[1, ], n, n)),
             S1 = t(matrix(locLin$S[2, ], n, n)),
             alpha = locLin$alpha))
}




SemiProf <- function(y, Z, coords, Q, 
                      nnIndx,nnIndxLU,
                      Kernel, h, 
                      GeomVariable,
                      nThreads){
  
 
  n <- length(y)
  storage.mode(y) <- "double" 
  storage.mode(Z) <- "double"
  storage.mode(n) <- "integer"
  
  storage.mode(coords) <- "double" 
  storage.mode(Q) <- "double" 
  
  storage.mode(nnIndx) <- "integer"
  storage.mode(nnIndxLU) <- "integer"
  storage.mode(Kernel) <- "integer"
  storage.mode(h) <- "double" 
  storage.mode(GeomVariable) <- "integer"
  storage.mode(nThreads) <- "integer"
  
  locPro = SemiAlphaPro(y, Z, n,
               coords, 
               Q,
               nnIndx, 
               nnIndxLU,
               Kernel,
               h,	
               GeomVariable,
               nThreads)

  return( list(S0 = t(matrix(locPro$S[1, ], n, n)),
               S1 = t(matrix(locPro$S[2, ], n, n)),
               alpha = locPro$alpha))
}
SemiPredPro <- function(y, Z, TestZ,
                        coords, 
                        TestCoords,
                        Q,
                        nnIndx, 
                        nnIndxLU,
                        Kernel,
                        h,	
                        GeomVariable,				 
                        nThreads){
  n <- length(y)
  nTest = length(TestZ)
  storage.mode(y) <- "double" 
  storage.mode(Z) <- "double"
  storage.mode(TestZ) <- "double"
  
  storage.mode(n) <- "integer"
  storage.mode(nTest) <- "integer"
  
  storage.mode(coords) <- "double" 
  storage.mode(TestCoords) <- "double" 
  
  storage.mode(Q) <- "double" 
  
  storage.mode(nnIndx) <- "integer"
  storage.mode(nnIndxLU) <- "integer"
  storage.mode(Kernel) <- "integer"
  storage.mode(h) <- "double" 
  storage.mode(GeomVariable) <- "integer"
  storage.mode(nThreads) <- "integer"
  

  return(semiProPred(y, Z, TestZ,
                      n, nTest, 
                      coords, 
                      TestCoords,
                      Q,
                      nnIndx, 
                      nnIndxLU,
                      Kernel,
                      h,	
                      GeomVariable,				 
                      nThreads))
}

profilePredict <- function(y, Z, TestZ, coords,
                           TestCoords, B, D, Q,
                           nnIndx, nnIndxLU, 
                           Kernel, h, 
                           GeomVariable,
                           nThreads)
{
  n = length(y)
  nTest = length(TestZ)
  
  storage.mode(y) <- "double"
  storage.mode(Z) <- "double"
  storage.mode(TestZ) <- "double"
  storage.mode(coords) <- "double"
  storage.mode(TestCoords) <- "double"
  
  storage.mode(B) <- "double"
  storage.mode(D) <- "double"
  storage.mode(Q) <- "double"
  
  storage.mode(nnIndx) <- "integer"
  storage.mode(nnIndxLU) <- "integer"
  
  storage.mode(n) <- "integer"
  storage.mode(nTest) <- "integer"
  storage.mode(Kernel) <- "integer"
  
  storage.mode(h) <- "double"
  storage.mode(GeomVariable) <- "integer"
  storage.mode(nThreads) <- "integer"
  
  
  return(semiPred(y, Z, TestZ, n, nTest,
                          coords, TestCoords, 
                          B, D, Q, nnIndx, nnIndxLU,
                          Kernel, h, GeomVariable, nThreads))
}










# semiPredict <- function(model, newZ, newCoords){
#   n = nrow(newZ)
#   N = nrow(model$data$order$data$coords)
#   p = ncol(newZ)
#   newZ <- as.matrix(newZ)
#   newCoords <- as.matrix(newCoords)
#   # if(p == model$data$input$Pz){}
#   semiPreDa <- matrix(0, nrow = n, ncol = model$data$input$Pz)
#   # L = chol(model$data$order$neighbors$Q)
#   if(model$model$profile){
#     for(i in 1:model$data$input$Pz) {
#       localPre <- SemiPredPro(model$data$order$data$semiY[, i],# -
#                               # rowSums(as.matrix(semiPreDa[, -i])),
#                               model$data$order$data$Z[, i],
#                               newZ[, i],
#                               model$data$order$data$coords,
#                               newCoords, 
#                               model$data$order$neighbors$Q,
#                               model$data$order$neighbors$nnIndx,
#                               model$data$order$neighbors$nnIndxLU,
#                               model$model$Kernel[1],
#                               model$model$SecMeanH[i], 
#                               model$model$GeomVariable,
#                               model$model$nThreads)
#       
#       # y.update <- L%*%model$data$order$data$semiY - 
#       #            (L - diag(N))%*% localPre$alpha[1,]
#       # localPre <- LocalPredict(y.update -
#       #                            rowSums(as.matrix(semiPreDa[, -i])), 
#       #                          model$data$order$data$Z[, i], 
#       #                          newZ[, i], 
#       #                          model$data$order$data$coords,
#       #                          newCoords, model$model$Kernel[1],
#       #                          model$model$meanH[i], 
#       #                          model$model$GeomVariable,
#       #                          model$model$nThreads)
#       
#       
#       # localPre <- profilePredict(model$data$order$data$semiY,
#       #                                  model$data$order$data$Z[, i],
#       #                                  newZ[, i], 
#       #                                  model$data$order$data$coords,
#       #                                  newCoords, 
#       #                                  model$data$order$neighbors$B, 
#       #                                  model$data$order$neighbors$varF, 
#       #                                  model$data$order$neighbors$Q,
#       #                                  model$data$order$neighbors$nnIndx,
#       #                                  model$data$order$neighbors$nnIndxLU,
#       #                                  model$model$Kernel[1],
#       #                                  model$model$meanH[i], 
#       #                                  model$model$GeomVariable,
#       #                                  model$model$nThreads)
#       semiPreDa[, i] <- localPre$alpha[1,]
#     }
#   }else{
#     for(i in 1:model$data$input$Pz) {
#       for(j in 1:1){
#       localPre <- LocalPredict(model$data$order$data$semiY[, i], #-
#                                       #rowSums(as.matrix(semiPreDa[, -i])), 
#                                      model$data$order$data$Z[, i], 
#                                      newZ[, i], 
#                                      model$data$order$data$coords,
#                                      newCoords, model$model$Kernel[1],
#                                      model$model$SecMeanH[i], 
#                                      model$model$GeomVariable,
#                                      model$model$nThreads)
#       # y.update <- L%*%model$data$order$data$semiY - 
#       #   (L - diag(N))%*% localPre$alpha[1,]
#       # localPre <- LocalPredict(y.update -
#       #                            rowSums(as.matrix(semiPreDa[, -i])), 
#       #                          model$data$order$data$Z[, i], 
#       #                          newZ[, i], 
#       #                          model$data$order$data$coords,
#       #                          newCoords, model$model$Kernel[1],
#       #                          model$model$meanH[i], 
#       #                          model$model$GeomVariable,
#       #                          model$model$nThreads)
#       semiPreDa[, i]<- localPre$alpha[1,]
#       }
#     }
#   }
# 
#   return(list(S0 = t(matrix(localPre$S[1, ], N, n)),
#               S1 = t(matrix(localPre$S[2, ], N, n)),
#               alpha = semiPreDa))
# }







# Predict <- function(model, newX, newZ, newCoords, newCovZ){
#   newX <- as.matrix(cbind(1, newX))
#   newZ <- as.matrix(newZ)
#   newCoords <- as.matrix(newCoords)
#   newCovX <- as.matrix(newCovX)
#   
#   # if(model$model$center){
#   #   for(i in 1:model$data$input$Pz){
#   #     Rz <- range(newZ[, i])
#   #     newZ[, i] <- (newZ[, i] - Rz[1])/(Rz[2]- Rz[1])
#   #   }
#   #   Rz <- range(newCovX)
#   #   newCovX <- (newCovX - Rz[1])/(Rz[2]- Rz[1])
#   # }
#   
#   
#   CovPre <- LocalPredict(model$data$order$residual^2, 
#                           model$data$order$data$covZ,
#                           newCovZ, model$data$order$data$coords,
#                           newCoords, model$model$Kernel[2],
#                           model$model$varH,
#                           model$model$nThreads)
#   
#   var.pre <- CovPre$alpha[1,]
#   # locpoly(newCovX, Fit$est$variance$sigmaSq, bandwidth = 0.25)
#   # var.pre = local_linear_kernel(pre.residual^2, newCovX, model$coords, 
#   #                               model$Kernel[2], 
#   #                               model$varH, model$GeomVariable, 
#   #                               model$n.omp.threads)$alpha[1, ]
#   
#   if(length(which(var.pre < 0))>= 1){
#     var.pre[which(var.pre < 0)] <- min(var.pre[- which(var.pre < 0)])
#   }
#   nn.indx.0 <- RANN::nn2(model$data$order$data$coords, newCoords,
#                          k = model$data$order$neighbors$n.neighbors)$nn.idx
#   n.0 <- nrow(newCoords)
#   n <- nrow(model$data$order$data$coords)
#   wPre <- vector()
#   for(i in 1:n.0)
#   {
#     Ci = sqrt(var.pre[i])*exp(-Rdist(newCoords[i,], 
#          model$data$order$data$coords[nn.indx.0[i,],])$Dist/model$model$phi)
#     
#     siGma <- exp(-Rdist(model$data$order$data$coords[nn.indx.0[i,], ],
#              model$data$order$data$coords[nn.indx.0[i,],])$Dist/model$model$phi)
#     
#     siGma <- diag(sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]])) %*% 
#             siGma %*% diag(sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]]))
#     # for(k1 in 1:model$data$order$neighbors$n.neighbors){
#     #   siGma[k1, ] <- siGma[k1, ] * sqrt(model$data$order$sigmaSq[nn.indx.0[i, k1]])
#     # }
#     # for(k2 in 1:model$data$order$neighbors$n.neighbors){
#     #   siGma[, k2] <- siGma[, k2] * sqrt(model$data$order$sigmaSq[nn.indx.0[i, k2]]) 
#     # }
#     
#     wPre[i] <- as.vector((Ci * sqrt(model$data$order$sigmaSq[nn.indx.0[i, ]])) %*% 
#                  solve(siGma) %*% model$data$order$residual[nn.indx.0[i, ]]) 
#   }
#   
#   #semiP <- semiPredict(model, newZ, newCoords)
#   semiPred <- rowSums(semiP$alpha)
#   FixEffect <- 0
#   if(!is.null(newX)){
#     FixEffect <- newX %*% rbind(model$est$mean$intercept, model$est$mean$beta)
#   } #+ wPre
#   pred <-  FixEffect + semiPred
#   # return(list(pred = pred))
#   return(list(FinalPred = pred, FixEffect = FixEffect,
#               semiPred = semiPred, RandmomPred = wPre, 
#               varPred = var.pre,
#               semiPreDa = semiP, CovPreDa = CovPre))
# }

semiPlot <- function(data, size){
  p <- ggplot(data, aes(x = x, y = gfun, col = group)) +
    # geom_point(size = 0.5) +
    geom_line(size = 1) +  theme_bw() +
    theme(axis.text = element_text(size = size[2], colour = "black")
          , axis.text.x  = element_text(angle = 0)
          , axis.title.x = element_text(size = size[1], colour = "black")
          , axis.title.y = element_text(size = size[1], colour = "black")
          , legend.title = element_text(size = size[1], colour = "black")
          , legend.text = element_text(size = size[2], colour = "black")
          # , legend.title = element_blank()
          # , legend.background = element_rect(fill="transparent")
          , legend.key.width = unit(5,"line")
          # , panel.grid.major = element_blank()
          # , panel.grid.minor = element_blank()
          , legend.position = "top"#c(0.6, 0.1)
          , legend.margin = margin(t = 1, unit='cm')
          , axis.ticks.length.y.right = unit(-0, "cm")
          , strip.text =  element_text(size = size[2], colour = "black")
          # , axis.text.y.right  = element_text(vjust = -2,
          #                                     hjust = -300,
          #                              margin = margin(l = 15, r = 2))
    )
  return(p)
}

Fitplot <- function(Fit, simDa, path = './figure/semiTemp.pdf'){
  ord <- Fit$data$order$ord
  residual <- Fit$data$order$residual[order(ord),, drop = F]
  
  pdf(file = path, width = 10, height = 5)
  par(mfrow = c(ncol(simDa$train.Z), 2))
  plot(simDa$train.VarX, residual, cex = 0.5, pch = 20,
       ylim = c(range(residual,  simDa$train.Var,
                      Fit$est$variance$sigmaSq)))
  points(simDa$train.VarX, simDa$train.Var, cex = 0.2) 
  points(simDa$train.VarX, Fit$est$variance$sigmaSq, cex = 0.2, col = "red")
  abline(h =0)
  # plot(simDa$train.VarX, residual, cex = 0.5, pch = 20)
  # points(simDa$train.VarX, simDa$train.Var, cex = 0.2) 
  # abline(h =0)
  
  col = c("black", "red", "green" ,"blue")
   pch = c(19, 19, 19, 19)
  # pch = c( ".",  ".",  ".", ".")
  for(nc in 1:ncol(simDa$train.Z))
  {
    y0 = simDa$train.h[, nc]
    y1 = Fit$est$mean$one.fx[, nc] + Fit$est$mean$one.fx.mean[nc]
    y2 = Fit$est$mean$fx[, nc] + Fit$est$mean$fx.mean[nc]
    y.w = simDa$train.h[, nc] + simDa$train.w
    if(nc == ncol(simDa$train.Z)){
      plot(simDa$train.Z[, nc], y0, 
           lwd = 2, cex = 0.2, 
           col = col[1], pch = pch[1],
           ylim = c(range(y0, y1 + Fit$est$mean$intercept,
                          y2 + Fit$est$mean$intercept,
                          y.w)))  #y1 + Fit$est$mean$intercept
    }else{
      plot(simDa$train.Z[, nc], y0, 
           lwd = 2, cex = 0.2, 
           col = col[1], pch = pch[1],
           ylim = c(range(y0, y1, y2, y.w))) 
    }
    if(nc == ncol(simDa$train.Z)){
      points(simDa$train.Z[, nc], y1 + Fit$est$mean$intercept
             ,lwd = 2, cex = 0.1, col = col[2], pch = pch[2]) 
      points(simDa$train.Z[, nc], y2 + Fit$est$mean$intercept
             ,lwd = 2, cex = 0.1, col = col[3], pch = pch[3])
    }else{
      points(simDa$train.Z[, nc], y1 
             #+ Fit$est$mean$intercept 
             , lwd = 2, cex = 0.1, col = col[2], pch = pch[2])
      points(simDa$train.Z[, nc], y2 
             #+ Fit$est$mean$intercept 
             , lwd = 2, cex = 0.1, col = col[3], pch = pch[3]) 
    }
    points(simDa$train.Z[, nc], y.w 
           #+ Fit$est$mean$intercept 
           , lwd = 1, cex = 0.1, 
           col = col[4], pch = pch[4]) 
   
    
    legend(quantile(simDa$train.Z[, nc], 0.65), 
           quantile(y.w, 0.95), 
           c("True", "the first", "the second", "+ noise"),
           col = col,
          # text.col = "green4", 
          # lty = c(2, -1, 1),
          pch = pch#,
          # merge = TRUE 
          # bg = 'gray90'
          )
    # da <- data.frame(x = simDa$train.bCov.xy[, 1],
    #                  y = simDa$train.bCov.xy[, 2],
    #                  bfun = simDa$train.bfun)
    # plot(da$x, da$y, cex = da$bfun)
    # da = reshape2::dcast(da, x ~ y, value.var = "bfun")
    # # da <- as.data.frame(da)
    # x <- da$x
    # y <- as.numeric(colnames(da)[-1])
    # bfun = as.matrix(da[, -1]) 
    # w[is.na(w)] = 0
    # contour(bfun, method = "edge", vfont = c("sans serif", "plain"))
    # image(x, y, bfun)
    # image.plot(x, y, bfun[, ncol(bfun):1])
    
  }
  dev.off()
}

spDensity <- function(xrounded, roundvalue = 1, burnin = 2, samples = 5, adaptive = FALSE, 
          gridsize = 200) 
{
  library(Kernelheaping)
  # gridx = seq(min(xrounded[, 1]) - 0.5 * roundvalue, max(xrounded[, 
  #                                                                 1]) + 0.5 * roundvalue, length = gridsize)
  # gridy = seq(min(xrounded[, 2]) - 0.5 * roundvalue, max(xrounded[, 
  #                                                                 2]) + 0.5 * roundvalue, length = gridsize)
  gridx = xrounded[, 1]
  gridy =  xrounded[, 2]
  
  Mestimates <- ks::kde(x = xrounded, 
                        H = diag(c(roundvalue,
                              roundvalue))^2,
                        gridsize = c(length(gridx), length(gridy)),
                        # eval.points = as.matrix(xrounded)
                        xmin = c(min(gridx), min(gridy)), xmax = c(max(gridx),
                                                                   max(gridy))
                        )
  resultDensity = array(dim = c(burnin + samples, length(gridx), 
                                length(gridy)))
  resultX = array(dim = c(samples + burnin, nrow(xrounded), 
                          2))
  rvalues = unique(xrounded)
  selectionGrid <- lapply(1:nrow(rvalues), function(k) {
    selectionX = which(Mestimates$eval.points[[1]] >= rvalues[k, 
                                                              1] - roundvalue * 0.5 & Mestimates$eval.points[[1]] < 
                         rvalues[k, 1] + roundvalue * 0.5)
    selectionY = which(Mestimates$eval.points[[2]] >= rvalues[k, 
                                                              2] - roundvalue * 0.5 & Mestimates$eval.points[[2]] < 
                         rvalues[k, 2] + roundvalue * 0.5)
    list(selectionX, selectionY)
  })
  delaigle = aggregate(list(length = rep(1, nrow(xrounded))), 
                       data.frame(xrounded), length)
  delaigle[, 3] = delaigle[, 3]/sum(delaigle[, 3])/roundvalue^2
  delaigleest = matrix(0, nrow = length(gridx), ncol = length(gridy))
  for (i in 1:nrow(delaigle)) {
    x = which(delaigle[i, 1] == rvalues[, 1] & delaigle[i, 
                                                        2] == rvalues[, 2])
    delaigleest[selectionGrid[[x]][[1]], selectionGrid[[x]][[2]]] = delaigle[i, 
                                                                             3]
  }
  for (j in 1:(burnin + samples)) {
    new = c()
    for (i in 1:nrow(rvalues)) {
      probs = as.vector(Mestimates$estimate[selectionGrid[[i]][[1]], 
                                            selectionGrid[[i]][[2]]])
      points = cbind(rep(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]], 
                         times = length(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]])), 
                     rep(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]], 
                         each = length(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]])))
      npoints = length(which(xrounded[, 1] == rvalues[i, 
                                                      1] & xrounded[, 2] == rvalues[i, 2]))
      new = rbind(new, points[sample(1:nrow(points), size = npoints, 
                                     replace = T, prob = probs), ])
    }
    if (adaptive == FALSE) {
      H <- ks::Hpi(x = new, binned = TRUE) * 2
    }
    if (adaptive == TRUE) {
      H <- ks::Hpi(x = new, binned = TRUE)
      H <- sqrt(sqrt(H[1, 1] * H[2, 2]))
    }
    if (adaptive == FALSE) {
      Mestimates <- ks::kde(x = new, H = H, gridsize = c(length(gridx), 
                                                         length(gridy)), bgridsize = c(length(gridx), 
                                                                                       length(gridy)), xmin = c(min(gridx), min(gridy)), 
                            xmax = c(max(gridx), max(gridy)), binned = TRUE)
    }
    if (adaptive == TRUE) {
      counts <- plyr::count(new)
      MestimatesAd <- sparr::bivariate.density(data = counts[, 
                                                             c(1:2)], pilotH = H, res = length(gridx), xrange = range(gridx), 
                                               yrange = range(gridy), adaptive = TRUE, comment = FALSE, 
                                               counts = counts[, 3])
      Mestimates$estimate = MestimatesAd$Zm
    }
    Mestimates$estimate[is.na(Mestimates$estimate)] = 1e-96
    resultDensity[j, , ] = Mestimates$estimate
    resultX[j, , ] = new
    print(paste("Iteration:", j, "of", burnin + 
                  samples))
  }
  Mestimates$estimate = apply(resultDensity[-c(1:burnin), , 
  ], c(2, 3), mean)
  est <- list(Mestimates = Mestimates, resultDensity = resultDensity, 
              resultX = resultX, xrounded = xrounded, gridx = gridx, 
              gridy = gridy, roundvalue = roundvalue, burnin = burnin, 
              samples = samples, adaptive = adaptive, delaigle = delaigleest)
  class(est) <- "bivrounding"
  return(est)
}

# bx <- Bessel::besselJ(0.2, nu=0.5)
Cfunc <- function(theta, d, m  = 5){
  lower = max(quantile(d, c(0.01)), 1e-2)
  upper = quantile(d, c(0.99))
  phi <- seq(lower, upper, , m )
  C <- 0
  theta <- exp(theta)
  for(i in 1:m){
    C <- C + theta[i]*Matern(d, range = phi[i]) 
  }
  return(C)
}

betaCov <- function(S, Q, X, R){
  SI <- (diag(nrow(X)) - S) %*% X
  S_D <- solve(t(SI) %*% Q %*% SI)
  V <- t(SI) %*% Q %*% (R) %*% t(Q) %*% SI 
  betaC <- S_D %*% V %*% S_D
  return(sum(log(diag(chol(betaC)))))
}





