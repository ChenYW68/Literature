library(RODBC)
DSN_01 <- odbcConnect(
                  "DSN_01",
                  uid = "myname",
                  pwd = "mypwd",
                  believeNRows = FALSE,
                  case = "toupper"
                )
source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
source("E:/Literature/semiBase/R/util.R")

n <- 100
Nt <- 20
r1 <- c(1e3, 1e3, 1e3)
r2 <- c(1e3, 1e3, 1e3)
str <- c("_TD_", "CG_5")  #c("I_", "01_05")
time <- seq(0, 1,, Nt)
seq <- c(1, 2, 4, 5)
c1 <- c(3)
theta.1 <- matrix(rep(f1(time), times = 1), 
                  nrow =  Nt, ncol = 1)
theta.2 <- matrix(rep(f2(time), times = 1), 
                  nrow =  Nt, ncol = 1)

iter <- 50
theta.1.1 <- matrix(rep(f1(time), times = iter), 
                    nrow =  Nt, ncol = iter)
theta.2.1 <- matrix(rep(f2(time), times = iter), 
                    nrow =  Nt, ncol = iter)


{
  #M <- c("WI", "WDt", "WDst", "WDstR")
  M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
  method <- M[seq[1]]
  tab <- paste0(method, str[1], n, "_", Nt, "_", str[2])
  # if(as.numeric(Phis) < 0.1){
  #   tab <- paste0(tab, substr(Phis, 4, 4))
  # }
  Result1 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
  load(paste0("./data/", tab, "_alpha_est.RData"))
  alpha.est1 <- alpha.est
  
  method <- M[seq[2]]
  # tab <- paste0(method, "1_", substr(sigma.sq.s, 3, 3),"_",
  #               n, "_", Nt, "_", substr(Phis, 3, 3))
  tab <- paste0(method, str[1], n, "_", Nt, "_", str[2])
  # if(as.numeric(Phis) < 0.1){
  #   tab <- paste0(tab, substr(Phis, 4, 4))
  # }
  Result2  <- sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
  load(paste0("./data/", tab, "_alpha_est.RData"))
  alpha.est2 <- alpha.est
  # 
  # method <- M[3]
  # tab <- paste0(method, "_", substr(sigma.sq.s, 3, 3),"_", 
  #               n, "_", Nt, "_", substr(Phis, 3, 3))  
  # Result3 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
  # load(paste0("./data/", tab, "_alpha_est.RData"))
  # alpha.est3 <- alpha.est
  
  method <- M[seq[3]]
  # tab <- paste0(method, "1_", substr(sigma.sq.s, 3, 3),"_",
  #               n, "_", Nt, "_", substr(Phis, 3, 3))
  tab <- paste0(method, str[1], n, "_", Nt, "_", str[2])
  # if(as.numeric(Phis) < 0.1){
  #   tab <- paste0(tab, substr(Phis, 4, 4))
  # }
  Result4  <- sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
  load(paste0("./data/", tab, "_alpha_est.RData"))
  alpha.est4 <- alpha.est
  
  method <- M[seq[4]]
  tab <- paste0(method, str[1], n, "_", Nt, "_", str[2])
  # if(as.numeric(Phis) < 0.1){
  #   # tab <- paste0(method, "1_", substr(sigma.sq.s, 3, 3),"_",
  #   #               n, "_", Nt, "_", substr(Phis, 3, 3),
  #   #               substr(Phis, 4, 4))
  #   tab <- paste0(tab, substr(Phis, 4, 4))
  # }
  
  Result5  <-sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
  load(paste0("./data/", tab, "_alpha_est.RData"))
  alpha.est5 <-alpha.est
 
  # method <- M[6]
  # tab <- paste0(method, "_", substr(sigma.sq.s, 3, 3),"_", 
  #               n, "_", Nt, "_", substr(Phis, 3, 3))  
  # Result6  <-sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
  # load(paste0("./data/", tab, "_alpha_est.RData"))
  # alpha.est6 <-alpha.est
  # (E(beta.hat) - beta) and sd(beta.hat) 
}

b1 <- 1
# b2 <- 5



{
  da <- rbind(data.frame(method = "WI", 
                         Bias.1 = r1[1]*(mean(Result1[1:iter, c1]) - b1), 
                         Sd.1 = r1[2]*sd((Result1[1:iter, c1])^1), 
                         MSE.1 = #r1[3]*mean((Result1[1:iter, c1] - b1)^2)
                           r1[3]*(mean(Result1[, c1]) - b1)^2 +
                           r1[3]*sd((Result1[, c1])^1)^2
                         ,
                         # Bias.2 = r[1]*(mean(Result1[, c1 + 1]) - b2), 
                         # Sd.2 = r[2]*sd((Result1[, c1 + 1])^1), 
                         # MSE.2 = r[3]*mean((Result1[, c1 + 1] - b2)^2),
                         Bias.3 = r2[1]*mean(rowMeans(alpha.est1[1:Nt, 1:iter])- theta.1), 
                         Sd.3 = r2[2]*mean(sqrt(rowVar(alpha.est1[1:Nt, 1:iter]))), 
                         MSE.3 = #r2[3]*mean((alpha.est1[1:Nt, 1:iter]- theta.1.1)^2)
                         r2[3]*(mean(rowMeans(alpha.est1[1:Nt, 1:iter])- theta.1)^2 +
                                  mean(sqrt(rowVar(alpha.est1[1:Nt, 1:iter])))^2)
                         ),
              data.frame(method = "WII", 
                         Bias.1 = r1[1]*(mean(Result2[1:iter, c1]) - b1), 
                         Sd.1 = r1[2]*sd((Result2[1:iter, c1])^1), 
                         MSE.1 = #r1[3]*mean((Result2[1:iter, c1] - b1)^2)
                         r1[3]*(mean(Result2[, c1]) - b1)^2 +
                         r1[3]*sd((Result2[, c1])^1)^2
                         ,
                         # Bias.2 = r[1]*(mean(Result1[, c1 + 1]) - b2), 
                         # Sd.2 = r[2]*sd((Result1[, c1 + 1])^1), 
                         # MSE.2 = r[3]*mean((Result1[, c1 + 1] - b2)^2),
                         Bias.3 = r2[1]*mean(rowMeans(alpha.est2[1:Nt, 1:iter])- theta.1), 
                         Sd.3 = r2[2]*mean(sqrt(rowVar(alpha.est2[1:Nt, 1:iter]))), 
                         MSE.3 = #r2[3]*mean((alpha.est2[1:Nt, 1:iter]- theta.1.1)^2)
                         r2[3]*(mean(rowMeans(alpha.est2[1:Nt, 1:iter])- theta.1)^2 +
                                  mean(sqrt(rowVar(alpha.est2[1:Nt, 1:iter])))^2)
              ),
              data.frame(method = "$\\text{WTC}$", 
                         Bias.1 = r1[1]*(mean(Result4[1:iter, c1]) - b1), 
                         Sd.1 = r1[2]*sd((Result4[1:iter, c1])^1), 
                         MSE.1 = #r1[3]*mean((Result4[1:iter, c1] - b1)^2)
                           r1[3]*(mean(Result4[, c1]) - b1)^2 +
                           r1[3]*sd((Result4[, c1])^1)^2
                         ,
                         # Bias.2 = r[1]*(mean(Result4[, c1 + 1]) - b2), 
                         # Sd.2 = r[2]*sd((Result4[, c1 + 1])^1), 
                         # MSE.2 = r[3]*mean((Result4[, c1 + 1] - b2)^2),
                         Bias.3 = r2[1]*mean(rowMeans(alpha.est4[1:Nt, 1:iter])- theta.1), 
                         Sd.3 = r2[2]*mean(sqrt(rowVar(alpha.est4[1:Nt, 1:iter]))), 
                         MSE.3 = #r2[3]*mean((alpha.est4[1:Nt, 1:iter]- theta.1.1)^2)
                         r2[3]*(mean(rowMeans(alpha.est4[1:Nt, 1:iter])- theta.1)^2 +
                                  mean(sqrt(rowVar(alpha.est4[1:Nt, 1:iter])))^2)
                         ),
              data.frame(method = "$\\text{WTC}_{w}$",
                         Bias.1 = r1[1]*(mean(Result5[1:iter, c1]) - b1), 
                         Sd.1 = r1[2]*sd((Result5[1:iter, c1])^1), 
                         MSE.1 = #r1[3]*mean((Result5[1:iter, c1] - b1)^2)
                           r1[3]*(mean(Result5[, c1]) - b1)^2 +
                           r1[3]*sd((Result5[, c1])^1)^2
                         ,
                         # Bias.2 = r[1]*(mean(Result5[, c1 + 1]) - b2), 
                         # Sd.2 = r[2]*sd((Result5[, c1 + 1])^1), 
                         # MSE.2 = r[3]*mean((Result5[, c1 + 1] - b2)^2),
                         Bias.3 = r2[1]*mean(rowMeans(alpha.est5[1:Nt, 1:iter])- theta.1), 
                         Sd.3 = r2[2]*mean(sqrt(rowVar(alpha.est5[1:Nt, 1:iter]))), 
                         MSE.3 = #r2[3]*mean((alpha.est5[1:Nt, 1:iter] - theta.1.1)^2)
                           r2[3]*(mean(rowMeans(alpha.est5[1:Nt, 1:iter])- theta.1)^2 +
                           mean(sqrt(rowVar(alpha.est5[1:Nt, 1:iter])))^2)
                         )#,
              # data.frame(method = "$\\text{WEC}_{wls}$", 
              #            Bias.1 = r[1]*(mean(Result6[, c1]) - b1), 
              #            Sd.1 = r[2]*sd((Re说sult6[, c1])^1), 
              #            MSE.1 = r[3]*mean((Result6[, c1] - b1)^2),
              #            Bias.2 = r[1]*(mean(Result6[, c1 + 1]) - b2), 
              #            Sd.2 = r[2]*sd((Result6[, c1 + 1])^1), 
              #            MSE.2 = r[3]*mean((Result6[, c1 + 1] - b2)^2),
              #            Bias.3 = mean(rowMeans(alpha.est6[1:Nt, ])- theta.1), 
              #            Sd.3 = mean(sqrt(rowVar(alpha.est6[1:Nt, ]))), 
              #            MSE.3 = mean((alpha.est6[1:Nt, ]- theta.1.1)^2))
  )
}
da[, -1] <- round(da[, -1], 4)
da
writexl::write_xlsx(da, path = "./data/Resum.xlsx")

# Result2[, 1]
# Result5[, 1]
mean(Result2[, 1])
mean(Result5[, 1])