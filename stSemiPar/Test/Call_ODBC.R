
# round(c(mean(Result1[, 4]) - 1,
# mean(Result2[, 4]) - 1, 
# mean(Result3[, 4]) - 1, 
# mean(Result4[, 4]) - 1, 
# sd((Result1[, 4] - 0)^1),
# sd((Result2[, 4] - 0)^1),
# sd((Result3[, 4] - 0)^1),
# sd((Result4[, 4] - 0)^1),
# mean((Result1[, 4] - 1)^2), 
# mean((Result2[, 4] - 1)^2), 
# mean((Result3[, 4] - 1)^2),
# mean((Result4[, 4] - 1)^2)
# ), 
# 8)

# round(c(mean(Result1[, 5]) - 5,
# mean(Result2[, 5]) - 5, 
# mean(Result3[, 5]) - 5,
# mean(Result4[, 5]) - 5, 
# sd((Result1[, 5] - 0)^1),
# sd((Result2[, 5] - 0)^1),
# sd((Result3[, 5] - 0)^1),
# sd((Result4[, 5] - 0)^1),
# mean((Result1[, 5] - 5)^2), 
# mean((Result2[, 5] - 5)^2), 
# mean((Result3[, 5] - 5)^2),
# mean((Result4[, 5] - 5)^2)
# ),
# 8)


# head(Result1[, c(3:5, 8, 11:15)])
# head(Result2[, c(3:5, 8, 11:15)])
# head(Result3[, c(3:5, 8, 11:15)])

# colMeans(Result1[,c(3:5, 8, 11:15)])
# colMeans(Result2[,c(3:5, 8, 11:15)])
# colMeans(Result3[,c(3:5, 8, 11:15)])

# theta
# round(c(mean(rowMeans(alpha.est1[1:Nt, ])- theta.1),
#         mean(rowMeans(alpha.est2[1:Nt, ])- theta.1),
#         mean(rowMeans(alpha.est3[1:Nt, ])- theta.1),
#         mean(rowMeans(alpha.est4[1:Nt, ])- theta.1),
#         mean(sqrt(rowVar(alpha.est1[1:Nt, ]))),
#         mean(sqrt(rowVar(alpha.est2[1:Nt, ]))),
#         mean(sqrt(rowVar(alpha.est3[1:Nt, ]))),
#         mean(sqrt(rowVar(alpha.est4[1:Nt, ]))),
#         mean((alpha.est1[1:Nt, ]- theta.1.1)^2),
#         mean((alpha.est2[1:Nt, ]- theta.1.1)^2),
#         mean((alpha.est3[1:Nt, ]- theta.1.1)^2),
#         mean((alpha.est4[1:Nt, ]- theta.1.1)^2)
#         ),
#         4)
# 
# round(c(mean(rowMeans(alpha.est1[21:(2*Nt), ])- theta.2),
#         mean(rowMeans(alpha.est2[21:(2*Nt), ])- theta.2),
#         mean(rowMeans(alpha.est3[21:(2*Nt), ])- theta.2),
#         mean(sqrt(rowVar(alpha.est1[11:(2*Nt), ]))),
#         mean(sqrt(rowVar(alpha.est2[21:(2*Nt), ]))),
#         mean(sqrt(rowVar(alpha.est3[21:(2*Nt), ]))),
# 		    mean((alpha.est1[21:(2*Nt), ]- theta.2.1)^2),
#         mean((alpha.est2[21:(2*Nt), ]- theta.2.1)^2),
#         mean((alpha.est3[21:(2*Nt), ]- theta.2.1)^2)), 4)
rm(list=ls())
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
source("./R/spT.validation.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
source("E:/Literature/semiBase/R/util.R")

n <- 50
Nt <- 20

time <- seq(0, 1,, Nt)
theta.1 <- matrix(rep(f1(time), times = 1), 
                  nrow =  Nt, ncol = 1)
theta.2 <- matrix(rep(f2(time), times = 1), 
                  nrow =  Nt, ncol = 1)

theta.1.1 <- matrix(rep(f1(time), times = 50), 
                    nrow =  Nt, ncol = 50)
theta.2.1 <- matrix(rep(f2(time), times = 50), 
                    nrow =  Nt, ncol = 50)

Phis <- "0.2"
sigma.sq.s <- 0.5


#M <- c("WI", "WDt", "WDst", "WDstR")
M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
method <- M[1]
tab <- paste0(method, "D_", substr(sigma.sq.s, 3, 3),"_", 
              n, "_", Nt, "_", substr(Phis, 3, 3)) 
Result1 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
load(paste0("./data/", tab, "_alpha_est.RData"))
alpha.est1 <- alpha.est

method <- M[2]
tab <- paste0(method, "I_", substr(sigma.sq.s, 3, 3),"_", 
              n, "_", Nt, "_", substr(Phis, 3, 3))  
Result2 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))

load(paste0("./data/", tab, "_alpha_est.RData"))
alpha.est2 <- alpha.est

method <- M[3]
tab <- paste0(method, "I_", substr(sigma.sq.s, 3, 3),"_", 
              n, "_", Nt, "_", substr(Phis, 3, 3))  
Result3 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
load(paste0("./data/", tab, "_alpha_est.RData"))
alpha.est3 <- alpha.est

method <- M[4]
tab <- paste0(method, "DD_", substr(sigma.sq.s, 3, 3),"_", 
              n, "_", Nt, "_", substr(Phis, 3, 3))  
Result4  <- sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
load(paste0("./data/", tab, "_alpha_est.RData"))
alpha.est4 <- alpha.est

method <- M[5]
tab <- paste0(method, "D_", substr(sigma.sq.s, 3, 3),"_", 
              n, "_", Nt, "_", substr(Phis, 3, 3))  
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
b1 <- 1
b2 <- 5
c1 <- c(5)
r <- c(1e3, 1e3, 1e3)
{
  da <- rbind(data.frame(method = "WI", 
                         Bias.1 = r[1]*(mean(Result1[, c1]) - b1), 
                         Sd.1 = r[2]*sd((Result1[, c1])^1), 
                         MSE.1 = r[3]*mean((Result1[, c1] - b1)^2),
                         # Bias.2 = r[1]*(mean(Result1[, c1 + 1]) - b2), 
                         # Sd.2 = r[2]*sd((Result1[, c1 + 1])^1), 
                         # MSE.2 = r[3]*mean((Result1[, c1 + 1] - b2)^2),
                         Bias.3 = mean(rowMeans(alpha.est1[1:Nt, ])- theta.1), 
                         Sd.3 = mean(sqrt(rowVar(alpha.est1[1:Nt, ]))), 
                         MSE.3 = mean((alpha.est1[1:Nt, ]- theta.1.1)^2)),
              
              # data.frame(method = "$\\text{WEC}_t$",
              #            Bias.1 = r[1]*(mean(Result2[, c1]) - b1),
              #            Sd.1 = r[2]*sd((Result2[, c1])^1),
              #            MSE.1 = r[3]*mean((Result2[, c1] - b1)^2),
              #            # Bias.2 = r[1]*(mean(Result2[, c1 + 1]) - b2),
              #            # Sd.2 = r[2]*sd((Result2[, c1 + 1])^1),
              #            # MSE.2 = r[3]*mean((Result2[, c1 + 1] - b2)^2),
              #            Bias.3 = mean(rowMeans(alpha.est2[1:Nt, ])- theta.1),
              #            Sd.3 = mean(sqrt(rowVar(alpha.est2[1:Nt, ]))),
              #            MSE.3 = mean((alpha.est2[1:Nt, ]- theta.1.1)^2)),
              # data.frame(method = "$\\text{WEC}_{tw}$",
              #            Bias.1 = r[1]*(mean(Result3[, c1]) - b1),
              #            Sd.1 = r[2]*sd((Result3[, c1])^1),
              #            MSE.1 = r[3]*mean((Result3[, c1] - b1)^2),
              #            # Bias.2 = r[1]*(mean(Result3[, c1 + 1]) - b2),
              #            # Sd.2 = r[2]*sd((Result3[, c1 + 1])^1),
              #            # MSE.2 = r[3]*mean((Result3[, c1 + 1] - b2)^2),
              #            Bias.3 = mean(rowMeans(alpha.est3[1:Nt, ])- theta.1),
              #            Sd.3 = mean(sqrt(rowVar(alpha.est3[1:Nt, ]))),
              #            MSE.3 = mean((alpha.est3[1:Nt, ]- theta.1.1)^2)),
              
              data.frame(method = "$\\text{WEC}_{st}$", 
                         Bias.1 = r[1]*(mean(Result4[, c1]) - b1), 
                         Sd.1 = r[2]*sd((Result4[, c1])^1), 
                         MSE.1 = r[3]*mean((Result4[, c1] - b1)^2),
                         # Bias.2 = r[1]*(mean(Result4[, c1 + 1]) - b2), 
                         # Sd.2 = r[2]*sd((Result4[, c1 + 1])^1), 
                         # MSE.2 = r[3]*mean((Result4[, c1 + 1] - b2)^2),
                         Bias.3 = mean(rowMeans(alpha.est4[1:Nt, ])- theta.1), 
                         Sd.3 = mean(sqrt(rowVar(alpha.est4[1:Nt, ]))), 
                         MSE.3 = mean((alpha.est4[1:Nt, ]- theta.1.1)^2)),
              data.frame(method = "$\\text{WEC}_{stw}$",
                         Bias.1 = r[1]*(mean(Result5[, c1]) - b1), 
                         Sd.1 = r[2]*sd((Result5[, c1])^1), 
                         MSE.1 = r[3]*mean((Result5[, c1] - b1)^2),
                         # Bias.2 = r[1]*(mean(Result5[, c1 + 1]) - b2), 
                         # Sd.2 = r[2]*sd((Result5[, c1 + 1])^1), 
                         # MSE.2 = r[3]*mean((Result5[, c1 + 1] - b2)^2),
                         Bias.3 = mean(rowMeans(alpha.est5[1:Nt, ])- theta.1), 
                         Sd.3 = mean(sqrt(rowVar(alpha.est5[1:Nt, ]))), 
                         MSE.3 = mean((alpha.est5[1:Nt, ]- theta.1.1)^2))#,
              # data.frame(method = "$\\text{WEC}_{wls}$", 
              #            Bias.1 = r[1]*(mean(Result6[, c1]) - b1), 
              #            Sd.1 = r[2]*sd((Result6[, c1])^1), 
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



