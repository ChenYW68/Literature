# # 1 loading packages ---------------------------------------
packages <- c("RandomFields", "fields", "Rcpp"
              , "data.table","plyr", "Matrix"
              , "tidyr", "dplyr", "MASS", "Hmisc"
              , "parallel", "fda"
              # , "spdep"
              # ,"sqldf","gpuR"
              # ,"reticulate","spBayes", "latex2exp"
              # ,"lubridate"
              # , "inlabru","INLA",
              , "progress", "RODBC"
              # , "ggplot2", "cowplot"
              # , "invgamma","DEoptim" , "MBA"
              , "SpecsVerification"
              , "scoringRules", "verification"
)  # , "geostatsp"
# ,'MASS'library(plyr); library(dplyr)
# 2  library
for(i in 1:length(packages))
{
  if(!lapply(packages[i], require,
             character.only = TRUE)[[1]])
  {
    install.packages(packages[i])
    # library(packages[i])
    lapply(packages[i], require,
           character.only = TRUE)
  }else{lapply(packages[i], require,
               character.only = TRUE)}
}
# x=lapply(packages, require, character.only = TRUE)
# rm(list=ls())
rm(i, packages)


