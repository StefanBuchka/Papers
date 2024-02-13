# Purpose: Build model formula for Bayesian joint models (estimated by brms package)
# 19.12.2023
# Author: Stefan Buchka

formula_maker <- function(lr_list = lr_list)
{
  if(length(lr_list[grep("ar1",lr_list$call),"call"]) > 0)
  {
    ar <- lr_list[grep("ar1",lr_list$call),"call"]
    ar <- strsplit(ar,"ar1")
    ar_temp <- lapply(ar, function(x){x[[1]]})
    ar <- ar[[1]][2]
    gr <- unlist(strsplit(split = '|',ar, fixed = T))[2]
    gr <- unlist(strsplit(split = ')',gr, fixed = T))[1]
    gr <- gsub(" ","",gr)
    
    ar <- paste0("ar(p = 1, gr = ",gr,")")
    
    ar <- paste0(ar_temp,ar)
    lr_list[grep("ar1",lr_list$call),"call"] <- ar
  }
  
  if(any(is.na(lr_list$logLik)))
  {
    lr_list <- lr_list[-which(is.na(lr_list$logLik)),]
  }
  
  return(lr_list)
}