# Purpose of this script: Creation of function to calculate the LRF (Information Theoretic approach)
# For Details see paper or appendix. 
# 25.11.2023  
# Author: Stefan Buchka

library(rlist)
r_2_function_IT <- function(model0 = NULL,
                            model1 = NULL,
                            ordinal = NULL,
                            TE = TE,
                            data = NULL)
{
  
if(!is.null(model0) & !is.null(model1))
{
  # If the two models two derive LRF are ordinal models
  if(ordinal == T)
  {
    # Extract log-likelihood
    L0 <- -2 * logLik(model0)
    L1 <- -2 * logLik(model1)
    
    # log-likelihood test statistics
    g2 <- -(L1 - L0)
    
    # Adjustments for ordinal endpoints (otherwise, the upper border of LRF is not 1)
    # For details: 
    # Alonso A, Bigirumurame T, Burzykowski T, Buyse M, Molenberghs G, Muchene L, et al. Applied surrogate endpoint evaluation methods with sas and r. CRC Press; 2016.
    l.ncpara <- qchisq(0.025, 1, g2)
    if (pchisq(g2, 1, 0) <= 0.95) {
      l.ncpara <- 0
    }
    u.ncpara <- qchisq(0.975, 1, g2)
    
    N.total <- length(unique(data$id))
    lrf <- 1 - exp(-g2/N.total)
    
    nact = NULL
    sformula = NULL
    data$TE <- as.numeric(unlist(data[,TE]))
    
    
    tryCatch({modRmax <<- clmm2(factor(TE) ~ 1, 
                                random = factor(id),
                                Hess = T,
                                data = data)},
             warning = function(w)
             {
               modRmax <<- clmm2(factor(TE) ~ 1, 
                                 Hess = T,
                                 data = data)
             })
    
    
    # Computation of LRF and its confidence interval
    R2h.max <- lrf/(1 - exp(2 * as.numeric(logLik(modRmax))/N.total))
    paral <- exp(-l.ncpara[!is.na(l.ncpara)]/N.total[!is.na(l.ncpara)])
    parau <- exp(-u.ncpara[!is.na(u.ncpara)]/N.total[!is.na(u.ncpara)])
    lower <- 1 - (sum(paral)/length(paral))
    upper <- 1 - (sum(parau)/length(parau))
    R2h.max.lb <- lower/(1 - exp(2 * as.numeric(logLik(modRmax))/N.total))
    R2h.max.ub <- min(1, upper/(1 - exp(2 * as.numeric(logLik(modRmax))/N.total)))
    R2h <- data.frame(cbind(R2h.max, R2h.max.lb, R2h.max.ub), stringsAsFactors = TRUE)
    colnames(R2h) <- c("R2h", "CI lower limit", "CI upper limit")
    rownames(R2h) <- c(" ")
    
    result <- data.frame(R_sq_ind = R2h$R2h,
                         lb = R2h$`CI lower limit`,
                         ub = R2h$`CI upper limit`)  
    
  }
  
  else
  {
    # Extract log-likelihood
    L0 <- -2 * logLik(model0)
    L1 <- -2 * logLik(model1)
    
    # log-likelihood ratio test statistics
    g2 <- -(L1 - L0)
    
    N.total <- length(unique(data$id))
    Alpha <- 0.05
    
    # Computation of LRF and its confidence interval
    R2h.single.value <- 1 - exp(-g2/N.total)
    k1 <- qchisq(Alpha, 1, g2)
    d1 <- qchisq((1 - Alpha), 1, g2)
    R2h.single.lb <- max(0, 1 - exp(-k1/N.total))
    R2h.single.ub <- min(1, 1 - exp(-d1/N.total))
    
    result <- data.frame(R_sq_ind = R2h.single.value,
                         lb = R2h.single.lb,
                         ub = R2h.single.ub)  
  }
  
}
 
  
else
{
  result <- data.frame(R_sq_ind = NA,
                       lb = NA,
                       ub = NA)      
}
  
  return(result)
}



