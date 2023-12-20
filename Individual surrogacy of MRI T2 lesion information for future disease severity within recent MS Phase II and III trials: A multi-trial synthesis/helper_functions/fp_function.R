# Purpose: Function to find fractional polynomials for time
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka

fp_func <- function (Outcome = "edss_score", # Outcome variable
                     Covariate = NULL, # Time variable
                     id = id, # id of patient
                     S = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3), # Possible powers for fractional polynomials
                     Max.M = 2, # Order of polynomials
                     fam = "gaussian", # Family of outcome variable
                     Dataset = data_loop) # Data 
{
  library(lme4)
  library(cAIC4)
  library(foreach)
  library(doParallel)
  library(rlist)
  
  if (Max.M > 2) {
    cat("The maximum allowed Max.M is 2. The requested Max.M is replaced by Max.M=5\n")
    Max.M <- c(3)
  }
  Covariate <- Dataset[,which(names(Dataset) == Covariate)]
  Outcome <- Dataset[, which(names(Dataset) == Outcome)]
  n = length(Outcome)
  Results_m1 <- Results_m2 <- NULL
  
  # Order 1
  for (j in 1:Max.M) {
    r <- length(S)
    n <- j
    aantal <- factorial(n + r - 1)/(factorial(n) * factorial(r - 1))
    if (aantal > 5000) {
      cat("Warning: A total of ", aantal, " models of degree ", 
          i, " can be fitted. Consider using \n", 
          sep = "")
      cat("a more restricted set S or a lower Max.M \n\n", 
          sep = "")
    }
  }
  
  Covariate[Covariate == 0] <- 1e-08
  if (Max.M >= 1) {
    
    # Find fractional polynomials for order one polynomials
    cl <- makeCluster(detectCores()-4)
    registerDoParallel(cl)
    
    Results_m1 <- foreach (i = 1:length(S)) %dopar% 
    {
      library(lme4)
      library(cAIC4)
      
      if (S[i] != 0) {
        term1 <- Covariate^S[i]
      }
      if (S[i] == 0) {
        term1 <- log(Covariate)
      }

      if(fam == "gaussian")
      {
        fit <<- glm(Outcome ~ term1, data = Dataset, family = fam)
        aic <<- AIC(fit)
      }   
      
      if(fam == "poisson")
      {
        fit <<- glm(Outcome ~ offset(log(Covariate)) + term1, data = Dataset, family = fam)
        aic <<- AIC(fit)
      } 
      
      mean_change1 <- mean(abs(term1-Covariate))
      LogLik <- logLik(fit)
      Results_m1 <- cbind(S[i], aic,LogLik,mean_change1)
      colnames(Results_m1) <- c("power1", "AIC","LogLik","mean_change1")
      Results_m1
    }
    
    stopCluster(cl)
  }
  
  Results_m1 <- list.rbind(Results_m1)
  Results_m2 <- NULL
  
  # Find fractional polynomials for order two polynomials
  if (Max.M >= 2) {
    
    cl <- makeCluster(detectCores()-4)
    registerDoParallel(cl)
    Results_m2 <- foreach (i = 1:length(S)) %dopar% 
    {
      library(lme4)
      library(cAIC4)
      library(tidyverse)
      
      modelsj <- vector(mode = "list", length = length(S))
      for (j in 1:length(S)) {
        if (S[i] != 0) {
          term1 <- Covariate^S[i]
        }
        if (S[i] == 0) {
          term1 <- log(Covariate)
        }
        if (S[j] != 0) {
          term2 <- Covariate^S[j]
        }
        if (S[j] == 0) {
          term2 <- log(Covariate)
        }
        if (S[i] == S[j] & S[i] != 0) {
          term2 <- (Covariate^S[j]) * log(Covariate)
        }
        if (S[i] == S[j] & S[i] == 0) {
          term2 <- term1^2
        }
        
           
           if(fam == "gaussian")
           {
             fit <<- glm(Outcome ~ term1 + term2, data = Dataset, family = fam)
             aic <<- AIC(fit)
           }   
           
           if(fam == "poisson")
           {
             fit <<- glm(Outcome ~ offset(log(Covariate)) + term1 + term2, data = Dataset, family = fam)
             aic <<- AIC(fit)
           } 
        
        mean_change1 <- mean(abs(term1-Covariate))
        mean_change2 <- mean(abs(term2-Covariate))
        
        LogLik <- logLik(fit)
        Results_m2 <- rbind(Results_m2, cbind(S[i], S[j],aic,LogLik,mean_change1,mean_change2))
      }
      colnames(Results_m2) <- c("power1", "power2","AIC","LogLik","mean_change1","mean_change2")
      Results_m2
    }
    
    stopCluster(cl)
  }
  
  Results_m2 <- list.rbind(Results_m2)
  
  
  AIC_tab <- bind_rows(as.data.frame(Results_m1),
                       as.data.frame(Results_m2))
  
  AIC_tab <- AIC_tab %>%
    mutate(lowest_AIC = min(AIC),
           AIC_delta = AIC - lowest_AIC) %>%
    arrange(AIC_delta)
  
  # Find the model with best AIC (Akaike Information Criterion)
  if(any(AIC_tab$AIC_delta >= 2))
  {
    # Arrangement strategy: 
    # 1.) Only models, that are considerable better are considered (AIC delta >= 2)
    # 2.) Arrangement of powers according to influence on time variable. The less the grade of transformation of the 
    #     time variable, the better. 
    
    AIC_tab <- AIC_tab %>%
      distinct() %>%
      mutate(mean_change2 = ifelse(is.na(mean_change2),0,mean_change2)) %>%
     filter(AIC_delta <= 4) %>%
      mutate(mean_change = mean_change1 + mean_change2) %>%
      arrange(mean_change, AIC_delta)
    
    fit <- AIC_tab[1,]
  }
  
  else
  {
    fit <- AIC_tab[which(AIC_tab$power1 == 1 & is.na(AIC_tab$power2)),]  
  }

  fit
}
