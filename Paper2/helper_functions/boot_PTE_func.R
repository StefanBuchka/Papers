

# Purpose of this script: Function to compute bootstrap intervals for PTE. 
# 19.12.2023
# Author: Stefan Buchka


library(MASS)
library(boot)

PET_boot <- function(data, indices,SE,TE)
{
  formula_2 <- paste0(TE," ~ treat")
  model_2 <- glm.nb(formula = as.formula(formula_2),
                    data = data[indices,])
  
  beta <- summary(model_2)$coefficients[2,1]
  
  
  formula_4 <- paste0(TE," ~ treat + ",SE)
  model_4 <- glm.nb(formula = as.formula(formula_4),
                    data = data[indices,])
  
  beta_s <- summary(model_4)$coefficients[2,1]
  
  PET <- 1-(beta_s/beta)
  
  return(c(beta,beta_s,PET))
}


