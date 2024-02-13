
# Purpose: Elucidate the behavior of R²_lambda for different settings.
# Details: see paper and appendix 
# 25.11.23
# Author: Stefan Buchka

# There are several settings:
# 1.) Set correlation (gamma_ST) between SEP and CEP (not shown in paper)

# 2.) Correlation induced by alpha (CEP(t) = SEP(t) * alpha + epsilon) (not shown in paper) 

# 3.) Correlation induced by alpha (CEP(t) = SEP(t) * alpha + epsilon)
# Note: The correlation between SEP and CEP starts after the halve of the measurement time points

library(MASS)
library(rlist)
library(data.table)
library(tidyverse)
library(ggplot2)



data_sim_func <- function(N_tp = c(2:10,20,50,100), # Number of measurement time points
                          sd_obs = 0.5, # standard deviation of observations
                          rho_auto = seq(0.05,0.95,0.05), # rho for auto-regressive auto-correlation
                          gamma_ST = seq(0.05,0.95,0.05), # set correlation between SEP and CEP (setting 1)
                          alpha = seq(0.1,5,length = 19), # alpha induces correlation between SEP and CEP (CEP = SEP * alpha + epsilon)
                          epsilon = seq(0.05,3,length = 19), # Error for CEP
                          seed = 82166)
{
  set.seed(seed)
  # function for autocorrelation (AR)
  ar1_cor <- function(n,
                      rho,
                      inverse = F)
  {
    exponent <- abs(matrix(1:n - 1, 
                           nrow = n,
                           ncol = n,
                           byrow = T) -
                      (1:n - 1))
    
    if(inverse == T)
    {
      exponent <- abs(matrix(n:1 - 1, 
                             nrow = n,
                             ncol = n,
                             byrow = T) -
                        (1:n - 1))
    }
    
    return <- rho^exponent
    return(return)
  }
  
  # Function to compute Lambda
  # Input: A covariance Matrix of joint model with correlated residuals (Sigma). 
  r_lambda_func <- function(M)
  {
    index <- nrow(M)/2
    n <- nrow(M)
    
    if(nrow(M) > 2)
    {
      lambda <- det(M)/(det(M[1:index,1:index])*det(M[(index+1):n,(index+1):n]))
    }
    
    else
    {
      lambda <- det(M)/(M[1,1]*M[2,2])    
    }
    
    return(lambda)
  }
  
  pos <- 1
  data_all <- NULL
  
  # loop for time points
  for(i_tp in 1:length(N_tp))
  { # loop for different correlation values between S and T
    for(i_gamma in 1:length(gamma_ST))
    {
      for(i_epsilon in 1:length(epsilon))
      {
        for(i_rho_auto in 1:length(rho_auto))
        {
          
          # auto-correlation within an outcome
          cor_con <- ar1_cor(n = N_tp[i_tp],
                             rho = rho_auto[i_rho_auto])
          # correlation between SEP and TEP (set by gamma_ST)
          cor_dis <- gamma_ST[i_gamma] * cor_con
          
          # Correlation induced by alpha (CEP = SEP * alpha + epsilon)
          cor_dis2 <- alpha[i_gamma] * cor_con
          cor_con_t <- (alpha[i_gamma]^2 + epsilon[i_epsilon]^2) * cor_con
          
          # When the correlation starts after halve of the time points
          cor_dis_dyn <- matrix(0,
                                ncol = N_tp[i_tp],
                                nrow = N_tp[i_tp])
          
                    
          tp_halve <- ceiling(N_tp[i_tp]/2)
          cor_dis_dyn[,ceiling(N_tp[i_tp]/2 + 0.1):ncol(cor_dis2)] <- cor_dis2[,1:tp_halve]
          
          
          mm_1 <- matrix(0,
                         nrow = floor(N_tp[i_tp]/2),
                         ncol = tp_halve)
          
          mm_2 <- cor_con[1:tp_halve,
                          1:tp_halve]
          mm_3 <- rbind(mm_1,mm_2)
          
          cor_con_t_dyn <- diag(rep(epsilon[i_epsilon]^2,N_tp[i_tp])) +
            alpha[i_gamma]^2 * cbind(matrix(0,
                                            ncol = floor(N_tp[i_tp]/2),
                                            nrow = N_tp[i_tp]),
                                     mm_3)
          
          
          # Creation of covariance matrix (Sigma in paper)
          
          # Setting 1.)
          cor <- cbind(cor_con,cor_dis)
          cor <- rbind(cor,cbind(cor_dis,cor_con))
          
          # Setting 2.)
          cor2 <- cbind(cor_con,cor_dis2)
          cor2 <- rbind(cor2,cbind(cor_dis2,cor_con_t))
          
          # setting 3.)
          cor_dyn <- cbind(cor_con,cor_dis_dyn)
          cor_dyn <- rbind(cor_dyn,cbind(t(cor_dis_dyn),cor_con_t_dyn))
          
          #creation of variance-covariance matrix
          D <- diag(sd_obs,N_tp[i_tp]*2)
          cov <- D%*%cor%*%D
          
          # computation of real R²_lambda values
          lambda <- r_lambda_func(cov)
          lambda2 <- r_lambda_func(cor2)
          lambda_dyn <- r_lambda_func(cor_dyn)
          R2_lambda <- 1-lambda
          R2_lambda2 <- 1-lambda2
          R2_lambda_dyn <- 1-lambda_dyn
          
          # Creation of data set
          data <- data.frame(study = pos,
                             R2_lambda = R2_lambda,
                             R2_lambda2 = R2_lambda2,
                             R2_lambda_dyn = R2_lambda_dyn,
                             rho_auto = rho_auto[i_rho_auto],
                             gamma_ST = gamma_ST[i_gamma],
                             alpha = alpha[i_gamma],
                             epsilon = epsilon[i_epsilon],
                             tp = N_tp[i_tp])
          
          data_all[[pos]] <- data
          pos <- pos + 1
          rm(list = c("data"))
          
          rm(list = c("R2_lambda","R2_lambda2","R2_lambda_dyn","cor_con","cor_dis",
                      "cor_dis_dyn","cor","cor2","cor_dyn","D","cov"))
          
          
        }
      }
    }
    print(i_tp)
  }
  
  # Create one data set with all studies included
  data_all <- list.rbind(data_all)
  data_all$id <- as.numeric(as.factor(paste0(data_all$study,"_",data_all$id)))
  
  return(data_all)
}

data <- data_sim_func(N_tp = c(2,3,4,5,6,7,8,10,25,50),
                      sd_obs = 0.5,
                      rho_auto = 0.8,
                      gamma_ST = seq(0.01,0.99,0.01),
                      alpha = seq(0,5,length = 99),
                      epsilon = c(0.05,0.3,0.6,1,2,3),
                      seed = 82166)

data$`Number time points` <- as.factor(data$tp)


# Plot for setting 1.) (not shown in paper)
p <- data %>%
  filter(!is.infinite(R2_lambda)) %>%
  dplyr::select(R2_lambda,gamma_ST,`Number time points`) %>%
  distinct() %>%
  ggplot(aes(x = gamma_ST, y = R2_lambda, color = `Number time points`)) +
  geom_path() +
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  labs(title = expression(R[Lambda]^{2} ~ "derived from data with different number of measurement time points")) +
  ylab(expression(R[Lambda]^{2})) +
  xlab("Set correlation between surrogate and true endpoint")

p



# Plot for setting 2.) (not shown in paper)
p2 <- data %>%
  filter(!is.infinite(R2_lambda2)) %>%
  dplyr::select(R2_lambda2,`Number time points`,epsilon,alpha,rho_auto) %>%
  mutate(epsilon = paste0("epsilon = ",epsilon),
         rho_auto = paste0("autocorrelation = ",rho_auto)) %>%
  distinct() %>%
  ggplot(aes(x = alpha, y = R2_lambda2, color = `Number time points`)) +
  geom_jitter(height = 0.02, size = 0.25) +
  geom_path(linewidth = 0.25) +
  facet_wrap(. ~ epsilon) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme(text = element_text(size = 7)) +
  labs(title = expression(R[Lambda]^{2} ~ "when correlation between outcomes arises after the have of measurement time points")) +
  ylab(expression(R[Lambda]^{2})) +
  xlab(expression(alpha ~ "(factor to derive CEP)")) +
  scale_color_brewer(palette = "Paired")

p2



# Plot for setting 3.) (Supplement Figure 10)
p3 <- data %>%
  filter(!is.infinite(R2_lambda_dyn)) %>%
  mutate(epsilon = paste0("epsilon = ",epsilon),
         rho_auto = paste0("AR1 = ",rho_auto)) %>%
  ggplot(aes(x = alpha, y = R2_lambda_dyn, color = `Number time points`)) +
  geom_jitter(height = 0.02, size = 0.25) +
  geom_path(linewidth = 0.25) +
  facet_wrap(. ~ epsilon) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme(text = element_text(size = 7)) +
  labs(title = expression(R[Lambda]^{2} ~ "when correlation between outcomes arises after the have of measurement time points")) +
  ylab(expression(R[Lambda]^{2})) +
  xlab(expression(alpha ~ "(factor to derive CEP)")) +
  scale_color_brewer(palette = "Paired")

p3







