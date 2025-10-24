
# Purpose of this script: simulation of data for the simulation study
# Procedure: for details see "methods" in the paper and "simulation" in appendix. 
# Author: Stefan Buchka, 25.11.2023



library(MCMCpack)
library(Matrix)
library(mvtnorm)
library(GenOrd)
library(tidyverse)
library(rlist)

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

# Function to create auto-correlated covariance structure (auto-regressive first order).
ar1_cor <- function(n,
                    rho)
{
  exponent <- abs(matrix(1:n - 1, 
                         nrow = n,
                         ncol = n,
                         byrow = T) -
                    (1:n - 1))
  return <- rho^exponent
  return(return)
}

# Vector of numbers of subjects in simulated data
N_pat = c(100,300,600)

# Vector of numbers of measurement time points in simulated data
N_tp = c(2,4)

# Factor for auto-correlation (auto-regressive structure first order)
rho_auto = 0.8

# Alpha to induce correlation between SEP (surrogate endpoint) and CEP (clinical endpoint).
alpha = c(0.1,0.5,2.25)

# Lambda for data to generate count data
lambda = 2

# Mean of generated Gaussian data
mu = 0

# Error term for data simulation of CEP
# CEP(t) = SEP(t) * alpha + epsilon
epsilon <- 1

# Gives position to store the generated data sets
pos <- 1

# Object to store generated data sets
data_all <- NULL


for(i_tp in 1:length(N_tp))
{
  for(i_alpha in 1:length(alpha))
  {
    
    # Generation of variance-covariance matrix Sigma with auto-regressive (first order) auto-correlation structure to
    # generate Gaussian, auto-correlated data for SEP. To generate Poisson distributed count data from correlated 
    # Gaussian data, the method of Barbiero and Ferrari was used. Details in publication: 
    # Barbiero A, Ferrari PA. Simulation of correlated Poisson variables. Applied Stochastic Models in Business and Industry. 2015;31:669–80.
    
    # number of correlated variables (in this case it is equivalent to the number
    # of time points, since a auto-correlation is simulated)
    k <- N_tp[i_tp]
    
    # configuration of distribution parameters
    lambda <- rep(lambda,N_tp[i_tp])
    
    # truncation error
    eps <- 0.0001
    
    # corresponding maximum value
    kmax <- qpois(1 - eps,lambda)
    
    
    # preliminary stage: TRUNCATION
    l <- NULL
    for(i in 1:k)
    {
      l[[i]] <- 0:kmax[i]
    }
    
    marg <- NULL
    for(i in 1:k)
    {
      marg[[i]] <- dpois(0:kmax[i],lambda[i])
      marg[[i]][kmax[i] + 1] <- 1 - sum(marg[[i]][1:(kmax[i])])
    }
    
    cm <- NULL
    for(i in 1:k)
    {
      cm[[i]] <- cumsum(marg[[i]])
      cm[[i]] <- cm[[i]][-(kmax[i] + 1)]
    }
    
    # Creation of auto-regressive (order 1) variance-covariance matrix
    cov_s <- ar1_cor(N_tp[i_tp],
                     rho = rho_auto)
    
    # Method of Barbiero and Ferrari 
    res <- ordcont(cm,cov_s,support = l)
    cov_s <- res[[1]]
    
    # Include error in Sigma
    D <- diag(epsilon,N_tp[i_tp])
    cov_s <- D%*%cov_s%*%D
    
    for(i_pat in 1:length(N_pat))
    {
      
      # Generate multivariate Gaussian data (SEP) with covariance matrix of SEP 
      data_s <- mvrnorm((N_pat[i_pat]),
                        c(rep(mu,k)),
                        cov_s)
      
      # Generate CEP from SEP
      # CEP(t) = SEP(t) * alpha + epsilon
      data_t <- apply(data_s,1,FUN = function(x){mvrnorm(1,
                                                         mu = (alpha[i_alpha] * x),
                                                         Sigma = cov_s)})
      data_t <- t(data_t)
      
      # Covariance matrix between SEP and CEP
      cov_st <- alpha[i_alpha] * cov_s
      
      # Covariance of CEP
      cov_t <- (alpha[i_alpha]^2 + epsilon^2) * cov_s
      
      # Covariance to calculate Lambda (Sigma)
      M_cov <- cbind(cov_s,cov_st)
      M_cov <- rbind(M_cov,cbind(cov_st,cov_t))
      
      # Generation of Poisson distributed SEP and CEP
      data_tp <- matrix(NA,ncol = ncol(data_s),nrow = nrow(data_s))
      data_sp <- matrix(NA,ncol = ncol(data_s),nrow = nrow(data_s))
       
      for(i in 1:k)
      {
        # Basis of this method is the Gaussian distributed SEP and CEP
        # Method according to Barbiero and Ferrari
        data_sp[,i] <- qpois(pnorm(data_s[,i],
                                  mean = mu,
                                  sd = epsilon),
                            lambda[i])
        
        data_tp[,i] <- qpois(pnorm(data_t[,i],
                                   mean = mu,
                                   sd = epsilon),
                             lambda[i])
        
        if(any(is.infinite(data_sp[,i])))
        {
          data_sp[is.infinite(data_sp[,i]),i] <- max(data_tp[is.finite(data_sp[,i]),i])
        }
        
        if(any(is.infinite(data_tp[,i])))
        {
          data_tp[is.infinite(data_tp[,i]),i] <- max(data_tp[is.finite(data_tp[,i]),i])
        }
        
      }
       
      data_s <- as.numeric(t(data_s))
      data_t <- as.numeric(t(data_t))
      data_sp <- as.numeric(t(data_sp))
      data_tp <- as.numeric(t(data_tp))
      
      
      # Transformation of Poisson distributed data to stabilize variance 
      # according to:
      # Huber W, Von Heydebreck A, Sültmann H, Poustka A, Vingron M. Variance stabilization applied to microarray data calibration and to the quantification of differential expression. Bioinformatics. 2002;18 suppl_1:S96–104. 
      
      source("scripts/trans_func.R")
      data_st <- trans_func(var = data_sp,
                            trans = "Huber",
                            bins = NULL)
      
      data_tt <- trans_func(var = data_tp,
                            trans = "Huber",
                            bins = NULL)
      
      
      # Creation of data set
      df <- data.frame(id = sort(rep(1:N_pat[i_pat],N_tp[i_tp])),
                       S = data_s,
                       `T` = data_t,
                       Sp = data_sp,
                       Tp = data_tp,
                       St = data_st,
                       Tt = data_tt,
                       treat = c(rep(-1,length(data_s)/2),
                                 rep(1,length(data_s)/2)))
      
      df$visit <- as.factor(rep(1:N_tp[i_tp],N_pat[i_pat]))
      df$R2_lambda <- 1 - r_lambda_func(M_cov)
      df$alpha <- alpha[i_alpha]
      df$study <- pos
      df$study_n <- N_pat[i_pat]
      df$tp <- N_tp[i_tp]
      
      # Add time variable (in dependency to visit)
      df$year_TE <- (as.numeric(df$visit) - 1)/2 + 0.1
      df$year_SE <- (as.numeric(df$visit) - 1)/2 + 0.1
      
      data_all[[pos]] <- df
      pos <- pos + 1
      
      rm(list = c("data_s","data_t","df"))
    }
  }
  rm(k)
}

# Combine all created data sets
data <- list.rbind(data_all)
data$id <- as.numeric(as.factor(paste0(data$study,"_",data$id)))

