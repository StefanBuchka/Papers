# Purpose of this script: Simulation study of Poisson SEP and Gaussian CEP (transformed Poisson outcomes)
# Note: High computational power necessary! Long run time of 1.) and 3.)
# All cores (except 4) of your machine will be used. But you can adjust the cl variable of "pblapply".
# 25.11.2023
# Author: Stefan Buchka

##### 1.) Individual surrogacy (longitudinal joint model considering correlated endpoints)
#####
if(!file.exists("poisson_gaus_t"))
{dir.create("poisson_gaus_t")}

if(!file.exists("poisson_gaus_t/MA"))
{dir.create("poisson_gaus_t/MA")}

set.seed(82166)
source("simulation/simulation_dataset.R")

# Individual surrogacy
bayes_model_func <- function(x = unique(data$study),
                             data = data,
                             TE = "T",
                             SE = "St",
                             time = "year_TE",
                             path = "poisson_gaus_t/MA")
{
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(rlist)
  
  # select study within function
  data <- data[data$study == x,]
  data$year <- data$year_SE
  
  data_loop <- data
  
  library(brms)
  data_loop$visit_com <- data_loop$visit
  
  
  formula_TE <- bf(formula = paste0(TE,"~ year_TE"),
                    family = gaussian())
  
  
  
  formula_SE <- bf(formula = paste0(SE,"~ year_SE"),
                    family = gaussian())
  
  
  
  # Estimation of the model
  tryCatch({
    model <<- brm(formula_SE + formula_TE + set_rescor(TRUE),
                  iter = 15000,
                  chains = 4,
                  cores = 4,
                  save_pars = save_pars(all = T),
                  control = list(adapt_delta = 0.95),
                  thin = 10,
                  data = data_loop)
  },
  error = function(e)
  {
    errors <<- e
    if(!file.exists(paste0(path,"/errors"))){dir.create(paste0(path,"/errors"))}
    save(errors,file = paste0(path,"/errors","/study_",x,".Rdata"))
  },
  warning = function(w)
  {
    # Catches warnings
    warn <<- w
    
    if(!file.exists(paste0(path,"/warnings"))){dir.create(paste0(path,"/warnings"))}
    save(warn,file = paste0(path,"/warnings/study_",x,".Rdata"))
  
    
    model <<- brm(formula_SE + formula_TE + set_rescor(TRUE),
                  iter = 20000,
                  chains = 4,
                  cores = 4,
                  save_pars = save_pars(all = T),
                  control = list(adapt_delta = 0.95),
                  thin = 10,
                  data = data_loop)
    
    warn <- warnings()
    save(warn,file = paste0(path,"/warnings/study_",x,".Rdata"))
  })
  
  if(!file.exists(paste0(path,"/models"))){dir.create(paste0(path,"/models"))}
  save(model,file = paste0(path,"/models/study_",x,".Rdata"))
  
  # Computation of lambda
  res <- residuals(model, summary = F)
  res_SE <- res[,,gsub("_","",SE)]
  res_TE <- res[,,gsub("_","",TE)]
  
  prepare_func <- function(x. = x,
                           res_SE. = res_SE,
                           res_TE. = res_TE,
                           i. = i,
                           data. = data_loop)
  {
    result_SE <- res_SE.[i.,data.$visit_com == x.]
    result_TE <- res_TE.[i.,data.$visit_com == x.]
    
    df <- data.frame(id = data.$id[which(data.$visit_com == x.)],
                     visit_com = data.$visit_com[which(data.$visit_com == x.)],
                     SE = result_SE,
                     TE = result_TE)
    
    return(df)
  }
  
  
  # Several opportunities to calculate lambda: 
  # - Complete observation
  # - complete paierwise observations (so only NAs for the two columns compared are removed and not all NAs)
  
  # - cov = Covariance
  # - cor = pearson correlation
  # - sp = spearman correlation
  # - ken = Kendall correlation 
  
  sigma_list_cov <- vector(mode = "list", length = nrow(res_SE))
  sigma_list_per <- vector(mode = "list", length = nrow(res_SE))
  sigma_list_sp <- vector(mode = "list", length = nrow(res_SE))
  sigma_list_ken <- vector(mode = "list", length = nrow(res_SE))
  
  sigma_list_cov_c <- vector(mode = "list", length = nrow(res_SE))
  sigma_list_per_c <- vector(mode = "list", length = nrow(res_SE))
  sigma_list_sp_c <- vector(mode = "list", length = nrow(res_SE))
  sigma_list_ken_c <- vector(mode = "list", length = nrow(res_SE))
  
  # How many missings in time point
  missings_perc <- vector(mode = "list", length = nrow(res_SE))
  
  for(i in 1:nrow(res_SE))
  {
    time_points <- lapply(unique(data_loop$visit_com),prepare_func)
    time_points <- list_rbind(time_points)
    
    time_points <- time_points %>%
      distinct(id,visit_com, .keep_all = T) %>%
      pivot_wider(names_from = visit_com,
                  values_from = c(SE,TE)) %>%
      select(-id)
    
    missings_perc[[i]] <- sapply(time_points,function(x){sum(is.na(x))})/nrow(time_points)
    time_points <- time_points[,which(missings_perc[[i]] <= 0.65)]
    
    sigma_list_cov[[i]] <- cov(time_points, 
                               use = "pairwise.complete.obs")
    
    sigma_list_per[[i]] <- cor(time_points, 
                               use = "pairwise.complete.obs")
    
    sigma_list_sp[[i]] <- cor(time_points, 
                              use = "pairwise.complete.obs",
                              method = "spearman")
    
    sigma_list_ken[[i]] <- cor(time_points, 
                               use = "pairwise.complete.obs",
                               method = "kendall")
    
    sigma_list_cov_c[[i]] <- cov(time_points, 
                                 use = "complete.obs")
    
    sigma_list_per_c[[i]] <- cor(time_points, 
                                 use = "complete.obs")
    
    sigma_list_sp_c[[i]] <- cor(time_points, 
                                use = "complete.obs",
                                method = "spearman")
    
    sigma_list_ken_c[[i]] <- cor(time_points, 
                                 use = "complete.obs",
                                 method = "kendall")
    
  }
  
  
  missings_perc <- unique(list.rbind(missings_perc))
  
  # Computes lambda of all matrices
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
  
  if(!file.exists(paste0(path,"/missings_perc"))){dir.create(paste0(path,"/missings_perc"))}
  save(missings_perc,file = paste0(path,"/missings_perc/study_",x,".Rdata"))
  
  lambda_cov <- lapply(sigma_list_cov,r_lambda_func)
  lambda_per <- lapply(sigma_list_per,r_lambda_func)
  lambda_sp <- lapply(sigma_list_sp,r_lambda_func)
  lambda_ken <- lapply(sigma_list_ken,r_lambda_func)
  
  lambda_cov_c <- lapply(sigma_list_cov_c,r_lambda_func)
  lambda_per_c <- lapply(sigma_list_per_c,r_lambda_func)
  lambda_sp_c <- lapply(sigma_list_sp_c,r_lambda_func)
  lambda_ken_c <- lapply(sigma_list_ken_c,r_lambda_func)
  
  lambda_cov <- unlist(lambda_cov)
  lambda_per <- unlist(lambda_per)
  lambda_sp <- unlist(lambda_sp)
  lambda_ken <- unlist(lambda_ken)
  lambda_cov_c <- unlist(lambda_cov_c)
  lambda_per_c <- unlist(lambda_per_c)
  lambda_sp_c <- unlist(lambda_sp_c)
  lambda_ken_c <- unlist(lambda_ken_c)
  
  lambda <- list(lambda_cov,lambda_per,lambda_sp,lambda_ken,lambda_cov_c,lambda_per_c,lambda_sp_c,lambda_ken_c)
  names(lambda) <- c("lambda_cov","lambda_per","lambda_sp","lambda_ken","lambda_cov_c","lambda_per_c","lambda_sp_c","lambda_ken_c")
  
  if(!file.exists(paste0(path,"/lambda"))){dir.create(paste0(path,"/lambda"))}
  save(lambda,file = paste0(path,"/lambda/study_",x,".Rdata"))
  
  if(!file.exists(paste0(path,"/data"))){dir.create(paste0(path,"/data"))}
  save(data,file = paste0(path,"/data/study_",x,".Rdata"))
}

library(pbapply)
library(parallel)
cl <- makeCluster(detectCores()-2)

pblapply(X = unique(data$study),
         FUN = bayes_model_func,
         data = data,
         TE = "T",
         SE = "St",
         time = "year_TE",
         path = "poisson_gaus_t/MA",
         cl = cl)

stopCluster(cl)

rm(list = ls())
#####

##### 2.) Individual Surrogacy: Information Theoretic approach
#####
if(!file.exists("poisson_gaus_t"))
{dir.create("poisson_gaus_t")}

if(!file.exists("poisson_gaus_t/IT"))
{dir.create("poisson_gaus_t/IT")}

set.seed(82166)
source("simulation/simulation_dataset.R")


# Individual Level Surrogacy
bayes_model_func <- function(x = unique(data$study),
                             data = data,
                             TE = "T",
                             SE = "St",
                             time = "year_TE",
                             path = "poisson_gaus_t/IT")
{
  library(lme4)
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(nlme)
  library(ggcorrplot)
  library(rlist)
  
  
  # To stabilise the Fractional polynomials, the following linear transformation will be done
  data$year <- data$year_TE
  
  # select study within function
  data <- data[data$study == x,]
  
  # Model estimation
  source("scripts/calculation_r_2.R")
  
  # zero inflated model
  library(glmmTMB)
  data_loop <- data
  data_loop$visit_com <- data_loop$visit
  
  tryCatch({
    model_0 <- glmmTMB(as.formula(paste0(TE, "~ year")),
                       REML = T,
                       control = glmmTMBControl(optimizer = optim,
                                                optArgs = list(method = "BFGS")),
                       data = data_loop)
    
    model_1 <- glmmTMB(as.formula(paste0(TE, "~",SE,"+ year")),
                       REML = T,
                       control = glmmTMBControl(optimizer = optim,
                                                optArgs = list(method = "BFGS")),
                       data = data_loop)
  },
  error = function(e)
  {
    errors <<- geterrmessage()
    if(!file.exists(paste0(path,"/errors"))){dir.create(paste0(path,"/errors"))}
    save(errors,file = paste0(path,"/errors/errors_0_",x,".Rdata"))
  },
  warning = function(w)
  {
    model_0 <<- glmmTMB(as.formula(paste0(TE, "~ year")),
                        REML = T,
                        control = glmmTMBControl(optimizer = optim,
                                                 optArgs = list(method = "BFGS")),
                        data = data_loop)
    
    model_1 <<- glmmTMB(as.formula(paste0(TE, "~",SE,"+ year")),
                        REML = T,
                        control = glmmTMBControl(optimizer = optim,
                                                 optArgs = list(method = "BFGS")),
                        data = data_loop)
    warn <<- w
    
    if(!file.exists(paste0(path,"/warn"))){dir.create(paste0(path,"/warn"))}
    save(warn,file = paste0(path,"/warn/warn_0_",x,".Rdata"))
  })
  
  if(!file.exists(paste0(path,"/models"))){dir.create(paste0(path,"/models"))}
  save(model_1,file = paste0(path,"/models/model_1_",x,".Rdata"))
  save(model_0,file = paste0(path,"/models/model_0_",x,".Rdata"))
  
  
  R2ht <- r_2_function_IT(model0 = model_0,
                          model1 = model_1,
                          ordinal = F,
                          data = data_loop)
  
  if(!file.exists(paste0(path,"/R2ht"))){dir.create(paste0(path,"/R2ht"))}
  save(R2ht,file = paste0(path,"/R2ht/",x,".Rdata"))
  
  rm(R2ht)
  
  
  # ordinal TE
  
  
  library(ordinal)
  
  data_loop <- data_loop %>%
    mutate(#`T` = abs(`T`),
      `T` = ifelse(`T` >= 4,4,`T`),
      `T` = floor(`T`*2)/2,
      year = year_TE)
  
  
  tryCatch({
    model_0_or <- clmm2(factor(`T`) ~ year ,
                        Hess = T,
                        data = data_loop)
    
    
    model_1_or <- clmm2(factor(`T`) ~ S + year ,
                        Hess = T,
                        data = data_loop)
  },
  error = function(e)
  {
    errors_or <<- geterrmessage()
    if(!file.exists(paste0(path,"/errors_or"))){dir.create(paste0(path,"/errors_or"))}
    save(errors_or,file = paste0(path,"/errors_or/errors_0_",x,".Rdata"))
  },
  warning = function(w)
  {
    model_0_or <<- clmm2(factor(`T`) ~ year ,
                         Hess = T,
                         data = data_loop)
    
    
    model_1_or <<- clmm2(factor(`T`) ~ S + year ,
                         Hess = T,
                         data = data_loop)
    warn_or <<- w
    
    if(!file.exists(paste0(path,"/warn_or"))){dir.create(paste0(path,"/warn_or"))}
    save(warn_or,file = paste0(path,"/warn_or/warn_0_",x,".Rdata"))
  })
  
  
  if(!file.exists(paste0(path,"/models_or"))){dir.create(paste0(path,"/models_or"))}
  save(model_1_or,file = paste0(path,"/models_or/model_1_",x,".Rdata"))
  save(model_0_or,file = paste0(path,"/models_or/model_0_",x,".Rdata"))
  
  
  R2ht_or <- r_2_function_IT(model0 = model_0_or,
                             model1 = model_1_or,
                             ordinal = F,
                             data = data_loop)
  
  if(!file.exists(paste0(path,"/R2ht_or"))){dir.create(paste0(path,"/R2ht_or"))}
  save(R2ht_or,file = paste0(path,"/R2ht_or/",x,".Rdata"))
  
  rm(R2ht_or)
  
  rm(data_loop,
     data_SE_loop,
     data_loop)
  
  
  
}

library(pbapply)
library(parallel)
cl <- makeCluster(detectCores()-2)

pblapply(X = unique(data$study),
         FUN = bayes_model_func,
         data = data,
         TE = "T",
         SE = "St",
         time = "year_TE",
         path = "poisson_gaus_t/IT",
         cl = cl)

stopCluster(cl)

rm(list = ls())
#####

##### 3.) Individual Surrogacy: Information Theoretic approach with iteration
#####

if(!file.exists("poisson_gaus_t"))
{dir.create("poisson_gaus_t")}

if(!file.exists("poisson_gaus_t/IT_iter"))
{dir.create("poisson_gaus_t/IT_iter")}



start <- Sys.time()
data_i <- NULL

for(iter in 1:1000)
{
  set.seed(82165 + iter)
  source("simulation/simulation_dataset.R")
  data_i[[iter]] <- data
  rm(data)
  print(iter)
}



# Individual Level Surrogacy
bayes_model_func <- function(x = unique(data_i[[1]]$study),
                             data = data_i,
                             TE = "T",
                             SE = "St",
                             time = "year_TE",
                             path = "poisson_gaus_t/IT_iter")
{
  library(lme4)
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(nlme)
  library(ggcorrplot)
  library(rlist)
  
  errors <- NULL
  errors_or <- NULL
  R2ht <- NULL
  
  warn <- NULL
  warn_or <- NULL
  R2ht_or <- NULL
  
  r2_true <- NULL
  
  for(i_data in 1:length(data))
  {
    
    data_loop <- data[[i_data]]
    
    # To stabilise the Fractional polynomials, the following linear transformation will be done
    data_loop$year <- data_loop$year_TE
    
    # select study within function
    data_loop <- data_loop[data_loop$study == x,]
    
    #Get True R^2
    r2_true[[i_data]] <- unique(data_loop$R2_lambda)
    
    # Model estimation
    source("scripts/calculation_r_2.R")
    
    # zero inflated model
    library(glmmTMB)
    data_loop$visit_com <- data_loop$visit
    
    tryCatch({
      model_0 <- glmmTMB(as.formula(paste0(TE, "~ year")),
                         REML = T,
                         control = glmmTMBControl(optimizer = optim,
                                                  optArgs = list(method = "BFGS")),
                         data = data_loop)
      
      model_1 <- glmmTMB(as.formula(paste0(TE, "~",SE,"+ year")),
                         REML = T,
                         control = glmmTMBControl(optimizer = optim,
                                                  optArgs = list(method = "BFGS")),
                         data = data_loop)
    },
    error = function(e)
    {
      errors[[i_data]] <<- e
    },
    warning = function(w)
    {
      model_0 <<- glmmTMB(as.formula(paste0(TE, "~ year")),
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
      
      model_1 <<- glmmTMB(as.formula(paste0(TE, "~",SE,"+ year")),
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
      warn[[i_data]] <<- w
    })
    
    
    
    R2ht[[i_data]] <- r_2_function_IT(model0 = model_0,
                                      model1 = model_1,
                                      ordinal = F,
                                      data = data_loop)
    
    
    
    # ordinal TE
    
    
    library(ordinal)
    
    data_loop <- data_loop %>%
      mutate(#`T` = abs(`T`),
        `T` = ifelse(`T` >= 4,4,`T`),
        `T` = floor(`T`*2)/2,
        year = year_TE)
    
    
    tryCatch({
      model_0_or <- clmm2(factor(`T`) ~ year ,
                          Hess = T,
                          data = data_loop)
      
      
      model_1_or <- clmm2(factor(`T`) ~ S + year ,
                          Hess = T,
                          data = data_loop)
    },
    error = function(e)
    {
      errors_or[[i_data]] <<- e
      
    },
    warning = function(w)
    {
      model_0_or <<- clmm2(factor(`T`) ~ year ,
                           Hess = T,
                           data = data_loop)
      
      
      model_1_or <<- clmm2(factor(`T`) ~ S + year ,
                           Hess = T,
                           data = data_loop)
      warn_or[[i_data]] <<- w
    })
    
    
    
    R2ht_or[[i_data]] <- r_2_function_IT(model0 = model_0_or,
                                         model1 = model_1_or,
                                         ordinal = F,
                                         data = data_loop)
    
    rm(list = c("data_loop","model_0","model_1","model_1_or","model_0_or"))
    #print(i_data)
  }
  
  
  # Save true R^2
  r2_true <- list.rbind(r2_true)
  if(!file.exists(paste0(path,"/r2_true"))){dir.create(paste0(path,"/r2_true"))}
  save(r2_true,file = paste0(path,"/r2_true/",x,".Rdata"))
  
  # Save output of Gaussian models
  R2ht <- list.rbind(R2ht)
  if(!file.exists(paste0(path,"/R2ht"))){dir.create(paste0(path,"/R2ht"))}
  save(R2ht,file = paste0(path,"/R2ht/",x,".Rdata"))
  
  warn <- list.rbind(warn)
  if(!file.exists(paste0(path,"/warn")) & length(warn) > 0){dir.create(paste0(path,"/warn"))}
  if(length(warn) > 0){save(warn,file = paste0(path,"/warn/warn_0_",x,".Rdata"))}
  
  errors <- list.rbind(errors)
  if(!file.exists(paste0(path,"/errors")) & length(errors) > 0){dir.create(paste0(path,"/errors"))}
  if( length(errors) > 0){  save(errors,file = paste0(path,"/errors/errors_0_",x,".Rdata"))}
  
  
  
  #  Save output of ordinal models
  R2ht_or <- list.rbind(R2ht_or)
  if(!file.exists(paste0(path,"/R2ht_or"))){dir.create(paste0(path,"/R2ht_or"))}
  save(R2ht_or,file = paste0(path,"/R2ht_or/",x,".Rdata"))
  
  warn_or <- list.rbind(warn_or)
  if(!file.exists(paste0(path,"/warn_or")) & length(warn_or) > 0){dir.create(paste0(path,"/warn_or"))}
  if(length(warn_or) > 0){save(warn_or,file = paste0(path,"/warn_or/warn_0_",x,".Rdata"))}
  
  errors_or <- list.rbind(errors_or)
  if(!file.exists(paste0(path,"/errors_or")) & length(errors_or) > 0){dir.create(paste0(path,"/errors_or"))}
  if(length(errors_or) > 0){  save(errors_or,file = paste0(path,"/errors_or/errors_0_",x,".Rdata"))}
  
}

library(pbapply)
library(parallel)
cl <- makeCluster(detectCores()-2)

pblapply(X = unique(data_i[[1]]$study),
         FUN = bayes_model_func,
         data = data_i,
         TE = "T",
         SE = "St",
         time = "year_TE",
         path = "poisson_gaus_t/IT_iter",
         cl = cl)

stopCluster(cl)

end <- Sys.time()
end-start

rm(list = ls())
#####