# Purpose of this script: Simulation study of Gaussian SEP and Poisson CEP
# Note: High computational power necessary! Long run time of 1.) and 3.)
# All cores (except 4) of your machine will be used. But you can adjust the cl variable of "pblapply".
# 25.11.2023
# Author: Stefan Buchka

##### 1.) Individual surrogacy (longitudinal joint model considering correlated endpoints)
#####

if(!file.exists("gaus_poisson"))
{dir.create("gaus_poisson")}

if(!file.exists("gaus_poisson/MA"))
{dir.create("gaus_poisson/MA")}

set.seed(82166)
source("simulation/simulation_dataset.R")

bayes_model_func <- function(x = unique(data$study),
                             data = data,
                             poisson = "Tp",
                             gaus = "S",
                             time = "year",
                             path = "gaus_poisson/MA")
{
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(rlist)
  
  # select study within function
  data <- data[data$study == x,]
  data$year <- data$year_TE
 
  # prepare data 
  data_loop <- data
  
 
  library(brms)
  
  # Creation poisson formulas for model
  data_loop$visit_com <- data_loop$visit
  data_loop$obs <- 1:nrow(data_loop)

  formula_poisson <- bf(formula = paste0(poisson,"~ offset(log(year_TE)) + year_TE + (1|c|obs)"),
                         family = poisson())
  
  
  formula_gaus <- bf(formula = paste0(gaus,"~ year_SE + (1|c|obs)"),
                      family = gaussian(),
                      sigma = 0.1)
  
  
  # Estimation of the model
  tryCatch({
    model <<- brm(formula_poisson + formula_gaus + set_rescor(F),
                  iter = 15000,
                  chains = 4,
                  cores = 4,
                  save_pars = save_pars(all = T),
                  control = list(adapt_delta = 0.95,
                                 max_treedepth = 10),
                  thin = 10,
                  data = data_loop)
  },
  error = function(e)
  {
    errors <<- e
    if(!file.exists(paste0(path,"/errors"))){dir.create(paste0(path,"/errors"))}
    save(errors,file = paste0(path,"/errors/study_",x,".Rdata"))
  },
  warning = function(w)
  {
    # Catches warnings
    warn <<- w
    
    if(!file.exists(paste0(path,"/warnings"))){dir.create(paste0(path,"/warnings"))}
    save(warn,file = paste0(path,"/warnings/study_",x,".Rdata"))
    
    
    model <<- brm(formula_poisson + formula_gaus + set_rescor(F),
                  iter = 20000,
                  chains = 4,
                  cores = 4,
                  save_pars = save_pars(all = T),
                  control = list(adapt_delta = 0.99),
                  thin = 10,
                  data = data_loop)
    
    
    warn <<- warnings()
    save(warn,file = paste0(path,"/warnings/study_",x,".Rdata"))
  })
  
  
  
  if(!file.exists(paste0(path,"/models"))){dir.create(paste0(path,"/models"))}
  save(model,file = paste0(path,"/models/study_",x,".Rdata"))
  
  # Computation of lambda
  res <- ranef(model, summary = F)$obs
  res_SE <- res[,,"S_Intercept"]
  res_TE <- res[,,"Tp_Intercept"]
  
  
  
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
  
  
  if(!file.exists(paste0(path,"/missings"))){dir.create(paste0(path,"/missings"))}
  save(missings_perc,file = paste0(path,"/missings/study_",x,".Rdata"))
  
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
  save(data_loop,file = paste0(path,"/data/study_",x,".Rdata"))
  
  rm(data_loop)
  
  
}

library(pbapply)
library(parallel)
cl <- makeCluster(detectCores()-4)

pblapply(X = unique(data$study),
         FUN = bayes_model_func,
         data = data,
         poisson = "Tp",
         gaus = "S",
         time = "year",
         path = "gaus_poisson/MA",
         cl = cl)

stopCluster(cl)

cl <- makeCluster(detectCores()-2)




rm(list = ls())
#####

##### 2.) Individual Surrogacy: Information Theoretic approach
#####

if(!file.exists("gaus_poisson"))
{dir.create("gaus_poisson")}

if(!file.exists("gaus_poisson/IT"))
{dir.create("gaus_poisson/IT")}

set.seed(82166)
source("simulation/simulation_dataset.R")

bayes_model_func <- function(x = unique(data_delta$study),
                             data = data,
                             SE = "S",
                             TE = "Tp",
                             time = "year_TE",
                             path = "gaus_poisson/IT")
{
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(glmmTMB)
  
  # select study within function
  data <- data[data$study == x,]
  
  # transform day variable to year
  data$year <- unlist(data[,time])
  
  data_loop <- data
  
  # Model estimation
  source("scripts/calculation_r_2.R")
  
  tryCatch({
    model_0_zi <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE")),
                          family = "poisson", 
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          ziformula = ~ 1,
                          data = data_loop)
    
    model_1_zi <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE")),
                          family = "poisson",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          ziformula = ~ 1,
                          data = data_loop)
  },
  error = function(e)
  {
    errors_zi <<- e
    if(!file.exists(paste0(path,"/errors_zi"))){dir.create(paste0(path,"/errors_zi"))}
    save(errors_zi,file = paste0(path,"/errors_zi/errors_0_zi_",x,".Rdata"))
  },
  warning = function(w)
  {
    model_0_zi <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           ziformula = ~ 1,
                           data = data_loop)
    
    model_1_zi <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           ziformula = ~ 1,
                           data = data_loop)
    warn_zi <<- w
    
    if(!file.exists(paste0(path,"/warn_zi"))){dir.create(paste0(path,"/warn_zi"))}
    save(warn_zi,file = paste0(path,"/warn_zi/warn_0_zi_",x,".Rdata"))
  })
  
  
  if(!file.exists(paste0(path,"/models_zi"))){dir.create(paste0(path,"/models_zi"))}
  save(model_1_zi,file = paste0(path,"/models_zi/model_1_zi_",x,".Rdata"))
  save(model_0_zi,file = paste0(path,"/models_zi/model_0_zi_",x,".Rdata"))
  
  
  R2ht_zi <- r_2_function_IT(model0 = model_0_zi,
                             model1 = model_1_zi,
                             ordinal = F,
                             data = data_loop)
  
  if(!file.exists(paste0(path,"/R2ht_zi"))){dir.create(paste0(path,"/R2ht_zi"))}
  save(R2ht_zi,file = paste0(path,"/R2ht_zi/",x,".Rdata"))
  
  rm(R2ht_zi)
  
  
  
  # negative binomial model
  
  tryCatch({
    model_0_nb <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE")),
                          family = "nbinom2",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
    
    model_1_nb <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE")),
                          family = "nbinom2",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
  },
  error = function(e)
  {
    errors_nb <<- geterrmessage()
    if(!file.exists(paste0(path,"/errors_nb"))){dir.create(paste0(path,"/errors_nb"))}
    save(errors_nb,file = paste0(path,"/errors_nb/errors_0_nb_",x,".Rdata"))
  },
  warning = function(w)
  {
    warn_nb <<- w
    
    model_0_nb <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE ")),
                           family = "nbinom2",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
    
    model_1_nb <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + ",SE,"+ year_TE ")),
                           family = "nbinom2",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
    
    if(!file.exists(paste0(path,"/warn_nb"))){dir.create(paste0(path,"/warn_nb"))}
    save(warn_nb,file = paste0(path,"/warn_nb/warn_0_nb_",x,".Rdata"))
  })
  
  
  if(!file.exists(paste0(path,"/models_nb"))){dir.create(paste0(path,"/models_nb"))}
  save(model_1_nb,file = paste0(path,"/models_nb/model_1_nb_",x,".Rdata"))
  save(model_0_nb,file = paste0(path,"/models_nb/model_0_nb_",x,".Rdata"))
  
  
  R2ht_nb <- r_2_function_IT(model0 = model_0_nb,
                             model1 = model_1_nb,
                             ordinal = F,
                             data = data_loop)
  
  if(!file.exists(paste0(path,"/R2ht_nb"))){dir.create(paste0(path,"/R2ht_nb"))}
  save(R2ht_nb,file = paste0(path,"/R2ht_nb/",x,".Rdata"))
  
  rm(R2ht_nb)
  
  
  
  
  
  # poisson model
  
  
  tryCatch({
    model_0_po <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE")),
                          family = "poisson",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
    
    model_1_po <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + ",SE,"+ year_TE")),
                          family = "poisson",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
  },
  error = function(e)
  {
    errors_po <<- geterrmessage()
    if(!file.exists(paste0(path,"/errors_po"))){dir.create(paste0(path,"/errors_po"))}
    save(errors_po,file = paste0(path,"/errors_po/errors_0_po_",x,".Rdata"))
  },
  warning = function(w)
  {
    warn_po <<- w
    
    model_0_po <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
    
    model_1_po <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
    
    if(!file.exists(paste0(path,"/warn_po"))){dir.create(paste0(path,"/warn_po"))}
    save(warn_po,file = paste0(path,"/warn_po/warn_0_po_",x,".Rdata"))
  })
  
  
  if(!file.exists(paste0(path,"/models_po"))){dir.create(paste0(path,"/models_po"))}
  save(model_1_po,file = paste0(path,"/models_po/model_1_po_",x,".Rdata"))
  save(model_0_po,file = paste0(path,"/models_po/model_0_po_",x,".Rdata"))
  
  
  R2ht_po <- r_2_function_IT(model0 = model_0_po,
                             model1 = model_1_po,
                             ordinal = F,
                             data = data_loop)
  
  if(!file.exists(paste0(path,"/R2ht_po"))){dir.create(paste0(path,"/R2ht_po"))}
  save(R2ht_po,file = paste0(path,"/R2ht_po/",x,".Rdata"))
  
  rm(R2ht_po)
  
  
  # Ordinal models
  library(ordinal)
  
  tryCatch({
    data_loop$Tp <- as.factor(as.numeric(cut_number(data_loop$Tp,2)))
  },
  error = function(e)
  {
    data_loop$Tp <<- as.factor(data_loop$Tp)
  })
  
  tryCatch({
    model_0_or <- clmm2(Tp ~ year_TE ,
                        Hess = T,
                        data = data_loop)
    
    
    model_1_or <- clmm2(Tp ~ S + year_TE ,
                        Hess = T,
                        data = data_loop)
  },
  error = function(e)
  {
    errors_or <<- geterrmessage()
    if(!file.exists(paste0(path,"/errors_or"))){dir.create(paste0(path,"/errors_or"))}
    save(errors_or,file = paste0(path,"/errors_or/errors_0_or_",x,".Rdata"))
  },
  warning = function(w)
  {
    warn_or <<- w
    
    model_0_or <<- clmm2(Tp ~ year_TE ,
                         Hess = T,
                         data = data_loop)
    
    
    model_1_or <<- clmm2(Tp ~ S + year_TE ,
                         Hess = T,
                         data = data_loop)
    
    if(!file.exists(paste0(path,"/warn_or"))){dir.create(paste0(path,"/warn_or"))}
    save(warn_or,file = paste0(path,"/warn_or/warn_0_or_",x,".Rdata"))
  })
  
  R2ht_or <- r_2_function_IT(model0 = model_0_or,
                             model1 = model_1_or,
                             ordinal = T,
                             TE = TE,
                             data = data_loop)
  
  if(!file.exists(paste0(path,"/R2ht_or"))){dir.create(paste0(path,"/R2ht_or"))}
  save(R2ht_or,file = paste0(path,"/R2ht_or/",x,".Rdata"))
  
  rm(R2ht_or)
  
  rm(data_loop)
  
  
  
}

library(pbapply)
library(parallel)
cl <- makeCluster(detectCores()-2)

pblapply(X = unique(data$study),
         FUN = bayes_model_func,
         data = data,
         SE = "S",
         TE = "Tp",
         time = "year_TE",
         path = "gaus_poisson/IT",
         cl = cl)

stopCluster(cl)

rm(list = ls())
#####

##### 3.) Individual Surrogacy: Information Theoretic approach with iteration
#####
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

if(!file.exists("gaus_poisson"))
{dir.create("gaus_poisson")}

if(!file.exists("gaus_poisson/IT_iter"))
{dir.create("gaus_poisson/IT_iter")}



bayes_model_func <- function(x = unique(data_i[[1]]$study),
                             data = data_i,
                             SE = "S",
                             TE = "Tp",
                             time = "year_TE",
                             path = "gaus_poisson/IT_iter")
{
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(glmmTMB)
  library(rlist)
  
  
  errors <- NULL
  warn <- NULL
  R2ht <- NULL
  
  errors_zi <- NULL
  warn_zi <- NULL
  R2ht_zi <- NULL
  
  errors_nb <- NULL
  warn_nb <- NULL
  R2ht_nb <- NULL
  
  errors_po <- NULL
  warn_po <- NULL
  R2ht_po <- NULL
  
  errors_or <- NULL
  warn_or <- NULL
  R2ht_or <- NULL
  
  r2_true <- NULL
  
for(i_data in 1:length(data))
{
  
  data_loop <- data[[i_data]]
  
  # select study within function
  data_loop <- data_loop[data_loop$study == x,]
  
  # transform day variable to year
  data_loop$year <- unlist(data_loop[,time])
  
  #Get True R^2
  r2_true[[i_data]] <- unique(data_loop$R2_lambda)
  
  # Model estimation
  source("scripts/calculation_r_2.R")
  
  model_0_zi <- NULL
  model_1_zi <- NULL
  
  tryCatch({
    model_0_zi <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE")),
                          family = "poisson", 
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          ziformula = ~ 1,
                          data = data_loop)
    
    model_1_zi <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE")),
                          family = "poisson",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          ziformula = ~ 1,
                          data = data_loop)
  },
  error = function(e)
  {
    errors_zi[[i_data]] <<- e
  },
  warning = function(w)
  {
    model_0_zi <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           ziformula = ~ 1,
                           data = data_loop)
    
    model_1_zi <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           ziformula = ~ 1,
                           data = data_loop)
    warn_zi[[i_data]] <<- w
    
  })
  
  
  
  
  R2ht_zi[[i_data]] <- r_2_function_IT(model0 = model_0_zi,
                             model1 = model_1_zi,
                             ordinal = F,
                             data = data_loop)

  
  # negative binomial model
  
  model_0_nb <- NULL
  model_1_nb <- NULL
  
  tryCatch({
    model_0_nb <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE")),
                          family = "nbinom2",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
    
    model_1_nb <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE")),
                          family = "nbinom2",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
  },
  error = function(e)
  {
    errors_nb[[i_data]] <<- geterrmessage()
  },
  warning = function(w)
  {
    warn_nb[[i_data]] <<- w
    
    model_0_nb <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE ")),
                           family = "nbinom2",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
    
    model_1_nb <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + ",SE,"+ year_TE ")),
                           family = "nbinom2",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
  })
  
  
  R2ht_nb[[i_data]] <- r_2_function_IT(model0 = model_0_nb,
                             model1 = model_1_nb,
                             ordinal = F,
                             data = data_loop)

  
  # poisson model
  
  model_0_po <- NULL
  model_1_po <- NULL
  
  tryCatch({
    model_0_po <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE")),
                          family = "poisson",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
    
    model_1_po <- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + ",SE,"+ year_TE")),
                          family = "poisson",
                          REML = T,
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS")),
                          data = data_loop)
  },
  error = function(e)
  {
    errors_po[[i_data]] <<- geterrmessage()
  },
  warning = function(w)
  {
    warn_po[[i_data]] <<- w
    
    model_0_po <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) + year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
    
    model_1_po <<- glmmTMB(as.formula(paste0(TE, "~ offset(log(year_TE)) +",SE,"+ year_TE ")),
                           family = "poisson",
                           REML = T,
                           control = glmmTMBControl(optimizer = optim,
                                                    optArgs = list(method = "BFGS")),
                           data = data_loop)
    
  })
  
  
  
  R2ht_po[[i_data]] <- r_2_function_IT(model0 = model_0_po,
                             model1 = model_1_po,
                             ordinal = F,
                             data = data_loop)
  
  
  # Ordinal models
  library(ordinal)
  
  
  model_0_or <- NULL
  model_1_or <- NULL
  
  tryCatch({
    data_loop$Tp <- as.factor(as.numeric(cut_number(data_loop$Tp,2)))
  },
  error = function(e)
  {
    data_loop$Tp <<- as.factor(data_loop$Tp)
  })
  
  tryCatch({
    model_0_or <- clmm2(Tp ~ year_TE ,
                        Hess = T,
                        data = data_loop)
    
    
    model_1_or <- clmm2(Tp ~ S + year_TE ,
                        Hess = T,
                        data = data_loop)
  },
  error = function(e)
  {
    errors_or[[i_data]] <<- geterrmessage()
  },
  warning = function(w)
  {
    warn_or[[i_data]] <<- w
    
    model_0_or <<- clmm2(Tp ~ year_TE ,
                         Hess = T,
                         data = data_loop)
    
    
    model_1_or <<- clmm2(Tp ~ S + year_TE ,
                         Hess = T,
                         data = data_loop)
  })
  
  R2ht_or[[i_data]] <- r_2_function_IT(model0 = model_0_or,
                             model1 = model_1_or,
                             ordinal = T,
                             TE = TE,
                             data = data_loop)
  
 
  rm(list = c("data_loop",grep("model",ls(),value = T)))
}
  
  # Save true R^2
  r2_true <- list.rbind(r2_true)
  if(!file.exists(paste0(path,"/r2_true"))){dir.create(paste0(path,"/r2_true"))}
  save(r2_true,file = paste0(path,"/r2_true/",x,".Rdata"))
  
  #  Save output of zero inflated models
  R2ht_zi <- list.rbind(R2ht_zi)
  if(!file.exists(paste0(path,"/R2ht_zi"))){dir.create(paste0(path,"/R2ht_zi"))}
  save(R2ht_zi,file = paste0(path,"/R2ht_zi/",x,".Rdata"))
  
  warn_zi <- list.rbind(warn_zi)
  if(!file.exists(paste0(path,"/warn_zi")) & length(warn_zi) > 0){dir.create(paste0(path,"/warn_zi"))}
  if(length(warn_zi) > 0){save(warn_zi,file = paste0(path,"/warn_zi/warn_0_",x,".Rdata"))}
  
  errors_zi <- list.rbind(errors_zi)
  if(!file.exists(paste0(path,"/errors_zi")) & length(errors_zi) > 0){dir.create(paste0(path,"/errors_zi"))}
  if(length(errors_zi) > 0){  save(errors_zi,file = paste0(path,"/errors_zi/errors_0_",x,".Rdata"))}
  
  
  #  Save output of poisson models
  R2ht_po <- list.rbind(R2ht_po)
  if(!file.exists(paste0(path,"/R2ht_po"))){dir.create(paste0(path,"/R2ht_po"))}
  save(R2ht_po,file = paste0(path,"/R2ht_po/",x,".Rdata"))
  
  warn_po <- list.rbind(warn_po)
  if(!file.exists(paste0(path,"/warn_po")) & length(warn_po) > 0){dir.create(paste0(path,"/warn_po"))}
  if(length(warn_po) > 0){save(warn_po,file = paste0(path,"/warn_po/warn_0_",x,".Rdata"))}
  
  errors_po <- list.rbind(errors_po)
  if(!file.exists(paste0(path,"/errors_po")) & length(errors_po) > 0){dir.create(paste0(path,"/errors_po"))}
  if(length(errors_po) > 0){  save(errors_po,file = paste0(path,"/errors_po/errors_0_",x,".Rdata"))}
  
  
  #  Save output of negative binomial models
  R2ht_nb <- list.rbind(R2ht_nb)
  if(!file.exists(paste0(path,"/R2ht_nb"))){dir.create(paste0(path,"/R2ht_nb"))}
  save(R2ht_nb,file = paste0(path,"/R2ht_nb/",x,".Rdata"))
  
  warn_nb <- list.rbind(warn_nb)
  if(!file.exists(paste0(path,"/warn_nb")) & length(warn_nb) > 0){dir.create(paste0(path,"/warn_nb"))}
  if(length(warn_nb) > 0){save(warn_nb,file = paste0(path,"/warn_nb/warn_0_",x,".Rdata"))}
  
  errors_nb <- list.rbind(errors_nb)
  if(!file.exists(paste0(path,"/errors_nb")) & length(errors_nb) > 0){dir.create(paste0(path,"/errors_nb"))}
  if(length(errors_nb) > 0){  save(errors_nb,file = paste0(path,"/errors_nb/errors_0_",x,".Rdata"))}
  
  
  
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
         SE = "S",
         TE = "Tp",
         time = "year_TE",
         path = "gaus_poisson/IT_iter",
         cl = cl)

stopCluster(cl)

rm(list = ls())
#####

