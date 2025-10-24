
# Purpose: Estimation of a Bayesian Gaussian-Gaussian joint model to derive R²_lambda
# Details: see paper and appendix 
# 25.11.23
# Author: Stefan Buchka

# Note: The clinical endpoint (CEP) is named TE and the surrogate endpoint (SEP) is named SE in the following 
# Note: long run time with high computational power. Parallelisation will be used. All cores except 4 will be in use. 
# You can change this via cl variable of the "pblapply" function

# Load data
library(tidyverse)
 load("mri.RData")
 load("edss.RData")


data_TE <- edss %>%
  dplyr::select(study,id,arm,edss_score,day,visit,visit_com) %>%
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(edss$day,edss$day))-1)*-1) %>%
  group_by(id,study) %>%
  mutate(delete = ifelse(!any(visit_com == 0),1,0)) %>% # remove patients without baseline measurement
  filter(delete == 0) %>%
  dplyr::select(-delete)

dim(data_TE)
dim(data_TE[complete.cases(data_TE),])

data_SE <- mri %>%
  dplyr::select(study,id,arm,t2_volume,day,visit,visit_com,unit) %>%
  filter(!is.na(t2_volume)) %>%
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(mri$day,mri$day))-1)*-1) %>%
  mutate(t2_volume = ifelse(unit == "mm3",t2_volume/1000,t2_volume),
         t2_volume = log(t2_volume + 0.1)) %>%
  dplyr::select(-unit) %>%
group_by(id,study) %>%
  mutate(delete = ifelse(!any(visit_com == 0),1,0)) %>% # remove patients without baseline measurement
  filter(delete == 0) %>%
  dplyr::select(-delete)

dim(data_SE)
dim(data_SE[complete.cases(data_SE),])

data_TE <- data_TE[which(data_TE$study %in% data_SE$study),]


# Individual surrogacy
bayes_model_func <- function(x = unique(data_TE$study),
                             data_TE = data_TE,
                             data_SE = data_SE,
                             forced_powers = NULL,
                             TE = "edss_score",
                             SE = "t2_volume",
                             time_SE = "day_mri",
                             time_TE = "day_edss")
{
  library(lme4)
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(nlme)
  library(ggcorrplot)
  library(rlist)
  
  # select study within function
  data_TE <- data_TE[data_TE$study == x,]
  data_SE <- data_SE[data_SE$study == x,]
  
  # transform day variable to year
  data_TE$year <- unlist(data_TE[,"day"]/365.25)
  data_SE$year <- unlist(data_SE[,"day"]/365.25)
  
  # To stabilize the Fractional polynomials, the following linear transformation will be done
  if(any(data_TE$year < 0.1))
  {
    data_TE$year <- data_TE$year + 0.1
  }
  
  if(any(data_SE$year < 0.1))
  {
    data_SE$year <- data_SE$year + 0.1
  }
  
  
  # figure out comparisons between control and intervention group
  control_gr <- c("SC IFNB-1a","Placebo","Interferon beta-1a","PLACEBO","INTERFERON BETA-1a","AVONEX")
  arms <- as.character(unique(data_TE$arm))
  control_gr <- arms[which(arms %in% control_gr)]
  test_gr <- arms[-which(arms %in% control_gr)]
  
  if("" %in% test_gr)
  {
    test_gr <- test_gr[-which(test_gr == "")]
  }
  
  comparisons <- crossing(control_gr,test_gr)
  
  # The analysis have to be carried out for each active drug vs. control group arm
  for(comp in 1:nrow(comparisons))
  {
    # prepare data 
    data_TE_loop <- data_TE %>%
      filter(arm %in% c(comparisons$control_gr[comp],comparisons$test_gr[comp]))
    
    data_SE_loop <- data_SE %>%
      filter(arm %in% c(comparisons$control_gr[comp],comparisons$test_gr[comp]))
    
    data_TE_loop$treat <- as.factor(ifelse(data_TE_loop$arm %in% control_gr,-1,1))
    data_SE_loop$treat <- as.factor(ifelse(data_SE_loop$arm %in% control_gr,-1,1))
    
    # Calculate Fractional polynomials
    source("scripts/fp_function.R")
    
    Covariate <- list(mode = "list", length = 2)
    Outcome <- list(mode = "list", length = 2)
    AIC_range <- list(mode = "list", length = 2)
    
    id <- data_TE_loop$id
    
    
    FP_TE <-   fp_func(Covariate = "year",
                       Outcome = TE, 
                       id = id,
                       Max.M = 2,
                       S = c(-2,-1,-0.5,0,0.5,1,2,3),
                       Dataset = as.data.frame(data_TE_loop))
    
    Covariate[[1]] <- data_TE_loop$year
    Outcome[[1]] <- eval(parse(text = paste0("data_TE_loop$",TE)))
  
    
    
    rm(id)
    id <- data_SE_loop$id
    FP_SE <-   fp_func(Covariate = "year",
                       Outcome = SE, 
                       id = id,
                       Max.M = 2,
                       S = c(-2,-1,-0.5,0,0.5,1,2,3),
                       Dataset = as.data.frame(data_SE_loop))
    
    rm(id)
    
    Covariate[[2]] <- data_SE_loop$year
    Outcome[[2]] <- eval(parse(text = paste0("data_SE_loop$",SE)))
    
    
    FP <- list(FP_TE,FP_SE)
    
    if(!file.exists("results/FP")){dir.create("results/FP")}
    save(FP,file = paste0("results/FP/",x,"_",comp,".Rdata"))
    
    term <- vector(mode = "list", length = 2)
    names(term) <- c("TE","SE")
    
    for(i in 1:length(FP))
    {
      
      #  FP degree == 1
      if(is.na(FP[[i]]$power2) & is.null(forced_powers))
      {
        power <- FP[[i]]$power1
        
        if(length(power) > 1)
        {
          power <- power[1]
        }
        
        if (power[1] != 0) {
          term1 <- Covariate[[i]]^power[1]
        }
        if (power[1] == 0) {
          term1 <- log(Covariate[[i]])
        }
        
        
        term[[i]] <- data.frame(term1 = term1)
        names(term[[i]]) <- paste0("term1_",names(term)[i])
        rm(term1)
      }
      
      
      # FP degree == 2
      else
      {
        power <- c(FP[[i]]$power1,FP[[i]]$power2)
        if(is.matrix(power))
        {
          power <- power[1,]
        }
        
        if(!is.null(forced_powers))
        {
          power <- forced_powers
        }
        
        if (power[1] != 0) {
          term1 <- Covariate[[i]]^power[1]
        }
        if (power[1] == 0) {
          term1 <- log(Covariate[[i]])
        }
        if (power[2] != 0) {
          term2 <- Covariate[[i]]^power[2]
        }
        if (power[2] == 0) {
          term2 <- log(Covariate[[i]])
        }
        if (power[1] == power[2] & power[1] != 0) {
          term2 <- (Covariate[[i]]^power[2]) * log(Covariate[[i]])
        }
        if (power[1] == power[2] & power[1] == 0) {
          term2 <- term1^2
        }
        
        term[[i]] <- as.data.frame(cbind(term1,term2))
        names(term[[i]]) <- paste0(c("term1","term2"),"_",names(term)[i])
        
        rm(term1)
        rm(term2)
        
      }
    }
    
    # The regression model only allows integer ids. 
    data_TE_loop$study_mod <- as.integer(as.factor(data_TE_loop$study))
    data_TE_loop$study_id <- paste0(data_TE_loop$study_mod,"_",data_TE_loop$id)
    
    data_SE_loop$study_mod <- as.integer(as.factor(data_SE_loop$study))
    data_SE_loop$study_id <- paste0(data_SE_loop$study_mod,"_",data_SE_loop$id)
    
    
    #data$study_id <- gsub("-",":",data$study_id)
    
    # bind Fractional polynomial terms to data set
    data_TE_loop <- cbind(data_TE_loop,term$TE)
    data_SE_loop <- cbind(data_SE_loop,term$SE)
    
    # combine two data sets for use for the bivariate longitodinal model
    data <- merge(data_TE_loop,
                  data_SE_loop,
                  by = c("study","id","arm","visit_com","treat","study_mod","study_id"),
                  suffixes = c("_TE","_SE"))
    
    library(brms)
    library(glmmTMB)
    library(mgcv)
    
    source("scripts/MA_model_sel.R")
    source("scripts/formula_maker.R")
    
    # Create fomrulars for model
    terms <- paste0(grep("term",names(data_TE_loop),value = T),collapse = "+")
    
    # model selection
    model_TE <- MA_model_sel_func(data_loop = data,
                              Outcome = TE,
                              family = "gaus",
                              year = "year_TE",
                              term = terms)
    
    model_TE$lr_list <- formula_maker(lr_list = model_TE$lr_list)
   
    formula_TE <- model_TE$lr_list$call[1]
    formula_TE <- bf(formula = formula_TE,
                     family = gaussian())
    
    
    terms <- paste0(grep("term",names(data_SE_loop),value = T),collapse = "+")
    model_SE <- MA_model_sel_func(data_loop = data,
                                  Outcome = SE,
                                  family = "gaus",
                                  year = "year_SE",
                                  term = terms)
    
    model_SE$lr_list <- formula_maker(lr_list = model_SE$lr_list)
    
    formula_SE <- model_SE$lr_list$call[1]
    formula_SE <- bf(formula = formula_SE,
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
                    data = data)
    },
    error = function(e)
    {
      errors <<- e
      if(!file.exists("results//errors")){dir.create("results/errors")}
      save(errors,file = paste0("results/errors/",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      # Catches warnings
      warn <<- w
      
      if(!file.exists("results/warnings/")){dir.create("results/warnings/")}
      save(warn,file = paste0("results/warnings/",x,"_",comp,".Rdata"))
      
      model <<- brm(formula_SE + formula_TE + set_rescor(TRUE),
                    iter = 20000,
                    chains = 4,
                    cores = 4,
                    save_pars = save_pars(all = T),
                    control = list(adapt_delta = 0.95),
                    thin = 10,
                    data = data)
      
      warn <- warnings()
      save(warn,file = paste0("results/warnings/",x,"_",comp,".Rdata"))
    })
    
    
    save(model,file = paste0("results/models/",x,"_",comp,".Rdata"))
    
    # Computation of lambda
    res <- residuals(model, summary = F)
    res_SE <- res[,,gsub("_","",SE)]
    res_TE <- res[,,gsub("_","",TE)]
    
    prepare_func <- function(x. = x,
                             res_SE. = res_SE,
                             res_TE. = res_TE,
                             i. = i,
                             data. = data)
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
      time_points <- lapply(unique(data$visit_com),prepare_func)
      time_points <- list_rbind(time_points)
      
      time_points <- time_points %>%
        distinct(id,visit_com, .keep_all = T) %>%
        arrange(visit_com) %>%
        pivot_wider(names_from = visit_com,
                    values_from = c(SE,TE)) %>%
        dplyr::select(-id)
      
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
    
    save(missings_perc,file = paste0("results/missings/",x,"_",comp,".Rdata"))
    
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
    
    save(lambda,file = paste0("results/lambda/",x,"_",comp,".Rdata"))
    
    if(!file.exists("results/data/")){dir.create("results/data/")}
    save(data,file = paste0("results/data/",x,"_",comp,".Rdata"))
    
    rm(data_TE_loop,
       data_SE_loop,
       data)
    
  }
}

library(pbapply)
library(parallel)
cl <- makeCluster(detectCores()-4)

 pblapply(X = unique(data_TE$study),
                           FUN = bayes_model_func,
                           data_TE = data_TE,
                           data_SE = data_SE,
                           forced_powers = NULL,
                           TE = "edss_score",
                           SE = "t2_volume",
                           time_SE = "day_mri",
                           time_TE = "day_edss",
                           cl = cl)

stopCluster(cl)
