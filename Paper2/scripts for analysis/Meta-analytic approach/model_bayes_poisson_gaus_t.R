# Purpose: Estimation of a Bayesian Gaussian-Gaussian joint model to derive R²_lambda 
# Note, the count variable is transformed
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka

# Note: The clinical endpoint (CEP) is named TE and the surrogate endpoint (SEP) is named SE in the following 
# Note: long run time with high computational power. Parallelisation will be used. All cores except 4 will be in use. 
# You can change this via cl variable of the "pblapply" function



library(tidyverse)
load("datasets/arr.RData")
load("datasets/mri.RData")
load("datasets/edss.RData")

# Data about EDSS
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

# Data about MRI
data_SE <- mri %>%
  dplyr::select(study,id,arm,t2_new_or_enlarged,day,visit,visit_com) %>%
  group_by(study,id) %>%
  filter(any(!is.na(t2_new_or_enlarged))) %>% #exclude studies, which did not measure new/enlarged T2 lesions
  ungroup() %>%
  mutate(t2_new_or_enlarged = ifelse(visit_com == 0,0,t2_new_or_enlarged)) %>% # replace baseline value by 0, since no new/enlarged lesions can be count at baseline
  filter(!is.na(t2_new_or_enlarged)) %>% 
  filter(visit != 800) %>% #visit 800 is accumulated t2 lesion count
  group_by(study,id,day) %>%
  mutate(n = n()) %>%
  filter(!(n == 2 & visit > 50)) %>% # filters unplanned visits
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(mri$day,mri$day))-1)*-1) %>%
  filter(!is.na(day))

dim(data_SE)
dim(data_SE[complete.cases(data_SE),])

data_TE <- data_TE[which(data_TE$study %in% data_SE$study),]

data <- merge(data_TE,
              data_SE,
              by = c("id","study","arm","visit_com"),
              suffixes = c("_TE","_SE"),
              all.x = T) %>%
  filter(!is.na(t2_new_or_enlarged) |
           visit_com == 0) %>% #remove all missing t2_lesion counts (except baseline visit. This is necessary to compute edss_delta from baseline)
  group_by(study,id) %>%
  mutate(edss_delta = edss_score - lag(edss_score)) %>% # replace baseline value by 0, since no new/enlarged lesions can be count at baseline
  filter(visit_com != 0) %>% # remove baseline count (since this is not the change in t2 lesions count)
  filter(!is.na(day_SE)) %>%
  arrange(study,id,arm)


dim(data_TE)
dim(data_SE)
dim(data)

# Transformation of poisson variable
source("scripts/trans_func.R")

data$t2_new_or_enlarged <- trans_func(var = data$t2_new_or_enlarged,
                                      trans = "Huber",
                                      bins = NULL)


bayes_model_func <- function(x = unique(data$study),
                             data = data,
                             forced_powers = NULL,
                             poisson = "t2_new_or_enlarged",
                             gaus = "edss_score",
                             time_TE = "day_TE",
                             time_SE = "day_SE",
                             path = "results_poisson_gaus_t")
{
  library(lme4)
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(nlme)
  library(ggcorrplot)
  library(rlist)
  
  # select study within function
  data <- data[data$study == x,]
  
  # transform day variable to year
  data$year_TE <- unlist(data[,time_TE]/365.25)
  data$year_SE <- unlist(data[,time_SE]/365.25)
  
  # To stabilise the Fractional polynomials, the following linear transformation will be done
  if(any(data$year_TE < 0.1))
  {
    data$year_TE <- data$year_TE + 0.1
  }
  
  if(any(data$year_SE < 0.1))
  {
    data$year_SE <- data$year_SE + 0.1
  }
  
  
  # figure out comparisons between control and intervention group
  control_gr <- c("SC IFNB-1a","Placebo","Interferon beta-1a","PLACEBO","INTERFERON BETA-1a","AVONEX")
  arms <- as.character(unique(data$arm))
  control_gr <- arms[which(arms %in% control_gr)]
  test_gr <- arms[-which(arms %in% control_gr)]
  
  if("" %in% test_gr)
  {
    test_gr <- test_gr[-which(test_gr == "")]
  }
  
  comparisons <- crossing(control_gr,test_gr)
  
  
  control_gr <- c("SC IFNB-1a","Placebo","Interferon beta-1a","PLACEBO","INTERFERON BETA-1a","AVONEX")
  arms <- as.character(unique(data$arm))
  control_gr <- arms[which(arms %in% control_gr)]
  test_gr <- arms[-which(arms %in% control_gr)]
  
  # The analysis have to be carried out for each active drug vs. control group arm
  for(comp in 1:nrow(comparisons))
  {
    # prepare data 
    data_loop <- data %>%
      filter(arm %in% c(comparisons$control_gr[comp],comparisons$test_gr[comp]))
    
    data_loop$treat <- as.factor(ifelse(data_loop$arm %in% control_gr,-1,1))
    
    
    # The regression model only allows integer ids. 
    data_loop$study_mod <- as.integer(as.factor(data_loop$study))
    data_loop$study_id <- paste0(data_loop$study_mod,"_",data_loop$id)
    
    
    # Calculate Fractional polynomials
    source("scripts/fp_function.R")
    
    id <- data_loop$id
    Covariate <- list(mode = "list", length = 2)
    Outcome <- list(mode = "list", length = 2)
    
    FP_gaus <-   fp_func(Covariate = "year_TE",
                         Outcome = gaus,
                         id = id,
                         Max.M = 2,
                         S = c(-2,-1,-0.5,0,0.5,1,2,3),
                         Dataset = as.data.frame(data_loop),
                         fam = "gaussian")
    
    Covariate[[1]] <- data_loop$year_TE
    Outcome[[1]] <- eval(parse(text = paste0("data_loop$",gaus)))
    
    FP_poisson <-   fp_func(Covariate = "year_SE",
                            Outcome = poisson,
                            id = id,
                            Max.M = 2,
                            S = c(-2,-1,-0.5,0,0.5,1,2,3),
                            Dataset = as.data.frame(data_loop),
                            fam = "gaussian")
    
    rm(id)
    
    Covariate[[2]] <- data_loop$year_SE
    Outcome[[2]] <- eval(parse(text = paste0("data_loop$",poisson)))
    
    
    FP <- list(FP_gaus,FP_poisson)
    
    term <- vector(mode = "list", length = 2)
    names(term) <- c("gaus","poisson")
    
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
    
    
    
    
    # bind Fractional polynomial terms to data set
    data_loop <- cbind(data_loop,term$gaus)
    data_loop <- cbind(data_loop,term$poisson)
    
    data_loop$obs <- 1:nrow(data_loop)
    
    library(brms)
    library(glmmTMB)
    source("scripts/MA_model_sel.R")
    source("scripts/formula_maker.R")
    
    # Creation poisson formulars for model
    
    terms <- paste0(grep("poisson",names(data_loop),value = T),collapse = "+")
    
    model_poisson <- MA_model_sel_func(data_loop = data_loop,
                                       Outcome = poisson,
                                       term = terms,
                                       family = "gaus")
    
    model_poisson$lr_list <- formula_maker(lr_list = model_poisson$lr_list)
    
    formula_poisson <- model_poisson$lr_list$call[1]
    
    formula_poisson <- bf(formula = formula_poisson,
                          family = gaussian())
    
    
    terms <- paste0(grep("gaus",names(data_loop),value = T),collapse = "+")
    
    # model selection
    model_gaus <- MA_model_sel_func(data_loop = data_loop,
                                    Outcome = gaus,
                                    family = "gaus",
                                    term = terms)
    
    model_gaus$lr_list <- formula_maker(lr_list = model_gaus$lr_list)
    
    formula_gaus <- model_gaus$lr_list$call[1]
    
    formula_gaus <- bf(formula = formula_gaus,
                       family = gaussian())
    
    
    # Estimation of the model
    tryCatch({
      model <<- brm(formula_poisson + formula_gaus + set_rescor(T),
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
      save(errors,file = paste0(path,"/errors/",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      # Catches warnings
      warn <<- w
      
      if(!file.exists(paste0(path,"/warnings"))){dir.create(paste0(path,"/warnings"))}
      save(warn,file = paste0(path,"/warnings/",x,"_",comp,".Rdata"))
      
      
      
      model <<- brm(formula_poisson + formula_gaus + set_rescor(T),
                    iter = 20000,
                    chains = 4,
                    cores = 4,
                    save_pars = save_pars(all = T),
                    control = list(adapt_delta = 0.99),
                    thin = 10,
                    data = data_loop)
      
      
      warn <<- warnings()
      save(warn,file = paste0(path,"/warnings/",x,"_",comp,".Rdata"))
    })
    
    
    
    if(!file.exists(paste0(path,"/models"))){dir.create(paste0(path,"/models"))}
    save(model,file = paste0(path,"/models/",x,"_",comp,".Rdata"))
    
    # Computation of lambda
    res <- residuals(model, summary = F)
    res_SE <- res[,,gsub("_","",poisson)]
    res_TE <- res[,,gsub("_","",gaus)]
    
    
    
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
    
    
    if(!file.exists(paste0(path,"/missings"))){dir.create(paste0(path,"/missings"))}
    save(missings_perc,file = paste0(path,"/missings/",x,"_",comp,".Rdata"))
    
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
    save(lambda,file = paste0(path,"/lambda/",x,"_",comp,".Rdata"))
    
    if(!file.exists(paste0(path,"/data"))){dir.create(paste0(path,"/data"))}
    save(data_loop,file = paste0(path,"/data/",x,"_",comp,".Rdata"))
    
    rm(data_loop)
    
  }
  
}

library(pbapply)
library(parallel)


cl <- makeCluster(detectCores()-4)

pblapply(X = unique(data_TE$study),
         FUN = bayes_model_func,
         data = data,
         forced_powers = NULL,
         poisson = "t2_new_or_enlarged",
         gaus = "edss_score",
         time_TE = "day_TE",
         time_SE = "day_SE",
         path = "results_poisson_gaus_t",
         cl = cl)

stopCluster(cl)




