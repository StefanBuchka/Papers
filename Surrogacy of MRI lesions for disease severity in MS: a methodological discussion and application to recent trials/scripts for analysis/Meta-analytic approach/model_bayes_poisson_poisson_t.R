# Purpose: Estimation of a Bayesian Gaussian-Gaussian joint model to derive R²_lambda 
# Note, the count variable is transformed
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka

# Note: The clinical endpoint (CEP) is named TE and the surrogate endpoint (SEP) is named SE in the following 
# Note: long run time with high computational power. Parallelisation will be used. All cores except 4 will be in use. 
# You can change this via cl variable of the "pblapply" function


library(tidyverse)
library(rlist)
load("datasets/arr.RData")
load("datasets/mri.RData")
load("datasets/edss.RData")

# Data about Relapses
data_TE <- arr %>%
  dplyr::select(study,id,arm,relapse,day) %>%
  arrange(study,arm,id,day) %>% 
  group_by(study,id) %>% 
  mutate(relapse = ifelse(sum(relapse != 0,na.rm = T) & #for some reason, there are lines of relapse persons, that contain NA. These are removed here
                            any(is.na(relapse)),1,relapse),
         day = ifelse(relapse == 1 & 
                        is.na(day),max(day,na.rm = T),day)) %>%
  distinct()

dim(data_TE)
dim(data_TE[complete.cases(data_TE),])

# Data about MRI
data_SE <- mri %>%
  dplyr::select(study,id,arm,t2_new_or_enlarged,day,visit,visit_com) %>%
  filter(!is.na(t2_new_or_enlarged)) %>% 
  filter(visit != 800) %>% #visit 800 is accumulated t2 lesion count
  group_by(study,id,day) %>%
  mutate(n = n()) %>%
  filter(!(n == 2 & visit > 50)) %>% 
  filter(visit_com != 0) %>% # remove baseline count (since this is not the change in t2 lesions count)
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(mri$day,mri$day))-1)*-1)

dim(data_SE)
dim(data_SE[complete.cases(data_SE),])

data_TE <- data_TE[which(data_TE$study %in% data_SE$study),]


# Data for interval approach. The visits were assigned to the difference of MRI measurement of two time points.

# The baseline t2 count were removed, since the base line count are the numbers of t2 lesions at base line. This is the the new or enlarged lesion count. 
# In the interval approach, lesions/relapses could occur between baseline and first measurement. Because of this, the baseline time point is included without 
# measurement of the lesion count. 
data_SE <- data_SE %>%
  group_by(study,id,arm) %>%
  group_modify(~ add_row(t2_new_or_enlarged = NA,day = 0,visit = NA, visit_com = 0, .x)) 


# Data for interval approach. 
# Relapses occur randomly. So no "visit" is available.
# The number of relapses, that occurred between two visits (T2 lesion count was measured), were assigned to the later visit of that two visits. 
data_delta <- merge(data_SE,
                    data_TE,
                    by = c("study","id","arm"),
                    suffixes = c("_SE","_TE"),all = T) %>% # merge two data sets
  filter(!is.na(day_SE)) %>%
  arrange(study,id,day_TE,day_SE) %>%
  mutate(rel_in = day_TE) %>%
  group_by(study,id) %>%
  mutate(day_TE = ifelse(day_TE > max(day_SE,na.rm = T) &
                           !is.na(relapse), # if the relapse is after the last visit, adjustments have to be done
                         max(day_TE,na.rm = T),day_TE),
         day_SE = ifelse(day_TE > max(c(day_SE,day_TE),na.rm = T) & # if the relapse is after the last visit, adjustments have to be done
                           !is.na(relapse),
                         max(day_SE,na.rm = T),
                         day_SE)) %>%
  group_by(study,id,rel_in) %>%
  mutate(day_SE = ifelse(day_SE == min(day_SE,na.rm = T) &
                           !is.na(relapse) &
                           visit_com == 0,0,day_SE), #replace the first day of visit by 0, to include relapses happens before the first MRI visit in the interval.
         day_SE_end = lead(day_SE,
                           default = max(c(day_SE,day_TE),na.rm = T)), # Add the last visit or relapse date (that one, which is later).
         t2_new_or_enlarged = lead(t2_new_or_enlarged)) %>% #assign the t2_lesion count to right time point  
  filter(day_SE != day_SE_end) %>%
  mutate(rel = ifelse(relapse == 1 &
                        (day_TE > day_SE & day_TE <= day_SE_end),1,0)) %>% # build the intervals
  group_by(study,id,day_SE) %>%
  mutate(rel = sum(rel)) %>%
  distinct(study,id,arm,t2_new_or_enlarged,day_SE,day_SE_end,visit,visit_com,rel, .keep_all = F) %>% # excludes duplicates
  mutate(rel = ifelse(is.na(rel),0,rel), # mark all patients without relapse with 0. 
         day_SE = day_SE + 1) %>% # add to all days one day to avoid 0
  group_by(study,id) %>%
  group_split()


# Excludes rows of observations, that only have one interval
func <- function(x)
{
  if(length(unique(x$day_SE)) == 1)
  {
    x <- x[nrow(x),]
  }  
  return(x)}

data_delta <- lapply(data_delta,func)
data_delta <- list.rbind(data_delta)

# Transformation of poisson variable

source("scripts/trans_func.R")

data_delta$t2_new_or_enlarged <- trans_func(var = data_delta$t2_new_or_enlarged,
                                            trans = "Huber",
                                            bins = NULL)

data_delta$rel <- trans_func(var = data_delta$rel,
                             trans = "Anscombe",
                             bins = NULL)


bayes_model_func <- function(x = unique(data_delta$study),
                             data = data_delta,
                             forced_powers = NULL,
                             SE = "t2_new_or_enlarged",
                             TE = "rel",
                             time = "day_SE",
                             path = "results_poisson_poisson_t/delta")
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
  data$year <- unlist(data[,time]/365.25)
  
  
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
    
    
    # create variable on observation level
    data_loop$obs <- 1:nrow(data_loop)
    
    data_loop <- data_loop %>%
      filter(!is.na(t2_new_or_enlarged))
    
    
    
    # Calculate Fractional polynomials
    source("scripts/fp_function.R")
    
    id <- data_loop$id
    Covariate <- list(mode = "list", length = 2)
    Outcome <- list(mode = "list", length = 2)
    
    FP_TE <-   fp_func(Covariate = "year",
                       Outcome = TE,
                       id = id,
                       Max.M = 2,
                       S = c(-2,-1,-0.5,0,0.5,1,2,3),
                       Dataset = as.data.frame(data_loop),
                       fam = "gaussian")
    
    Covariate[[1]] <- data_loop$year
    Outcome[[1]] <- eval(parse(text = paste0("data_loop$",TE)))
    
    FP_SE <-   fp_func(Covariate = "year",
                       Outcome = SE,
                       id = id,
                       Max.M = 2,
                       S = c(-2,-1,-0.5,0,0.5,1,2,3),
                       Dataset = as.data.frame(data_loop),
                       fam = "gaussian")
    
    rm(id)
    
    Covariate[[2]] <- data_loop$year
    Outcome[[2]] <- eval(parse(text = paste0("data_loop$",SE)))
    
    
    FP <- list(FP_TE,FP_SE)
    
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
    
    # bind Fractional polynomial terms to data set
    data_loop <- cbind(data_loop,term$TE)
    data_loop <- cbind(data_loop,term$SE)
    
    
    library(brms)
    library(glmmTMB)
    source("scripts/MA_model_sel.R")
    source("scripts/formula_maker.R")
    # Create poisson fomrulas for model
    
    
    terms <- paste0(grep("_TE",names(data_loop),value = T),collapse = "+")
    
    model_TE <- MA_model_sel_func(data_loop = data_loop,
                                  Outcome = TE,
                                  term = terms,
                                  family = "gaus")
    
    model_TE$lr_list <- formula_maker(lr_list = model_TE$lr_list)
    
    formula_TE <- model_TE$lr_list$call[1]
    
    formula_TE <- bf(formula = formula_TE,
                     family = gaussian())
    
    
    
    terms <- paste0(grep("term1_SE|term2_SE",names(data_loop),value = T),collapse = "+")
    
    model_SE <- MA_model_sel_func(data_loop = data_loop,
                                  Outcome = SE,
                                  term = terms,
                                  family = "gaus")
    
    model_SE$lr_list <- formula_maker(lr_list = model_SE$lr_list)
    
    formula_SE <- model_SE$lr_list$call[1]
    
    formula_SE <- bf(formula = formula_SE,
                     family = gaussian())
    
    
    
    # Estimation of the model
    tryCatch({
      model <<- brm(formula_SE + formula_TE + set_rescor(T),
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
      save(errors,file = paste0(path,"/errors/",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      # Catches warnings
      warn <<- w
      
      if(!file.exists(paste0(path,"/warnings"))){dir.create(paste0(path,"/warnings"))}
      save(warn,file = paste0(path,"/warnings/",x,"_",comp,".Rdata"))
      

      model <<- brm(formula_SE + formula_TE + set_rescor(T),
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
    res_SE <- res[,,gsub("_","",SE)]
    res_TE <- res[,,gsub("_","",TE)]
    
    # res <- VarCorr(model, summary = F)$obs$cor
    # res_SE <- res[,,"t2volume_Intercept"]
    # res_TE <- res[,,"rel_Intercept"]
    
    
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

pblapply(X = unique(data_delta$study),
         FUN = bayes_model_func,
         data = data_delta,
         forced_powers = NULL,
         SE = "t2_new_or_enlarged",
         TE = "rel",
         time = "day_SE",
         path = "results_poisson_poisson_t/delta",
         cl = cl)

stopCluster(cl)











