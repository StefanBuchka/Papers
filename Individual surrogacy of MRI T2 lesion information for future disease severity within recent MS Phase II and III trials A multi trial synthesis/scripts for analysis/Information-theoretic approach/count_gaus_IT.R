# Purpose: Estimation of LRF using the IT methodology
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka

# Note: The clinical endpoint (CEP) is named TE and the surrogate endpoint (SEP) is named SE in the following 


library(tidyverse)
 load("arr.RData")
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
  dplyr::select(study,id,arm,t2_new_or_enlarged,day,visit,visit_com) %>%
  filter(!is.na(t2_new_or_enlarged)) %>% 
  filter(visit != 800) %>% #visit 800 is cumulated t2 lesion count
  group_by(study,id,day) %>%
  mutate(n = n()) %>%
  filter(!(n == 2 & visit > 50)) %>% 
  filter(visit_com != 0) %>% # remove baseline count (since this is not the change in t2 lesions count)
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(mri$day,mri$day))-1)*-1)

dim(data_SE)
dim(data_SE[complete.cases(data_SE),])

data_TE <- data_TE[which(data_TE$study %in% data_SE$study),]

data <- merge(data_TE,
              data_SE,
              by = c("id","study","arm","visit","visit_com"),
              suffixes = c("_TE","_SE"))

dim(data_TE)
dim(data_SE)
dim(data)

bayes_model_func <- function(x = unique(data$study),
                             data = data,
                             forced_powers = NULL,
                             poisson = "t2_new_or_enlarged",
                             gaus = "edss_score",
                             time = "day_SE",
                             path = "results_count_gaus_IT")
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
    
    
    # The regression model only allows integer ids. 
    data_loop$study_mod <- as.integer(as.factor(data_loop$study))
    data_loop$study_id <- paste0(data_loop$study_mod,"_",data_loop$id)
    
    
    # Calculate Fractional polynomials
    source("scripts/fp_function.R")
    
    id <- data_loop$id
    
    FP_gaus <-   fp_func(Covariate = "year",
                         Outcome = gaus,
                         id = id,
                         Max.M = 2,
                         S = c(-2,-1,-0.5,0,0.5,1,2,3),
                         Dataset = as.data.frame(data_loop))
    
    rm(id)
    
    Covariate <- data_loop$year
    Outcome <- eval(parse(text = paste0("data_loop$",gaus)))
    
    #  FP degree == 1
    if(is.na(FP_gaus$power2) & is.null(forced_powers))
    {
      power <- FP_gaus$power1
      
      if(length(power) > 1)
      {
        power <- power[1]
      }
      
      if (power[1] != 0) {
        term1 <- Covariate^power[1]
      }
      if (power[1] == 0) {
        term1 <- log(Covariate)
      }
      
      
      term <- data.frame(term1 = term1)
      rm(term1)
    }
    
    
    # FP degree == 2
    else
    {
      power <- c(FP_gaus$power1,FP_gaus$power2)
      if(is.matrix(power))
      {
        power <- power[1,]
      }
      
      if(!is.null(forced_powers))
      {
        power <- forced_powers
      }
      
      if (power[1] != 0) {
        term1 <- Covariate^power[1]
      }
      if (power[1] == 0) {
        term1 <- log(Covariate)
      }
      if (power[2] != 0) {
        term2 <- Covariate^power[2]
      }
      if (power[2] == 0) {
        term2 <- log(Covariate)
      }
      if (power[1] == power[2] & power[1] != 0) {
        term2 <- (Covariate^power[2]) * log(Covariate)
      }
      if (power[1] == power[2] & power[1] == 0) {
        term2 <- term1^2
      }
      
      term <- as.data.frame(cbind(term1,term2))
      
      rm(term1)
      rm(term2)
      
    }
    
    
    
    
    #data$study_id <- gsub("-",":",data$study_id)
    
    # bind Fractional polynomial terms to data set
    data_loop <- cbind(data_loop,term)
    
    data_loop$obs <- 1:nrow(data_loop)
    
    term_TE <- paste0(grep("term",names(data_loop),value = T),collapse = "+")
    
    # Estimation of the model
    # Model estimation
    source("scripts/IT_mode_sel.R") # model selection function
    source("scripts/calculation_r_2.R")
    
    # zero inflated model
    library(glmmTMB)
    
    models <- IT_model_sel_func(data_loop = data_loop,
                                family = "gaus",
                                TE = gaus,
                                SE = poisson,
                                term_TE = term_TE)
    
    tryCatch({
      model_0 <- eval(models[[1]]$call)
    },
    error = function(e)
    {
      errors <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors"))){dir.create(paste0(path,"/errors"))}
      save(errors,file = paste0(path,"/errors/errors_0_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_0 <<- eval(models[[1]]$call)
      warn <<- w
      
      if(!file.exists(paste0(path,"/warn"))){dir.create(paste0(path,"/warn"))}
      save(warn,file = paste0(path,"/warn/warn_0_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1 <- eval(models[[2]]$call)
      
    },
    error = function(e)
    {
      errors <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors"))){dir.create(paste0(path,"/errors"))}
      save(errors,file = paste0(path,"/errors/errors_1_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn <<- w
      
      model_1 <<- eval(models[[2]]$call)
      
      if(!file.exists(paste0(path,"/warn"))){dir.create(paste0(path,"/warn"))}
      save(warn,file = paste0(path,"/warn/warn_1_",x,"_",comp,".Rdata"))
    })
    
    if(!file.exists(paste0(path,"/models"))){dir.create(paste0(path,"/models"))}
    save(model_1,file = paste0(path,"/models/model_1_",x,"_",comp,".Rdata"))
    save(model_0,file = paste0(path,"/models/model_0_",x,"_",comp,".Rdata"))
    
    
    R2ht <- r_2_function_IT(model0 = model_0,
                            model1 = model_1,
                            ordinal = F,
                            data = data_loop)
    
    if(!file.exists(paste0(path,"/R2ht"))){dir.create(paste0(path,"/R2ht"))}
    save(R2ht,file = paste0(path,"/R2ht/",x,"_",comp,".Rdata"))
    
    rm(R2ht)
    
    
    
    
    
    # ordinal TE
    
    
    library(ordinal)
    data_loop <- data_loop %>%
      mutate(edss_score = ifelse(edss_score >= 4,4,edss_score))
    
    models_or <- IT_model_sel_func(data_loop = data_loop,
                                   family = "ordinal",
                                   TE = gaus,
                                   SE = poisson,
                                   term_TE = term_TE)
    
    tryCatch({
      model_0_or <<- eval(models_or[[1]]$call)
    },
    error = function(e)
    {
      errors_or <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_or"))){dir.create(paste0(path,"/errors_or"))}
      save(errors_or,file = paste0(path,"/errors_or/errors_0_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_0_or <<- eval(models_or[[1]]$call)
      warn_or <<- w
      
      if(!file.exists(paste0(path,"/warn_or"))){dir.create(paste0(path,"/warn_or"))}
      save(warn_or,file = paste0(path,"/warn_or/warn_0_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_or <<- eval(models_or[[2]]$call)
      
    },
    error = function(e)
    {
      errors_or <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_or"))){dir.create(paste0(path,"/errors_or"))}
      save(errors_or,file = paste0(path,"/errors_or/errors_1_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_or <<- w
      
      model_1_or <<- eval(models_or[[2]]$call)
      
      if(!file.exists(paste0(path,"/warn_or"))){dir.create(paste0(path,"/warn_or"))}
      save(warn_or,file = paste0(path,"/warn_or/warn_1_",x,"_",comp,".Rdata"))
    })
    
    if(!file.exists(paste0(path,"/models_or"))){dir.create(paste0(path,"/models_or"))}
    save(model_1_or,file = paste0(path,"/models_or/model_1_",x,"_",comp,".Rdata"))
    save(model_0_or,file = paste0(path,"/models_or/model_0_",x,"_",comp,".Rdata"))
    
    
    R2ht_or <- r_2_function_IT(model0 = model_0_or,
                               model1 = model_1_or,
                               ordinal = F,
                               data = data_loop)
    
    if(!file.exists(paste0(path,"/R2ht_or"))){dir.create(paste0(path,"/R2ht_or"))}
    save(R2ht_or,file = paste0(path,"/R2ht_or/",x,"_",comp,".Rdata"))
    
    rm(R2ht_or)
    
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
         time = "day_SE",
         path = "results_count_gaus_IT",
         cl = cl)

stopCluster(cl)

cl <- makeCluster(detectCores()-4)

