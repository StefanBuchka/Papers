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
  select(study,id,arm,edss_score,day,visit,visit_com) %>%
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(edss$day,edss$day))-1)*-1) %>%
  group_by(id,study) %>%
  mutate(delete = ifelse(!any(visit_com == 0),1,0)) %>% # remove patients without baseline measurement
  filter(delete == 0) %>%
  select(-delete)

dim(data_TE)
dim(data_TE[complete.cases(data_TE),])

data_SE <- mri %>%
  select(study,id,arm,t2_volume,day,visit,visit_com,unit) %>%
  filter(!is.na(t2_volume)) %>%
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(mri$day,mri$day))-1)*-1) %>%
  mutate(t2_volume = ifelse(unit == "mm3",t2_volume/1000,t2_volume),
         t2_volume = log(t2_volume + 0.1)) %>%
  select(-unit) %>%
  group_by(id,study) %>%
  mutate(delete = ifelse(!any(visit_com == 0),1,0)) %>% # remove patients without baseline measurement
  filter(delete == 0) %>%
  select(-delete)

dim(data_SE)
dim(data_SE[complete.cases(data_SE),])

data_TE <- data_TE[which(data_TE$study %in% data_SE$study),]


# Individual Level Surrogacy
bayes_model_func <- function(x = unique(data_TE$study),
                             data_TE = data_TE,
                             data_SE = data_SE,
                             forced_powers = NULL,
                             TE = "edss_score",
                             SE = "t2_volume",
                             time_SE = "day_mri",
                             time_TE = "day_edss",
                             path = "results_gaus_gaus_IT")
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
  
  # To stabilise the Fractional polynomials, the following linear transformation will be done
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
  
  
  control_gr <- c("SC IFNB-1a","Placebo","Interferon beta-1a","PLACEBO","INTERFERON BETA-1a","AVONEX")
  arms <- as.character(unique(data_TE$arm))
  control_gr <- arms[which(arms %in% control_gr)]
  test_gr <- arms[-which(arms %in% control_gr)]
  
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
    
    data_loop <- merge(data_TE_loop,
                       data_SE_loop,
                       by = c("study","id","arm","visit_com","treat","study_mod","study_id"),
                       suffixes = c("_TE","_SE"))
    
    
    # Create fomrulars for model
    term_TE <- paste0(grep("term",names(data_TE_loop),value = T),collapse = "+")
    term_SE <- paste0(grep("term",names(data_SE_loop),value = T),collapse = "+")

    
    # Model estimation
    source("scripts/IT_mode_sel.R") # model selection function
    source("scripts/calculation_r_2.R")
    
    # zero inflated model
    library(glmmTMB)
    
    models <- IT_model_sel_func(data_loop = data_loop,
                                family = "gaus",
                                TE = TE,
                                SE = SE,
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
      mutate(edss_score = ifelse(edss_score >= 4,4,edss_score),
             year = year_TE)
  
    models_or <- IT_model_sel_func(data_loop = data_loop,
                                   family = "ordinal",
                                   TE = TE,
                                   SE = SE,
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
    
    rm(data_TE_loop,
       data_SE_loop,
       data_loop)
    
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
         path = "results_gaus_gaus_IT",
         cl = cl)

stopCluster(cl)

