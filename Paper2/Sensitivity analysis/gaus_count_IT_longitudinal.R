# Purpose of this script: Sensitivity analyses using the IT methodology
# Longitudinal approach
# New/newly enlarged T2 volume (SEP)/ number of relapses (CEP)
# For details see publication and appendix
# 19.12.2023
# Author: Stefan Buchka


library(tidyverse)
library(rlist)
load("arr.RData")
load("mri.RData")
load("edss.RData")


data_TE <- arr %>%
  dplyr::select(study,id,arm,relapse,day) %>%
  arrange(study,arm,id,day) %>% 
  group_by(study,id) %>% 
  mutate(relapse = ifelse(sum(relapse != 0,na.rm = T) & #for some reason, there are lines of relapse persons, that contain NA. These are removed here
                            any(is.na(relapse)),1,relapse),
         day = ifelse(relapse == 1 & 
                        is.na(day),max(day,na.rm = T),day)) %>%
  filter(study %in% c("CAMMS223","CAMMS323","CAMMS324","CFTY720D2301","CFTY720D2309","WA21092","WA21093")) %>%
  distinct() %>%
  mutate(relapse = ifelse(is.na(relapse),0,relapse),
         relapse = ifelse(day >= 900,0,relapse)) %>% #exclude all relapses raised after two years (+ tolerance)
  mutate(ind = ifelse(day > 365,1,0)) %>%
  group_by(study,id,ind) %>%
  summarise(relapse = sum(relapse,na.rm = T),
            arm = unique(arm)) %>%
  ungroup()


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
  dplyr::select(-delete) %>%
  filter(study %in% c("CAMMS223","CAMMS323","CAMMS324","CFTY720D2301","CFTY720D2309","WA21092","WA21093")) %>%
  filter(visit_com %in% c(0,11,12))

dim(data_SE)
dim(data_SE[complete.cases(data_SE),])


data_TE <- data_TE[which(data_TE$study %in% data_SE$study),]

data <- merge(data_SE,
              data_TE,
              by = c("study","id","arm"),
              all.x = T,
              suffixes = c("_SE","_TE")) %>%
  group_by(study,id) %>%
  mutate(relapse = ifelse(is.na(relapse),0,relapse)) %>%
  filter((visit_com == 0 & ind == 0) |
           (visit_com %in% c(11,12) & ind == 1) |
           is.na(ind)) %>%
  arrange(study,id,visit_com) %>%
  mutate(relapse = cumsum(relapse),
         n = n()) %>% 
  as.data.frame()



bayes_model_func <- function(x = unique(data$study),
                             data = data,
                             forced_powers = NULL,
                             SE = "t2_volume",
                             TE = "relapse",
                             time = "day",
                             path = "results_gaus_count_IT_longitudinal")
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
    
    FP <-   fp_func(Covariate = "year",
                    Outcome = TE,
                    id = id,
                    Max.M = 2,
                    S = c(-2,-1,-0.5,0,0.5,1,2,3),
                    Dataset = as.data.frame(data_loop),
                    fam = "poisson")
    
    rm(id)
    
    Covariate <- data_loop$year
    Outcome <- eval(parse(text = paste0("data_loop$",TE)))
    
    #  FP degree == 1
    if(is.na(FP$power2) & is.null(forced_powers))
    {
      power <- FP$power1
      
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
      power <- c(FP$power1,FP$power2)
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
    
    
    
    # Model estimation
    source("scripts/IT_mode_sel.R") # model selection function
    source("scripts/calculation_r_2.R")
    
    # zero inflated model
    library(glmmTMB)
    
    models_zi <- IT_model_sel_func(data_loop = data_loop,
                                   family = "zero_inflated",
                                   TE = TE,
                                   SE = SE,
                                   term_TE = term_TE)
    
    tryCatch({
      model_0_zi <- eval(models_zi[[1]]$call)
    },
    error = function(e)
    {
      errors_zi <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_zi"))){dir.create(paste0(path,"/errors_zi"))}
      save(errors_zi,file = paste0(path,"/errors_zi/errors_0_zi_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_0_zi <<- eval(models_zi[[1]]$call)
      warn_zi <<- w
      
      if(!file.exists(paste0(path,"/warn_zi"))){dir.create(paste0(path,"/warn_zi"))}
      save(warn_zi,file = paste0(path,"/warn_zi/warn_0_zi_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_zi <- eval(models_zi[[2]]$call)
      
    },
    error = function(e)
    {
      errors_zi <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_zi"))){dir.create(paste0(path,"/errors_zi"))}
      save(errors_zi,file = paste0(path,"/errors_zi/errors_1_zi_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_zi <<- w
      
      model_1_zi <<- eval(models_zi[[2]]$call)
      
      if(!file.exists(paste0(path,"/warn_zi"))){dir.create(paste0(path,"/warn_zi"))}
      save(warn_zi,file = paste0(path,"/warn_zi/warn_1_zi_",x,"_",comp,".Rdata"))
    })
    
    if(!file.exists(paste0(path,"/models_zi"))){dir.create(paste0(path,"/models_zi"))}
    save(model_1_zi,file = paste0(path,"/models_zi/model_1_zi_",x,"_",comp,".Rdata"))
    save(model_0_zi,file = paste0(path,"/models_zi/model_0_zi_",x,"_",comp,".Rdata"))
    
    
    R2ht_zi <- r_2_function_IT(model0 = model_0_zi,
                               model1 = model_1_zi,
                               ordinal = F,
                               data = data_loop)
    
    if(!file.exists(paste0(path,"/R2ht_zi"))){dir.create(paste0(path,"/R2ht_zi"))}
    save(R2ht_zi,file = paste0(path,"/R2ht_zi/",x,"_",comp,".Rdata"))
    
    rm(R2ht_zi)
    
    
    
    
    # negative binomial model
    
    models_nb <- IT_model_sel_func(data_loop = data_loop,
                                   family = "nbinomial",
                                   TE = TE,
                                   SE = SE,
                                   term_TE = term_TE)
    
    tryCatch({
      model_0_nb <- eval(models_nb[[1]]$call)
    },
    error = function(e)
    {
      errors_nb <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_nb"))){dir.create(paste0(path,"/errors_nb"))}
      save(errors_nb,file = paste0(path,"/errors_nb/errors_0_nb_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_nb <<- w
      model_0_nb <<- eval(models_nb[[1]]$call)
      if(!file.exists(paste0(path,"/warn_nb"))){dir.create(paste0(path,"/warn_nb"))}
      save(warn_nb,file = paste0(path,"/warn_nb/warn_0_nb_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_nb <- eval(models_nb[[2]]$call)
    },
    error = function(e)
    {
      errors_nb <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_nb"))){dir.create(paste0(path,"/errors_nb"))}
      save(errors_nb,file = paste0(path,"/errors_nb/errors_1_nb_/",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_nb <<- w
      model_1_nb <<- eval(models_nb[[2]]$call)
      if(!file.exists(paste0(path,"/warn_nb"))){dir.create(paste0(path,"/warn_nb"))}
      save(warn_nb,file = paste0(path,"/warn_nb/warn_1_nb_",x,"_",comp,".Rdata"))
    })
    
    if(!file.exists(paste0(path,"/models_nb"))){dir.create(paste0(path,"/models_nb"))}
    save(model_1_nb,file = paste0(path,"/models_nb/model_1_nb_",x,"_",comp,".Rdata"))
    save(model_0_nb,file = paste0(path,"/models_nb/model_0_nb_",x,"_",comp,".Rdata"))
    
    
    R2ht_nb <- r_2_function_IT(model0 = model_0_nb,
                               model1 = model_1_nb,
                               ordinal = F,
                               data = data_loop)
    
    if(!file.exists(paste0(path,"/R2ht_nb"))){dir.create(paste0(path,"/R2ht_nb"))}
    save(R2ht_nb,file = paste0(path,"/R2ht_nb/",x,"_",comp,".Rdata"))
    
    rm(R2ht_nb)
    
    
    
    
    
    # poisson model
    
    models_po <- IT_model_sel_func(data_loop = data_loop,
                                   family = "poisson",
                                   TE = TE,
                                   SE = SE,
                                   term_TE = term_TE)
    
    tryCatch({
      model_0_po <- eval(models_po[[1]]$call)
    },
    error = function(e)
    {
      errors_po <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_po"))){dir.create(paste0(path,"/errors_po"))}
      save(errors_po,file = paste0(path,"/errors_po/errors_0_po_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_po <<- w
      
      model_0_po <<- eval(models_po[[1]]$call)
      
      if(!file.exists(paste0(path,"/warn_po"))){dir.create(paste0(path,"/warn_po"))}
      save(warn_po,file = paste0(path,"/warn_po/warn_0_po_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_po <- eval(models_po[[2]]$call)
    },
    error = function(e)
    {
      errors_po <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_po"))){dir.create(paste0(path,"/errors_po"))}
      save(errors_po,file = paste0(path,"/errors_po/errors_1_po_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_po <<- w
      
      model_1_po <<- eval(models_po[[2]]$call)
      
      if(!file.exists(paste0(path,"/warn_po"))){dir.create(paste0(path,"/warn_po"))}
      save(warn_po,file = paste0(path,"/warn_po/warn_1_po_",x,"_",comp,".Rdata"))
    })
    
    if(!file.exists(paste0(path,"/models_po"))){dir.create(paste0(path,"/models_po"))}
    save(model_1_po,file = paste0(path,"/models_po/model_1_po_",x,"_",comp,".Rdata"))
    save(model_0_po,file = paste0(path,"/models_po/model_0_po_",x,"_",comp,".Rdata"))
    
    
    R2ht_po <- r_2_function_IT(model0 = model_0_po,
                               model1 = model_1_po,
                               ordinal = F,
                               data = data_loop)
    
    if(!file.exists(paste0(path,"/R2ht_po"))){dir.create(paste0(path,"/R2ht_po"))}
    save(R2ht_po,file = paste0(path,"/R2ht_po/",x,"_",comp,".Rdata"))
    
    rm(R2ht_po)
    
    
    # Ordinal models
    library(ordinal)
    
    data_loop$rel_cat <- ifelse(data_loop$relapse >=2,2,data_loop$relapse)
    data_loop$rel_cat <- data_loop$rel_cat + 1
    
    #debugonce("IT_model_sel_func")
    models_or <- IT_model_sel_func(data_loop = data_loop,
                                   family = "ordinal",
                                   TE = "rel_cat",
                                   SE = SE,
                                   term_TE = term_TE)
    
    
    
    tryCatch({
      model_0_or <<- eval(models_or[[1]]$call)
    },
    error = function(e)
    {
      errors_or <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_or"))){dir.create(paste0(path,"/errors_or"))}
      save(errors_or,file = paste0(path,"/errors_or/errors_0_or_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_or <<- w
      
      model_0_or <<- eval(models_or[[1]]$call)
      
      if(!file.exists(paste0(path,"/warn_or"))){dir.create(paste0(path,"/warn_or"))}
      save(warn_or,file = paste0(path,"/warn_or/warn_0_or_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_or <<- eval(models_or[[2]]$call)
      
    },
    error = function(e)
    {
      errors_or <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_or"))){dir.create(paste0(path,"/errors_or"))}
      save(errors_or,file = paste0(path,"/errors_or/errors_1_or_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      warn_or <<- w
      
      model_1_or <<- eval(models_or[[2]]$call)
      
      if(!file.exists(paste0(path,"/warn_or"))){dir.create(paste0(path,"/warn_or"))}
      save(warn_or,file = paste0(path,"/warn_or/warn_1_or_",x,"_",comp,".Rdata"))
    })
    
    
    if(!file.exists(paste0(path,"/models_or"))){dir.create(paste0(path,"/models_or"))}
    save(model_1_or,file = paste0(path,"/models_or/model_1_or_",x,"_",comp,".Rdata"))
    save(model_0_or,file = paste0(path,"/models_or/model_0_or_",x,"_",comp,".Rdata"))
    
    
    R2ht_or <- r_2_function_IT(model0 = model_0_or,
                               model1 = model_1_or,
                               ordinal = T,
                               TE = TE,
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
         SE = "t2_volume",
         TE = "relapse",
         time = "day",
         path = "results_gaus_count_IT_longitudinal",
         cl = cl)

stopCluster(cl)

