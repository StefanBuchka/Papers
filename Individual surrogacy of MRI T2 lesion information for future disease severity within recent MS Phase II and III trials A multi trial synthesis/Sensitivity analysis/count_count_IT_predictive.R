
# Purpose of this script: Sensitivity analyses using the IT methodology
# Predictive approach
# New/newly enlarged T2 lesions (SEP)/ number of relapses (CEP)
# For details see publication and appendix
# 19.12.2023
# Author: Stefan Buchka


library(tidyverse)
library(rlist)
library(MASS)

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
                        is.na(day),max(day,na.rm = T),day),
         indi = ifelse(day >= 365,T,F),
         indi = ifelse(is.na(indi),F,indi)) %>%
  distinct() %>%
  group_by(study,id) %>%
  mutate(relapse = ifelse(is.na(relapse) | indi == F,0,relapse)) %>%
  summarise(study = unique(study),
            id = unique(id),
            arm = unique(arm),
            relapse = sum(relapse, na.rm = T)) %>%
  group_by(study,id) %>%
  summarise(arm = unique(arm),
            relapse = sum(relapse)) %>%
  filter(study %in% c("CAMMS323","CAMMS324","CFTY720D2301","CFTY720D2309","WA21092","WA21093")) %>%
  as.data.frame()

dim(data_TE)
dim(data_TE[complete.cases(data_TE),])

data_SE <- mri %>%
  filter(study %in% c("CAMMS323","CAMMS324","CFTY720D2301","CFTY720D2309","WA21092","WA21093")) %>%
  filter(visit_com %in% c(6,11,12)) %>%
  dplyr::select(study,id,arm,t2_new_or_enlarged,day,visit,visit_com) %>%
  filter(!is.na(t2_new_or_enlarged)) %>% 
  group_by(study,id,day) %>%
  mutate(n = n()) %>%
  filter(!(n == 2 & visit > 50)) %>% 
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(mri$day,mri$day))-1)*-1) %>%
  group_by(study,id) %>%
  mutate(n_obs = n()) %>%
  group_by(study) %>%
  mutate(n_visit = length(unique(visit_com))) %>%
  group_by(study,id) %>%
  summarise(arm = unique(arm),
            t2_new_or_enlarged = sum(t2_new_or_enlarged),
            n_obs = unique(n_obs),
            n_visit = unique(n_visit)) %>%
  as.data.frame()


dim(data_SE)
dim(data_SE %>% filter(n_obs == n_visit))
data_SE <- data_SE %>% filter(n_obs == n_visit)
dim(data_SE[complete.cases(data_SE),])


data_TE <- data_TE[which(data_TE$study %in% data_SE$study),]

data <- merge(data_SE,
              data_TE,
              by = c("study","id","arm"),
              all.x = T) %>%
  mutate(relapse = ifelse(is.na(relapse),0,relapse),
         t2_new_or_enlarged = ifelse(t2_new_or_enlarged >=8,8,t2_new_or_enlarged))










bayes_model_func <- function(x = unique(data$study),
                             data = data,
                             SE = "t2_new_or_enlarged",
                             TE = "relapse",
                             path = "results_count_count_IT_predictive")
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
    
    
    # zero inflated model
    library(glmmTMB)
    
    formula_0 <- paste0(TE," ~ treat")
    formula_1 <- paste0(TE," ~ treat + ",SE)
    
    tryCatch({
      model_0_zi <- glmmTMB(formula = as.formula(formula_0),
                            ziformula = ~ 1,
                            family = "poisson",
                            data = data_loop)
    },
    error = function(e)
    {
      errors_zi <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_zi"))){dir.create(paste0(path,"/errors_zi"))}
      save(errors_zi,file = paste0(path,"/errors_zi/errors_0_zi_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_0_zi <<- glmmTMB(formula = as.formula(formula_0),
                             ziformula = ~ 1,
                             family = "poisson",
                             data = data_loop)
      warn_zi <<- w
      
      if(!file.exists(paste0(path,"/warn_zi"))){dir.create(paste0(path,"/warn_zi"))}
      save(warn_zi,file = paste0(path,"/warn_zi/warn_0_zi_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_zi <- glmmTMB(formula = as.formula(formula_1),
                            ziformula = ~ 1,
                            family = "poisson",
                            data = data_loop)
      
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
      
      model_1_zi <<- glmmTMB(formula = as.formula(formula_1),
                             ziformula = ~ 1,
                             family = "poisson",
                             data = data_loop)
      
      
      if(!file.exists(paste0(path,"/warn_zi"))){dir.create(paste0(path,"/warn_zi"))}
      save(warn_zi,file = paste0(path,"/warn_zi/warn_1_zi_",x,"_",comp,".Rdata"))
    })
    
    if(!file.exists(paste0(path,"/models_zi"))){dir.create(paste0(path,"/models_zi"))}
    save(model_1_zi,file = paste0(path,"/models_zi/model_1_zi_",x,"_",comp,".Rdata"))
    save(model_0_zi,file = paste0(path,"/models_zi/model_0_zi_",x,"_",comp,".Rdata"))
    
    source("scripts/calculation_r_2.R")
    R2ht_zi <- r_2_function_IT(model0 = model_0_zi,
                               model1 = model_1_zi,
                               ordinal = F,
                               data = data_loop)
    
    if(!file.exists(paste0(path,"/R2ht_zi"))){dir.create(paste0(path,"/R2ht_zi"))}
    save(R2ht_zi,file = paste0(path,"/R2ht_zi/",x,"_",comp,".Rdata"))
    
    rm(R2ht_zi)
    
    
    
    
    # negative binomial model
    
    tryCatch({
      model_0_nb <- glmmTMB(as.formula(formula_0), 
                            family = "nbinom2",
                            data = data_loop)
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
      
      model_0_nb <<- glmmTMB(as.formula(formula_0), 
                             family = "nbinom2",
                             data = data_loop)
      
      if(!file.exists(paste0(path,"/warn_nb"))){dir.create(paste0(path,"/warn_nb"))}
      save(warn_nb,file = paste0(path,"/warn_nb/warn_0_nb_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_nb <- glmmTMB(as.formula(formula_1), 
                            family = "nbinom2",
                            data = data_loop)
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
      
      model_1_nb <<- glmmTMB(as.formula(formula_1), 
                             family = "nbinom2",
                             data = data_loop)
      
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
    
    tryCatch({
      model_0_po <- glmmTMB(as.formula(formula_0), 
                            family = "poisson",
                            data = data_loop)
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
      
      model_0_po <<- glmmTMB(as.formula(formula_0), 
                             family = "poisson",
                             data = data_loop)
      
      if(!file.exists(paste0(path,"/warn_po"))){dir.create(paste0(path,"/warn_po"))}
      save(warn_po,file = paste0(path,"/warn_po/warn_0_po_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_po <- glmmTMB(as.formula(formula_1), 
                            family = "poisson",
                            data = data_loop)
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
      
      model_1_po <<- glmmTMB(as.formula(formula_1), 
                             family = "poisson",
                             data = data_loop)
      
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
    
    
    data_loop$t2_cat <- ifelse(data_loop$t2_new_or_enlarged >=8,8,data_loop$t2_new_or_enlarged)
    
    data_loop$rel <- data_loop$rel + 1
    data_loop$t2_cat <- data_loop$t2_cat + 1
    
    
    formula_0 <- (paste0("factor(",TE,") ~ treat"))
    text_0 <- paste0("clmm2(",formula_0,",Hess = T,data = data_loop)")
    
    formula_1 <- (paste0("factor(",TE,") ~ as.numeric(", SE, ") + treat"))
    text_1 <- paste0("clmm2(",formula_1,",Hess = T,data = data_loop)")
    
    tryCatch({
      model_0_or <<- eval(parse(text = text_0))
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
      
      model_0_or <<- eval(parse(text = text_0))
      
      if(!file.exists(paste0(path,"/warn_or"))){dir.create(paste0(path,"/warn_or"))}
      save(warn_or,file = paste0(path,"/warn_or/warn_0_or_",x,"_",comp,".Rdata"))
    })
    
    
    tryCatch({
      model_1_or <<- eval(parse(text = text_1))
      
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
      
      model_1_or <<- eval(parse(text = text_1))
      
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

pblapply(X = unique(data$study),
         FUN = bayes_model_func,
         data = data,
         SE = "t2_new_or_enlarged",
         TE = "relapse",
         path = "results_count_count_IT_predictive",
         cl = cl)

stopCluster(cl)