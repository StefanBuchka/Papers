# Purpose: Estimation of LRF using the IT methodology
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka

# Note: The clinical endpoint (CEP) is named TE and the surrogate endpoint (SEP) is named SE in the following 


library(tidyverse)
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
  distinct()

dim(data_TE)
dim(data_TE[complete.cases(data_TE),])

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

bayes_model_func <- function(x = unique(data_delta$study),
                             data = data_delta,
                             forced_powers = NULL,
                             SE = "t2_new_or_enlarged",
                             TE = "rel",
                             time = "day_SE",
                             path = "results_count_count_IT")
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
    
    data_loop <- data_loop %>%
      filter(!is.na(t2_new_or_enlarged))
    
    
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
    
    
    data_loop$t2_cat <- ifelse(data_loop$t2_new_or_enlarged >=8,8,data_loop$t2_new_or_enlarged)
    
    data_loop$rel <- data_loop$rel + 1
    data_loop$t2_cat <- data_loop$t2_cat + 1
    
    #debugonce("IT_model_sel_func")
    models_or <- IT_model_sel_func(data_loop = data_loop,
                                   family = "ordinal",
                                   TE = TE,
                                   SE = "t2_cat",
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
         data = data_delta,
         forced_powers = NULL,
         SE = "t2_new_or_enlarged",
         TE = "rel",
         time = "day_SE",
         path = "results_count_count_IT",
         cl = cl)

stopCluster(cl)
