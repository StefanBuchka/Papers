# Purpose: Estimation of LRF using the IT methodology
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka

# Note: The clinical endpoint (CEP) is named TE and the surrogate endpoint (SEP) is named SE in the following 
# Note: count data was transformed

library(tidyverse)
load("datasets/arr.RData")
load("datasets/mri.RData")
load("datasets/edss.RData")


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

# The baseline t2 count were removed, since the base line count are the numbers of t2 lesions at base line. This is the the new or enlarged lesoin count. 
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
                           visit_com == 0,0,day_SE), #replace the first day of visit by 0, to inclued relapses happend before the first MRI visit in the interval.
         day_SE_end = lead(day_SE,
                           default = max(c(day_SE,day_TE),na.rm = T)), # Add the last visit or relapse date (that one, which is later).
         t2_new_or_enlarged = lead(t2_new_or_enlarged)) %>% #assigne the t2_lesion count to right time point  
  filter(day_SE != day_SE_end) %>%
  mutate(rel = ifelse(relapse == 1 &
                        (day_TE > day_SE & day_TE <= day_SE_end),1,0)) %>% # bild the intervals
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
                             path = "results_count_count_IT_t")
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
                    Dataset = as.data.frame(data_loop))
    
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
    
    
    
    
    
    # bind Fractional polynomial terms to data set
    data_loop <- cbind(data_loop,term)
    
    data_loop$obs <- 1:nrow(data_loop)
    
    term_TE <- paste0(grep("term",names(data_loop),value = T),collapse = "+")
    
    data_loop <- data_loop %>%
      filter(!is.na(t2_new_or_enlarged))
    
    
    # Model estimation
    source("scripts/IT_mode_sel.R") # model selection function
    source("scripts/calculation_r_2.R")
    
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
      mutate(rel = ifelse(rel >= 1,2,1),
             year = year)
    
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
         path = "results_count_count_IT_t",
         cl = cl)

stopCluster(cl)









