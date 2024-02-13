
# Purpose of this script: Sensitivity analyses using the Prentice criteria to validate surrogate. 
# PTE (proportion explained) is surrogacy metric. (Replication of paper of Sormani about individual surrogacy)
# Associative approach
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
                        is.na(day),max(day,na.rm = T),day)) %>%
  distinct() %>%
  group_by(study,id) %>%
  mutate(relapse = ifelse(is.na(relapse),0,relapse)) %>%
  summarise(study = unique(study),
            id = unique(id),
            arm = unique(arm),
            relapse = sum(relapse, na.rm = T)) %>%
  group_by(study,id) %>%
  summarise(arm = unique(arm),
            relapse = sum(relapse)) %>%
  as.data.frame()

dim(data_TE)
dim(data_TE[complete.cases(data_TE),])

data_SE <- mri %>%
  filter(!(study == "WA21493")) %>%
  dplyr::select(study,id,arm,t2_new_or_enlarged,day,visit,visit_com) %>%
  filter(!is.na(t2_new_or_enlarged)) %>% 
  filter(visit != 800) %>% #visit 800 is cumulated t2 lesion count
  group_by(study,id,day) %>%
  mutate(n = n()) %>%
  filter(!(n == 2 & visit > 50)) %>% 
  filter(visit_com != 0) %>% # remove baseline count (since this is not the change in t2 lesions count)
  arrange(study,arm,id,day) %>%
  dplyr::mutate(day = day + (min(c(mri$day,mri$day))-1)*-1) %>%
  filter(visit_com != 999) %>%
  filter(!(study == "CFTY720D1201" & visit_com == 777)) %>%
  filter(!(study == "CFTY720D2301" & visit_com == 777)) %>%
  filter(!(study == "CFTY720D2309" & visit_com == 24)) %>%
  filter(!(study == "CFTY720D2309" & visit_com == 501)) %>%
  filter(!(study == "WA21493" & visit_com == 22)) %>%
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

  


#####
prentice <- function(x = unique(data$study),
                             data = data,
                             SE = "t2_new_or_enlarged",
                             TE = "relapse",
                             path = "results_replication_associative")
{
  library(lme4)
  library(pbapply)
  library(parallel)
  library(tidyverse)
  library(nlme)
  library(ggcorrplot)
  library(rlist)
  library(boot)
  
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
  
  data$treat <- ifelse(data$arm %in% control_gr,-1,1)
  
  
  # The analysis have to be carried out for each active drug vs. control group arm
  for(comp in 1:nrow(comparisons))
  {
    # prepare data 
    data_loop <- data %>%
      filter(arm %in% c(comparisons$control_gr[comp],comparisons$test_gr[comp]))
  
    library(MASS)
    # Prentice criteria
    ## 1. criterion
    
    formula_1 <- paste0(SE," ~ treat")
    
    tryCatch({
      model_1 <- glm.nb(formula = as.formula(formula_1),
                            data = data_loop)
    },
    error = function(e)
    {
      errors_1 <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_1"))){dir.create(paste0(path,"/errors_1"))}
      save(errors_1,file = paste0(path,"/errors_1/errors_1_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_1 <<- glm.nb(formula = as.formula(formula_1),
                        data = data_loop)
      warn_1 <<- w
      
      if(!file.exists(paste0(path,"/warn_1"))){dir.create(paste0(path,"/warn_1"))}
      save(warn_1,file = paste0(path,"/warn_1/warn_1_",x,"_",comp,".Rdata"))
    })
    
    
    if(!file.exists(paste0(path,"/model_1"))){dir.create(paste0(path,"/model_1"))}
    save(model_1,file = paste0(path,"/model_1/model_1_",x,"_",comp,".Rdata"))
    
    
    
    
    
    ## 2. criterion
    formula_2 <- paste0(TE," ~ treat")
    
    tryCatch({
      model_2 <- glm.nb(formula = as.formula(formula_2),
                        data = data_loop)
    },
    error = function(e)
    {
      errors_2 <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_2"))){dir.create(paste0(path,"/errors_2"))}
      save(errors_2,file = paste0(path,"/errors_2/errors_2_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_1 <<- glm.nb(formula = as.formula(formula_1),
                         data = data_loop)
      warn_2 <<- w
      
      if(!file.exists(paste0(path,"/warn_2"))){dir.create(paste0(path,"/warn_2"))}
      save(warn_2,file = paste0(path,"/warn_2/warn_2_",x,"_",comp,".Rdata"))
    })
    
    
    if(!file.exists(paste0(path,"/model_2"))){dir.create(paste0(path,"/model_2"))}
    save(model_2,file = paste0(path,"/model_2/model_2_",x,"_",comp,".Rdata"))
    
    
    
    
    ## 3. criterion
    formula_3 <- paste0(TE," ~ ",SE)
    
    tryCatch({
      model_3 <- glm.nb(formula = as.formula(formula_3),
                        data = data_loop)
    },
    error = function(e)
    {
      errors_3 <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_3"))){dir.create(paste0(path,"/errors_3"))}
      save(errors_3,file = paste0(path,"/errors_3/errors_3_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_3 <<- glm.nb(formula = as.formula(formula_3),
                         data = data_loop)
      warn_3 <<- w
      
      if(!file.exists(paste0(path,"/warn_3"))){dir.create(paste0(path,"/warn_3"))}
      save(warn_3,file = paste0(path,"/warn_3/warn_3_",x,"_",comp,".Rdata"))
    })
    
    
    if(!file.exists(paste0(path,"/model_3"))){dir.create(paste0(path,"/model_3"))}
    save(model_3,file = paste0(path,"/model_3/model_3_",x,"_",comp,".Rdata"))
    
    
    
    
    
    ## 4. Criterion
    formula_4 <- paste0(TE," ~ treat + ",SE)
    
    tryCatch({
      model_4 <- glm.nb(formula = as.formula(formula_4),
                        data = data_loop)
    },
    error = function(e)
    {
      errors_4 <<- geterrmessage()
      if(!file.exists(paste0(path,"/errors_4"))){dir.create(paste0(path,"/errors_4"))}
      save(errors_4,file = paste0(path,"/errors_4/errors_4_",x,"_",comp,".Rdata"))
    },
    warning = function(w)
    {
      model_4 <<- glm.nb(formula = as.formula(formula_4),
                         data = data_loop)
      warn_4 <<- w
      
      if(!file.exists(paste0(path,"/warn_4"))){dir.create(paste0(path,"/warn_4"))}
      save(warn_4,file = paste0(path,"/warn_4/warn_4_",x,"_",comp,".Rdata"))
    })
    
    
    if(!file.exists(paste0(path,"/model_4"))){dir.create(paste0(path,"/model_4"))}
    save(model_4,file = paste0(path,"/model_4/model_4_",x,"_",comp,".Rdata"))
    
    
    P_model1 <- summary(model_1)$coefficients
    rownames(P_model1) <- c("Intercept", "Treatment")
    
    P_model2 <- summary(model_2)$coefficients
    rownames(P_model2) <- c("Intercept", "Treatment")
    
    P_model3 <- summary(model_3)$coefficients
    rownames(P_model3) <- c("Intercept", "Surrogate")
    
    P_model4 <- summary(model_4)$coefficients
    rownames(P_model4) <- c("Intercept", "Treatment", 
                            "Surrogate")
    
    P_crit1 <- summary(model_1)$coefficients[2,4]
    P_crit2 <- summary(model_2)$coefficients[2,4]
    P_crit3 <- summary(model_3)$coefficients[2,4]
    P_crit4 <- data.frame(summary(model_4)$coefficients, stringsAsFactors = TRUE)[2, 4]
    
    Alpha <- 0.05
    
    if ((P_crit1 < Alpha & P_crit2 < Alpha & P_crit3 < Alpha & 
         P_crit4 > Alpha) == TRUE) {
      Prentice.Passed <- TRUE
    }
    if ((P_crit1 < Alpha & P_crit2 < Alpha & P_crit3 < Alpha & 
         P_crit4 > Alpha) == FALSE) {
      Prentice.Passed <- FALSE
    }
    fit <- list(Prentice.Model.1 = P_model1, Prentice.Model.2 = P_model2, 
                Prentice.Model.3 = P_model3, Prentice.Model.4 = P_model4, 
                Prentice.Passed = Prentice.Passed, Call = match.call())
    class(fit) <- "Prentice"
    
    
    if(!file.exists(paste0(path,"/results"))){dir.create(paste0(path,"/results"))}
    save(fit,file = paste0(path,"/results/results_",x,"_",comp,".Rdata"))
    
    
    # Estimation of confidence intervals for PTE using bootstrap
    source("scripts/boot_PTE_func.R")
    
    boot_results <- boot(data = data_loop,
                         statistic = PET_boot,
                         SE = SE,
                         TE = TE,
                         R = 1000)
    
    ci_beta <- boot.ci(boot_results,
                       type = "bca",
                       index = 1)
    
    ci_beta_s <- boot.ci(boot_results,
                         type = "bca",
                         index = 2)
    
    ci_PTE <- boot.ci(boot_results,
                      type = "bca",
                      index = 3)
    
    
    
    df <- data.frame(alpha = round(fit$Prentice.Model.1[2,1],2),
                     p_alpha = ifelse(P_crit1 < 0.001,"< 0.001",round(P_crit1,3)),
                     beta = round(fit$Prentice.Model.2[2,1],2),
                     p_beta = ifelse(P_crit2 < 0.001,"< 0.001",round(P_crit2,3)),
                     gamma = round(fit$Prentice.Model.3[2,1],2),
                     p_gamma = ifelse(P_crit3 < 0.001,"< 0.001",round(P_crit3,3)),
                     beta_s = round(fit$Prentice.Model.4[2,1],2),
                     p_beta_s = ifelse(P_crit4 < 0.001,"< 0.001",round(P_crit4,3)),
                     Prentice_Passed = Prentice.Passed,
                     PTE = round(1-(fit$Prentice.Model.4[2,1]/fit$Prentice.Model.2[2,1]),2),
                     PTE_lb = round(ci_PTE$bca[1,4],2),
                     PTE_ub = round(ci_PTE$bca[1,5],2))
    
    if(!file.exists(paste0(path,"/table"))){dir.create(paste0(path,"/table"))}
    save(df,file = paste0(path,"/table/table_",x,"_",comp,".Rdata"))
    
   rm(list = c("model_1","model_2","model_3","model_4",
               "formula_1","formula_2","formula_3","formula_4",
               "fit","df","data_loop","ci_PTE"))
    
  }
  
}

library(pbapply)
library(parallel)
cl <- makeCluster(detectCores()-4)

pblapply(X = unique(data$study),
         FUN = prentice,
         data = data,
         SE = "t2_new_or_enlarged",
         TE = "relapse",
         path = "results_replication_associative",
         cl = cl)

stopCluster(cl)
#####



