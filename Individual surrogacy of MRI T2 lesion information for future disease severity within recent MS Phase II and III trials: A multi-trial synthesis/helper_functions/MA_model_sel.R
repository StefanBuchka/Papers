# Purpose: Find the best model for the MA approach. In general, three models per outcome model of the MA approach were considered: 
# The model selection is conducted for each outome separatley. Afterwards, the joint model with both model formulas will be estimated.
# 1.) Interaction between treatment and time (non-linear time trajectories were considered by fractional polynomials (variable name: term_TE))
# 2.) Model with treatment and time variable (non-linear time trajectories were considered by fractional polynomials (variable name: term_TE))
# 3.) Model only adjusted for treatment
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka


MA_model_sel_func <- function(data_loop = data_loop,
                              Outcome = TE,
                              family = "gaus",
                              year = "year",
                              trial = F,
                              term = NULL)
{
  
 
  data_loop$visit_com <- as.factor(data_loop$visit_com)
  data_loop$visit_com2 <- numFactor(data_loop$visit_com)
  data_loop$id_fac <- as.factor(data_loop$id)
  model <- NULL
  

  if(family == "gaus")
  {

      # formula_0_2 <- as.formula(paste0(Outcome,"~","treat*(",term,") + (1|id)"))
      # formula_0_3 <- as.formula(paste0(Outcome,"~","treat+",term,"+(1|id)"))
      # formula_0_4 <- as.formula(paste0(Outcome,"~ treat + (1|id)"))
      formula_0_5 <- as.formula(paste0(Outcome,"~","treat*(",term,")")) # Interaction treatment and time
      formula_0_6 <- as.formula(paste0(Outcome,"~","treat+",term)) # No Interaction between treatment and time
      formula_0_7 <- as.formula(paste0(Outcome,"~ treat")) # model adjusted for treatment only
      
      # Autocorrelation
      # formula_0_8 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ (1|id) + ar1(visit_com + 0|id)"))
      # formula_0_9 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ (1|id) + ar1(visit_com + 0|id)"))
      # formula_0_10 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ ar1(visit_com + 0|id)"))
      # formula_0_11 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ ar1(visit_com + 0|id)"))
      # formula_0_12 <- as.formula(paste0(Outcome,"~","treat + (1|id) + ar1(visit_com + 0|id)"))
      # formula_0_13 <- as.formula(paste0(Outcome,"~","treat + ar1(visit_com + 0|id)"))
      # 
      # formula_0_14 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ (1|id) + exp(visit_com2 + 0|id)"))
      # formula_0_15 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ (1|id) + exp(visit_com2 + 0|id)"))
      # formula_0_16 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ exp(visit_com2 + 0|id)"))
      # formula_0_17 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ exp(visit_com2 + 0|id)"))
      # formula_0_18 <- as.formula(paste0(Outcome,"~","treat + (1|id) + exp(visit_com2 + 0|id)"))
      # formula_0_19 <- as.formula(paste0(Outcome,"~","treat + exp(visit_com2 + 0|id)"))
      # 
      # formula_0_20 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ (1|id) + mat(visit_com2 + 0|id)"))
      # formula_0_21 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ (1|id) + mat(visit_com2 + 0|id)"))
      # formula_0_22 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ mat(visit_com2 + 0|id)"))
      # formula_0_23 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ mat(visit_com2 + 0|id)"))
      # formula_0_24 <- as.formula(paste0(Outcome,"~","treat + (1|id) + mat(visit_com2 + 0|id)"))
      # formula_0_25 <- as.formula(paste0(Outcome,"~","treat + mat(visit_com2 + 0|id)"))
      # 
      # formula_0_26 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ (1|id) + gau(visit_com2 + 0|id)"))
      # formula_0_27 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ (1|id) + gau(visit_com2 + 0|id)"))
      # formula_0_28 <- as.formula(paste0(Outcome,"~","treat*(",term,")","+ gau(visit_com2 + 0|id)"))
      # formula_0_29 <- as.formula(paste0(Outcome,"~","treat + (",term,")","+ gau(visit_com2 + 0|id)"))
      # formula_0_30 <- as.formula(paste0(Outcome,"~","treat + (1|id) + gau(visit_com2 + 0|id)"))
      # formula_0_31 <- as.formula(paste0(Outcome,"~","treat + gau(visit_com2 + 0|id)"))
      
      formula_text <- paste0("list(",paste0("formula_0_",c(5:7),collapse = ","),")")
      formula <- eval(parse(text = formula_text))
      
      library(foreach)
      library(doParallel)
      
      cl <- makeCluster(detectCores()-4)
      registerDoParallel(cl)
      
      # Fit considered models
      model <- foreach(i = 1:length(formula)) %dopar%
      {
        library(glmmTMB)
        warn <- NA
        
        tryCatch({
          model <- glmmTMB(as.formula(formula[[i]]), 
                            family = gaussian(),
                            data = data_loop)
          list(model,warn)
        },
        warning = function(w)
        {
          warn <- w # collect warnings during model fitting
          model <- glmmTMB(as.formula(formula[[i]]), 
                           family = gaussian(),
                           data = data_loop)
          list(model,warn)
        },
        error = function(e)
        {
          # If a model is to complex to fit, then the simplest model (only treatment variable included) will bit fit.
          model <- glmmTMB(as.formula(paste0(Outcome,"~ treat")), 
                            family = gaussian(),
                            data = data_loop)
          list(model,warn)
        })
    }   
    
      stopCluster(cl)
      
     # Select best model by AIC
      lr <- data.frame(model_name = 1:length(model),
                       AIC = sapply(model,AIC),
                       logLik = sapply(model,logLik),
                       call = sapply(model,function(x){as.character(x$call)[2]}),
                       warnings = list.rbind(warn))
    
 
  }
  
  
  if(family == "poisson")
  {

      # formula_0_2 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat*(",term,") + (1|id)"))
      # formula_0_3 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat+",term,"+(1|id)"))
      # formula_0_4 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) +  treat + (1|id)"))
      formula_0_5 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat*(",term,")")) # Interaction treatment and time
      formula_0_6 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat+",term)) # No Interaction between treatment and time
      formula_0_7 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat")) # model adjusted for treatment only
      
      # Autocorrelation
      # formula_0_8 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat*(",term,")","+ (1|id) + ar1(visit_com + 0|id)"))
      # formula_0_9 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat + (",term,")","+ (1|id) + ar1(visit_com + 0|id)"))
      # formula_0_10 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat*(",term,")","+ ar1(visit_com + 0|id)"))
      # formula_0_11 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat + (",term,")","+ ar1(visit_com + 0|id)"))
      # formula_0_12 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat + (1|id) + ar1(visit_com + 0|id)"))
      # formula_0_13 <- as.formula(paste0(Outcome,"~ offset(log(",year,")) + treat + ar1(visit_com + 0|id)"))
    
      
      formula_text <- paste0("list(",paste0("formula_0_",5:7,collapse = ","),")")
      formula <- eval(parse(text = formula_text))

      
      
      library(foreach)
      library(doParallel)
      
      cl <- makeCluster(detectCores()-4)
      registerDoParallel(cl)
      
      # Fit considered models
      model <- foreach(i = 1:length(formula)) %dopar%
      {
        library(glmmTMB)
        warn <- NA
        
        tryCatch({
          model <- glmmTMB(as.formula(formula[[i]]), 
                                family = poisson(),
                                data = data_loop)
          list(model,warn)
        },
        warning = function(w)
        {
          warn <- w # collect warnings during model fitting
          model <- glmmTMB(as.formula(formula[[i]]), 
                           family = poisson(),
                           data = data_loop)
          list(model,warn)
        },
        error = function(e)
        { # If a model is to complex to fit, then the simplest model (only treatment variable included) will bit fit.
          model <- glmmTMB(as.formula(paste0(Outcome,"~ treat")), 
                           family = poisson(),
                           data = data_loop)
          list(model,warn)
        })
      }   
      
      stopCluster(cl)
    
    
      
      warn <- lapply(model,function(x){x[[2]]})
      model <- lapply(model,function(x){x[[1]]})
      
      # Select best model by AIC
      lr <- data.frame(model_name = 1:length(model),
                       AIC = sapply(model,AIC),
                       logLik = sapply(model,logLik),
                       call = sapply(model,function(x){as.character(x$call)[2]}),
                       warnings = list.rbind(warn))
    
  }
  
 

  # Select best model by AIC
    lr <- lr %>%
      arrange(AIC)
  
    
  
    model <- list(model,lr)
    names(model) <- c("model","lr_list")
  
  return(model)
  
}

