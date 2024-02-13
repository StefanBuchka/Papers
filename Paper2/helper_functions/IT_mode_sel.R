# Purpose: Find the best model for the IT approach. In general, three models per for each model of the IT approach were considered:  
# 1.) Interaction between treatment and time (non-linear time trajectories were considered by fractional polynomials (variable name: term_TE))
# 2.) Model with treatment and time variable (non-linear time trajectories were considered by fractional polynomials (variable name: term_TE))
# 3.) Model only adjusted for treatment
# Details: see paper and appendix 
# 19.12.2023
# Author: Stefan Buchka



library(ordinal)
IT_model_sel_func <- function(data_loop = data_loop,
                              family = "zero_inflated",
                              TE = TE,
                              SE = SE,
                              term_TE = NULL,
                              trial = NULL)
{

  if(family != "gaus") # If the clinical endpoint (CEP) is number of relapses, an offset was included 
  {
    
    # Model for IT approach without the surrogate 
    
    # formula_0_1 <- as.formula(paste0(TE,"~ offset(log(year)) + treat*(",term_TE,") + (1|id)"))
    # formula_0_2 <- as.formula(paste0(TE,"~ offset(log(year)) +treat*(",term_TE,") + (1|id)"))
    # formula_0_3 <- as.formula(paste0(TE,"~ offset(log(year)) +treat+",term_TE,"+(1|id)"))
    # formula_0_4 <- as.formula(paste0(TE,"~ offset(log(year)) + treat + (1|id)"))
    formula_0_5 <- as.formula(paste0(TE,"~ offset(log(year)) +treat*(",term_TE,")")) # Interaction treatment and time
    formula_0_6 <- as.formula(paste0(TE,"~ offset(log(year)) +treat+",term_TE)) # No Interaction between treatment and time
    formula_0_7 <- as.formula(paste0(TE,"~ offset(log(year)) + treat")) # model adjusted for treatment only
    formula_0 <- list(formula_0_5,formula_0_6,formula_0_7)
    
    # Model for IT approach wit the surrogate
    
    # formula_1_1 <- as.formula(paste0(TE,"~ offset(log(year)) + ",SE,"+treat*(",term_TE,") + (1|id)"))
    # formula_1_2 <- as.formula(paste0(TE,"~ offset(log(year)) + ",SE,"+treat*(",term_TE,") + (1|id)"))
    # formula_1_3 <- as.formula(paste0(TE,"~ offset(log(year)) + ",SE,"+treat+",term_TE,"+(1|id)"))
    # formula_1_4 <- as.formula(paste0(TE,"~ offset(log(year)) + ",SE,"+treat + (1|id)"))
    formula_1_5 <- as.formula(paste0(TE,"~ offset(log(year)) + ",SE,"+treat*(",term_TE,")")) # Interaction treatment and time
    formula_1_6 <- as.formula(paste0(TE,"~ offset(log(year)) + ",SE,"+treat+",term_TE)) # No Interaction between treatment and time
    formula_1_7 <- as.formula(paste0(TE,"~ offset(log(year)) + ",SE,"+treat")) # model adjusted for treatment only
    formula_1 <- list(formula_1_5,formula_1_6,formula_1_7)
  }
  
  if(family == "gaus")
  {
    # Model for IT approach without the surrogate
    
    # formula_0_1 <- as.formula(paste0(TE,"~  treat*(",term_TE,") + (1|id)"))
    # formula_0_2 <- as.formula(paste0(TE,"~ treat*(",term_TE,") + (1|id)"))
    # formula_0_3 <- as.formula(paste0(TE,"~ treat+",term_TE,"+(1|id)"))
    # formula_0_4 <- as.formula(paste0(TE,"~  treat + (1|id)"))
    formula_0_5 <- as.formula(paste0(TE,"~ treat*(",term_TE,")")) # Interaction treatment and time
    formula_0_6 <- as.formula(paste0(TE,"~ treat+",term_TE)) # No Interaction between treatment and time
    formula_0_7 <- as.formula(paste0(TE,"~  treat")) # model adjusted for treatment only
    formula_0 <- list(formula_0_5,formula_0_6,formula_0_7)
    
    # Model for IT approach with the surrogate
    
    # formula_1_1 <- as.formula(paste0(TE,"~  ",SE,"+treat*(",term_TE,") + (1|id)"))
    # formula_1_2 <- as.formula(paste0(TE,"~  ",SE,"+treat*(",term_TE,") + (1|id)"))
    # formula_1_3 <- as.formula(paste0(TE,"~  ",SE,"+treat+",term_TE,"+(1|id)"))
    # formula_1_4 <- as.formula(paste0(TE,"~  ",SE,"+treat + (1|id)"))
    formula_1_5 <- as.formula(paste0(TE,"~  ",SE,"+treat*(",term_TE,")")) # Interaction treatment and time
    formula_1_6 <- as.formula(paste0(TE,"~  ",SE,"+treat+",term_TE)) # No Interaction between treatment and time
    formula_1_7 <- as.formula(paste0(TE,"~  ",SE,"+treat")) # model adjusted for treatment only
    formula_1 <- list(formula_1_5,formula_1_6,formula_1_7)
  }  
  
  model_0 <- vector(mode = "list", length = length(formula_0))
  model_1 <- vector(mode = "list", length = length(formula_1))
  model_temp <- NULL
  
  i <- 1
  
  # Estimate the modles with different model families
  while(i <= length(model_0))
  {    
    if(family == "zero_inflated")
    {
    
      model_temp <- NULL
      tryCatch(
        {
          model_0[[i]] <- glmmTMB(as.formula(formula_0[[i]]), 
                                  ziformula = ~ 1,
                                  family = "poisson",
                                  data = data_loop)
          
        },
        error = function(e)
        {
          model_temp <<- glmmTMB(as.formula(formula_0[[7]]), 
                                  ziformula = ~ 1,
                                  family = "poisson",
                                  data = data_loop)
        })
      
      if(!is.null(model_temp))
      {
        model_0[[i]] <- model_temp
        model_temp <- NULL
      }
      
      
      tryCatch(
        {
          model_1[[i]] <- glmmTMB(as.formula(formula_1[[i]]), 
                                  ziformula = ~ 1,
                                  family = "poisson",
                                  data = data_loop)
          
        },
        error = function(e)
        {
          model_temp <<- glmmTMB(as.formula(formula_1[[7]]), 
                                 ziformula = ~ 1,
                                 family = "poisson",
                                 data = data_loop)
        })
      
      if(!is.null(model_temp))
      {
        model_1[[i]] <- model_temp
        model_temp <- NULL
      }
      
    }
    
    
    
    if(family == "gaus")
    {
      
      model_temp <- NULL
      
      model_0[[i]] <- glmmTMB(as.formula(formula_0[[i]]), 
                              family = gaussian(),
                              data = data_loop)
      
      
      model_1[[i]] <- glmmTMB(as.formula(formula_1[[i]]), 
                              family = gaussian(),
                              data = data_loop)
      
    }
    
    
    
    if(family == "nbinomial")
    {
      
      model_temp <- NULL
      tryCatch(
        {
          model_0[[i]] <- glmmTMB(as.formula(formula_0[[i]]), 
                                 family = "nbinom2",
                                 data = data_loop)
          
        },
        error = function(e)
        {
          model_temp <<- glmmTMB(as.formula(formula_0[[7]]), 
                                 family = "nbinom2",
                                 data = data_loop)
        })
      
      if(!is.null(model_temp))
      {
        model_0[[i]] <- model_temp
        model_temp <- NULL
      }
      
    
      
      model_temp <- NULL
      tryCatch(
        {
          model_1[[i]] <- glmmTMB(as.formula(formula_1[[i]]), 
                                  family = "nbinom2",
                                  data = data_loop)
          
        },
        error = function(e)
        {
          model_temp <<- glmmTMB(as.formula(formula_1[[7]]), 
                                 family = "nbinom2",
                                 data = data_loop)
        })
      
      if(!is.null(model_temp))
      {
        model_1[[i]] <- model_temp
      }
        
    }
    
    
    if(family == "poisson")
    {
      model_temp <- NULL
      
      tryCatch(
        {
          model_0[[i]] <- glmmTMB(as.formula(formula_0[[i]]), 
                                  family = "poisson",
                                  data = data_loop)
          
        },
        error = function(e)
        {
          model_temp <<- glmmTMB(as.formula(formula_0[[7]]), 
                                  family = "poisson",
                                  data = data_loop)
        })
      
      if(!is.null(model_temp))
      {
        model_0[[i]] <- model_temp
        model_temp <- NULL
      }
      
      
      tryCatch(
        {
          model_1[[i]] <- glmmTMB(as.formula(formula_1[[i]]), 
                                  family = "poisson",
                                  data = data_loop)
          
        },
        error = function(e)
        {
          model_temp <<- glmmTMB(as.formula(formula_1[[7]]), 
                                  family = "poisson",
                                  data = data_loop)
        })
      
      if(!is.null(model_temp))
      {
        model_1[[i]] <- model_temp
        model_temp <- NULL
      }
    }
    
    i <- i + 1
    
    if(family == "ordinal")
    {
      
      # Model for IT approach without the surrogate
      formula_0_1 <- (paste0("factor(",TE,")~  treat*(",term_TE,")")) # Interaction treatment and time
      formula_0_2 <- (paste0("factor(",TE,")~  treat + (",term_TE,")")) # No Interaction between treatment and time
      formula_0_3 <- (paste0("factor(",TE,")~  treat")) # model adjusted for treatment only
      formula_0 <- list(formula_0_1,formula_0_2,formula_0_3)
      
      # Model for IT approach with the surrogate
      formula_1_2 <- (paste0("factor(",TE,")~  as.numeric(", SE, ") + treat*(",term_TE,")")) # Interaction treatment and time
      formula_1_3 <- (paste0("factor(",TE,")~  as.numeric(", SE, ") + treat + (",term_TE,")")) # No Interaction between treatment and time
      formula_1_4 <- (paste0("factor(",TE,")~  as.numeric(", SE, ") + treat")) # model adjusted for treatment only
      formula_1 <- list(formula_1_2,formula_1_3,formula_1_4)
      
      model_0 <- vector(mode = "list", length = length(formula_0))
      model_1 <- vector(mode = "list", length = length(formula_1))
      
      model_temp <- NULL
      
      for(q in 1:3)
      {
        text <- paste0("clmm2(",formula_0[[q]],",Hess = T,data = data_loop)")
        tryCatch(
          {
            model_0[[q]] <- eval(parse(text = text))
            
          },
          error = function(e)
          {
            
            text <<- paste0("clmm2(",formula_0[[3]],",Hess = T,data = data_loop)")
            model_temp <<- eval(parse(text = text))
          })
        
        if(!is.null(model_temp))
        {
          model_0[[q]] <- model_temp
          model_temp <- NULL
        }
        
        rm(text)
        # text <- paste0("clmm2(",formula_0[[q]],",Hess = T,data = data_loop)")
        # tryCatch(
        #   {
        #     model_0[[q+3]] <- eval(parse(text = text))
        #     
        #   },
        #   error = function(e)
        #   {
        #     text <<- paste0("clmm2(",formula_0[[3]],",Hess = T,data = data_loop)")
        #     model_temp <<- eval(parse(text = text))
        #   })
        # 
        # 
        # if(!is.null(model_temp))
        # {
        #   model_0[[q+3]] <- model_temp
        #   model_temp <- NULL
        # }
        

        text <- paste0("clmm2(",formula_1[[q]],",Hess = T,data = data_loop)")
        tryCatch(
          {
            model_1[[q]] <- eval(parse(text = text))
            
          },
          error = function(e)
          {
            text <<- paste0("clmm2(",formula_1[[3]],",Hess = T,data = data_loop)")
            model_temp <<- eval(parse(text = text))
          })
        
        if(!is.null(model_temp))
        {
          model_1[[q]] <- model_temp
          model_temp <- NULL
        }
        
        
        rm(text)
        # text <- paste0("clmm2(",formula_1[[q]],",Hess = T,data = data_loop)")
        # tryCatch(
        #   {
        #     model_1[[q+3]] <- eval(parse(text = text))
        #     
        #   },
        #   error = function(e)
        #   {
        #     
        #     text <<- paste0("clmm2(",formula_1[[3]],",Hess = T,data = data_loop)")
        #     model_temp <<- eval(parse(text = text))
        #   })
        # 
        # if(!is.null(model_temp))
        # {
        #   model_1[[q+3]] <- model_temp
        #   model_temp <- NULL
        # }
        # 
        # rm(text)
        
        i <- length(model_0) + i
        
      }
      
      
    }
    
  }  
  
  # Extract AIC and log likelihood
  lr_0 <- data.frame(model_name = 1:length(model_0),
                     AIC = sapply(model_0,AIC),
                     logLik = sapply(model_0,logLik))
  
  lr_1 <- data.frame(model_name = 1:length(model_1),
                     AIC = sapply(model_1,AIC),
                     logLik = sapply(model_1,logLik))
  
  
  
  lr_0 <- lr_0 %>%
    filter(!is.na(logLik)) %>%
    arrange(AIC)
  
  lr_1 <- lr_1 %>%
    filter(!is.na(logLik)) %>%
    arrange(AIC)
  
  # Get best model selected by AIC.
  model_0 <- model_0[[lr_1[1,"model_name"]]] # This takes the best model with surrogate, so that both models are comparable
  model_1 <- model_1[[lr_1[1,"model_name"]]]
  
  if(family == "ordinal")
  {
    model_0 <- deparse(model_1$call)
    model_0 <- gsub(paste0("as.numeric\\(",SE,"\\) \\+"),"",model_0) #removes Surrogate from formula
    model_0 <- eval(parse(text = model_0))
  }
  
  
  return(list(model_0,model_1))
  
}
