library(boot)
library(tidyverse)
library(rlist)
library(ggrepel)

set.seed(82166)

# Bootstrap-function for confidence intervals
boot_function <- function(data, indices) {
  
  d <- data[indices, ]  # Resample mit Indizes
  
  # Linear regression models according to Prentice criteria to derive PTE.
  # Note, this models are the same as for deriving R²hindiv based on the information theoretic approach based on ideas of Alonso.
  model_CEP <- model_CE_on_Treatment <- lm(CEP ~ `Study arm`, data = d)
  model_CEP_SEP <- model_CE_on_Treatment_SE <- lm(CEP ~ `Study arm` + SEP, data = d)
  
  # Calculation of PTE
  treatment_effect_only <- coef(model_CE_on_Treatment)["`Study arm`Untreated"]
  adjusted_treatment_effect <- coef(model_CE_on_Treatment_SE)["`Study arm`Untreated"]
  PTE <- (treatment_effect_only - adjusted_treatment_effect) / treatment_effect_only
  
  # Test, if Deviances differ
  D_CEP <- sum(residuals(model_CEP, "deviance")^2)
  D_CEP_SEP <- sum(residuals(model_CEP_SEP, "deviance")^2)
  D_delta <- D_CEP - D_CEP_SEP
  df_diff <- df.residual(model_CEP) - df.residual(model_CEP_SEP)
  p_val_dev <- pchisq(D_delta, df = df_diff, lower.tail = FALSE)
  
  # Likelihood Ratio Test statistic
  LL_CEP <- logLik(model_CEP)
  LL_CEP_SEP <- logLik(model_CEP_SEP)
  LR_delta <- -2 * (LL_CEP - LL_CEP_SEP)
  p_val_lr <- pchisq(LR_delta, df = 1, lower.tail = FALSE) # df = 1, because full model (model_CEP_SEP) has one more parameter than reduced model (model_CEP)
  
  # Calculation of R2hindiv
  
  H_CEP <- -LL_CEP / nrow(d)
  H_CEP_SEP <- -LL_CEP_SEP / nrow(d)
  MI <- H_CEP - H_CEP_SEP
  R2h <- 1 - exp(-2 * MI)
  
  # Calcualte correlation
  r <- cor.test(d$SEP,d$CEP)$estimate
  
  return(c(PTE = PTE, R2h = R2h, p_val_dev = p_val_dev,p_val_lr = p_val_lr,r = r))
}

IT_vs_dev_LR <- function(N = N,
                         var_SEP = var_SEP,
                         var_CEP = var_CEP,
                         boot_function = boot_function,
                         n_boot = 1000)
{
  # Data simulation
  SEP_untreated <- rnorm(N, 6, var_SEP)
  CEP_untreated <- rnorm(N, 0, var_CEP) + 1.5 * SEP_untreated
  
  SEP_treated <- rnorm(N, 4, var_SEP)
  CEP_treated <- rnorm(N, 0, var_CEP) + 1.5 * SEP_treated
  
  SEP <- c(SEP_untreated, SEP_treated)
  CEP <- c(CEP_untreated, CEP_treated)
  group <- factor(rep(c("Untreated", "Treated"), each = N))
  
  data <- tibble(SEP = SEP,
                 CEP = CEP,
                 `Study arm` = group)
  
  # Bootstrap for CIs
  boot_res <- boot(data = data, statistic = boot_function, R = n_boot)
  
  # Result tibble
  boot_estimates <- tibble(
    #Estimate = boot_res$t0,
    #Bias = colMeans(boot_res$t) - boot_res$t0,
    `Value (Bias corrected)` = boot_res$t0 - (colMeans(boot_res$t) - boot_res$t0),
    `2.5%` = apply(boot_res$t, 2, quantile, 0.025),
    `97.5%` = apply(boot_res$t, 2, quantile, 0.975)
  )
  boot_estimates$Metric <- c("PTE", "R2hindiv","p-value Deviance","p-value LR","r")
  boot_estimates$N <- 2*N # untreated and treated
  boot_estimates$n_boot <- n_boot 
  boot_estimates$var_SEP <- var_SEP 
  boot_estimates$var_CEP <- var_CEP 
  
  boot_estimates <- boot_estimates %>%
    mutate(across(where(is.numeric), ~ format(round(.x, 2), scientific = FALSE)))
  
  
  return(boot_estimates)
  
}




var_CEP <- seq(0.5,2,0.5)
var_SEP <- c(0.5)

results <- NULL
pos <- 1


for(ii in 1:length(var_SEP))
{
  for(i in 1:length(var_CEP))
  {
    # Calculation of PTE; R2hindiv, and p-values for both Likelihood ratios and deviance difference.
    results[[pos]] <-  IT_vs_dev_LR(N = 1000,
                                    var_SEP = var_SEP[ii],
                                    var_CEP = var_CEP[i],
                                    boot_function,
                                    n_boot = 1000)
    
    pos <- pos + 1
  }
  
  print(paste0(ii,"/",length(var_SEP)))
}

results <- list.rbind(results)


results <- results %>%
  mutate(across(-Metric, as.numeric),
         var_SEP = paste0("Variance of surrogate endpoint: ",var_SEP))




# Creation of a plot
results %>%
  filter(Metric %in% c("PTE","R2hindiv","p-value LR")) %>%
  ggplot() +
  geom_errorbar(aes(x = var_CEP, ymin = `2.5%`, ymax = `97.5%`, color = Metric), width = 0.2) +
  geom_point(aes(x = var_CEP, y = `Value (Bias corrected)`, color = Metric)) +
  labs(title = expression("Effect of SEP/CEP variance on PTE and" ~ R[hindiv]^2),
       x = "Variance of the clinical enpoint",
       y = expression("PTE or"~R[hindiv]^2),
       color = "Metric", fill = "Metric") +
  coord_cartesian(ylim = c(0, 1.2), clip = "on") +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.1)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  facet_wrap(. ~ var_SEP, scales = "free_x")
