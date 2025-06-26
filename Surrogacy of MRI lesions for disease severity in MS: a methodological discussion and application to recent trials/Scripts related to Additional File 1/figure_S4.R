# Required packages
library(ggplot2)
library(dplyr)

# Wilson confidence interval function (suitable for proportions)
wilson_ci <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level)/2)
  p <- x / n
  denom <- 1 + z^2/n
  centre <- (p + z^2/(2*n)) / denom
  halfwidth <- z * sqrt((p*(1 - p) + z^2/(4*n)) / n) / denom
  lower <- pmax(0, centre - halfwidth)
  upper <- pmin(1, centre + halfwidth)
  return(list(lower = lower, centre = centre, upper = upper))
}

# Input values for sensitivity and specificity (reported in Table 1 of  Wattjes MP, Rovira À, Miller D, Yousry TA, Sormani MP, De Stefano N, et al. 
# MAGNIMS consensus guidelines on the use of MRI in multiple sclerosis--establishing disease prognosis and monitoring patients. 
# Nature Reviews Neurology. 2015;11(10):597–607) 

sens_list <- c(0.71, 0.24, 0.34, 0.68, 0.61)
spec_list <- c(0.77, 0.97, 0.90, 0.80, 0.83)
prev <- seq(0.01, 0.99, by = 0.01) # Potential T2 lesion prevalence
n <- 500  # assumed total sample size (impprtant for confidence bands)

# Empty dataframe
df_all <- data.frame()

# Loop over all combinations
for (i in 1:length(sens_list)) {
      
      # Define specificity and specificity
      sens <- sens_list[i]
      spec <- spec_list[i]

      # Contingency tables
      TP <- sens * prev * n
      FP <- (1 - spec) * (1 - prev) * n
      FN <- (1 - sens) * prev * n
      TN <- spec * (1 - prev) * n
      
      # PPV & NPV
      ppv_val <- TP / (TP + FP)
      npv_val <- TN / (TN + FN)
      
      # ppv_val = (sens * prev) / ((sens * prev) + ((1 - spec) * (1 - prev)))
      # npv_val = (spec * (1 - prev)) / ((spec * (1 - prev)) + ((1 - sens) * prev))
      
      # CIs
      ppv_ci <- wilson_ci(TP, TP + FP)
      npv_ci <- wilson_ci(TN, TN + FN)
      
      # Label for this combination
      label <- paste0("Sens=", sens, " | Spec=", spec)
      
      # Merge into temporary dataframe
      df_tmp <- data.frame(
        Prevalence = prev,
        Value = c(ppv_val, npv_val),
        Lower = c(ppv_ci$lower, npv_ci$lower),
        Upper = c(ppv_ci$upper, npv_ci$upper),
        Type = rep(c("PPV", "NPV"), each = length(prev)),
        Combo = label
      )
      
      df_all <- bind_rows(df_all, df_tmp)
      rm(df_tmp,sens,spec)
}

# Plot
df_all %>%
  ggplot(aes(x = Prevalence, y = Value, group = interaction(Type, Combo))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Type), alpha = 0.15, color = NA) +
  geom_line(aes(color = Type), size = 1) +
  facet_wrap(~ Combo) +
  scale_color_manual(values = c("PPV" = "red", "NPV" = "blue")) +
  scale_fill_manual(values = c("PPV" = "red", "NPV" = "blue")) +
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  labs(title = paste0("PPV & NPV over Prevalence (n = ",n,")"),
       x = "Prevalence", y = "Predictive Value", color = "", fill = "") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"))
