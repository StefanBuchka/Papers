# Purpose of this script: Function to transform Poisson distributed data to normality
# 25.11.2023
# Author: Stefan Buchka


library(MASS)

trans_func <- function(var = NULL,
                       bins = NULL,
                       trans = NULL)
{

# Method accoding to:
#  Huber W, Von Heydebreck A, Sültmann H, Poustka A, Vingron M. Variance stabilization applied to microarray data calibration and to the quantification of differential expression. Bioinformatics. 2002;18 suppl_1:S96–104.  
if(trans == "Huber")
{
  if(is.null(bins))
  {
    bins <- 2*IQR(var[which(!is.na(var) == T)])/length(var[which(!is.na(var) == T)])^(1/3)
    bins <- floor((max(var[which(!is.na(var) == T)]) - min(var[which(!is.na(var) == T)]))/bins)
  }
  
  df <- data.frame(var = var[which(!is.na(var) == T)])
  df <- df %>%
    mutate(bins = ntile(var, n = bins)) %>%
    group_by(bins) %>%
    summarise(mean = mean(var),
              sd = sd(var),
              n = n()) %>%
    as.data.frame()
  
  m <- lm(mean ~ sd, data = df)
  
  coef <- coef(m)
  var[which(!is.na(var) == T)] <- asinh(coef[1] + coef[2] * var[which(!is.na(var) == T)])
  
}
# Method according to:   
# Anscombe FJ. The transformation of Poisson, binomial and negative-binomial data. Biometrika. 1948;35:246–54.  
  if(trans == "Anscombe")
  {
    estimates <- fitdistr(var[which(!is.na(var) == T)],"Negative Binomial")
    
    r <- var[which(!is.na(var) == T)]
    c <- 0.09
    k <- estimates$estimate[1]
    var[which(!is.na(var) == T)] <- asinh(sqrt((r + c)/(k - 2*c)))
  }
  
  return(var)
}

