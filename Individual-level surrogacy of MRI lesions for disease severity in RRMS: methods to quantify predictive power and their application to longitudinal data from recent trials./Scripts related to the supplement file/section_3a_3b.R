
# Load required library
library(dplyr)


####### Example 2x2 contingency table:
# Rows: SEP (0 or 1)
# Columns: CEP (0 or 1)

# Deginition of 2 by 2 table:
contingency_table <- matrix(c(1000, 500, 500, 1000),  # order: SEP=0&CEP=0, SEP=0&CEP=1, SEP=1&CEP=0, SEP=1&CEP=1
                            nrow = 2,
                            byrow = TRUE)
rownames(contingency_table) <- c("SEP = 0", "SEP = 1")
colnames(contingency_table) <- c("CEP = 0", "CEP = 1")

# Show the contingency table
print(contingency_table)

# Calculate total number of samples
N <- sum(contingency_table)

# Marginal probabilities
P_CEP <- colSums(contingency_table) / N  # P(CEP = 0), P(CEP = 1)
P_SEP <- rowSums(contingency_table) / N  # P(SEP = 0), P(SEP = 1)

# a) Shannon Entropy H(CEP)
H_CEP <- -sum(P_CEP * log(P_CEP))
cat("Shannon Entropy H(CEP):", round(H_CEP,2), "Nats\n")

# b) Conditional Entropy H(CEP | SEP)

# Calculate P(CEP = c | SEP = s) and then H(CEP | SEP = s) for each SEP
H_CEP_given_SEP <- 0
for (s in 1:2) {  # SEP=0 (row 1) and SEP=1 (row 2)
  P_CEP_given_SEP <- contingency_table[s, ] / sum(contingency_table[s, ])
  H_CEP_given_SEP_s <- -sum(P_CEP_given_SEP * log(P_CEP_given_SEP), na.rm = TRUE)  # handle log(0)
  H_CEP_given_SEP <- H_CEP_given_SEP + P_SEP[s] * H_CEP_given_SEP_s
}

cat("Conditional Entropy H(CEP|SEP):", round(H_CEP_given_SEP,2), "Nats\n")

# c) Mutual Information MI(CEP, SEP)
MI <- H_CEP - H_CEP_given_SEP
cat("Mutual Information MI(CEP, SEP):", MI, "Nats\n")

# d) R²bindiv (binary version of LRF - the likelihood reduction factor)
H_SEP <- -sum(P_SEP * log(P_SEP))
H_min <- min(H_SEP, H_CEP)
R2bindiv <- MI/H_min
cat("R²bindiv:", R2bindiv, "\n")

rm(list = ls()[which(ls() != "contingency_table")])


####### High and significant odds ratios, but low LRF

# MI and a logistic regression model
# Load required packages
if (!require("data.table")) install.packages("data.table")
library(data.table)
#
# ---- 1. Simulate Binary Risk Factors and Outcome ----
set.seed(123)
n <- 1000
X1 <- rbinom(n, 1, 0.5)
X2 <- rbinom(n, 1, 0.4)
X3 <- rbinom(n, 1, 0.3)

# Simulate binary outcome using a logistic model
logit <- -1 + 1.2 * X1 + 0.8 * X2 + 1.5 * X3
p <- 1 / (1 + exp(-logit))
Y <- rbinom(n, 1, p)

# ---- 2. Fit Logistic Regression Model ----
data <- data.table(Y, X1, X2, X3)
model <- glm(Y ~ X1 + X2 + X3, data = data, family = binomial())
summary(model)
data[, prob := predict(model, type = "response")]

# ---- 3. Prevalence of Risk Combinations ----
data[, risk_group := paste0("X1_", X1, "_X2_", X2, "_X3_", X3)]
risk_summary <- data[, .(count = .N,prevalence = mean(Y),
                         mean_pred_prob = mean(prob)), by = risk_group]
#
print(risk_summary)
#
# calculate MI for the whole setting
#
str(risk_summary)
#
P.x<-risk_summary$count/sum(risk_summary$count)
P.y<-c(sum(P.x*(1-risk_summary$mean_pred_prob)),sum(P.x*risk_summary$mean_pred_prob))
P.xy<-c(P.x*(1-risk_summary$mean_pred_prob),P.x*(risk_summary$mean_pred_prob))
P.xy.below<-c(P.x*P.y[1],P.x*P.y[2])
MI.log.reg<-sum(P.xy*log2(P.xy/P.xy.below))
#
1-exp(-2*MI.log.reg)
#






####### Estimate mutual information by logistic regression

# Simulate data

# 1.) concordant pairs
SEP <- CEP <- rep(0,contingency_table[1,1])

SEP <- c(SEP,rep(1,contingency_table[2,2]))
CEP <- c(CEP,rep(1,contingency_table[2,2]))

# 2) disconcordant pairs
SEP <- c(SEP,rep(1,contingency_table[2,1]))
SEP <- c(SEP,rep(0,contingency_table[1,2]))

CEP <- c(CEP,rep(0,contingency_table[2,1]))
CEP <- c(CEP,rep(1,contingency_table[1,2]))


# 1.) concordant pairs
SEP <- CEP <- rep(0,contingency_table[1,1])

SEP <- c(SEP,rep(1,contingency_table[2,2]))
CEP <- c(CEP,rep(1,contingency_table[2,2]))

# 2) disconcordant pairs
SEP <- c(SEP,rep(1,contingency_table[2,1]))
SEP <- c(SEP,rep(0,contingency_table[1,2]))

CEP <- c(CEP,rep(0,contingency_table[2,1]))
CEP <- c(CEP,rep(1,contingency_table[1,2]))


# Assumption: CEP and SEP are binary variables (0/1)

# Step 1: Fit logistic regression models
# Null model: only intercept (used to estimate H(CEP))
model_CEP <- glm(CEP ~ 1, family = binomial)

# Full model: CEP predicted by SEP (used to estimate H(CEP | SEP))
model_CEP_SEP <- glm(CEP ~ SEP, family = binomial)

# SEP model: for H(SEP) computation
model_SEP <- glm(SEP ~ 1, family = binomial)

N <- length(SEP)
H_CEP <- -logLik(model_CEP)[1]/N
H_SEP <- -logLik(model_SEP)[1]/N
H_CEP_given_SEP <- -logLik(model_CEP_SEP)[1]/N

MI <- H_CEP - H_CEP_given_SEP
R2bindiv <- MI/min(H_CEP,H_SEP)


# Output results
cat("Mutual Information (MI):", MI, "\n")
cat("R²bindiv (individual-level surrogacy):", R2bindiv, "\n")
