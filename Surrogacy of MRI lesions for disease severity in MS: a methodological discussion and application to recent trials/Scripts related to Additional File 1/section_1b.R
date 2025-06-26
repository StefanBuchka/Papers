

PPV.1 <- 0.3 * 0.71 / (0.3 * 0.71 + 0.7 * 0.29)
PPV.2 <- 0.2 * 0.61 / (0.2 * 0.61 + 0.8 * 0.17)
#
#
nn<-100
#
# Not treated
#
S.0 <- rnorm(nn,6,0.5)
C.0 <- rnorm(nn,0,0.5) + 1.5*S.0
#
# Treated
#
S.1 <- rnorm(nn,4,0.5)
C.1 <- rnorm(nn,0,0.5) + 1.5*S.1
#
SS<-c(S.0,S.1)
CC<-c(C.0,C.1)
#
library(ggplot2)

dat.1 <- data.frame(
  CE = CC,
  SE = SS,
  treat = factor(rep(c("Untreated", "Treated"), each = nn))
)

ggplot(dat.1, aes(x = SE, y = CE, color = treat)) +
  geom_point() +
  labs(x = "Surrogate", y = "Clinical Endpoint", color = "Group") +
  annotate("text", x = 6.5, y = 6, label = "Untreated", hjust = 0) +
  annotate("text", x = 3.5, y = 8, label = "Treated", hjust = 0) +
  theme_minimal()
#
#
#
dat.1<-data.frame(CE=CC,SE=SS,treat=rep(c(0,1),c(nn,nn)))
lm.0<-lm(SE~treat,data=dat.1)
cc.0<-coefficients(lm.0)
res.1.SE<-residuals(lm.0)

lm.1<-lm(CE~treat,data=dat.1)
summary(lm.1)
cc.1<-coefficients(lm.1)
res.1.CE<-residuals(lm.1)

lm.2<-lm(CE~SE,data=dat.1)     
summary(lm.2)

lm.3<-lm(CE~treat+SE,data=dat.1) 
summary(lm.3)
cc.3<-coefficients(lm.3)
PTE<-(cc.1[2]-cc.3[2])/cc.1[2]
PTE
#
# Meta-analytisch
#
cov(cbind(res.1.CE,res.1.SE))
cor(cbind(res.1.CE,res.1.SE))
#
# informations-theoretisch
#
glm.1<-glm(CE~treat,data=dat.1)
summary(glm.1)$aic
ll.1<-(-0.5)*(summary(glm.1)$aic+2)
glm.2<-glm(CE~treat+SE,data=dat.1)
ll.2<-(-0.5)*(summary(glm.2)$aic+4)
DDD<-(ll.2-ll.1)/(2*nn)
LRF<-1-exp(-2*DDD)
LRF
#
#
# Situation mit mehr Streuung
#
nn<-100
#
# Not treated
#
S.0 <- rnorm(nn,6,0.5)
C.0 <- rnorm(nn,0,2) + 1.5*S.0
#
# Treated
#
S.1 <- rnorm(nn,4,0.5)
C.1 <- rnorm(nn,0,2) + 1.5*S.1
#
SS<-c(S.0,S.1)
CC<-c(C.0,C.1)
#
dat.2 <- data.frame(
  CE = CC,
  SE = SS,
  treat = factor(rep(c("Untreated", "Treated"), each = nn))
)

ggplot(dat.2, aes(x = SE, y = CE, color = treat)) +
  geom_point() +
  labs(x = "Surrogate", y = "Clinical Endpoint", color = "Group") +
  annotate("text", x = 6.5, y = 2, label = "Untreated", hjust = 0) +
  annotate("text", x = 3.5, y = 12, label = "Treated", hjust = 0) +
  theme_minimal()
#
#
#
dat.2<-data.frame(CE=CC,SE=SS,treat=rep(c(0,1),c(nn,nn)))
lm.0<-lm(SE~treat,data=dat.2)
cc.0<-coefficients(lm.0)
res.1.SE<-residuals(lm.0)
lm.1<-lm(CE~treat,data=dat.2)
summary(lm.1)
cc.1<-coefficients(lm.1)
res.1.CE<-residuals(lm.1)
lm.2<-lm(CE~SE,data=dat.2)     
summary(lm.2)
lm.3<-lm(CE~treat+SE,data=dat.2) 
summary(lm.3)
cc.3<-coefficients(lm.3)
#
# Prentice
#
PTE<-(cc.1[2]-cc.3[2])/cc.1[2]
PTE
#
# meta-analytisch
#
cov(cbind(res.1.CE,res.1.SE))
cor(cbind(res.1.CE,res.1.SE))
#
# informations-theoretisch
#
glm.1<-glm(CE~treat,data=dat.2)
summary(glm.1)$aic
ll.1<-(-0.5)*(summary(glm.1)$aic+2)
glm.2<-glm(CE~treat+SE,data=dat.2)
ll.2<-(-0.5)*(summary(glm.2)$aic+4)
DDD<-(ll.2-ll.1)/(2*nn)
LRF<-1-exp(-2*DDD)
LRF

