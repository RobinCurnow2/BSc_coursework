#Checks that all required packaged have been installed
packages <- c("GGally")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


envData <- read.table(file = "Environmental_data.txt", header = TRUE)
Y <- envData$ozone
X1 <- envData$radiation
X2 <- envData$temperature
X3 <- envData$wind

n <- length(Y)

#(A)
pairs(envData)

library(GGally)
ggpairs(envData)

#(B), (ii)
model <- lm(Y ~ X1 + X2 + X3)
X <- cbind(rep(1,n),X1, X2, X3)
XTXInv <- solve(t(X) %*% X)
XTY <- t(X) %*% Y
YTY <- t(Y) %*% Y

#(iii)
betaHat <- XTXInv %*% XTY
q = length(betaHat)

#(iv)
M = diag(n)-X%*%(XTXInv)%*%t(X)
SSR <- t(betaHat) %*% t(X) %*% Y
SST <- t(Y) %*% Y
SSE <- t(Y)%*%M%*%Y
SSTC <- SST-(sum(Y)^2)/n

SS <- c(SSE, SSR, SSTC)
DF <- c(q-1, n-q, n-1)
varSum <- data.frame(Source = c("SSE", "SSR", "SSTC"), 
                            SS = SS, DF = DF, 
                            meanSq= c(SS[1]/DF[1], SS[2]/DF[2], NA))
library(xtable)
anova_table <- xtable(varSum,caption = "Anova sumarry for air pollution data", align=c(rep("|c",length(names(varSum))),"|c|"))

#(v)
sigmaSq <- SSE/(n-q)
SE <- rep(0,q-1)
for (i in seq(2:q)){
  SE[i] <- sqrt(sigmaSq*XTXInv[i+1,i+1])
  }

tStat <- 0
for(i in seq(1:length(SE))){
  tStat[i] <- betaHat[i+1]/SE[i] 
}

pVals <- (1 - pt(abs(tStat), n-q))*2
  
FStat <- (varSum[3,2]-varSum[1,2])*varSum[2,3]/(varSum[1,2]*varSum[1,3])
pVal <- 1- pf(FStat, q-1, n-q)

#(vi)
Rsq <- (SSTC-SSE)/SSTC

#(vii)
#calculates the predicted values of ozone concentration fo the two days 
x0 <- c(1, 100, 70, 10)
x1 <- c(1, 50, 80, 10)
y0Hat <- t(x0) %*% betaHat
y1Hat <- t(x1) %*% betaHat
yDiffHat <- abs(y0Hat-y1Hat)
conf <- 0.975

epsilon_0 <- qt(conf, n-q)*sqrt((sigmaSq)*(1+t(x0) %*% XTXInv %*% x0))
epsilon_1 <- qt(conf, n-q)*sqrt((sigmaSq)*(1+t(x1) %*% XTXInv %*% x1))

CI_0 <- c(y0Hat - epsilon_0,y0Hat + epsilon_0)
CI_1 <- c(y1Hat - epsilon_1, y0Hat + epsilon_1)


#(viii))
residuals(model)
ggpairs(data.frame(res, Y, X1, X2, X3))

#(C)
#(ii)
BetaX2_int <- qt(conf, n-q-1)*sqrt((SSE*XTXInv[3,3])/(n-q))
X2_CI <- c(betaHat[3] - BetaX2_int, betaHat[3] + BetaX2_int)
X2_CI

X_True <- c(confint(model, level = 0.95)[3],confint(model, level = 0.95)[7])

#D
#(i)
model_int <- lm(Y~X1+ X2+ X3 + X1*X2 + X1*X3 + X2*X3)
model_int[["coefficients"]]

#(ii)
summary(model_int)
summary(model_int)$r.squared

#(iii)
confint(model_int, level = 0.95)
CI_table <- xtable(confint(model_int, level = 0.95), caption = 
                  "Confidence intervals for parameters beta_0-beta_6 in 
                  the updated regression model", align = c("|c", "|c", "|c|"))

#(v)
SSE_int <- sum(residuals(model_int)^2)
F <- (SSE- SSE_int)*(n-q)/(SSE_int*3)
pVal_comp <- 1-pf(F, 3,n-q)

