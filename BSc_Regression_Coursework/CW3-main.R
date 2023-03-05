rm(list=ls())   #clears any pre-existing variables and data
install.packages("readr", "xtable","bivariate")   #ensures all neccesarry packages are install
data <- sweden_ins_data
x <- sweden_ins_data$claims
y <- sweden_ins_data$payment

#1
plot(x, y, xlab = "claims", ylab="payment", main="claims vs payment for insurance claims in Sweden",cex.main=0.8)

#2
OLS = function(xvals, yvals){
  #defines the vectors required to complete the question as zero vectors
  a<-0
  b<-0
  yhat<-0
  residual<-0
  #uses a loop to find the sums required to calculate estimates of alpha and beta from the data
  for (i in 1:length(xvals)) {
    a <- a+((xvals[i]-mean(xvals))*yvals[i])
    b <- b+(xvals[i]-mean(xvals))^2
  }
  betahat <- a/b
  alphahat <- mean(y)-betahat*mean(x)
  #uses values of alpha and betahat to calculate the fitted values and errors 
  for (i in 1:length(xvals)) {
    yhat[i]<-alphahat+betahat*x[i]
    residual[i] <- yvals[i]-yhat[i]
  }
  #uses the number of points in the data to calculate the degrees of freedom
  dof <- length(xvals)-2
  #puts the answers into a list, assigns names to each item in the list and returns these to the console
  answers = list(c(alphahat,betahat),yhat, residual, dof)
  names(answers)<-(c("alpha/betahat","fitted values","residuals","degrees of Freedom"))
  return(answers)
}

#calls the OLS function above with the insurance data as input. Stores the returned values as a list named "ans" 
ans<-OLS(x,y)

#3
xx <- seq(min(x),max(x),length.out = 1000)
yy <- (ans[["alpha/betahat"]][1]+ans[["alpha/betahat"]][2]*xx)
lines(xx,yy,col="blue")

SSE <- sum(ans[["residuals"]]^2)
SSR <- sum((ans[["fitted values"]]-mean(y))^2)
sigmahatsq <- SSE/ans[["degrees of Freedom"]]
F=SSR/(SSE/(length(x)-2))
F-qf(0.95,1,length(x))

alphahat <- ans[["alpha/betahat"]][1]
betahat <- ans[["alpha/betahat"]][2]

x1=80
est <- alphahat+betahat*x1
Sxx <-sum((x-mean(x))^2)
error <- qt(0.975,length(x)-2)*sqrt(sigmahatsq*((1/length(x))+((mean(x)-x1)^2)/(Sxx)))
ci_upper <- est+error
ci_lower <- est-error
CI <- c(ci_lower,ci_upper)
CI

#4
#a graphics window that fits 2 plots side by side
par(mfrow=c(1,2))
#vector containing the standardized residuals of the data
sdres <- rstandard(lm(y~x))
#plots of fitted values and claims against the standardised residuals
plot(ans[["fitted values"]], sdres, xlab = "fitted values", ylab = "standardised residuals",
     main = "Fitted values against standardised residuals")
abline(h=0)
abline(h=c(2,-2),col="red",lty=2)
plot(x, sdres, xlab = "claims", ylab = "standardised residuals"
     ,main = "claims against standardised residuals")
abline(h=0)
abline(h=c(2,-2),col="red",lty=2)

#5
#sets plot window to contain only one plot
par(mfrow=c(1,1))
#orders the standard residuals
Q1 <-sort(sdres)
#finds the quantile estimates corresponding to the order statistics
k<-seq(1:length(x))
p_k = (k-3/8)/(length(x)+1/4)
Q2 <- qnorm(p_k)

#plots Q1 against Q2 
plot(Q2, Q1, pch=20, cex=0.5,xlab="standard normal quantiles",
ylab = "standardised residuals", main="Normal quantile-quantile plot of the standardised residuals",
cex.main=0.75, cex.lab=0.6)
abline(a=mean(sdres), b=sd(sdres), col="red")

#6
#kstest <- ks.test(sdres, y=pnorm, alternative = "two.sided")

kstestx.ecdf<-(1:length(x))/length(x)
#x.ecdf

yyy<-pnorm(sort(sdres))
#y

diff1=x.ecdf-yyy
md1=max(diff1)
dn.plus=max(md1,0)
x.KS1=sort(sdres)[dn.plus==diff1]

diff2=yyy-x.ecdf
md2=max(diff2)
md2=md2+(1/length(x))
dn.minus=max(md2,0)
x.KS2=sort(sdres)[dn.minus==diff2+(1/length(x))]

KSstat=max(dn.minus, dn.plus)
KSstat

if(dn.minus < dn.plus) x.KSstat=x.KS1
if(dn.minus > dn.plus) x.KSstat=x.KS2
x.KSstat

#7
#plots empirical and expected cdfs  
plot(sort(sdres),x.ecdf, type="s",cex.lab=0., xlab="standard residuals",ylab = 
       "F(x)")
    segments(x0=min(sdres),y0=0, x1=min(sdres), y1=x.ecdf[1])
    segments(x0=max(sdres),y0=1, x1=1000, y1=1)
    segments(-1000,y0=0, x1=min(sdres), y1=0)

xxx <- seq(from = -3, to = 3, length.out = 1000)    
lines(xxx, pnorm(xxx), col="red")    
abline(v=x.KSstat, col="blue", lty=2)
title(cex.main=0.75,main="Empirical and true cdf for standardised residuals from insurance data",
      sub=paste("KS-Stat = ", as.character(round(KSstat, digits=3)),sep=""),adj=1)

#8
#Function to simulate the sampling distribution of the Kolmogorov-Smirnov test 
#statistic from a specified number random N(0,1) observations

simKS = function(nsim,n){
  KS.stat = 0
  for (i in 1:nsim) {
    xdat = rnorm(n)
    KS <- ks.test(xdat, y=pnorm, alternative = "two.sided")
    KS.stat[i] = KS$statistic
  }
  return(KS.stat)
}


#Histogram showing the distribution and kernal density estimate of this
KS.stat = simKS(10000,length(x))
hist(KS.stat, freq=FALSE, breaks = 21, main = paste("Histogram of 
     estimated sampling distribution of K-S test statistic from n= ",length(x),
     "N(0,1) observations"), xlab = "K-S test statistic", 
     xlim=c(0,0.3), ylim=c(0,14), col="lightblue", cex.main=0.9)
lines(density(KS.stat), col="magenta", lwd=1.2)

#D and p-value for the standardised residuals of Swedish insurance claims data
ks.test(sdres, y=pnorm, alternative = "two.sided")
#5% critical value for K-S statistic
quantile(KS.stat, prob=0.95)

