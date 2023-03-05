#MATH38172 - GLM Coursework
#Robin Curnow, 10171119
#13/04/2022
install.packages("texreg")
library(texreg)

#a)
#sets working directory and reads excel file into R
setwd("~/GLM")
birthweight <- read.csv("birthweight.csv")

#shows the names and types of the variables in the data
str(birthweight)
#lists the different categories for Race used in the data set
levels(factor(birthweight$Race))
#ensures "White" is used as the baseline in the model
race1 <- factor(birthweight$Race, levels = c("White", "Black", "Other"))
#fits a logistic regression model to the data and outputs a summary of it to the console
bw1 <- glm(Low ~ Age+MotherWeight+race1+Smoke+Hypertension+Premature+UI+Visits,
           family=binomial, data = birthweight)
summary(bw1)
#creates LATEX syntax for a table showing the information in the summary table
texreg(bw1, digits=5, single.row = TRUE)

#c)
#fits a logistic regression model to the data and outputs a summary of it to the console
bw2 <- glm(Low ~ MotherWeight+race+Smoke+Hypertension, family=binomial, data = birthweight)
summary(bw2)
#creates LATEX syntax for a table showing the information in the summary table
texreg(bw2, digits=5, single.row = TRUE)

#e)
a <- coef(bw2)
eta <- c(1, 130, 1, 0, 0, 0)%*%a
mu <- exp(eta)/(1+exp(eta))
mu

#f)
p <- 0.2
#provides a point estimate of the weight using the MLEs
eta1 <- log(p/(1-p))
x_1 <- (eta1 - a[1])/a[2]
x_1
#defines the gradient function of the function of the regression parameters
gradh <- c((-1)/a[2],(a[1]-eta1)/a[2]^2, 0, 0, 0, 0)
x_1CI <- c(x_1) + c(-1,1)*qnorm(0.975)*c(sqrt(t(gradh)%*%(vcov(bw2)%*%gradh)))
x_1CI

#g
#defines a dummy variable for ’history of premature labour’ called prematureHist
prematureHist <- rep(0,nrow(birthweight))
#sets the data for a mother who has had one or more premature birth to 1 in PrematureHist and other cells to 0
for(i in c(1:nrow(birthweight))){
  if(birthweight$Premature[i] == 0) {
    prematureHist[i] <- 0
  } else {
    prematureHist[i] <- 1
  }
}
#checks that the prematureHist data matches the original data.
table(PrematureHist, birthweight$Premature)
#sets visits data a categorical variable
visits1 <- factor(birthweight$Visits)
#checks the different number of times that people have had physician visits (there are 6)
levels(visits1)

#defines a third glm with these new categorical factors, summaries it and creates LATEX
output for the summary
bw3 <- glm(Low ~ Age+MotherWeight+race1+Smoke+Hypertension+prematureHist+UI+visits1,
           family=binomial, data = birthweight)
summary(bw3)
texreg(bw3, digits=5, single.row = TRUE)
table(birthweight$Smoke, prematureHist)