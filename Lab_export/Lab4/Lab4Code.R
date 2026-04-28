library(dplyr)
library(ggplot2)
library(GGally)
library(ggfortify)
library(Hmisc)
library(lmtest)
library(car)
library(pROC)
library(boot)

kidney = read.csv("kidney_simple.csv")
names(kidney)
dim(kidney)

# No missing data
sum(!complete.cases(kidney))/nrow(kidney)

# Set black and female as factors
kidney$FEMALE = as.factor(kidney$FEMALE)
kidney$BLACK = as.factor(kidney$BLACK)
kidney$DISEASE = as.factor(kidney$DISEASE)
attach(kidney)

# Basic EDA
ggcorr(kidney,label=TRUE)
plot(DISEASE,SODIUM)
plot(DISEASE,WEIGHT)
plot(DISEASE,BMI)
plot(DISEASE,AGE)
plot(DISEASE,CYS)
table(DISEASE,FEMALE)
table(DISEASE,BLACK)

# New way to plot categorical vs continuous
ggplot(data=kidney) + geom_density(aes(x=SODIUM,fill=DISEASE),alpha = 0.2)
ggplot(data=kidney) + geom_density(aes(x=WEIGHT,fill=DISEASE),alpha = 0.2)
ggplot(data=kidney) + geom_density(aes(x=BMI,fill=DISEASE),alpha = 0.2)
ggplot(data=kidney) + geom_density(aes(x=AGE,fill=DISEASE),alpha = 0.2)
ggplot(data=kidney) + geom_density(aes(x=CYS,fill=DISEASE),alpha = 0.2)

# Fit a logistic regression model with all variables
glm.all = glm(DISEASE~.,data=kidney,family="binomial")
summary(glm.all)

# Comparing nested models
glm.1 = glm(DISEASE~BLACK+FEMALE+AGE+CYS,data=kidney,family="binomial")
lrtest(glm.1,glm.all)

# Variance inflation factors
vif(glm.all)

# Keep weight this time
glm.2 = glm(DISEASE~BLACK+FEMALE+AGE+CYS+WEIGHT,data=kidney,family="binomial")
lrtest(glm.2,glm.all)

# Plot using autoplot
autoplot(glm.2)

# We can predict probabilities or odds
glm.probs = predict(glm.2,type="response")
glm.odds = exp(predict(glm.2,type="link"))
glm.pred = (glm.probs > 0.5)
table(glm.pred,DISEASE)

# ROC curve and AUC
roc.2 = roc(DISEASE,glm.probs)
auc(roc.2)
ggroc(data=roc.2)

# Calibration
ordering = order(glm.probs)
avg.probs = numeric(10)
avg.disease = numeric(10)
for (i in 1:10){
  start = (i-1)*155+1
  end = min(1547,start+154)
  avg.probs[i] = mean(glm.probs[ordering[start:end]])
  avg.disease[i] = mean(DISEASE[ordering[start:end]]==TRUE)
}
ggplot()+geom_point(aes(x=avg.probs,y=avg.disease))+geom_abline(aes(slope=1,intercept=0),col="red")

# Pearson and deviance residuals - and standardized
pear.res = residuals(glm.2,type="pearson")
std.pear = pear.res/sqrt(1 - hatvalues(glm.2))
dev.res = residuals(glm.2,type="deviance")
std.dev = dev.res/sqrt(1-hatvalues(glm.2))

# Can still find outliers and high leverage points
c.dist = cooks.distance(glm.2)
summary(c.dist)

summary(std.pear)
kidney[std.pear > 2,]

# Use 10-fold cross-validation to compare model 1 and model 2
folds = sample(1:10, nrow(kidney), replace = TRUE)
cv.mse = matrix(NA, nrow = 10, ncol = 2)

for (i in 1:10){
  glm.1=glm(DISEASE~BLACK+FEMALE+AGE+CYS,family="binomial",data=kidney[folds!=i,])
  glm.2=glm(DISEASE~BLACK+FEMALE+AGE+CYS+WEIGHT,family="binomial",data=kidney[folds!=i,])
  y = as.numeric(DISEASE[folds == i])-1
  pred.1 = predict(glm.1,type="response",newdata=kidney[folds==i,])
  pred.2 = predict(glm.2,type="response",newdata=kidney[folds==i,])
  cv.mse[i,1] = mean((y-pred.1)^2)
  cv.mse[i,2] = mean((y-pred.2)^2)
}
avg.cv.mse = apply(cv.mse,2,mean)
sd.cv.mse = apply(cv.mse,2,sd)

# Bootstrap - two ways
intercept = function(data,indices){
  d = data[indices,]
  return(coef(glm(DISEASE~BLACK+FEMALE+AGE+CYS+WEIGHT,data=d,family="binomial"))[1])
}

boot(data=kidney,statistic=intercept,R=1000)

intercept = function(formula,data,indices){
  d = data[indices,]
  return(coef(glm(formula,data=d,family="binomial"))[1])
}

boot(data=kidney,statistic=intercept,R=1000,formula=DISEASE~BLACK+FEMALE+AGE+CYS+WEIGHT)
