library(dplyr)
library(ggplot2)
library(ggfortify)
library(GGally)
library(Hmisc)
library(stats)
library(glmnet)

kidney = read.csv("kidney.csv")
dim(kidney)
names(kidney)

# Missing data
perc.missing = apply(apply(kidney, 2, is.na),2, mean)
kidney = kidney %>% select(-hbpstatus, -UPHO, -TRIG, -LDL, -HDL, -HBA1C)
kidney = na.omit(kidney)
dim(kidney)
attach(kidney)

# Convert to factors
kidney$BLACK = as.factor(BLACK) 
kidney$FEMALE = as.factor(FEMALE) 
kidney$Diabetes = as.factor(Diabetes)

# Look at variable relationships
ggcorr(kidney,label=TRUE)
ggplot(data=kidney) + geom_density(aes(x=GFR,fill=BLACK),alpha = 0.2)
ggplot(data=kidney) + geom_density(aes(x=GFR,fill=FEMALE),alpha = 0.2)
ggplot(data=kidney) + geom_density(aes(x=GFR,fill=Diabetes),alpha = 0.2)

# Split into test and train set
split = sample(1:nrow(kidney), size = floor(0.75*nrow(kidney)))
kidney.train = kidney[split,]
kidney.test = kidney[-split,]

# Step function for forward selection
lm.null = lm(GFR~1,data=kidney.train)
lm.full = lm(GFR~.,data=kidney.train)
forward.lm = step(lm.null,scope=list(lower=lm.null,upper=lm.full),direction="forward")
forward.lm$anova
summary(forward.lm)
forward.mse = mean((kidney.test$GFR - predict(forward.lm, kidney.test))^2)

ggplot()+geom_line(aes(x=1:12, y=forward.lm$anova$AIC))+ylab("AIC")+xlab("Number of Variables")

# Try out for backward selection and compare models and test MSE. 
backward.lm = step(lm.full,scope=list(lower=lm.null,upper=lm.full),direction="backward")
backward.lm$anova
summary(backward.lm)
backward.mse = mean((kidney.test$GFR - predict(backward.lm, kidney.test))^2)


# Lasso and ridge require data in matrix format - use model.matrix which transforms categorical variables
X=model.matrix(GFR~.,kidney)[,-1]
Y=kidney$GFR

# First look at output for ridge regression
grid = 10^seq(10, -2, length = 100)
ridge.lm = glmnet(X,Y,alpha=0,lambda=grid)
plot(ridge.lm)
ridge.lm$lambda[20]
coef(ridge.lm)[,20]
ridge.lm$lambda[50]
coef(ridge.lm)[,50]

# 10-fold cross-validation for ridge
folds = sample(1:10, nrow(kidney), replace = TRUE)
cv.mse = matrix(NA, nrow = 10, ncol = 100)

for (i in 1:10){
  ridge.lm=glmnet(X[folds!=i,],Y[folds!=i],alpha=0,lambda=grid)
  for (j in 1:100){
    cv.mse[i,j] = mean((predict(ridge.lm,s=grid[j],newx=X[folds==i,])-Y[folds==i])^2)
  }
}
avg.cv.mse = apply(cv.mse,2,mean)
sd.cv.mse = apply(cv.mse,2,sd)
ridge.lambda.min = grid[which.min(avg.cv.mse)]
ggplot()+geom_line(aes(x=log(ridge.lm$lambda),y=avg.cv.mse)) +
  geom_errorbar(aes(x=log(ridge.lm$lambda),ymin=avg.cv.mse-sd.cv.mse, ymax=avg.cv.mse+sd.cv.mse),col="gray") +
  xlab("Log(lambda)")+ylab("Avg MSE")+
  geom_vline(aes(xintercept=ridge.lambda.min),col="blue")

# Get final model chosen by k-fold cross-validation
# The coefficients are always returned on the original scale.
coef(ridge.lm)[,which.min(avg.cv.mse)]

# Look at lasso coefficients - difference to ridge plot
lasso.lm = glmnet(X,Y,alpha=1,lambda=grid)
plot(lasso.lm)

# Can use predict to get coefficients for particular lambda 
predict(lasso.lm,type="coefficients",s=15)[1:20,]

# Now run with lasso using cv.glmnet, plot and find lambda$min
cv.out=cv.glmnet(X,Y,alpha=1,nfolds=10,type.measure="mse")
plot(cv.out)
lasso.lambda.min=cv.out$lambda.min

predict(lasso.lm,type="coefficients",s=lasso.lambda.min)[1:20,]

