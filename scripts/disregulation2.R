library('data.table')
library('ggplot2')
library('glmnet')

library('plyr')
library('readr')
library('caret')
library('repr')
library('dplyr')
# Disregulation score for individauls

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
#chroms=1:10
env="ALL"
pheno="grain_yield_15"

totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)

beta_z=totalrares%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]


x=data.frame(id=unique(beta_z$ID),stringsAsFactors=F)
for(time1 in times){
	tbetaz=beta_z[beta_z$time==time1,]
	di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
	di=as.data.frame(di,stringsAsFactors=F)
	x[,time1]=di[match(x$id,di$ID),]$di_score
}
inds=x$id
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code


y=phenotypes[inds,pheno,drop=F]
y=y[!is.na(y[,pheno]),,drop=F]


newinds=rownames(y)

rownames(x)=x$id
x=x[,-1]
x=x[newinds,]

#dat=cbind(y,x)
impx=makeX(x,na.impute=TRUE)
dat=cbind(y,impx)
y=unlist(y)

fit=glmnet(impx,y,alpha=0)

png(sprintf('QTT/%s_%s_disregulation_glmnet.png',pheno,env))
print(plot(fit,label=TRUE))
dev.off()

### cross-validation

set.seed(100) 

# 70% of data

index = sample(1:nrow(dat), 0.7*nrow(dat)) 

train = dat[index,] # Create the training data 
test = dat[-index,] # Create the test data

dim(train)
dim(test)

cols = times

pre_proc_val <- preProcess(train[,cols], method = c("center", "scale"))

train[,cols] = predict(pre_proc_val, train[,cols])
test[,cols] = predict(pre_proc_val, test[,cols])

summary(train)

### Regularization
cols_reg = c(times,'grain_yield_15')

dummies <- dummyVars(grain_yield_15 ~ ., data = dat[,cols_reg])

train_dummies = predict(dummies, newdata = train[,cols_reg])

test_dummies = predict(dummies, newdata = test[,cols_reg])

print(dim(train_dummies)); print(dim(test_dummies))

### Ridge Regression

x = as.matrix(train_dummies)
y_train = train$grain_yield_15

x_test = as.matrix(test_dummies)
y_test = test$grain_yield_15

lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)

summary(ridge_reg)


### Optimal lambda value

cv_ridge <- cv.glmnet(x, y_train, alpha = 0, lambda = lambdas)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

# 3.981072


#### Predict values #####

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))

  
  # Model performance metrics
data.frame(
  RMSE = RMSE,
  Rsquare = R_square
)
  
}

# Prediction and evaluation on train data
predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = x)
eval_results(y_train, predictions_train, train)

# Prediction and evaluation on test data
predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
eval_results(y_test, predictions_test, test)

#cvfit <- cv.glmnet(impx, y,nfolds=10)

#pdf(sprintf('QTT/%s_%s_disregulation_glmnet_CV.pdf',pheno,env))
#print(plot(cvfit))
#dev.off()