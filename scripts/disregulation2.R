library('data.table')
library('ggplot2')
library('glmnet')

library('lme4')
library('lmerTest')
library('plyr')
library('readr')
library('caret')
library('repr')
library('dplyr')
# Disregulation score for individauls
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes %>% group_by(Loc.Year.Treat) %>% reframe(r=cor(tkw_15,grain_yield_15,use="complete.obs"))
# A tibble: 8 × 2
#  Loc.Year.Treat         r
#  <chr>              <dbl>
#1 ALL                0.240
#2 BLOIS_2014_OPT     0.327
#3 BLOIS_2017_OPT     0.187
#4 EXP_STPAUL_2017_WD 0.325
#5 GRANEROS_2015_OPT  0.190
#6 NERAC_2016_WD      0.111
#7 STPAUL_2017_WD     0.332
#8 SZEGED_2017_OPT    0.187


times=c("WD_0712","WD_0718","WD_0720","WD_0727")
#chroms=1:10



log_inverse=function(x){
  	return(2^x)
}

eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - (SSE / SST)
  RMSE = sqrt(SSE/nrow(df))

  
  # Model performance metrics
data.frame(
  RMSE = RMSE,
  Rsquare = R_square
)
  
}

########## 
env="EXP_STPAUL_2017_WD"
env="ALL"
# thousand kernel weight
pheno="tkw_15"

####### Di froma avg total expression ######
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code
inds=rownames(phenotypes)
y=phenotypes[inds,pheno,drop=F]
y=y[!is.na(y[,pheno]),,drop=F]
inds=rownames(y)


dat=data.frame(y=y,stringsAsFactors=F)
names(dat)=pheno
for(time1 in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	rownames(unlog)=rownames(exp)
	avg_exp = apply(unlog,MARGIN=2,mean)
	avg_exp[avg_exp<1]=0
	names(avg_exp)=names(exp)
	avg_exp=sort(avg_exp,decreasing=TRUE)
	avg_exp=avg_exp[1:5000] # grap top 5000 highest expressed genes
	
	unlog=unlog[,names(avg_exp)]
	# re-log the input data
	#avg_logexp=log2(avg_exp)
	#avg_logexp[is.infinite(avg_logexp)]=0
	# For each individual, get the average squared diff from the avg expression
	# Unlog expression for both, get the average difference, and then relog
	dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
	rownames(dev)=rownames(unlog)
	colnames(dev)=colnames(unlog)
	# relog log2cpm
	di=apply(dev,MARGIN=1,function(x) log10(mean(x,na.rm=T)))
	dat[,time1]=di[match(rownames(dat),names(di))]
}


###### Linear Model #######



cor(dat$tkw_15,dat$WD_0712,use="complete.obs")
# 0.08152518 tkw_15 ALL
# 0.05109589  tkw_15 EXP_STPAUL_2017_WD
cor(dat$tkw_15,dat$WD_0718,use="complete.obs")
# 0.09604271 tkw_15 ALL
# 0.04683859 tkw_15 EXP_STPAUL_2017_WD
cor(dat$tkw_15,dat$WD_0720,use="complete.obs")
# 0.08378463 tkw_15 ALL
# 0.09644231 tkw_15 EXP_STPAUL_2017_WD
cor(dat$tkw_15,dat$WD_0727,use="complete.obs")
# 0.1317172 tkw_15 ALL
# 0.1200785 tkw_15 EXP_STPAUL_2017_WD

mod12=lm(tkw_15 ~ WD_0712,dat)
anova(mod12)
#Analysis of Variance Table

#Response: tkw_15 ALL
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1   123.3  123.28   0.542 0.4637
#Residuals 81 18424.4  227.46 

#Response: tkw_15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1    76.4   76.38  0.1885 0.6655
#Residuals 72 29177.6  405.24 

mod18=lm(tkw_15 ~ WD_0718,dat)
anova(mod18)
#Response: tkw_15 ALL
#           Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1   278.2  278.19  1.3127 0.2538
#Residuals 141 29881.1  211.92

#Response: tkw_15 EXP_STPAUL_2017_WD
#            Df Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1    115  115.33  0.2792 0.5981
#Residuals 127  52454  413.03 

mod20=lm(tkw_15 ~ WD_0720,dat)
anova(mod20)
#Response: tkw_15 ALL
#           Df Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1    328  327.81  1.5411 0.2158
#Residuals 218  46370  212.71


#Response: tkw_15 EXP_STPAUL_2017_WD
#            Df Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1    693  692.60  1.8026  0.181
#Residuals 192  73771  384.23

mod27=lm(tkw_15 ~ WD_0727,dat)
anova(mod27)
#Response: tkw_15 ALL
#           Df Sum Sq Mean Sq F value  Pr(>F)  
#WD_0727     1    671  670.83  3.3899 0.06714 .
#Residuals 192  37995  197.89   


#Response: tkw_15 EXP_STPAUL_2017_WD
#            Df Sum Sq Mean Sq F value Pr(>F)
#WD_0727     1    974  974.32  2.5602 0.1114
#Residuals 175  66598  380.56
dat$ID=rownames(dat)
datmelt=reshape2::melt(dat[,c('tkw_15','ID',times)],c('ID','tkw_15'))
datmelt=datmelt[complete.cases(datmelt),]
datmelt$ID_time=paste0(datmelt$ID,'-',datmelt$variable)

datmelt2 = datmelt %>% group_by(ID) %>% reframe(mean_di=mean(value,na.rm=T),tkw_15=unique(tkw_15))

full_model=lm(tkw_15 ~ variable + value,datmelt)

full_model2=lmer(tkw_15 ~ (1|ID) + variable + value,datmelt)
anova(full_model2,ddf='K')
#Type III Analysis of Variance Table with Kenward-Roger's method
#      Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#value 1238.9  1238.9     1 377.9  5.9526 0.01515 *
#
p1=ggplot(aes(x=value,y=tkw_15),data=datmelt) + geom_point(aes(color=variable)) 

png(sprintf('QTT/%s_%s_disregulation_score_avg_totexp_plot.png',pheno,env))
print(p1)
dev.off()



###### glmnet #######
x=dat[,times]
# drop individuals with fewer than 2 timepoints
xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
x=x[xn<3,]
dim(x)
impx=makeX(x,na.impute=TRUE)
y=y[rownames(x),,drop=F]
dat=cbind(y,impx)
y=unlist(y)

fit=glmnet(impx,y,alpha=0)

png(sprintf('QTT/%s_%s_avg_totexp_disregulation_glmnet.png',pheno,env))
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
cols_reg = c(times,'tkw_15')

dummies <- dummyVars(tkw_15 ~ ., data = dat[,cols_reg])

train_dummies = predict(dummies, newdata = train[,cols_reg])

test_dummies = predict(dummies, newdata = test[,cols_reg])

print(dim(train_dummies)); print(dim(test_dummies))

### Ridge Regression

x = as.matrix(train_dummies)
y_train = train$tkw_15

x_test = as.matrix(test_dummies)
y_test = test$tkw_15

#lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian')

summary(ridge_reg)

### Optimal lambda value

cv_ridge <- cv.glmnet(x, y_train, alpha = 0)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

# 63.24908 tkw_15 ALL
#17.69236 tkw_15 EXP_STPAUL_2017_WD

# Prediction and evaluation on train data
predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = x)
eval_results(y_train, predictions_train, train)

# Prediction and evaluation on test data
predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
eval_results(y_test, predictions_test, test)

# tkw_15 ALL
#      RMSE    Rsquare
#1      RMSE    Rsquare
#1 14.26912 0.01102374
#      RMSE    Rsquare
#1 14.54346 0.00870857

  
# tkw_15 EXP_STPAUL_2017_WD
#       RMSE    Rsquare
#1 20.02741 0.04046428
#      RMSE    Rsquare
#1 19.53892 -0.0687139

###### Di from average zscore (from top 5k genes) #######
totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	#allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]

# switch from summarize to reframe
totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]

# Add beta_z score to totalrares
totalrares$beta_z=beta_z[match(totalrares$gene_time_founder,beta_z$gene_time_founder),]$beta_z
# 
x=data.frame(id=unique(totalrares$ID),stringsAsFactors=F)
for(time1 in times){
	tbetaz=totalrares[totalrares$time==time1,]
	# di is average disregulation
	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
	# di is sum of disregulation **2
	di=tbetaz %>% group_by(ID) %>% reframe(di_score=mean(beta_z^2,na.rm=T))

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

dat=cbind(y,x)


#datmelt=reshape2::melt(dat[,c(pheno,times)],pheno)
dat$ID=rownames(dat)
datmelt=reshape2::melt(dat[,c('tkw_15','ID',times)],c('ID','tkw_15'))
datmelt=datmelt[complete.cases(datmelt),]
datmelt$ID_time=paste0(datmelt$ID,'-',datmelt$variable)

datmelt2 = datmelt %>% group_by(ID) %>% reframe(mean_di=mean(value,na.rm=T),tkw_15=unique(tkw_15))

full_model=lm(tkw_15 ~ mean_di,datmelt2)

full_model3=lm(tkw_15 ~ variable + value,datmelt)
anova(full_model3)

full_model2=lmer(tkw_15 ~ (1|variable) + value,datmelt)
anova(full_model2,ddf='K')




p1=ggplot(aes(x=value,y=tkw_15),data=datmelt) + geom_point(aes(color=variable)) 

png(sprintf('QTT/%s_%s_disregulation_avg_zscore_plot.png',pheno,env))
print(p1)
dev.off()


cor(dat$tkw_15,dat$WD_0712,use="complete.obs")
# -0.00908917tkw_15 ALL
# -0.05465077 tkw_15 EXP_STPAUL_2017_WD
cor(dat$tkw_15,dat$WD_0718,use="complete.obs")
# -0.02008408 tkw_15 ALL
# -0.07213502 tkw_15 EXP_STPAUL_2017_WD
cor(dat$tkw_15,dat$WD_0720,use="complete.obs")
# -0.1286326tkw_15 ALL
# -0.1044709 tkw_15 EXP_STPAUL_2017_WD
cor(dat$tkw_15,dat$WD_0727,use="complete.obs")
# -0.02280601 tkw_15 ALL
# -0.02885656 tkw_15 EXP_STPAUL_2017_WD

###### Linear Model #######
# WD_0712
mod12=lm(tkw_15 ~ WD_0712,dat)
anova(mod12)
#Analysis of Variance Table

#Response: tkw_15 ALL
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1     1.5   1.532  0.0067  0.935
#Residuals 81 18546.2 228.965

#Analysis of Variance Table

#Response: tkw_15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1    87.4   87.37  0.2157 0.6437
#Residuals 72 29166.6  405.09 

mod18=lm(tkw_15 ~ WD_0718,dat)
anova(mod18)
#Response: tkw_15 ALL
#           Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1    12.2  12.165  0.0569 0.8118
#Residuals 141 30147.1 213.809

#Response: tkw_15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1    274  273.54  0.6643 0.4166
#Residuals 127  52296  411.78  

mod20=lm(tkw_15 ~ WD_0720,dat)
anova(mod20)
#Response: tkw_15 ALL
#           Df Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1    773  772.68  3.6678 0.05678 .
#Residuals 218  45925  210.67

#Response: tkw_15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1    813  812.71  2.1186 0.1471
#Residuals 192  73651  383.60

mod27=lm(tkw_15 ~ WD_0727,dat)
anova(mod27)
#Response: tkw_15 ALL
#           Df Sum Sq Mean Sq F value Pr(>F)
#WD_0727     1     20  20.111  0.0999 0.7523
#Residuals 192  38646 201.281  

#Response: tkw_15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0727     1     56   56.27  0.1458  0.703
#Residuals 175  67516  385.81

###### glmnet #######

xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
x=x[xn<3,]
dim(x)

impx=makeX(x,na.impute=TRUE)
y=y[rownames(x),,drop=F]
dat=cbind(y,impx)
y=unlist(y)

fit=glmnet(impx,y,alpha=0)

png(sprintf('QTT/%s_%s_avg_zscore_disregulation_glmnet.png',pheno,env))
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
cols_reg = c(times,'tkw_15')

dummies <- dummyVars(tkw_15 ~ ., data = dat[,cols_reg])

train_dummies = predict(dummies, newdata = train[,cols_reg])

test_dummies = predict(dummies, newdata = test[,cols_reg])

print(dim(train_dummies)); print(dim(test_dummies))

### Ridge Regression

x = as.matrix(train_dummies)
y_train = train$tkw_15

x_test = as.matrix(test_dummies)
y_test = test$tkw_15

#lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian')

summary(ridge_reg)


### Optimal lambda value

cv_ridge <- cv.glmnet(x, y_train, alpha = 0)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

# 408.3829tkw_15 ALL
# 1210.581 tkw_15 EXP_STPAUL_2017_WD




# Prediction and evaluation on train data
predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = x)
eval_results(y_train, predictions_train, train)

# Prediction and evaluation on test data
predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
eval_results(y_test, predictions_test, test)

# tkw_15 ALL
#       RMSE     Rsquare
#1 14.46558 0.002165917
#      RMSE    Rsquare
#1 14.35169 -0.0165037

# tkw_15 EXP_STPAUL_2017_WD
#        RMSE      Rsquare
#1 19.53116 0.0004676685
#      RMSE      Rsquare
#1 21.05419 9.822137e-05


#################### #################### 
#################### #################### 
#################### #################### 
#################### #################### 
#################### #################### 
env="EXP_STPAUL_2017_WD"
env="ALL"
# grain yield
pheno="grain_yield_15"


log_inverse=function(x){
  	return(2^x)
}

eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - (SSE / SST)
  RMSE = sqrt(SSE/nrow(df))

  
  # Model performance metrics
data.frame(
  RMSE = RMSE,
  Rsquare = R_square
)
  
}

####### Di froma avg total expression ######
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code
inds=rownames(phenotypes)
y=phenotypes[inds,pheno,drop=F]
y=y[!is.na(y[,pheno]),,drop=F]
inds=rownames(y)


dat=data.frame(y=y,stringsAsFactors=F)
names(dat)=pheno
for(time1 in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	rownames(unlog)=rownames(exp)
	avg_exp = apply(unlog,MARGIN=2,mean)
	avg_exp[avg_exp<1]=0
	names(avg_exp)=names(exp)
	avg_exp=sort(avg_exp,decreasing=TRUE)
	avg_exp=avg_exp[1:5000] # grap top 5000 highest expressed genes
	
	unlog=unlog[,names(avg_exp)]
	# re-log the input data
	#avg_logexp=log2(avg_exp)
	#avg_logexp[is.infinite(avg_logexp)]=0
	# For each individual, get the average squared diff from the avg expression
	# Unlog expression for both, get the average difference, and then relog
	dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
	rownames(dev)=rownames(unlog)
	colnames(dev)=colnames(unlog)
	# relog log2cpm
	di=apply(dev,MARGIN=1,function(x) log10(mean(x,na.rm=T)))
	dat[,time1]=di[match(rownames(dat),names(di))]
}


###### Linear Model #######



cor(dat$grain_yield_15,dat$WD_0712,use="complete.obs")
# 0.06166244 grain_yield_15 ALL
# 0.05956349  grain_yield_15 EXP_STPAUL_2017_WD
cor(dat$grain_yield_15,dat$WD_0718,use="complete.obs")
# 0.04842427 grain_yield_15 ALL
# 0.05305196 grain_yield_15 EXP_STPAUL_2017_WD
cor(dat$grain_yield_15,dat$WD_0720,use="complete.obs")
# -0.01901306 grain_yield_15 ALL
# -0.05234582 grain_yield_15 EXP_STPAUL_2017_WD
cor(dat$grain_yield_15,dat$WD_0727,use="complete.obs")
# 0.03756198 grain_yield_15 ALL
#  -0.004692406 grain_yield_15 EXP_STPAUL_2017_WD

mod12=lm(grain_yield_15 ~ WD_0712,dat)
anova(mod12)
#Analysis of Variance Table

#Response: grain_yield_15 ALL
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1   13.4  13.380  0.3092 0.5797
#Residuals 81 3505.5  43.278

#Response: grain_yield_15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1     37  37.019  0.2564 0.6142
#Residuals 72  10397 144.405 

mod18=lm(grain_yield_15 ~ WD_0718,dat)
anova(mod18)
#Response: grain_yield_15 ALL
#           Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1     10  10.001  0.3314 0.5657
#Residuals 141   4255  30.178 

#Response: grain_yield_15 EXP_STPAUL_2017_WD
#            Df Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1    37.4  37.412  0.3556  0.552
#Residuals 126 13255.0 105.198

mod20=lm(grain_yield_15 ~ WD_0720,dat)
anova(mod20)
#Response: grain_yield_15 ALL
#           Df Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1    2.6   2.612  0.0788 0.7791
#Residuals 218 7222.8  33.132


#Response: grain_yield_15 EXP_STPAUL_2017_WD
#            Df Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1     62  61.999  0.5248 0.4697
#Residuals 191  22565 118.139

mod27=lm(grain_yield_15 ~ WD_0727,dat)
anova(mod27)
#Response: grain_yield_15 ALL
#           Df Sum Sq Mean Sq F value  Pr(>F)  
#WWD_0727     1    9.1   9.082  0.2713 0.6031
#Residuals 192 6427.9  33.478 


#Response: grain_yield_15 EXP_STPAUL_2017_WD
#            Df Sum Sq Mean Sq F value Pr(>F)
#WD_0727     1     0.5   0.473  0.0038 0.9507
#Residuals 174 21459.1 123.328
datmelt=reshape2::melt(dat[,c('grain_yield_15',times)],'grain_yield_15')

p1=ggplot(aes(x=value,y=grain_yield_15),data=datmelt) + geom_point(aes(color=variable)) 

png(sprintf('QTT/%s_%s_disregulation_score_avg_totexp_plot.png',pheno,env))
print(p1)
dev.off()

###### glmnet #######


x=dat[,times]
xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
x=x[xn<3,]
dim(x)

impx=makeX(x,na.impute=TRUE)
y=y[rownames(x),,drop=F]
dat=cbind(y,impx)
y=unlist(y)

fit=glmnet(impx,y,alpha=0)

png(sprintf('QTT/%s_%s_avg_totexp_disregulation_glmnet.png',pheno,env))
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

#lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian')

summary(ridge_reg)

### Optimal lambda value

cv_ridge <- cv.glmnet(x, y_train, alpha = 0)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

# 213.1109 grain_yield_15 ALL
#1216.665 grain_yield_15 EXP_STPAUL_2017_WD

# Prediction and evaluation on train data
predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = x)
eval_results(y_train, predictions_train, train)

# Prediction and evaluation on test data
predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
eval_results(y_test, predictions_test, test)

# grain_yield_15 ALL
#      RMSE Rsquare
#1 5.993931       0
#      RMSE     Rsquare
#1 5.103489 -0.09136803

  
# grain_yield_15 EXP_STPAUL_2017_WD
#     RMSE      Rsquare
#1 10.69589 0.0002061998
#      RMSE       Rsquare
#1 11.26273 -0.0003417997

###### Di from average zscore (from top 5k genes) #######
totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	#allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]

# switch from summarize to reframe
totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),gene_time=unique(gene_time),max_f=unique(max_f))

beta_z=totf%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
beta_z=as.data.frame(beta_z,stringsAsFactors=F)
beta_z=beta_z[!is.na(beta_z$beta_z),]

# Add beta_z score to totalrares
totalrares$beta_z=beta_z[match(totalrares$gene_time_founder,beta_z$gene_time_founder),]$beta_z
# 
x=data.frame(id=unique(totalrares$ID),stringsAsFactors=F)
for(time1 in times){
	tbetaz=totalrares[totalrares$time==time1,]
	# di is average disregulation
	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
	# di is sum of disregulation **2
	di=tbetaz %>% group_by(ID) %>% reframe(di_score=mean(beta_z^2,na.rm=T))

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

dat=cbind(y,x)


datmelt=reshape2::melt(dat[,c(pheno,times)],pheno)

p1=ggplot(aes(x=value,y=grain_yield_15),data=datmelt) + geom_point(aes(color=variable)) 

png(sprintf('QTT/%s_%s_disregulation_avg_zscore_plot.png',pheno,env))
print(p1)
dev.off()


cor(dat$grain_yield_15,dat$WD_0712,use="complete.obs")
# 0.1262777 grain_yield_15 ALL
# 0.005983731 grain_yield_15 EXP_STPAUL_2017_WD
cor(dat$grain_yield_15,dat$WD_0718,use="complete.obs")
# -0.01878229 grain_yield_15 ALL
# 0.03289901 grain_yield_15 EXP_STPAUL_2017_WD
cor(dat$grain_yield_15,dat$WD_0720,use="complete.obs")
# 0.05169034grain_yield__15 ALL
# 0.05476768 grain_yield__15 EXP_STPAUL_2017_WD
cor(dat$grain_yield_15,dat$WD_0727,use="complete.obs")
# 0.1679906 grain_yield__15 ALL
# 0.04760901 grain_yield__15 EXP_STPAUL_2017_WD

###### Linear Model #######
# WD_0712
mod12=lm(grain_yield_15 ~ WD_0712,dat)
anova(mod12)
#Analysis of Variance Table

#Response: grain_yield__15 ALL
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1   56.1  56.113  1.3126 0.2553
#Residuals 81 3462.8  42.750 

#Analysis of Variance Table

#Response: grain_yield__15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712    1     0.4   0.374  0.0026 0.9596
#Residuals 72 10433.8 144.914   

mod18=lm(grain_yield_15 ~ WD_0718,dat)
anova(mod18)
#Response: grain_yield__15 ALL
#           Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1    1.5  1.5046  0.0498 0.8238
#Residuals 141 4263.5 30.2378

#Response: grain_yield__15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0718     1    14.4  14.387  0.1365 0.7124
#Residuals 126 13278.0 105.381

mod20=lm(grain_yield_15 ~ WD_0720,dat)
anova(mod20)
#Response: grain_yield__15 ALL
#           Df Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1   19.3  19.306   0.584 0.4456
#Residuals 218 7206.1  33.056

#Response: grain_yield__15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0720     1    67.9  67.869  0.5746 0.4494
#Residuals 191 22558.7 118.109 

mod27=lm(grain_yield_15 ~ WD_0727,dat)
anova(mod27)
#Response: grain_yield__15 ALL
#           Df Sum Sq Mean Sq F value Pr(>F)
#WD_0727     1  181.7  181.66  5.5758 0.01921 *
#Residuals 192 6255.3   32.58  

#Response: grain_yield_15 EXP_STPAUL_2017_WD
#          Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0727     1    48.6  48.641  0.3953 0.5304
#Residuals 174 21411.0 123.052  

###### glmnet #######

xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
x=x[xn<3,]
dim(x)

impx=makeX(x,na.impute=TRUE)
y=y[rownames(x),,drop=F]
dat=cbind(y,impx)
y=unlist(y)

fit=glmnet(impx,y,alpha=0)

png(sprintf('QTT/%s_%s_avg_zscore_disregulation_glmnet.png',pheno,env))
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

#lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian')

summary(ridge_reg)


### Optimal lambda value

cv_ridge <- cv.glmnet(x, y_train, alpha = 0)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

#  4.421875 grain_yield__15 ALL
# 746.9854 grain_yield__15 EXP_STPAUL_2017_WD




# Prediction and evaluation on train data
predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = x)
eval_results(y_train, predictions_train, train)

# Prediction and evaluation on test data
predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
eval_results(y_test, predictions_test, test)

# grain_yield_15 ALL
#       RMSE    Rsquare
#1 5.458528 0.06083844
#      RMSE     Rsquare
#1 6.041712 -0.04071393

# grain_yield_15 EXP_STPAUL_2017_WD
#      RMSE      Rsquare
#1 10.63319 0.0005849253
#      RMSE      Rsquare
#1 11.41057 -0.008483968