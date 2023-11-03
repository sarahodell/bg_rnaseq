#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])

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
#env="EXP_STPAUL_2017_WD"
#env="ALL"
# thousand kernel weight
#pheno="tkw_15"


###### Di from average zscore (from top 5k genes) #######
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code
inds=rownames(phenotypes)
phenotypes=phenotypes[inds,pheno,drop=F]
phenotypes=phenotypes[!is.na(phenotypes[,pheno]),,drop=F]
inds=rownames(phenotypes)

#plot_list=list()
#count=1
all_cv=c()

predictions=data.frame(ID=inds,stringsAsFactors=F)
rownames(predictions)=predictions$ID
predictions=cbind(predictions,as.data.frame(matrix(nrow=nrow(predictions),ncol=4)))
names(predictions)=c('ID',"WD_0712","WD_0718","WD_0720","WD_0727")

for(time1 in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	rownames(exp)=exp$V1
	exp=exp[,-1]
	inter=intersect(inds,rownames(exp))
	exp=exp[inter,]
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	rownames(unlog)=rownames(exp)
	colnames(unlog)=colnames(exp)
	avg_exp = apply(unlog,MARGIN=2,mean)
	names(avg_exp)=names(unlog)
	avg_exp=sort(avg_exp,decreasing=TRUE)
	avg_exp=avg_exp[1:5000] # grap top 5000 highest expressed genes
	#unlog=unlog[,names(avg)]
	unlog=unlog[inter,names(avg_exp)]
	avg_logexp=log2(avg_exp)
	absdev=sapply(seq(1,length(avg_exp)),function(x) log2(abs(unlog[,x]-avg_exp[x])))
	#direction=sapply(seq(1,length(avg_exp)),function(x) ifelse((unlog[,x]-avg_exp[x])>0,1,-1))
	dev=absdev #*direction
	#rownames(unlog)=rownames(exp)
	#avg_exp = apply(unlog,MARGIN=2,mean)
	#avg_exp[avg_exp<1]=0
	#dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x]))
	rownames(dev)=rownames(unlog)
	colnames(dev)=colnames(unlog)
	y=phenotypes[rownames(dev),pheno,drop=F]
	n_inds=nrow(dev)
	n_test=20
	n_folds=floor(n_inds/n_test) # at least 20 in test
	print(sprintf("%s %.0f Fold Validation",time1,n_folds))
	breaks=round(seq(1,n_inds,length.out=(n_folds+1)))
	breakn=breaks[2]-breaks[1]# 45 in test, 180 in train
	seq_end=n_inds-breakn
	draw = sample(1:n_inds,n_inds)
	blup_predictions=data.frame(ID=rownames(dev),stringsAsFactors=F)
	rownames(blup_predictions)=blup_predictions$ID
	blup_predictions=cbind(blup_predictions,as.data.frame(matrix(nrow=nrow(dev),ncol=n_folds)))
	names(blup_predictions)=c('ID',breaks[1:n_folds])
	cv_results=c()
	for(i in breaks[1:n_folds]){
		r_blup1=NA
		while(is.na(r_blup1)){
			i_end=i+(breakn-1)
			#set.seed(i) 
			# randomly order samples
			index=seq(i,i_end)
			subdraw=draw[index]
			x_test = as.matrix(dev[subdraw,]) # Create the training data 
			y_test = y[subdraw,]
			x_train = as.matrix(dev[-subdraw,])
			y_train = y[-subdraw,]
			ridge_reg = glmnet(x_train, y_train,alpha = 0,family='gaussian')
			#summary(ridge_reg)
			### Optimal lambda value from validation set
			png(sprintf('QTT/images/%s_%s_%s_totexp_score_cv_%.0f.png',time1,pheno,env,i))
			cv_ridge <- cv.glmnet(x_train, y_train, nfolds=n_folds,alpha = 0,family='gaussian')
			plot(cv_ridge)
			dev.off()
			optimal_lambda1 <- cv_ridge$lambda.min
			optimal_lambda1
			optimal_lambda2 <- cv_ridge$lambda.1se
			optimal_lambda2
			# Prediction and evaluation on train data
			predictions_train <- predict(ridge_reg, s = c(optimal_lambda1,optimal_lambda2), newx = x_train,type = "response",exact=FALSE)
			train_eval1=eval_results(y_train, predictions_train[,1], x_train)
			train_eval2=eval_results(y_train, predictions_train[,2], x_train)
			# Prediction and evaluation on test data
			predictions_test <- predict(ridge_reg, s = c(optimal_lambda1,optimal_lambda2), newx = x_test,type = "response",exact=FALSE)
			test_eval1=eval_results(y_test, predictions_test[,1], x_test)
			test_eval2=eval_results(y_test, predictions_test[,2], x_test)
			test_blups=y[rownames(predictions_test),]
			r_blup1=cor(predictions_test[,1],test_blups,use="complete.obs")
			r_blup2=cor(predictions_test[,2],test_blups,,use="complete.obs")
			# Save predicted BLUPs
		}
		test_ids=rownames(x_test)
		blup_predictions[test_ids,c(paste0(i))]=predictions_test
		line1=data.frame(time=time1,i=i,lambda.min=optimal_lambda1,train_RMSE_1=train_eval1[1],train_R2_1=train_eval1[2],test_RMSE_1=test_eval1[1],test_R2_1=test_eval1[2],test_r_1=r_blup1,lambda.1se=optimal_lambda2,train_RMSE_2=train_eval2[1],train_R2_2=train_eval2[2],test_RMSE_2=test_eval2[1],test_R2_2=test_eval2[2],test_r_2=r_blup2)
		#line=data.frame(time=time1,i=i,lambda=optimal_lambda1,train_RMSE=train_eval[1],train_R2=train_eval[2],test_RMSE=test_eval[1],test_R2=test_eval[2],test_r=r_blup)
		cv_results=rbind(cv_results,line1)
	}
	blup_predictions$mean_pred=apply(blup_predictions[,2:(n_folds+1)],MARGIN=1,function(x) mean(x,na.rm=T))
	predictions[,time1]=blup_predictions[match(rownames(predictions),rownames(blup_predictions)),]$mean_pred
	cv_results=as.data.frame(cv_results,stringsAsFactors=F)
	names(cv_results)=c('time','i','lambda.min','train_RMSE_1','train_R2_1','test_RMSE_1','test_R2_1','test_cor_1','lambda.1se','train_RMSE_2','train_R2_2','test_RMSE_2','test_R2_2','test_cor_2')
	all_cv=rbind(all_cv,cv_results)
}

all_cv=as.data.frame(all_cv,stringsAsFactors=F)
fwrite(all_cv,sprintf('QTT/%s_%s_totexp_cross_validation_results2.txt',pheno,env),row.names=F,quote=F,sep='\t')

predictions=as.data.frame(predictions,stringsAsFactors=F)
predictions$pheno=phenotypes[match(rownames(predictions),rownames(phenotypes)),pheno]
fwrite(predictions,sprintf('QTT/%s_%s_totexp_fitted_values2.txt',pheno,env),row.names=F,quote=F,sep='\t')

print(cor(predictions[,times],predictions$pheno,use="complete.obs"))


summary=all_cv %>% group_by(time) %>% reframe(mean_test_R2=mean(test_R2_1,na.rm=T),mean_train_R2=mean(train_R2_1,na.rm=T),mean_test_cor=mean(test_cor_1,na.rm=T))
print(summary)

#tkw_15 ALL
# A tibble: 4 × 4
#   time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.120         0.0488       -0.217 
#2 WD_0718      -0.211         0.214         0.0313
#3 WD_0720      -0.0718        0.0755       -0.143 
#4 WD_0727      -0.0911        0.203        -0.0300


#tkw_15 EXP_STPAUL_2017_WD
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.316          0           NaN    
#2 WD_0718      -0.0631         0.687         0.228
#3 WD_0720      -0.0489         0.172         0.109
#4 WD_0727      -0.132          0.147        -0.263

#grain_yield_15 ALL
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.299         0.0252        -0.487
#2 WD_0718      -0.156         0.0982        -0.508
#3 WD_0720      -0.0235        0.542          0.190
#4 WD_0727       0.0662        0.574          0.357

#grain_yield_15 EXP_STPAUL_2017_WD
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.230          0          NaN     
#2 WD_0718      -0.109          0.256        0.0986
#3 WD_0720       0.0891         0.745        0.399 
#4 WD_0727       0.0943         0.831        0.440 


########### 
#env="EXP_STPAUL_2017_WD"
##env="ALL"
## thousand kernel weight
#pheno="grain_yield_15"
#
#
####### Di from average zscore (from top 5k genes) #######
#phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
#phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
#rownames(phenotypes)=phenotypes$Genotype_code
#inds=rownames(phenotypes)
#phenotypes=phenotypes[inds,pheno,drop=F]
#phenotypes=phenotypes[!is.na(phenotypes[,pheno]),,drop=F]
#inds=rownames(phenotypes)
#
##plot_list=list()
##count=1
#all_cv=c()
#
##dat=data.frame(y=phenotypes,stringsAsFactors=F)
##names(dat)=pheno
#for(time1 in times){
#
#	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
#	rownames(exp)=exp$V1
#	exp=exp[,-1]
#	inter=intersect(inds,rownames(exp))
#	avg_exp = apply(exp,MARGIN=2,mean)
#	names(avg_exp)=names(exp)
#	avg_exp=sort(avg_exp,decreasing=TRUE)
#	avg_exp=avg_exp[1:5000] # grap top 5000 highest expressed genes
#	#unlog=unlog[,names(avg
#	exp=exp[inter,avg_exp]
#	dev=sapply(seq(1,length(avg_exp)),function(x) (exp[,x]-avg_exp[x]))
#
#	#unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
#	#rownames(unlog)=rownames(exp)
#	#avg_exp = apply(unlog,MARGIN=2,mean)
#	#avg_exp[avg_exp<1]=0
#	#dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x]))
#	rownames(dev)=rownames(exp)
#	colnames(dev)=colnames(exp)
#	
#	cv_results=c()
#	#set.seed(100)
#	y=phenotypes[rownames(dev),pheno,drop=F]
#	#y$tkw_15_scaled=(y$tkw_15-mean(y$tkw_15,na.rm=T))/sd(y$tkw_15,na.rm=T)
#	n_inds=nrow(dev)
#	n_folds=10
#	breaks=round(seq(1,n_inds,length.out=n_folds))
#	breakn=breaks[2]-breaks[1]# 45 in test, 180 in train
#	seq_end=n_inds-breakn
#	draw = sample(1:n_inds,n_inds)
#	for(i in breaks[1:(n_folds-1)]){
#		i_end=i+(breakn-1)
#		#set.seed(i) 
#		# randomly order samples
#		index=seq(i,i_end)
#		subdraw=draw[index]
#		x_test = as.matrix(dev[subdraw,]) # Create the training data 
#		y_test = y[subdraw,]
#	
#		x_train = as.matrix(dev[-subdraw,])
#		y_train = y[-subdraw,]
#		ridge_reg = glmnet(x_train, y_train,alpha = 0,family='gaussian')
#		#summary(ridge_reg)
#		### Optimal lambda value from validation set
#		png(sprintf('QTT/images/%s_%s_%s_totexp_score_cv.png',time1,pheno,env))
#		cv_ridge <- cv.glmnet(x_train, y_train, nfolds=10,alpha = 0,family='gaussian')
#		plot(cv_ridge)
#		dev.off()
#		optimal_lambda1 <- cv_ridge$lambda.min
#		optimal_lambda1
#		optimal_lambda2 <- cv_ridge$lambda.1se
#		optimal_lambda2
#		# Prediction and evaluation on train data
#		predictions_train <- predict(ridge_reg, s = c(optimal_lambda1,optimal_lambda2), newx = x_train,type = "response",exact=FALSE)
#		train_eval1=eval_results(y_train, predictions_train[,1], x_train)
#		train_eval2=eval_results(y_train, predictions_train[,2], x_train)
#		# Prediction and evaluation on test data
#		predictions_test <- predict(ridge_reg, s = c(optimal_lambda1,optimal_lambda2), newx = x_test,type = "response")
#		test_eval1=eval_results(y_test, predictions_test[,1], x_test)
#		test_eval2=eval_results(y_test, predictions_test[,2], x_test)
#		test_blups=y[rownames(predictions_test),]
#		r_blup1=cor(predictions_test[,1],test_blups)
#		r_blup2=cor(predictions_test[,2],test_blups)
#		# Save predicted BLUPs
#		test_ids=rownames(x_test)
#		#blup_predictions[test_ids,c(paste0(i))]=predictions_test
#		line1=data.frame(time=time1,i=i,lambda.min=optimal_lambda1,train_RMSE_1=train_eval1[1],train_R2_1=train_eval1[2],test_RMSE_1=test_eval1[1],test_R2_1=test_eval1[2],test_r_1=r_blup1,lambda.1se=optimal_lambda2,train_RMSE_2=train_eval2[1],train_R2_2=train_eval2[2],test_RMSE_2=test_eval1[1],test_R2_2=test_eval2[2],test_r_2=r_blup2)
#		#line=data.frame(time=time1,i=i,lambda=optimal_lambda1,train_RMSE=train_eval[1],train_R2=train_eval[2],test_RMSE=test_eval[1],test_R2=test_eval[2],test_r=r_blup)
#		cv_results=rbind(cv_results,line1)
#	}
#	cv_results=as.data.frame(cv_results,stringsAsFactors=F)
#	names(cv_results)=c('time','i','lambda.min','train_RMSE_1','train_R2_1','test_RMSE_1','test_R2_1','test_cor_1','lambda.1se','train_RMSE_2','train_R2_2','test_RMSE_2','test_R2_2','test_cor_2')
#	all_cv=rbind(all_cv,cv_results)
#}
#
#all_cv=as.data.frame(all_cv,stringsAsFactors=F)
#fwrite(all_cv,sprintf('QTT/%s_%s_totexp_cross_validation_results.txt',pheno,env),row.names=F,quote=F,sep='\t')
#
#
#summary=all_cv %>% group_by(time) %>% reframe(mean_test_R2=mean(test_R2_1,na.rm=T),mean_train_R2=mean(train_R2_1,na.rm=T),mean_test_cor=mean(test_cor_1,na.rm=T))
#print(summary)

# Grain_yield_15 ALL
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.328         0.0226       -0.314 
#2 WD_0718      -0.0688        0.0119        0.0301
#3 WD_0720      -0.0839        0.0321        0.0238
#4 WD_0727      -0.0944        0.0368        0.0428


# grain_yield_!5 EXP_STPAUL_2017_EXP
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.262         0.0736       -0.138 
#2 WD_0718      -0.199         0.0243       -0.128 
#3 WD_0720      -0.0423        0.0457        0.0867
#4 WD_0727      -0.0523        0.0747        0.189

#x=dat[,times]
## drop individuals with fewer than 2 timepoints
#xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
#x=x[xn<3,]
#dim(x)
#
#n_na=apply(x,MARGIN=2,function(x) sum(is.na(x)))
#xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
#
####### glmnet #######
#times2=c("WD_0720","WD_0727")
## filter our individuals with 
#x=x[,times2]
#xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
#n_na=apply(x,MARGIN=2,function(x) sum(is.na(x)))
#
#x=x[xn<2,]
#dim(x)
#
## 225 samples - 15-fold cross-validation
#### cross-validation
#x=x[complete.cases(x),]
#y=y[rownames(x),,drop=F]
#dat=cbind(y,x)
#
#
#n_inds=nrow(x)
#n_folds=10
#
#
#breaks=round(seq(1,n_inds,length.out=n_folds))
#breakn=breaks[2]-breaks[1]# 45 in test, 180 in train
#seq_end=n_inds-breakn
##n_train=n_inds-breakn
#
#
#blup_predictions=data.frame(ID=rownames(dat),stringsAsFactors=F)
#rownames(blup_predictions)=blup_predictions$ID
#
#
#
#cv_results=c()
##set.seed(100)
#blup_predictions=cbind(blup_predictions,as.data.frame(matrix(nrow=nrow(dat),ncol=n_folds)))
#names(blup_predictions)=c('ID',breaks)
#
#y$tkw_15_scaled=(y$tkw_15-mean(y$tkw_15,na.rm=T))/sd(y$tkw_15,na.rm=T)
#
#draw = sample(1:n_inds, n_inds)
#for(i in breaks[1:(n_folds-1)]){
#	i_end=i+(breakn-1)
#	#set.seed(i) 
#	# randomly order samples
#	index=seq(i,i_end)
#	subdraw=draw[index]
#	x_test = x[subdraw,] # Create the training data 
#	y_test = y[subdraw,'tkw_15_scaled']
#	
#	train=dat[-subdraw,]
#	test=dat[subdraw,]
#	x_train = x[-subdraw,] 
#	y_train = y[-subdraw,'tkw_15_scaled']
#	
#	# Inner -fold, for this training set, what is the best lambda?
#	#for(j in seq(1,n_train,breakn)){
#	#	j_end=j+(breakn-1)
#	#	v_index=seq(j,j_end)
#		#validation=train[v_index,]
#		#intrain=train[-v_index,]
#		
#		#### Run the model
#	cols_reg = c(times2,'tkw_15')
#	dummies <- dummyVars(tkw_15 ~ ., data = dat[,cols_reg])
#	train_dummies = predict(dummies, newdata = train[,cols_reg])
#	test_dummies = predict(dummies, newdata = test[,cols_reg])
#	### Ridge Regression
#	impx=makeX(x_train,x_test,na.impute=TRUE)
#	#x_train = as.matrix(train_dummies)
#	#y_train = train$tkw_15
#
#	#x_test = as.matrix(test_dummies)
#	#y_test = test$tkw_15
#
#	#lambdas <- 10^seq(2, -3, by = -.1)
#	ridge_reg = glmnet(impx$x, y_train,alpha = 0,family='gaussian')
#	#summary(ridge_reg)
#
#	### Optimal lambda value from validation set
#	#png(sprintf('QTT/images/%s_%s_z_score_cv.png',pheno,env))
#	cv_ridge <- cv.glmnet(impx$x, y_train, nfolds=10,alpha = 0,family='gaussian')
#	#print(cv_ridge)
#	#dev.off()
#	optimal_lambda <- cv_ridge$lambda.min
#	optimal_lambda
#	# Prediction and evaluation on train data
#	predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = as.matrix(impx$x),type = "response")
#	train_eval=eval_results(y_train, predictions_train, impx$x)
#
#	# Prediction and evaluation on test data
#	predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = impx$xtest,type = "response")
#	test_eval=eval_results(y_test, predictions_test, impx$xtest)
#	
#	test_blups=y[rownames(predictions_test),]$tkw_15_scaled
#	r_blup=cor(predictions_test,test_blups)
#	
#	# Save predicted BLUPs
#	test_ids=rownames(test)
#	blup_predictions[test_ids,c(paste0(i))]=predictions_test
#	
#	line=data.frame(i=i,lambda=optimal_lambda,train_RMSE=train_eval[1],train_R2=train_eval[2],test_RMSE=test_eval[1],test_R2=test_eval[2],test_r=r_blup)
#	cv_results=rbind(cv_results,line)
#	#}
#
#}
#cv_results=as.data.frame(cv_results,stringsAsFactors=F)
#names(cv_results)=c('i','lambda','train_RMSE','train_R2','test_RMSE','test_R2','test_cor')
#
#mean(cv_results$train_R2)
##[1] 0.01289348  tkw_15 EXP_STPAUL_2017_WD
#
##0.02224652  tkw_15 ALL
#
#mean(cv_results$test_R2)
##[1] -0.08159384tkw_15 EXP_STPAUL_2017_WD
##-0.08247933 tkw_15 ALL
#
#mean(cv_results$test_cor,na.rm=T)
##[1] 0.1049881 tkw_15 EXP_STPAUL_2017_WD
## 0.1071496 tkw_15 ALL
#
#fwrite(cv_results,'QTT/%s_%s_totexp_cross_validation_results.txt',row.names=F,quote=F,sep='\t')
#
#blup_predictions$mean_pred=apply(blup_predictions[,2:(n_folds+1)],MARGIN=1,function(x) mean(x,na.rm=T))
#
#cor(blup_predictions$mean_pred,y,use="complete.obs")
##         tkw_15 ALL
##[1,]  -0.02944679
#
##         tkw_15 EXP_STPAUL_2017_WD
##[1,]-0.06290345
#
##fit=glmnet(impx,y,alpha=0)
#
##png(sprintf('QTT/%s_%s_avg_zscore_disregulation_glmnet_2.png',pheno,env))
##print(plot(fit,label=TRUE))
##dev.off()
#
##predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
#
#
##################### #################### 
##################### #################### 
##################### #################### 
##################### #################### 
##################### #################### 
#env="EXP_STPAUL_2017_WD"
#env="ALL"
## grain yield
#pheno="grain_yield_15"
#
#
#log_inverse=function(x){
#  	return(2^x)
#}
#
#eval_results <- function(true, predicted, df) {
#  SSE <- sum((predicted - true)^2)
#  SST <- sum((true - mean(true))^2)
#  R_square <- 1 - (SSE / SST)
#  RMSE = sqrt(SSE/nrow(df))
#
#  
#  # Model performance metrics
#data.frame(
#  RMSE = RMSE,
#  Rsquare = R_square
#)
#  
#}
#
#phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
#phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
#rownames(phenotypes)=phenotypes$Genotype_code
#inds=rownames(phenotypes)
#y=phenotypes[inds,pheno,drop=F]
#y=y[!is.na(y[,pheno]),,drop=F]
#inds=rownames(y)
#
#
#dat=data.frame(y=y,stringsAsFactors=F)
#names(dat)=pheno
#for(time1 in times){
#	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
#	rownames(exp)=exp$V1
#	exp=exp[,-1]
#	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
#	rownames(unlog)=rownames(exp)
#	avg_exp = apply(unlog,MARGIN=2,mean)
#	avg_exp[avg_exp<1]=0
#	names(avg_exp)=names(exp)
#	avg_exp=sort(avg_exp,decreasing=TRUE)
#	avg_exp=avg_exp[1:5000] # grap top 5000 highest expressed genes
#	
#	unlog=unlog[,names(avg_exp)]
#	# re-log the input data
#	#avg_logexp=log2(avg_exp)
#	#avg_logexp[is.infinite(avg_logexp)]=0
#	# For each individual, get the average squared diff from the avg expression
#	# Unlog expression for both, get the average difference, and then relog
#	dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
#	rownames(dev)=rownames(unlog)
#	colnames(dev)=colnames(unlog)
#	# relog log2cpm
#	di=apply(dev,MARGIN=1,function(x) log10(mean(x,na.rm=T)))
#	dat[,time1]=di[match(rownames(dat),names(di))]
#}
#
#x=dat[,times]
## drop individuals with fewer than 2 timepoints
#xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
#x=x[xn<3,]
#dim(x)
#
#n_na=apply(x,MARGIN=2,function(x) sum(is.na(x)))
#xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
#
####### glmnet #######
#times2=c("WD_0720","WD_0727")
## filter our individuals with 
#x=x[,times2]
#xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
#n_na=apply(x,MARGIN=2,function(x) sum(is.na(x)))
#
#x=x[xn<2,]
#dim(x)
#
## 225 samples - 15-fold cross-validation
#### cross-validation
#x=x[complete.cases(x),]
#y=y[rownames(x),,drop=F]
#dat=cbind(y,x)
#
#n_inds=nrow(x)
#n_folds=10
#
#
#breaks=round(seq(1,n_inds,length.out=n_folds))
#breakn=breaks[2]-breaks[1]# 45 in test, 180 in train
#seq_end=n_inds-breakn
##n_train=n_inds-breakn
#
#
#blup_predictions=data.frame(ID=rownames(dat),stringsAsFactors=F)
#rownames(blup_predictions)=blup_predictions$ID
#
#
#
#cv_results=c()
##set.seed(100)
#blup_predictions=cbind(blup_predictions,as.data.frame(matrix(nrow=nrow(dat),ncol=n_folds)))
#names(blup_predictions)=c('ID',breaks)
#
#y$grain_yield_15_scaled=(y$grain_yield_15-mean(y$grain_yield_15,na.rm=T))/sd(y$grain_yield_15,na.rm=T)
#
#draw = sample(1:n_inds, n_inds)
#for(i in breaks[1:(n_folds-1)]){
#	i_end=i+(breakn-1)
#	#set.seed(i) 
#	# randomly order samples
#	index=seq(i,i_end)
#	subdraw=draw[index]
#	x_test = as.matrix(x[subdraw,]) # Create the training data 
#	y_test = y[subdraw,'grain_yield_15']
#	
#	#train=dat[-subdraw,]
#	#test=dat[subdraw,]
#	x_train = as.matrix(x[-subdraw,])
#	y_train = y[-subdraw,'grain_yield_15']
#	
#	# Inner -fold, for this training set, what is the best lambda?
#	#for(j in seq(1,n_train,breakn)){
#	#	j_end=j+(breakn-1)
#	#	v_index=seq(j,j_end)
#		#validation=train[v_index,]
#		#intrain=train[-v_index,]
#		
#		#### Run the model
#	#cols_reg = c(times2,'grain_yield_15')
#	#dummies <- dummyVars(grain_yield_15 ~ ., data = dat[,cols_reg])
#	#train_dummies = predict(dummies, newdata = train[,cols_reg])
#	#test_dummies = predict(dummies, newdata = test[,cols_reg])
#	### Ridge Regression
#	#impx=makeX(x_train,x_test,na.impute=TRUE)
#	#x_train = as.matrix(train_dummies)
#	#y_train = train$tkw_15
#
#	#x_test = as.matrix(test_dummies)
#	#y_test = test$tkw_15
#
#	#lambdas <- 10^seq(2, -3, by = -.1)
#	ridge_reg = glmnet(x_train, y_train,alpha = 0,family='gaussian')
#	#summary(ridge_reg)
#
#	### Optimal lambda value from validation set
#	#png(sprintf('QTT/images/%s_%s_z_score_cv.png',pheno,env))
#	cv_ridge <- cv.glmnet(x_train, y_train, nfolds=10,alpha = 0,family='gaussian')
#	#print(cv_ridge)
#	#dev.off()
#	optimal_lambda <- cv_ridge$lambda.min
#	optimal_lambda
#	# Prediction and evaluation on train data
#	predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = x_train,type = "response")
#	train_eval=eval_results(y_train, predictions_train, impx$x)
#
#
#	# Prediction and evaluation on test data
#	predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = ix_test,type = "response")
#	test_eval=eval_results(y_test, predictions_test, impx$xtest)
#	
#	test_blups=y[rownames(predictions_test),]$grain_yield_15_scaled
#	r_blup=cor(predictions_test,test_blups)
#	
#	# Save predicted BLUPs
#	test_ids=rownames(test)
#	blup_predictions[test_ids,c(paste0(i))]=predictions_test
#	
#	line=data.frame(i=i,lambda=optimal_lambda,train_RMSE=train_eval[1],train_R2=train_eval[2],test_RMSE=test_eval[1],test_R2=test_eval[2],test_r=r_blup)
#	cv_results=rbind(cv_results,line)
#	#}
#
#}
#cv_results=as.data.frame(cv_results,stringsAsFactors=F)
#names(cv_results)=c('i','lambda','train_RMSE','train_R2','test_RMSE','test_R2','test_cor')
#
#mean(cv_results$train_R2)
##[1]0.0001641364  grain_yield_15 EXP_STPAUL_2017_WD
#
## 3.903799e-05  grain_yield_15 ALL
#
#mean(cv_results$test_R2)
##[1] -0.09957765 grain_yield_15 EXP_STPAUL_2017_WD
##-0.01858153 grain_yield_15 ALL
#
#mean(cv_results$test_cor,na.rm=T)
##[1] -0.4699129 grain_yield_15 EXP_STPAUL_2017_WD
## -0.3403155 grain_yield_15 ALL
#
#fwrite(cv_results,'QTT/%s_%s_cross_validation_results.txt',row.names=F,quote=F,sep='\t')
#
#blup_predictions$mean_pred=apply(blup_predictions[,2:(n_folds+1)],MARGIN=1,function(x) mean(x,na.rm=T))
#
#cor(blup_predictions$mean_pred,y,use="complete.obs")
#
##[1] -0.2933744 grain_yield_15 EXP_STPAUL_2017_WD
## -0.1385591 grain_yield_15 ALL