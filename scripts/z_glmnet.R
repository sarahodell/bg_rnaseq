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

#totalrares=c()
#for(time1 in times){
#	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
#	#allrares=allrares[allrares$max_f!="B73_inra",]
#	totalrares=rbind(totalrares,allrares)
#}
#totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares=fread('eqtl/results/beta_merge_all_v3.txt',data.table=F)
#totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)

# first I need to break up by gene_time_founder
#totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

#totalrares=totalrares[!is.na(totalrares$max_f),]
#totalrares=totalrares[totalrares$max_f!="",]

# switch from summarize to reframe
#totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),gene_time=unique(gene_time),max_f=unique(max_f))

#beta_z=totf%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
#beta_z=as.data.frame(beta_z,stringsAsFactors=F)
#beta_z=beta_z[!is.na(beta_z$beta_z),]

# Add beta_z score to totalrares
#totalrares$beta_z=beta_z[match(totalrares$gene_time_founder,beta_z$gene_time_founder),]$beta_z


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

#dat=data.frame(y=phenotypes,stringsAsFactors=F)
#names(dat)=pheno
for(time1 in times){
	tbetaz=totalrares[totalrares$time==time1,]
	tbetaz$abs_betaz=abs(tbetaz$beta_z)
	genes=unique(tbetaz$Gene_ID)
	# di is average disregulation
	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
	# di is sum of disregulation **2
	inter=intersect(rownames(phenotypes),unique(totalrares$ID))
	dev=data.frame(id=inter,stringsAsFactors=F)
	dev=cbind(dev,as.data.frame(matrix(nrow=nrow(dev),ncol=length(genes))))
	names(dev)=c('id',genes)
	for(i in 1:nrow(dev)){
		subt=tbetaz[tbetaz$ID==dev$id[i],]
		dev[i,genes]=subt[match(names(dev)[-1],subt$Gene_ID),]$abs_betaz
	}
	dev[is.na(dev)]=0
	rownames(dev)=dev$id
	dev=dev[,-1]
	nullinds=colSums(apply(dev,MARGIN=1,function(x) x==0))
	nullinds=names(nullinds[nullinds==ncol(dev)])
	dev=dev[!(rownames(dev) %in% nullinds),]
	inter=intersect(rownames(phenotypes),rownames(dev))
	y=phenotypes[inter,pheno,drop=F]
	n_inds=nrow(dev)
	n_test=20
	n_folds=floor(n_inds/n_test) # at least 20 in test
	print(sprintf("%s %.0f Fold Validation",time1,n_folds))
	breaks=round(seq(1,n_inds,length.out=(n_folds+1)))
	breakn=breaks[2]-breaks[1]# 45 in test, 180 in train
	seq_end=n_inds-breakn
	blup_predictions=data.frame(ID=rownames(dev),stringsAsFactors=F)
	rownames(blup_predictions)=blup_predictions$ID
	blup_predictions=cbind(blup_predictions,as.data.frame(matrix(nrow=nrow(dev),ncol=n_folds)))
	names(blup_predictions)=c('ID',breaks[1:n_folds])
	dat=cbind(y,dev)
	cv_results=c()
	draw = sample(1:n_inds,n_inds)
	for(i in breaks[1:n_folds]){
		r_blup1=NA
		while(is.na(r_blup1)){
			i_end=i+(breakn-1)
			#set.seed(i) 
			# randomly order samples
			index=seq(i,i_end)
			subdraw=draw[index]
			cols_reg = colnames(dat)
			test = dat[subdraw,] # Create the training data 
			train = dat[-subdraw,]
			#y_test = y[subdraw,]
			#test_data=cbind(y_test,x_test)
			#colnames(test_data)=cols_reg
			#x_train = as.matrix(dev[-subdraw,])
			#y_train = y[-subdraw,]
			#train_data=cbind(y_train,x_train)
			#colnames(train_data)=cols_reg
			f=as.formula(paste0(pheno,' ~ .'))
			dummies <- dummyVars(f, data = dat[,cols_reg])
			train_dummies = predict(dummies, newdata = train[,cols_reg])
			test_dummies = predict(dummies, newdata = test[,cols_reg])
			#print(dim(train_dummies)); print(dim(test_dummies))
			x = as.matrix(train_dummies)
			y_train = train[,pheno]
			x_test = as.matrix(test_dummies)
			y_test = test[,pheno]
			ridge_reg = glmnet(x, y_train,alpha = 0,family='gaussian')
			#summary(ridge_reg)
			### Optimal lambda value from validation set
			#png(sprintf('QTT/images/%s_%s_%s_z_score_cv_%.0f.png',time1,pheno,env,i))
			cv_ridge <- cv.glmnet(x, y_train, nfolds=n_folds,alpha = 0,family='gaussian')
			#plot(cv_ridge)
			#dev.off()
			optimal_lambda1 <- cv_ridge$lambda.min
			optimal_lambda1
			optimal_lambda2 <- cv_ridge$lambda.1se
			optimal_lambda2
			# Prediction and evaluation on train data
			predictions_train <- predict(ridge_reg, s = c(optimal_lambda1,optimal_lambda2), newx = x,type = "response",exact=FALSE)
			train_eval1=eval_results(y_train, predictions_train[,1], train)
			train_eval2=eval_results(y_train, predictions_train[,2], train)
			# Prediction and evaluation on test data
			predictions_test <- predict(ridge_reg, s = c(optimal_lambda1,optimal_lambda2), newx = x_test,type = "response",exact=FALSE)
			test_eval1=eval_results(y_test, predictions_test[,1], test)
			test_eval2=eval_results(y_test, predictions_test[,2], test)
			test_blups=y[rownames(predictions_test),]
			r_blup1=cor(predictions_test[,1],test_blups)
			r_blup2=cor(predictions_test[,2],test_blups)
			# Save predicted BLUPs
			test_ids=rownames(x_test)
		}
		blup_predictions[test_ids,c(paste0(i))]=predictions_test[,1]
		line1=data.frame(time=time1,i=i,lambda.min=optimal_lambda1,train_RMSE_1=train_eval1[1],train_R2_1=train_eval1[2],test_RMSE_1=test_eval1[1],test_R2_1=test_eval1[2],test_r_1=r_blup1,lambda.1se=optimal_lambda2,train_RMSE_2=train_eval2[1],train_R2_2=train_eval2[2],test_RMSE_2=test_eval2[1],test_R2_2=test_eval2[2],test_r_2=r_blup2)
		#line=data.frame(time=time1,i=i,lambda=optimal_lambda1,train_RMSE=train_eval[1],train_R2=train_eval[2],test_RMSE=test_eval[1],test_R2=test_eval[2],test_r=r_blup)
		cv_results=rbind(cv_results,line1)	
	}
	cv_results=as.data.frame(cv_results,stringsAsFactors=F)
	names(cv_results)=c('time','i','lambda.min','train_RMSE_1','train_R2_1','test_RMSE_1','test_R2_1','test_cor_1','lambda.1se','train_RMSE_2','train_R2_2','test_RMSE_2','test_R2_2','test_cor_2')
	blup_predictions$mean_pred=apply(blup_predictions[,2:(n_folds+1)],MARGIN=1,function(x) mean(x,na.rm=T))
	predictions[,time1]=blup_predictions[match(rownames(predictions),rownames(blup_predictions)),]$mean_pred
	all_cv=rbind(all_cv,cv_results)
}

all_cv=as.data.frame(all_cv,stringsAsFactors=F)
fwrite(all_cv,sprintf('QTT/%s_%s_z_score_cross_validation_results_all.txt',pheno,env),row.names=F,quote=F,sep='\t')

phenos=c("tkw_15","grain_yield_15")
envs=c("ALL","EXP_STPAUL_2017_WD")

for(env in envs){
	for(pheno in phenos){
		print(pheno)
		print(env)
		all_cv=fread(sprintf('QTT/%s_%s_z_score_cross_validation_results2.txt',pheno,env),data.table=F)
		summary=all_cv %>% group_by(time) %>% reframe(mean_test_R2=mean(test_R2_1,na.rm=T),sd_test_R2=sd(test_R2_1,na.rm=T),mean_train_R2=mean(train_R2_1,na.rm=T),sd_train_R2=sd(train_R2_1,na.rm=T),mean_test_cor=mean(test_cor_1,na.rm=T))
		print(summary)
	}
}

all_cv=fread(sprintf('QTT/%s_%s_z_score_cross_validation_results_all.txt',pheno,env),data.table=F)
summary=all_cv %>% group_by(time) %>% reframe(mean_test_R2=mean(test_R2_1,na.rm=T),sd_test_R2=sd(test_R2_1,na.rm=T),mean_train_R2=mean(train_R2_1,na.rm=T),sd_train_R2=sd(train_R2_1,na.rm=T),mean_test_cor=mean(test_cor_1,na.rm=T))
print(summary)

# tkw_15 ALL
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712     -0.0129          0.332        0.146 
#2 WD_0718     -0.0442          0.106       -0.0331
#3 WD_0720     -0.0467          0.600        0.145 
#4 WD_0727     -0.00920         0.324        0.118
predictions=as.data.frame(predictions,stringsAsFactors=F)
predictions$pheno=phenotypes[match(rownames(predictions),rownames(phenotypes)),pheno]
fwrite(predictions,sprintf('QTT/%s_%s_zscore_fitted_values_all.txt',pheno,env),row.names=F,quote=F,sep='\t')

cor(predictions[,times],predictions$pheno,use="complete.obs")
#WD_0712  0.2229664
#WD_0718 -0.4976741
#WD_0720  0.3374909
#WD_0727  0.3035828

# tkw_15 EXP_STPAUL_2017_WD
## A tibble: 4 × 4
#   time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.0420         0.208        0.0778
#2 WD_0718      -0.0567         0.201        0.0473
#3 WD_0720      -0.0515         0.499        0.166 
#4 WD_0727      -0.0324         0.280        0.101 

# grain_yield_15 ALL
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712    -0.0886           0.260        0.0312
#2 WD_0718    -0.0406           0.261        0.0465
#3 WD_0720    -0.00380          0.340        0.105 
#4 WD_0727     0.000150         0.613        0.203

# grain_yield_15 EXP_STPAUL_2017_WD
# time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.101          0.174       0.0237 
#2 WD_0718      -0.0354         0.139      -0.0497 
#3 WD_0720      -0.0335         0.200      -0.00426
#4 WD_0727      -0.0735         0.452       0.134

########## 
#env="EXP_STPAUL_2017_WD"
##env="ALL"
## grain_yield
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
#	tbetaz=totalrares[totalrares$time==time1,]
#	genes=unique(tbetaz$Gene_ID)
#	# di is average disregulation
#	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
#	# di is sum of disregulation **2
#	inter=intersect(rownames(phenotypes),unique(totalrares$ID))
#	dev=data.frame(id=inter,stringsAsFactors=F)
#	#dev=data.frame(id=unique(totalrares$ID),stringsAsFactors=F)
#	dev=cbind(dev,as.data.frame(matrix(nrow=nrow(dev),ncol=length(genes))))
#	names(dev)=c('id',genes)
#	for(i in 1:nrow(dev)){
#		subt=tbetaz[tbetaz$ID==dev$id[i],]
#		dev[i,genes]=subt[match(names(dev)[-1],subt$Gene_ID),]$beta_z
#	}
#	dev[is.na(dev)]=0
#	
#
#	rownames(dev)=dev$id
#	dev=dev[,-1]
#	
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
#		png(sprintf('QTT/images/%s_%s_%s_z_score_score_cv.png',time1,pheno,env))
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
#fwrite(all_cv,sprintf('QTT/%s_%s_z_score_cross_validation_results.txt',pheno,env),row.names=F,quote=F,sep='\t')
#
#
#summary=all_cv %>% group_by(time) %>% reframe(mean_test_R2=mean(test_R2_1,na.rm=T),mean_train_R2=mean(train_R2_1,na.rm=T),mean_test_cor=mean(test_cor_1,na.rm=T))
#print(summary)

# grain_yield_15 ALL
# A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712      -0.103          0.268       -0.0295
#2 WD_0718      -0.0546         0.239       -0.0560
#3 WD_0720      -0.0196         0.346        0.184 
#4 WD_0727       0.0126         0.623        0.226

# grain_yield_15 EXP_STPAUL_2017_WD
# A tibble: 4 × 4
#  time    mean_test_R2 mean_train_R2 mean_test_cor
#  <chr>          <dbl>         <dbl>         <dbl>
#1 WD_0712     -0.0489         0.109        -0.0519
#2 WD_0718     -0.0509         0.0430       -0.0287
#3 WD_0720     -0.122          0.174        -0.0251
#4 WD_0727     -0.00958        0.460         0.169

#totalrares=c()
#for(time1 in times){
#	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
#	#allrares=allrares[allrares$max_f!="B73_inra",]
#	totalrares=rbind(totalrares,allrares)
#}
#totalrares=as.data.frame(totalrares,stringsAsFactors=F)
#totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
#
## first I need to break up by gene_time_founder
#totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)
#
#totalrares=totalrares[!is.na(totalrares$max_f),]
#totalrares=totalrares[totalrares$max_f!="",]
#
## switch from summarize to reframe
#totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),gene_time=unique(gene_time),max_f=unique(max_f))
#
#beta_z=totf%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
#beta_z=as.data.frame(beta_z,stringsAsFactors=F)
#beta_z=beta_z[!is.na(beta_z$beta_z),]
#
## Add beta_z score to totalrares
#totalrares$beta_z=beta_z[match(totalrares$gene_time_founder,beta_z$gene_time_founder),]$beta_z
## 
#x=data.frame(id=unique(totalrares$ID),stringsAsFactors=F)
#for(time1 in times){
#	tbetaz=totalrares[totalrares$time==time1,]
#	# di is average disregulation
#	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
#	# di is sum of disregulation **2
#	di=tbetaz %>% group_by(ID) %>% reframe(di_score=mean(beta_z^2,na.rm=T))
#
#	di=as.data.frame(di,stringsAsFactors=F)
#	x[,time1]=di[match(x$id,di$ID),]$di_score
#}
#inds=x$id
#phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
#
#phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
#rownames(phenotypes)=phenotypes$Genotype_code
#
#
#y=phenotypes[inds,pheno,drop=F]
#y=y[!is.na(y[,pheno]),,drop=F]
#
#
#newinds=rownames(y)
#
#rownames(x)=x$id
#x=x[,-1]
#x=x[newinds,]
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
#
#y=y[rownames(x),,drop=F]
#dat=cbind(y,x)
##y=unlist(y)
#
#
#
#
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
##[1] 0.004876447  tkw_15 EXP_STPAUL_2017_WD
#
## 0.01617424  tkw_15 ALL
#
#mean(cv_results$test_R2)
##[1] -0.1248758 tkw_15 EXP_STPAUL_2017_WD
##-0.1356803 tkw_15 ALL
#
#mean(cv_results$test_cor,na.rm=T)
##[1] 0.07783118 tkw_15 EXP_STPAUL_2017_WD
## 0.04820048 tkw_15 ALL
#
#fwrite(cv_results,'QTT/%s_%s_cross_validation_results.txt',row.names=F,quote=F,sep='\t')
#
#blup_predictions$mean_pred=apply(blup_predictions[,2:(n_folds+1)],MARGIN=1,function(x) mean(x,na.rm=T))
#
#cor(blup_predictions$mean_pred,y,use="complete.obs")
##         tkw_15 ALL
##[1,]  -0.143375
#
##         tkw_15 EXP_STPAUL_2017_WD
##[1,]-0.231752
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
####### Di from average zscore (from top 5k genes) #######
#totalrares=c()
#for(time1 in times){
#	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
#	#allrares=allrares[allrares$max_f!="B73_inra",]
#	totalrares=rbind(totalrares,allrares)
#}
#totalrares=as.data.frame(totalrares,stringsAsFactors=F)
#totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)
#
## first I need to break up by gene_time_founder
#totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)
#
#totalrares=totalrares[!is.na(totalrares$max_f),]
#totalrares=totalrares[totalrares$max_f!="",]
#
## switch from summarize to reframe
#totf=totalrares %>% group_by(gene_time_founder) %>% reframe(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),gene_time=unique(gene_time),max_f=unique(max_f))
#
#beta_z=totf%>% group_by(Gene_ID,time) %>% mutate(beta_z=(beta-mean(beta,na.rm=T))/sd(beta,na.rm=T))
#beta_z=as.data.frame(beta_z,stringsAsFactors=F)
#beta_z=beta_z[!is.na(beta_z$beta_z),]
#
## Add beta_z score to totalrares
#totalrares$beta_z=beta_z[match(totalrares$gene_time_founder,beta_z$gene_time_founder),]$beta_z
## 
#x=data.frame(id=unique(totalrares$ID),stringsAsFactors=F)
#for(time1 in times){
#	tbetaz=totalrares[totalrares$time==time1,]
#	# di is average disregulation
#	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
#	# di is sum of disregulation **2
#	di=tbetaz %>% group_by(ID) %>% reframe(di_score=mean(beta_z^2,na.rm=T))
#
#	di=as.data.frame(di,stringsAsFactors=F)
#	x[,time1]=di[match(x$id,di$ID),]$di_score
#}
#inds=x$id
#phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
#
#phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
#rownames(phenotypes)=phenotypes$Genotype_code
#
#
#y=phenotypes[inds,pheno,drop=F]
#y=y[!is.na(y[,pheno]),,drop=F]
#
#
#newinds=rownames(y)
#
#rownames(x)=x$id
#x=x[,-1]
#x=x[newinds,]
#
#times2=c("WD_0720","WD_0727")
## filter our individuals with 
#x=x[,times2]
#xn=apply(x,MARGIN=1,function(x) sum(is.na(x)))
#n_na=apply(x,MARGIN=2,function(x) sum(is.na(x)))
#
#x=x[xn<2,]
#dim(x)
#
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
#y$grain_yield_15_scaled=(y$grain_yield_15-mean(y$grain_yield_15,na.rm=T))/sd(y$grain_yield_15,na.rm=T)
#
#draw = sample(1:n_inds, n_inds)
#for(i in breaks[1:(n_folds-1)]){
#	i_end=i+(breakn-1)
#	#set.seed(i) 
#	# randomly order samples
#	index=seq(i,i_end)
#	subdraw=draw[index]
#	x_test = x[subdraw,] # Create the training data 
#	y_test = y[subdraw,'grain_yield_15_scaled']
#	
#	train=dat[-subdraw,]
#	test=dat[subdraw,]
#	x_train = x[-subdraw,] 
#	y_train = y[-subdraw,'grain_yield_15_scaled']
#	
#	# Inner -fold, for this training set, what is the best lambda?
#	#for(j in seq(1,n_train,breakn)){
#	#	j_end=j+(breakn-1)
#	#	v_index=seq(j,j_end)
#		#validation=train[v_index,]
#		#intrain=train[-v_index,]
#		
#		#### Run the model
#	cols_reg = c(times2,'grain_yield_15')
#	dummies <- dummyVars(grain_yield_15 ~ ., data = dat[,cols_reg])
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
##[1]0.0002296072  grain_yield_15 EXP_STPAUL_2017_WD
#
## 0.02262143  grain_yield_15 ALL
#
#mean(cv_results$test_R2)
##[1] -0.1261684 grain_yield_15 EXP_STPAUL_2017_WD
##-0.1517769 grain_yield_15 ALL
#
#mean(cv_results$test_cor,na.rm=T)
##[1] 0.007688066 grain_yield_15 EXP_STPAUL_2017_WD
## 0.1518709 grain_yield_15 ALL
#
#fwrite(cv_results,'QTT/%s_%s_cross_validation_results.txt',row.names=F,quote=F,sep='\t')
#
#blup_predictions$mean_pred=apply(blup_predictions[,2:(n_folds+1)],MARGIN=1,function(x) mean(x,na.rm=T))
#
#cor(blup_predictions$mean_pred,y,use="complete.obs")
#
##[1] -0.2808178 grain_yield_15 EXP_STPAUL_2017_WD
## -0.1140374 grain_yield_15 ALL