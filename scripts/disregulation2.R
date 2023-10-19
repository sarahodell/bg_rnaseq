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
pheno="tkw_15"

ft_df=fread('eqtl/data/FT_genelist.txt',data.table=F)
##### High GERP
totalrares=fread('eqtl/results/all_5kb_rare_counts_high_GERP.txt',data.table=F)

genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

log_inverse=function(x){
  	return(2^x)
}

phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code


#y=phenotypes[inds,pheno,drop=F]
y=phenotypes[,pheno,drop=F]
y=y[!is.na(y[,pheno]),,drop=F]


# Add beta_z score to totalrares'
totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
#totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

totf=totalrares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(total),gene_time=unique(gene_time),max_f=unique(max_f))

#totf=totalrares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID.x),time=unique(time.x),chr=unique(chr.x),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(total),gene_time=unique(gene_time.x),max_f=unique(max_f))

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
	di=tbetaz %>% group_by(ID) %>% reframe(di_score=sum(beta_z^2,na.rm=T))
	
	di=as.data.frame(di,stringsAsFactors=F)
	di$di_score=log10(di$di_score)
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
#dat=cbind(y,impx)
dat=cbind(y,x)
# WD_0712
mod12=lm(tkw_15 ~ WD_0712,dat)
anova(mod12)

mod18=lm(tkw_15 ~ WD_0718,dat)
anova(mod18)

mod20=lm(tkw_15 ~ WD_0720,dat)
anova(mod20)

datmelt=melt(dat[,c('tkw_15',times)],'tkw_15')

p1=ggplot(aes(x=value,y=tkw_15),data=datmelt) + geom_point(aes(color=variable)) +
xlab("Dysregulation Score (log10(Di))") + ylab("BLUP Scores Thousand-Kernel Weight (kgs)") +
theme_classic()

png(sprintf('QTT/%s_%s_high_gerp_disregulation_score_Zsum_plot.png',pheno,env))
print(p1)
dev.off()



dat=data.frame(grain_yield_15=y,stringsAsFactors=F)
for(time1 in times){
	gerpgenes=unique(totalrares[totalrares$time==time1,]$Gene_ID)
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	#exp=exp[,c('V1',top5k)]
	rownames(exp)=exp$V1
	exp=exp[,-1]
	exp=exp[,gerpgenes]
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	avg_exp = apply(unlog,2,mean)
	names(avg_exp)=names(exp)
	avg_exp[avg_exp<1]=0
	# re-log the input data
	avg_logexp=log2(avg_exp)
	avg_logexp[is.infinite(avg_logexp)]=0
	# For each individual, get the average squared diff from the avg expression
	# Unlog expression for both, get the average difference, and then relog
	dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
	rownames(dev)=rownames(exp)
	colnames(dev)=colnames(exp)
	# relog log2cpm
	di=apply(dev,MARGIN=1,function(x) log10(mean(x,na.rm=T)))
	
	dat[,time1]=di[match(rownames(dat),names(di))]
}



cor(dat$grain_yield_15,dat$WD_0712,use="complete.obs")

#0.1074302

cor(dat$grain_yield_15,dat$WD_0712,use="complete.obs")
#0.06185197

cor(dat$grain_yield_15,dat$WD_0720,use="complete.obs")
#-0.0432486

cor(dat$grain_yield_15,dat$WD_0727,use="complete.obs")
#0.01262483


#####
dat$tkw_15 = phenotypes[match(rownames(dat),phenotypes$Genotype_code),]$tkw_15

cor(dat$tkw_15,dat$WD_0712,use="complete.obs")
# 0.07398719
cor(dat$tkw_15,dat$WD_0718,use="complete.obs")
# 0.05249369
cor(dat$tkw_15,dat$WD_0720,use="complete.obs")
# 0.07572549
cor(dat$tkw_15,dat$WD_0727,use="complete.obs")
# 0.1120877

datmelt=melt(dat[,c('tkw_15',times)],'tkw_15')

p1=ggplot(aes(x=value,y=tkw_15),data=datmelt) + geom_point(aes(color=variable)) 

png(sprintf('QTT/%s_%s_high_gerp_disregulation_score_totexp_plot.png',pheno,env))
print(p1)
dev.off()

#########

totalrares=c()
for(time1 in times){
	allrares=fread(sprintf('eqtl/results/rare_counts_%s_max_f.txt',time1),data.table=F)
	allrares=allrares[allrares$max_f!="B73_inra",]
	totalrares=rbind(totalrares,allrares)
}
totalrares=as.data.frame(totalrares,stringsAsFactors=F)
totalrares$gene_time=paste0(totalrares$Gene_ID,'-',totalrares$time)

# first I need to break up by gene_time_founder
totalrares$gene_time_founder=paste0(totalrares$Gene_ID,'-',totalrares$time,'-',totalrares$max_f)

totalrares=totalrares[!is.na(totalrares$max_f),]
totalrares=totalrares[totalrares$max_f!="",]
totf=totalrares %>% group_by(gene_time_founder) %>% summarize(Gene_ID=unique(Gene_ID),time=unique(time),chr=unique(chr),beta=unique(beta),beta_rank=unique(beta_rank),rare_count=unique(rare_count),gene_time=unique(gene_time),max_f=unique(max_f))

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
	di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2,na.rm=T))

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
#dat=cbind(y,impx)
dat=cbind(y,x)
# WD_0712
mod12=lm(grain_yield_15 ~ WD_0712,dat)
anova(mod12)

#Analysis of Variance Table

#Response: grain_yield_15
#           Df  Sum Sq Mean Sq F value Pr(>F)
#WD_0712     1     2.8   2.764  0.0234 0.8785
#Residuals 248 29295.5 118.127 

mod18=lm(grain_yield_15 ~ WD_0718,dat)
anova(mod18)

mod20=lm(grain_yield_15 ~ WD_0720,dat)
anova(mod20)

mod27=lm(grain_yield_15 ~ WD_0727,dat)
anova(mod27)

# Make plot
datmelt=melt(dat,'grain_yield_15')

p1=ggplot(aes(x=value,y=grain_yield_15),data=datmelt) + geom_point(aes(color=variable)) 

png(sprintf('QTT/%s_%s_disregulation_score_plot.png',pheno,env))
print(p1)
dev.off()

##### Using tkw_15 #####
pheno="tkw_15"
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code


y=phenotypes[inds,pheno,drop=F]
y=y[!is.na(y[,pheno]),,drop=F]

x=data.frame(id=unique(totalrares$ID),stringsAsFactors=F)
for(time1 in times){
	tbetaz=totalrares[totalrares$time==time1,]
	# di is average disregulation
	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2)/length(unique(Gene_ID)))
	# di is sum of disregulation **2
	di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2,na.rm=T))

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

dat=cbind(y,x)
# WD_0712
mod12=lm(tkw_15 ~ WD_0712,dat)
anova(mod12)

mod18=lm(tkw_15 ~ WD_0718,dat)
anova(mod18)

mod20=lm(tkw_15 ~ WD_0720,dat)
anova(mod20)

mod27=lm(tkw_15 ~ WD_0727,dat)
anova(mod27)

#none of them are significant

##### Using average disregulation instead of sum
x=data.frame(id=unique(totalrares$ID),stringsAsFactors=F)
for(time1 in times){
	tbetaz=totalrares[totalrares$time==time1,]
	# di is average disregulation
	di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2,na.rm=T)/length(unique(Gene_ID)))
	# di is sum of disregulation **2
	#di=tbetaz %>% group_by(ID) %>% summarize(di_score=sum(beta_z^2,na.rm=T))

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
dat=cbind(y,x)
# WD_0712
mod12=lm(tkw_15 ~ WD_0712,dat)
anova(mod12)


# none of them are significant


##### Using the measure of disregulation that Kremling used (total expression, not local additive effect)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

log_inverse=function(x){
  	return(2^x)
}

phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

phenotypes=phenotypes[phenotypes$Loc.Year.Treat==env,]
rownames(phenotypes)=phenotypes$Genotype_code


#y=phenotypes[inds,pheno,drop=F]
y=phenotypes[,pheno,drop=F]
y=y[!is.na(y[,pheno]),,drop=F]


dat=data.frame(tkw_15=y,stringsAsFactors=F)
for(time1 in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	#exp=exp[,c('V1',top5k)]
	rownames(exp)=exp$V1
	exp=exp[,-1]
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	avg_exp = apply(unlog,2,mean)
	names(avg_exp)=names(exp)
	avg_exp[avg_exp<1]=0
	# re-log the input data
	avg_logexp=log2(avg_exp)
	avg_logexp[is.infinite(avg_logexp)]=0
	# For each individual, get the average squared diff from the avg expression
	# Unlog expression for both, get the average difference, and then relog
	dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
	rownames(dev)=rownames(exp)
	colnames(dev)=colnames(exp)
	# relog log2cpm
	di=apply(dev,MARGIN=1,function(x) log10(mean(x,na.rm=T)))
	
	dat[,time1]=di[match(rownames(dat),names(di))]
}

cor(dat$tkw_15,dat$WD_0712,use="complete.obs")
# 0.05176112
cor(dat$tkw_15,dat$WD_0718,use="complete.obs")
# 0.04683624
cor(dat$tkw_15,dat$WD_0720,use="complete.obs")
# 0.09817568
cor(dat$tkw_15,dat$WD_0727,use="complete.obs")
# 0.1199425

dat$grain_yield_15 = phenotypes[match(rownames(dat),phenotypes$Genotype_code),]$grain_yield_15

cor(dat$grain_yield_15,dat$WD_0712,use="complete.obs")
#0.05809461

cor(dat$grain_yield_15,dat$WD_0718,use="complete.obs")
#0.05352837

cor(dat$grain_yield_15,dat$WD_0720,use="complete.obs")
#-0.05156357

cor(dat$grain_yield_15,dat$WD_0727,use="complete.obs")
# -0.004765058


# Not correlated

##### What about for just the top 5k expressed genes (excluding FT genes)
ftgenes=ft_df$Gene_ID
dat=data.frame(tkw_15=y,stringsAsFactors=F)
for(time1 in times){
	exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted_FIXED.txt',time1),data.table=F)
	#exp=exp[,c('V1',top5k)]
	rownames(exp)=exp$V1
	exp=exp[,-1]
	tmpdf=fread(sprintf('eqtl/results/eQTL_%s_freq_chi_data.txt',time1),data.table=F)
	tmpdf$gene_time_snp=paste0(tmpdf$Trait,'-',tmpdf$time,'-',tmpdf$X_ID)
	#top5k=unique(tmpdf[tmpdf$rank<=5000,]$Trait)
	top5k=unique(tmpdf[tmpdf$rank<=5000,]$Trait)
	exp=exp[,top5k]
	exp=exp[,!(names(exp) %in% ftgenes)]
	unlog=data.frame(lapply(exp,log_inverse),stringsAsFactors=F)
	avg_exp = apply(unlog,2,mean)
	names(avg_exp)=names(exp)
	avg_exp[avg_exp<1]=0
	# re-log the input data
	avg_logexp=log2(avg_exp)
	avg_logexp[is.infinite(avg_logexp)]=0
	# For each individual, get the average squared diff from the avg expression
	# Unlog expression for both, get the average difference, and then relog
	dev=sapply(seq(1,length(avg_exp)),function(x) (unlog[,x]-avg_exp[x])**2)
	rownames(dev)=rownames(exp)
	colnames(dev)=colnames(exp)
	# relog log2cpm
	di=apply(dev,MARGIN=1,function(x) log10(mean(x,na.rm=T)))
	
	dat[,time1]=di[match(rownames(dat),names(di))]
}

cor(dat$tkw_15,dat$WD_0712,use="complete.obs")
# 0.05054033
cor(dat$tkw_15,dat$WD_0718,use="complete.obs")
# 0.04681356
cor(dat$tkw_15,dat$WD_0720,use="complete.obs")
# 0.09486369
cor(dat$tkw_15,dat$WD_0727,use="complete.obs")
#0.1201153

dat$grain_yield_15 = phenotypes[match(rownames(dat),phenotypes$Genotype_code),]$grain_yield_15

cor(dat$grain_yield_15,dat$WD_0712,use="complete.obs")
#0.05851975

cor(dat$grain_yield_15,dat$WD_0718,use="complete.obs")
#0.05301685

cor(dat$grain_yield_15,dat$WD_0720,use="complete.obs")
#-0.04998232

cor(dat$grain_yield_15,dat$WD_0727,use="complete.obs")
#-0.004644701


datmelt=melt(dat[,c('tkw_15',times)],'tkw_15')

p1=ggplot(aes(x=value,y=tkw_15),data=datmelt) + geom_point(aes(color=variable)) 

png(sprintf('QTT/%s_%s_disregulation_score_totexp_plot.png',pheno,env))
print(p1)
dev.off()

#### What if I exclude FT genes? Does this change anything?




##### glmnet ######
y=unlist(y)

fit=glmnet(impx,y,alpha=0)

png(sprintf('QTT/%s_%s_high_gerp_disregulation_glmnet.png',pheno,env))
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

lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)

summary(ridge_reg)


### Optimal lambda value

cv_ridge <- cv.glmnet(x, y_train, alpha = 0, lambda = lambdas)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda

# 63.09573


#### Predict values #####

# Compute R^2 from true and predicted values
 
#cvfit <- cv.glmnet(impx, y,nfolds=10)

#pdf(sprintf('QTT/%s_%s_disregulation_glmnet_CV.pdf',pheno,env))
#print(plot(cvfit))
#dev.off()

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

