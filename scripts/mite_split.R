#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('stringr')

pheno="male_flowering_d6"
env="ALL"
chr="8"


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
mite=c("NO_MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","NO_MITE")
mite_late=c("FV252_inra","C103_inra","F492","A632_usa")
# Read in Kinship Matrix
K=fread(sprintf('../GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)

phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
data=data[data$Loc.Year.Treat==env,]
data=data[!is.na(data$y),]
data$y = (data$y - mean(data$y))/sd(data$y)


mite_founder=fread('metadata/vgt1_max_prob_founder.txt',data.table=F)

mite_founder$mite=mite[match(mite_founder$founder,founders)]


mite_prob=fread('../GridLMM/mite_probabilities.txt',data.table=F)
mite_founder$final=mite_prob[match(mite_founder$ID,mite_prob$ID),]$final
rownames(mite_founder)=mite_founder$ID


mite_prob=mite_prob[data$ID,]

has_mite=mite_prob[mite_prob$final>=0.9,]$ID

adj_chroms=c("5","9")
if(chr %in% adj_chroms){
		X_list=readRDS(sprintf('phenotypes/bg%s_adjusted_genoprobs.rds',chr))
}else{
	X_list=readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))

}
inds=rownames(X_list[[1]])
data=data[data$ID %in% has_mite,]

inter=intersect(data$ID,inds)
X_list=lapply(X_list,function(x) x[inter,])
#founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


#Make B73 the first in the list so that it is the one that is dropped
#names(X_list)=founders
#X_list=X_list[new_founders]
K=K[inter,inter]
rownames(data)=data$ID
data=data[inter,]



pheno="male_flowering_d6"
chr="8"


library('data.table')
library('dplyr')
library('lme4')
library('ggplot2')
library('reshape2')
library('tibble')
library('dplyr')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
mite=c("NO_MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","MITE","MITE","MITE","NO_MITE","MITE","MITE","MITE","NO_MITE")

mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)
has_mite=which(mite==T,mite)
no_mite=which(mite==F,mite)

bg=readRDS('../genotypes/probabilities/geno_probs/raw/bg8_genoprobs.rds')
dimnames(bg[[1]])[[2]]=founders
pmap=fread('../genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)


founder_probs = readRDS(sprintf('../genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
#names(founder_probs)=founders
mite_start=135.947816*1e6
mite_end=135.946644*1e6


#PZE-108076585   8 135943387
#28984   AX-91102970   8 135947736
#28985   AX-91102985   8 135947876
#28986   AX-91102982   8 135948460
m2=c('PZE-108076585','AX-91102970','AX-91102985','AX-91102982')

X3=do.call(cbind,lapply(bg,function(x) x[,,'AX-91102985']))
mite_prob3=rowSums(X3[,has_mite])
f3=unlist(apply(X3,MARGIN=1,function(x) names(which.max(x))))

#'AX-91769748'
X4=do.call(cbind,lapply(bg,function(x) x[,,'AX-91769748']))
mite_prob4=rowSums(X4[,has_mite])
f4=unlist(apply(X4,MARGIN=1,function(x) names(which.max(x))))


region=pmap[pmap$pos>135800000 & pmap$pos<136100000,]$marker
markers=dimnames(founder_probs[[1]])[[2]]
markers=markers[markers %in% region]
markers=c("AX-91102912","AX-91103124")
#snp="AX-91102858"
sub8=lapply(founder_probs,function(x) x[,markers])
l=dim(sub8[[1]])[1]
mitefounders=founder_probs[c(has_mite)]
mite_prob = rowSums(do.call(cbind,lapply(mitefounders,function(x) x[,markers[1]])))
#mite_prob=ifelse(mite_prob>=0.85,1,0)
mite_prob2 = rowSums(do.call(cbind,lapply(mitefounders,function(x) x[,markers[2]])))

t=data.frame(ID=names(mite_prob),m1=mite_prob,m2=mite_prob2,stringsAsFactors=F)
t$m3=mite_prob3[match(t$ID,names(mite_prob3))]
t$m4=mite_prob4[match(t$ID,names(mite_prob4))]

t$final=rowMeans(t[,c('m3','m4')])

#'EB.09S.H.00072','EB.09S.H.00052'

t=t[,c('ID','m3','m4','final')]
names(t)=c('ID','AX-91102985','AX-91769748','final')
t$founder_AX91102985=f3[match(t$ID,names(f3))]
t$founder_AX91769748=f4[match(t$ID,names(f4))]

#t=t[,c('ID','final')]
fwrite(t,'phenotypes/mite_probabilities.txt',quote=F,row.names=F,sep='\t')

poi=c('EB.09S.H.00052','EB.09S.H.00072','EB.10H.H.00060','EB.09S.H.00494')

#mite_prob2=ifelse(mite_prob2>=0.85,1,0)
#mite_prob=rowSums(lapply(mitefounders,function(x) x))
#Probability that each individual has the MITE
#mite_prob=sapply(seq(1,length(region)), function(x) sapply(seq(1,l), function(i) sum(sub8[i,has_mite,x])))
#rownames(mite_prob)=dimnames(founder_probs[[1]])[[1]]
#colnames(mite_prob)=region

marker='AX-91102970'
mprob=rowSums(bg[[1]][,has_mite,marker])
t2=data.frame(ID=names(mprob),final=unlist(mprob),stringsAsFactors=F)
fwrite(t2,'GridLMM/mite_probabilities.txt',quote=F,row.names=F,sep='\t')

#at_mite=as.data.frame(mite_prob[,marker])
#at_mite$ID=rownames(at_mite)
#names(at_mite)=c(marker,'ID')
#at_mite=at_mite[,c('ID',marker)]

#fwrite(at_mite,'GridLMM/mite_probabilities.txt',quote=F,row.names=F,sep='\t')
mitebeta=
#qDTA8_1
# STPAUL
results=fread('QTL/MITE_only/results/Biogemma_chr8_male_flowering_d6_x_STPAUL_2017_WD_MITE_only_founderprobs.txt',data.table=F)
snp="AX-91069527"
result=results[results$X_ID==snp,]    
betas=unlist(result[founders])
betas[-1]=betas[1]+betas[-1]

# OH43 allele is early
# VA85 allele is early
# B73 allele is ~0  
# B96 allele is very late

# F252 allele is late
# A632 allele is slightly early
# F492 allele is slightly early
# C103 allele is very early

# CO255 allele is early
# A654 allele is slightly early
# F2 allele is early
# EP1 allele is ~0
# D105 allele is ~0
# W117 allele is ~0
# DK63 allele is late
# ND245 allele is early

# Is there an interaction between alleles at MITE and this SNP? 
# Can I only test this with SNPs, or can I use founders?

new_betas=c()
new8=qtlmite[qtlmite$CHR==8,]
pmap=fread('../genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)
for(i in 1:nrow(new8)){
	row=new8[i,]
	snp=row$SNP
	chr=row$CHR
	pheno=row$phenotype
	env=row$environment
	pos=pmap[pmap$marker==snp,]$pos
	results=fread(sprintf('QTL/MITE_only/results/Biogemma_chr%s_%s_x_%s_MITE_only_founderprobs.txt',chr,pheno,env),data.table=F)
	result=results[results$X_ID==snp,]   
	
	betas=unlist(result[founders])
	wn=which(!is.na(betas))[1]
	betas[-wn]=betas[wn]+betas[-wn] 
	#betas[-1]=betas[1]+betas[-1]
	line=c(pheno,env,snp,pos,betas)
	new_betas=rbind(new_betas,line)
}

new_betas=as.data.frame(new_betas,stringsAsFactors=F)
names(new_betas)=c('phenotype','environment','SNP','BP',founders)
row.names(new_betas)=seq(1,nrow(new_betas))


ibd=fread('../ibd_segments/refinedibd/600K/bg8_refined_ibd_blocks.txt',data.table=F)
pmap=fread('../genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)
pos=139651619
gene="Zm00001d011152"
snp="AX-91104035"
base=7
chr="8"
ibd[ibd$start<=pos & ibd$end>=pos,]
#h =15
h=15
xfile=sprintf('../genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h)
hapfile=readRDS(xfile)
markers=dimnames(hapfile[[1]])[[2]]
region=pmap[pmap$marker %in% markers,]
which.min(abs(region$pos-pos))
hsnp="AX-91103979"

# F252 and A632 are the same haplotype here


