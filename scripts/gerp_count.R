#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('ggplot2')
library('stringr')


chr=10

allgerp=c()
for(chr in 1:10){
	gerp=fread(sprintf('datasets/gerp_v4/Zea_mays.v4.chr%s.gerpcolRatesPlusOGs.txt.gz',chr),data.table=F)
	split1=sapply(seq(1,nrow(gerp)),function(x) as.numeric(strsplit(gerp$V1[x],' ')[[1]][2]))
	gerp$POS=as.numeric(split1)
	gerp$CHR=chr
	names(gerp)=c('SITE','NR','RS','Sor','Set','POS','CHR')
	bbreaks=c(min(gerp$RS)-0.1,0,2,4,max(gerp$RS)+0.1)
	gerp$bin=cut(gerp$RS, breaks=bbreaks,labels=F)
	founder=fread(sprintf('datasets/hapmap/chr%s_founder_rare_alleles.txt',chr),data.table=F)
	names(founder)=c("SNP","CHR","POS","REF","ALT","A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra","DK63","EP1_inra","F492", "FV252_inra", "FV2_inra","MBS847","ND245","OH43_inra","VA85","W117_inra")
	#overlap=merge(rares,hapsnps,by='SNP')
	founders=c("MBS847","B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
	founder=founder[,c("SNP","CHR","POS","REF","ALT",founders)]
	founder[,founders]=apply(founder[,founders],MARGIN=2,function(x) ifelse(x=="0/0",0,ifelse(x=="1/1",1,ifelse(x=="2/2",2,NA))))
	# remove multiallelic sites?
	founder$ALT2=NA
	mult=which(grepl(',',founder$ALT))
	alts=founder$ALT[mult]
	alt1=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[1]])
	alt2=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[2]])
	founder[mult,]$ALT=alt1
	founder[mult,]$ALT2=alt2
	founder=founder[,c("SNP","CHR","POS","REF","ALT","ALT2",founders)]
	founder=merge(founder,gerp,by='POS')
	allgerp=rbind(allgerp,founder)
}

allgerp=as.data.frame(allgerp)
fwrite(allgerp,'datasets/all_founder_rare_alleles_GERP_scores.txt',row.names=F,quote=F,sep='\t')
#### Classify GERP scores based on ?
#−2< GERP ≤0, nearly neutral; 0< GERP ≤2, slightly deleterious; 2< GERP ≤4, moderately deleterious; GERP >4, strongly deleterious

#Rejected substitutions (RS) defined as the number of substitutions expected under neutrality minus the number of substitutions “observed” at the position.
# Thus, positive scores represent a substitution deficit (which would be expected for sites under selective constraint),
# while negative scores represent a substitution surplus.

# How many of these are in the founders?



table(allgerp$bin)

#     1      2      3      4 
# 68640  53605  31878 112248

# 266371 total

# 42% predicted to be strongly deleterious

### Look at disregulation broken up by GERP?


### Does dysregulation at genes with more strongly deleterious genes better predict yield?


### For each gene - how many sites does it have within a 10kb window that are strongly or moderately deleterious?