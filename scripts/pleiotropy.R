#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

trans=fread('eqtl/results/all_trans_fdr_peaks.txt',data.table=F)
snp_info=trans %>% group_by(X_ID) %>% count
msnp=snp_info[snp_info$n>2,]$X_ID
mtrans=trans[trans$X_ID %in% msnp,]

genecors=as.data.frame(matrix(nrow=0,ncol=8))

cor_plots=list()
count=1

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


for(m in msnp){
	subdf=trans[trans$X_ID==m,]
	sgenes=unique(subdf$Trait)
	for(i in 1:(length(sgenes)-1)){
		gene1=sgenes[i]
		row1=subdf[subdf$Trait==gene1,]
		time1=row1$time
		chr1=row1$CHR
		results=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results.rds',time1,chr1))
		w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene1)))
		results=results[[w]]
		results=results[results$X_ID==m,]
		betas1=unlist(results[,c(6,10:24)])
		wn=which(!is.na(betas1))[1]
		betas1[-wn]=betas1[-wn]+betas1[wn]
		for(j in (i+1):length(sgenes)){
			gene2=sgenes[j]
			row2=subdf[subdf$Trait==gene2,]
			time2=row2$time
			chr2=row2$CHR
			results=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results.rds',time2,chr2))
			w=which(unlist(lapply(results,function(x) unique(x$Trait)==gene2)))
			results=results[[w]]
			results=results[results$X_ID==m,]
			betas2=unlist(results[,c(6,10:24)])
			wn=which(!is.na(betas2))[1]
			betas2[-wn]=betas2[-wn]+betas2[wn]
			
			test=cor.test(betas1,betas2,use="complete.obs")
			r=test$estimate
			p=test$p.value
			line=data.frame(gene1=gene1,gene2=gene2,time1=time1,time2=time2,chr=chr1,snp=m,r=r,pvalue=p,stringsAsFactors=F)
			genecors=rbind(genecors,line)
			
			pldf=data.frame(beta1=betas1,beta2=betas2,founder=founders)
			p1=ggplot(pldf,aes(x=beta1,y=beta2)) + geom_point(aes(color=founder)) +
			xlab(sprintf("%s Effect Size",gene1)) + ylab(sprintf('%s Effect Size',gene2)) +
			ggtitle(sprintf('%s Effect Size, r=%.3f',m,r))
			cor_plots[[count]]=p1
			count=count+1
		}		
	}
}

fwrite(genecors,'eqtl/results/trans_eQTL_effect_size_corr.txt',row.names=F,quote=F,sep='\t')

pdf('eqtl/images/trans_eQTL_effect_size_corr.pdf')
for(i in 1:length(cor_plots)){
	print(cor_plots[[i]])
}
dev.off()