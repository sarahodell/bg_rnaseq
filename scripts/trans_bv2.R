#!/usr/bin/env Rscript

#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
time=as.character(args[[1]])
t_chr=as.character(args[[2]])
cores=as.numeric(args[[3]])

library('data.table')
library('dplyr')
library('parallel')
library('MASS')
library('stringr')
# Get cis fitted values
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

bounds=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)
bounds=bounds[bounds$time==time,]
bounds$gene_snp=paste0(bounds$gene,'-',bounds$SNP)

bounds=merge(bounds,genetable,by.x='gene',by.y='Gene_ID')
bounds=bounds[bounds$CHROM==t_chr,]
#trans=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_weights_results_FIXED.rds',time,t_chr))
#tgenes=unlist(lapply(trans,function(x) unique(x$Trait)))
#genetable=genetable[genetable$Gene_ID %in% tgenes,]
#genetable=genetable[genetable$CHROM==t_chr,]
#cgenes=unique(genetable$Gene_ID)

samples=fread(sprintf('eqtl/data/%s_samples_FIXED.txt',time),data.table=F)
inter=samples$id

#transbv=readRDS(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_bvs_FIXED.rds',time,t_chr))
cisbv=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_bv_FIXED.txt',time,t_chr),data.table=F)
rownames(cisbv)=cisbv$V1
cisbv=cisbv[,-1]
cisbv=cisbv[inter,]

tgenes=unique(bounds$gene)


get_bvs=function(rep){
	gene=tgenes[rep]
	subbounds=bounds[bounds$gene==gene,]
	bvs=cisbv[,gene,drop=F]
	tchroms=unique(subbounds$CHR)
	for(chr in tchroms){
		transbv=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%s_bv_FIXED.txt',time,chr),data.table=F)
		pos=which(grepl(gene,colnames(transbv)))
		transbv=transbv[,pos,drop=F]
		bvs=cbind(bvs,transbv)
	}
	df=c()
	for(i in 1:(ncol(bvs)-1)){
		name1=colnames(bvs)[i]
		for(j in (i+1):ncol(bvs)){
			name2=colnames(bvs)[j]
			r=cor(bvs[,i],bvs[,j])
			line=data.frame(gene=gene,name1=name1,name2=name2,r=r,stringsAsFactors=F)
			df=rbind(df,line)
		}	
	}
	return(df)
}

n_reps=1:length(tgenes)
#n_reps=1:5

print(system.time({
full_results=mclapply(n_reps,get_bvs,mc.cores=cores)
}))

all_cors=rbindlist(full_results)
all_cors=as.data.frame(all_cors,stringsAsFactors=F)

cutoff=0.5
all_cors=all_cors[all_cors$r>=cutoff,]

all_cors=merge(all_cors,bounds,by.x='name1',by.y='gene_snp',all.x=TRUE)
all_cors=merge(all_cors,bounds,by.x='name2',by.y='gene_snp',all.x=TRUE)

if(nrow(all_cors)!=0){
	fwrite(all_cors,sprintf('eqtl/trans/results/eQTL_%s_c%s_correlations_FIXED.txt',time,t_chr),row.names=F,quote=F,sep='\t')
}
