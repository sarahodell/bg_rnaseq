#!/usr/bin/env Rscript

library('data.table')
library('dplyr')

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

bounds=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)
bounds$gene_snp=paste0(bounds$gene,'-',bounds$SNP)

bounds=merge(bounds,genetable,by.x='gene',by.y='Gene_ID')

eqtl=fread('eqtl/results/all_cis_eQTL_weights_fdr_hits_FIXED.txt',data.table=F)
eqtl$gene_time=paste0(eqtl$Trait,'-',eqtl$time)
# Grab only the highest cis SNP
eqtl2= eqtl %>% group_by(gene_time) %>% slice(which.max(value))
eqtl=as.data.frame(eqtl2)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}
eqtl$block_start=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$start
eqtl$block_end=all_founder_blocks[match(eqtl$X_ID,all_founder_blocks$focal_snp),]$end




all_cors=c()
for(time in times){
	for(chr in chroms){
		cors=fread(sprintf('eqtl/trans/results/eQTL_%s_c%s_correlations_FIXED.txt',time,chr),data.table=F)
		all_cors=rbind(all_cors,cors)
	}
}
all_cors=as.data.frame(all_cors,stringsAsFactors=F)
# 8,849 rows
# 3,577 genes

tsus=all_cors[complete.cases(all_cors),]

tsus$chr_pair=interaction(tsus$CHR.x,tsus$CHR.y)
tsus$r2=0


chrpairs=unique(interaction(tsus$CHR.x,tsus$CHR.y))

newt=c()
for(cp in chrpairs){
	subtsus=tsus[tsus$chr_pair==cp,]
	c1=unique(subtsus$CHR.x)
	c2=unique(subtsus$CHR.y)
	ldfile=sprintf('ld/snp/SNP_ld_c%s_c%s_matrix.txt',c1,c2)
	if(file.exists(ldfile)){
		ld=fread(ldfile,data.table=F)
		r2=sapply(seq(1,nrow(subtsus)),function(x) get_r2(subtsus$SNP.x[x],subtsus$SNP.y[x]))
	}else{
		tmp=c1
		c1=c2
		c2=tmp
		ldfile=sprintf('ld/snp/SNP_ld_c%s_c%s_matrix.txt',c1,c2)
		ld=fread(ldfile,data.table=F)
		r2=sapply(seq(1,nrow(subtsus)),function(x) get_r2(subtsus$SNP.y[x],subtsus$SNP.x[x]))
	}
	subtsus$r2=r2
	newt=rbind(newt,subtsus)
}

fwrite(newt,'eqtl/results/trans_pair_high_cor.txt',row.names=F,quote=F,sep='\t')


# 8,504 rows
# 3,372 genes

same=tsus[tsus$CHR.x==tsus$CHR.y,]
same$distance=abs(same$BP.y-same$BP.x)
summary(same$distance/1e6)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  25.00   28.90   40.76   48.46   61.79  261.70

# 7820 rows on same chr

# 684 on different chr

diff=tsus[tsus$CHR.x!=tsus$CHR.y,]
difft=as.data.frame(table(interaction(diff$CHR.x,diff$CHR.y)))
summary(difft$Freq)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00    5.00    7.00    6.84    9.00   16.00

sus=all_cors[!complete.cases(all_cors),]

sus=sus[,c('gene','time.y','CHROM.y','START.y','END.y','name2','r','CHR.y','SNP.y','BP.y','left_bound_bp.y','right_bound_bp.y')]
names(sus)=c('gene','time','CHROM','START','END','transeqtl','r','CHR','SNP','BP','left_bound_bp','right_bound_bp')

sus$gene_time=paste0(sus$gene,'-',sus$time)
sus$cis_SNP=eqtl[match(sus$gene_time,eqtl$gene_time),]$X_ID
sus$cis_BP=eqtl[match(sus$gene_time,eqtl$gene_time),]$BP

fwrite(sus,'eqtl/results/cis_trans_high_cor.txt',row.names=F,quote=F,sep='\t')
sus$chr_pair=interaction(sus$CHROM,sus$CHR)
sus$r2=0


chrpairs=unique(interaction(sus$CHROM,sus$CHR))

get_r2=function(x1,x2){
	if(x1 %in% ld$V1 & x2 %in% colnames(ld)){
		r2=ld[ld$V1==x1,x2]
		if(is.na(r2)){
			r2=0
		}
	}else{
		r2=NA
	}
	return(r2)
}

get_fr2=function(x1,x2){
	r2=fld[fld$SNP_A==x1 & fld$SNP_B==x2,]$r2
	return(r2)
}

fld=fread('eqtl/data/all_founder_ld.txt',data.table=F)
chrpairs=unique(interaction(sus$CHROM,sus$CHR))

newsus=c()
for(cp in chrpairs){
	subsus=sus[sus$chr_pair==cp,]
	c1=unique(subsus$CHROM)
	c2=unique(subsus$CHR)
	ldfile=sprintf('ld/snp/SNP_ld_c%s_c%s_matrix.txt',c1,c2)
	if(file.exists(ldfile)){
		ld=fread(ldfile,data.table=F)
		r2=sapply(seq(1,nrow(subsus)),function(x) get_r2(subsus$cis_SNP[x],subsus$SNP[x]))
	}else{
		tmp=c1
		c1=c2
		c2=tmp
		ldfile=sprintf('ld/snp/SNP_ld_c%s_c%s_matrix.txt',c1,c2)
		ld=fread(ldfile,data.table=F)
		r2=sapply(seq(1,nrow(subsus)),function(x) get_r2(subsus$SNP[x],subsus$cis_SNP[x]))
	}
	fr2=sapply(seq(1,nrow(subsus)),function(x) get_fr2(subsus$cis_SNP[x],subsus$SNP[x]))
	subsus$r2=r2
	subsus$fr2=fr2
	newsus=rbind(newsus,subsus)
}

fwrite(newsus,'eqtl/results/cis_trans_high_cor.txt',row.names=F,quote=F,sep='\t')

# 345 rows
# 282 genes

# 54 instances of high correlation between 4 and 5
# 14 instances of high correlation between 9 and 10


