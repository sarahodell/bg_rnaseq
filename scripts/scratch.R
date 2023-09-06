#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')




# What Factors match up for whole and residual MegaLMM?
times=c("WD_0712","WD_0718","WD_0720","WD_0727")

df=c()
for(time in times){
	lambda_r=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_Lambda_means_FIXED.txt',time),data.table=F)
	rownames(lambda_r)=lambda_r$V1
	lambda=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means_FIXED.txt',time),data.table=F)
	rownames(lambda)=lambda$V1
	ginter=intersect(lambda$V1,lambda_r$V1)
	lambda=lambda[ginter,]
	lambda_r=lambda_r[ginter,]
	factors=names(lambda)[-1]
	factors_r=names(lambda_r)[-1]
	#factorsm=which.max(c(length(factors),length(factors_r)))
	#if(factorsm==1){
	#	factorlist=factors
	#}else{
	#	factorlist=factors_r
	#}
	for(factor in factors){
		r1=lambda[,c('V1',factor)]
		rownames(r1)=r1$V1
		fmax_cor=0
		fmax=0
		for(i in factors_r){
			column=lambda_r[,i]
			r=cor(r1[,factor],column)
			if(abs(r)>abs(fmax_cor)){
				fmax_cor=r
				fmax=i
			}
		}
		line=data.frame(time=time,factor=factor,factor_r=fmax,fmax_cor=fmax_cor,stringsAsFactors=F)
		df=rbind(df,line)
	}
}


fwrite(df,'MegaLMM/MegaLMM_whole_residual_factor_pairs.txt',row.names=F,quote=F,sep='\t')

df2=c()
# What factors match up across time points?
for(t1 in 1:3){
	time1=times[t1]
	for(t2 in (t1+1):4){
		time2=times[t2]
		lambda1=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means_FIXED.txt',time1),data.table=F)
		rownames(lambda1)=lambda1$V1
		lambda2=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means_FIXED.txt',time2),data.table=F)
		rownames(lambda2)=lambda2$V1
		ginter=intersect(lambda1$V1,lambda2$V1)
		lambda1=lambda1[ginter,]
		lambda2=lambda2[ginter,]
		
		factors1=names(lambda1)[-1]
		factors2=names(lambda2)[-1]

		for(factor in factors1){
			r1=lambda1[,c('V1',factor)]
			rownames(r1)=r1$V1
			fmax_cor=0
			fmax=0
			for(i in factors2){
				column=lambda2[,i]
				r=cor(r1[,factor],column)
				if(abs(r)>abs(fmax_cor)){
					fmax_cor=r
					fmax=i
				}
			}
			line=data.frame(time1=time1,factor1=factor,time2=time2,factor2=fmax,fmax_cor=fmax_cor,stringsAsFactors=F)
			df2=rbind(df2,line)
		}
	}
}

fwrite(df2,'MegaLMM/MegaLMM_timepoint_factor_pairs.txt',row.names=F,quote=F,sep='\t')


# Across timepoints for residuals
df3=c()
# What factors match up across time points?
for(t1 in 1:3){
	time1=times[t1]
	for(t2 in (t1+1):4){
		time2=times[t2]
		lambda1=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_Lambda_means_FIXED.txt',time1),data.table=F)
		rownames(lambda1)=lambda1$V1
		lambda2=fread(sprintf('MegaLMM/MegaLMM_%s_residuals_all_Lambda_means_FIXED.txt',time2),data.table=F)
		rownames(lambda2)=lambda2$V1
		ginter=intersect(lambda1$V1,lambda2$V1)
		lambda1=lambda1[ginter,]
		lambda2=lambda2[ginter,]
		
		factors1=names(lambda1)[-1]
		factors2=names(lambda2)[-1]

		for(factor in factors1){
			r1=lambda1[,c('V1',factor)]
			rownames(r1)=r1$V1
			fmax_cor=0
			fmax=0
			for(i in factors2){
				column=lambda2[,i]
				r=cor(r1[,factor],column)
				if(abs(r)>abs(fmax_cor)){
					fmax_cor=r
					fmax=i
				}
			}
			line=data.frame(time1=time1,factor1=factor,time2=time2,factor2=fmax,fmax_cor=fmax_cor,stringsAsFactors=F)
			df3=rbind(df3,line)
		}
	}
}

fwrite(df3,'MegaLMM/MegaLMM_timepoint_residuals_factor_pairs.txt',row.names=F,quote=F,sep='\t')


chr="8"
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
genetable=genetable[genetable$CHROM==chr,]
founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
pmap=fread('../genotypes/qtl2/startfiles/Biogemma_pmap_c8.csv',data.table=F)
ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
ft_genes=unique(ft_genelist$V4)

prop_var=fread('MegaLMM/MegaLMM_WD_0727_prop_variance.txt',data.table=F)
fgenes=prop_var[prop_var$Factor2>=0.1,'V1']

#### 
f2_27=fread('eqtl/results/Factor2_trans_WD_0727_eQTL_fkeep_hits.txt',data.table=F)
fsnp="AX-91107495"
f2_27$block_start=founder_blocks[match(f2_27$SNP,founder_blocks$focal_snp),]$start
f2_27$block_end=founder_blocks[match(f2_27$SNP,founder_blocks$focal_snp),]$end

env1=f2_27
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(genetable)
#env2$end=env2$end-1
setkey(env2,START,END)
comparison=foverlaps(env1,env2,by.x=c('block_start','block_end'),by.y=c('START','END'),nomatch=NULL)
which(comparison$Gene_ID %in% fgenes)



# 1) How are factor effect sizes correlated with qDTA8 effect sizes?
qtl=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_hits.txt',data.table=F)
qtl$block_start=founder_blocks[match(qtl$SNP,founder_blocks$focal_snp),]$start
qtl$block_end=founder_blocks[match(qtl$SNP,founder_blocks$focal_snp),]$end

line=data.frame(start=145758585,end=150866243)
# 8 145758585 150866243

env1=line
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(genetable)
#env2$end=env2$end-1
setkey(env2,START,END)
comparison2=foverlaps(env1,env2,by.x=c('start','end'),by.y=c('START','END'),nomatch=NULL)
which(comparison2$Gene_ID %in% fgenes)

comparison2[c(3,75),]
#          Gene_ID CHROM     START       END     start       end
#1: Zm00001d011285     8 145825712 145826506 145758585 150866243
#2: Zm00001d011377     8 148947109 148951347 145758585 150866243

# highest SNP
line2=data.frame(start=146188918,end=146530752)
env1=line2
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(genetable)
#env2$end=env2$end-1
setkey(env2,START,END)
comparison3=foverlaps(env1,env2,by.x=c('start','end'),by.y=c('START','END'),nomatch=NULL)
which(comparison3$Gene_ID %in% fgenes)
#          Gene_ID CHROM     START       END     start       end
#1: Zm00001d011308     8 146317656 146322089 146188918 146530752
prop_var[prop_var$V1=="Zm00001d011308",'Factor2']
#[1] 0.5839009

env1=qtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(genetable)
#env2$end=env2$end-1
setkey(env2,START,END)
comparison4=foverlaps(env1,env2,by.x=c('block_start','block_end'),by.y=c('START','END'),nomatch=NULL)
which(comparison4$Gene_ID %in% fgenes)

which(unique(comparison4$Gene_ID) %in% fgenes)
#  19  33  58 116
which(unique(comparison4$Gene_ID) %in% ft_genes)
# 3 73

 
results=fread('QTL/Biogemma_chr8_male_flowering_d6_x_EXP_STPAUL_2017_WD_vst_founderprobs.txt',data.table=F)
results$value=-log10(results$p_value_ML)
results$pos=pmap[match(results$X_ID,pmap$marker),]$pos
results=results[order(results$pos),]
which(results$X_ID==snp)
which.max(results$value)
betas=unlist(unname(results[results$X_ID==fsnp,c(6:21)]))
betas[-1]=betas[1]+betas[-1]

results[191:203,]

chr="8"
#snp 8:152630937-153132010
# QTL highest SNP 8:146188918-146530752

fresults=fread('eqtl/trans/results/WD_0727_c8_Factor2_trans_results.txt',data.table=F)
fresults$value=-log10(fresults$p_value_ML)
fbetas=unlist(unname(fresults[fresults$X_ID==fsnp,c(6,8:22)]))
fbetas[-1]=fbetas[1]+fbetas[-1]

# highest QTL SNP correlation with highest Factor2 trans-eQTL effect size is 0.1677
# highest Factor2 trans-eQTL effect size correlation with flowering time effet size for that
# SNP is 0.2540641



#### Overlap of trans-eQTL and interLD regions
all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}



ld=fread('../stats/ld_decay/circos/ld_bundled_links_filtered.txt',data.table=F)
ld$V1=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld$V1[x],'chr')[[1]][2]))
ld$V4=sapply(seq(1,nrow(ld)),function(x) as.numeric(strsplit(ld$V4[x],'chr')[[1]][2]))
names(ld)=c('CHR_A','START_A','END_A','CHR_B','START_B','END_B','SIZE')
eqtl=fread('eqtl/results/cis_trans_eQTL_hits.txt',data.table=F)
#eqtl$BP_end=all_founder_blocks[match(eqtl$SNP,all_founder_blocks$focal_snp),]$end

# Overlap of LD A with eQTL variants
env1=eqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(ld)
#env2$end=env2$end-1
setkey(env2,CHR_A,START_A,END_A)
acomp=foverlaps(env1,env2,by.x=c('CHR','BP','BP_end'),by.y=c('CHR_A','START_A','END_A'),nomatch=NULL)

# Overlap of LD B with eQTL variants
env1=eqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(ld)
#env2$end=env2$end-1
setkey(env2,CHR_B,START_B,END_B)
bcomp=foverlaps(env1,env2,by.x=c('CHR','BP','BP_end'),by.y=c('CHR_B','START_B','END_B'),nomatch=NULL)

# Overlap of LD A with eQTL genes
env1=eqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(ld)
#env2$end=env2$end-1
setkey(env2,CHR_A,START_A,END_A)
agcomp=foverlaps(env1,env2,by.x=c('gene_chr','gene_start','gene_end'),by.y=c('CHR_A','START_A','END_A'),nomatch=NULL)

# Overlap of LD B with eQTL genes
env1=eqtl
#env1$BP_start=env1$BP-5000
#env1$BP_end=env1$BP+5000
env1=as.data.table(env1)
env2=as.data.table(ld)
#env2$end=env2$end-1
setkey(env2,CHR_B,START_B,END_B)
bgcomp=foverlaps(env1,env2,by.x=c('gene_chr','gene_start','gene_end'),by.y=c('CHR_B','START_B','END_B'),nomatch=NULL)


# "Zm00001d035217" "Zm00001d047592" variants are in inter LD regions -

# Variant on chr 5 is in LD with regions on chr 1,2,3,4 ,6, 7, 8, 9, and 10 - ALL chromosome!
# for Zm00001d035217, there is not high LD in the location of the gene, but on either side of it
# This region in chr 5 is in LD with region on chr 8 that are near zmrap2.7 and the factor 2 eqtL variant
# ZmRap2.7 is : (Chr8: 136007716..136013584) - in v5 it is collapsed into 



######## Enrichment of FT genes in Factors

ft_genelist=fread('../selection/FT_gene_list_AGPv4.bed',data.table=F)
# Enrichment of flowering time genes?
cutoff=0.1
test=subvar[subvar$Factor2>=cutoff,]$V1
gtable=genetable[genetable$Gene_ID %in% test,]


find_nearest_snp=function(row){
    index=which.min(abs(row$START-pmap[pmap$chr==row$CHROM,]$pos))
    return(pmap[index,]$marker)
}

nearest_snps=sapply(seq(1,nrow(gtable)),function(x) find_nearest_snp(gtable[x,]))
gtable$SNP=nearest_snps

snplist=gtable[,'SNP',drop=F]
fwrite(snplist,'WD_0727_Factor2_snplist.txt',row.names=F,col.names=F,sep='\t',quote=F)

true_nft=length(intersect(ft_genelist$V4,test))

ngenes=length(test)
allgenes=length(genes)
nfts=c()
for(i in 1:1000){
	draw=sample(seq(1,allgenes),ngenes)
	nft=intersect(ft_genelist$V4,genes[draw])
	nfts=c(nfts,length(nft))
}

# 0.1 cutoff 60 genes
#summary(nfts)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  39.00   53.75   58.00   57.36   61.00   74.00 
# quantile(nfts,0.95)
#95% 
# 67 

# 0.5 cutoff 5 genes
summary(nfts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   7.000   9.000   9.166  11.000  23.000 
#quantile(nfts,0.05)
#5% 
# 4

################## Pair factors across time points
library('data.table')


times=c("WD_0718","WD_0712","WD_0720","WD_0727")

all_genes=c()
for(time in times){
	#exp=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
	geneh2s=fread(sprintf('eqtl/data/lme4qtl_%s_h2s.txt',time),data.table=F)
	kept_genes=geneh2s[geneh2s$h2>0 ,]$gene
	#exp=exp[,c(kept_genes)]
	print(length(all_genes))
	#current=kept_genes
	if(time==times[1]){
		all_genes=kept_genes
	}else{
		all_genes=intersect(all_genes,kept_genes)
	}
}

coln=c()
for(i in 1:3){
	for(j in (i+1):4){
		comb=paste0(times[i],'-',times[j])
		combr=paste0(times[i],'-',times[j],'_r2')
		coln=c(coln,c(comb,combr))
	}
}


# 582 genes shared across the 4 timepoints

data=data.frame(matrix(ncol = 13,nrow=22))
names(data)=coln
data$WD_0718_factor=c( "Factor1" , "Factor2" , "Factor3" , "Factor4" , "Factor5"  ,"Factor6" ,
"Factor7" , "Factor8" , "Factor10", "Factor11", "Factor12", "Factor13",
"Factor14", "Factor16", "Factor17" ,"Factor19", "Factor20" ,"Factor21",
"Factor22" ,"Factor23", "Factor24" ,"Factor26")
rownames(data)=data$WD_0718_factor

for(i in 1:3){
	time1=times[i]
	lambda1=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means.txt',time1),data.table=F)
	fdf1=fread(sprintf('MegaLMM/MegaLMM_%s_factor_correlations.txt',time1),data.table=F)
	count=c()
	for(k in 1:nrow(fdf1)){
  		count=c(count,sum(abs(fdf1[k,c(4:6)])>=0.9,na.rm=T))
	}
	strong=which(count==3)
	factors1=fdf1[strong,]$r1_factor_1
	rownames(lambda1)=lambda1$V1
	lambda1=lambda1[,-1]
	lambda1=lambda1[all_genes,factors1]
	
	data=data.frame(matrix(ncol = 13,nrow=length(factors1)))
	facol=paste0(time1,'_factor')
	names(data)=c(facol,coln)
	data[,facol]=factors1
	rownames(data)=factors1

	for(j in (i+1):4){
		time2=times[j]
		comb=paste0(time1,'-',time2)
		combr=paste0(time1,'-',time2,'_r2')
		
		
		lambda2=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means.txt',time2),data.table=F)
		fdf2=fread(sprintf('MegaLMM/MegaLMM_%s_factor_correlations.txt',time2),data.table=F)
		count=c()
		for(k in 1:nrow(fdf2)){
  			count=c(count,sum(abs(fdf2[k,c(4:6)])>=0.9,na.rm=T))
		}
		strong=which(count==3)
		factors2=fdf2[strong,]$r1_factor_1
		rownames(lambda2)=lambda2$V1
		lambda2=lambda2[,-1]
		lambda2=lambda2[all_genes,factors2]
		
		for(l in factors1){
  			max_factor=NA
  			max_cor=0
  			for(m in factors2){
    			cur_cor=cor(lambda1[,l],lambda2[,m],use="complete.obs")
    			if(abs(cur_cor)>abs(max_cor)){
      				max_cor=cur_cor
      				max_factor=m
    			}
    			
  			}
  			data[l,combr]=max_cor
  			data[l,comb]=max_factor
		}
		fwrite(data,sprintf('MegaLMM/MegaLMM_%s_Lambda_correlations.txt',time1),row.names=F,quote=F,sep='\t')
	}
	
}


### What factors do eQTL load on in their timepoints?


allhits=fread('eqtl/results/all_eQTL_hits.txt',data.table=F)
eqtl=allhits[allhits$class!='factor',]

for(time in unique(eqtl$time)){
	prop_var=fread(sprintf('MegaLMM/MegaLMM_%s_prop_variance.txt',time),data.table=F)
	subqtl=eqtl[eqtl$time==time,]
	subgenes=unique(subqtl$Trait)
	print(time)
	for(gene in subgenes){
		#gene=subgenes[i]
		infactor=names(prop_var[,-1])[which(prop_var[prop_var$V1==gene,-1]>0.1)]
		print(gene)
		print(infactor)
	}
}

# "WD_0712"

#[1] "Zm00001d031961"
#[1] "Factor34" "Factor38"

#[1] "Zm00001d012882"
#[1] "Factor3"  "Factor38"

#[1] "Zm00001d025548"
#[1] "Factor1"

#[1] "Zm00001d039372"
#[1] "Factor1"  "Factor13"

#[1] "Zm00001d039421"
#[1] "Factor1"  "Factor3"  "Factor19" "Factor38" "Factor40"

#[1] "Zm00001d046047"
#[1] "Factor1"  "Factor38"


# "WD_0720"

#[1] "Zm00001d025017"
#[1] "Factor14"

#[1] "Zm00001d039315"
#[1] "Factor8"  "Factor12" "Factor15"


# "WD_0727"

#[1] "Zm00001d006483"
#[1] "Factor7"  "Factor9"  "Factor12"

#[1] "Zm00001d035217"
#[1] "Factor14"

#[1] "Zm00001d047592"
#[1] "Factor7"  "Factor14"


#"WD_0718"

#[1] "Zm00001d022126"
#[1] "Factor3"

#[1] "Zm00001d042747"
#[1] "Factor3"  "Factor4"  "Factor11"

#[1] "Zm00001d051801"
#[1] "Factor3" "Factor4" "Factor7"


