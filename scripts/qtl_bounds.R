#!/usr/bin/env Rscript

library('data.table')
library('dplyr')



#qtl=fread('Biogemma_QTL.csv',data.table=F)
#base_list=c(7,7,8,6,7,6,7,7,6,8)
#qtl=fread('QTL/all_adjusted_QTL_peaks_trimmed.txt',data.table=F)
qtl=fread('QTL/MITE_only/all_MITE_only_QTL_peaks.txt',data.table=F)

#qtl=fread('Biogemma_10p_QTL.csv',data.table=F)
alt_right_bound_bp=c()
alt_right_bound_snp=c()

for (i in seq(1,dim(qtl)[1])){
	chr=qtl[i,]$CHR
    pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
    snp=qtl[i,]$SNP
    m=pmap[pmap$marker==snp,]$pos
    fdropped=readRDS(sprintf('../genotypes/probabilities/geno_probs/dropped/bg%.0f_dropped_markers_genoprobs.rds',chr))
    dropped_markers=fdropped[[which(unlist(lapply(fdropped,function(x) x$marker==snp)))]]$linked
    sub=pmap[pmap$marker %in% dropped_markers,]
    rownames(sub)=seq(1,dim(sub)[1])
    right=which.max(sub$pos)
    right_bound=sub[right,]$pos
    right_snp=sub[right,]$marker
    #print(c(i,length(right_bound)))
    alt_right_bound_bp=c(alt_right_bound_bp,right_bound)
    alt_right_bound_snp=c(alt_right_bound_snp,right_snp)
 	 print(c(i,length(alt_right_bound_bp)))
}


# Once you have the highest SNP of a QTL, you can calculate the bounds


all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}

bounds=c()
for(i in 1:nrow(qtl)){
	row=qtl[i,]
	chr=row$CHR
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
	snp=row$SNP
	pheno=row$phenotype
	env=row$environment
	#pos=pmap[pmap$marker==snp,]$pos
	max_pos=pmap[pmap$marker==snp,]$pos
	#results=fread(sprintf('QTL/adjusted/Biogemma_chr%s_%s_x_%s_adjusted_founderprobs.txt',chr,pheno,env),data.table=F)
	results=fread(sprintf('QTL/MITE_only/results/Biogemma_chr%s_%s_x_%s_MITE_only_founderprobs.txt',chr,pheno,env),data.table=F)

  	# Set window size to a 100MB region around the highest SNP
 	window_range=c(max_pos-5e7, max_pos + 5e7)
  	max_logp=-log10(results[results$X_ID==snp,'p_value_ML'])
  	window=results[-log10(results$p_value_ML)>=(max_logp-2),c('X_ID','p_value_ML')]
  	window=merge(window,pmap,by.x='X_ID',by.y='marker')
  	# Grab only SNPs within that window
  	window=window[window$pos>=window_range[1] & window$pos<=window_range[2],]
  	rownames(window)=seq(1,dim(window)[1])
  	left_bound_snp=window[which.min(window$pos),]$X_ID
  	left_bound_bp=min(window$pos)
  	right_bound_snp=window[which.max(window$pos),]$X_ID
  	right_bound_bp=as.numeric(max(window$pos))
  	bounds=rbind(bounds,c(pheno,env,chr,snp,max_pos,left_bound_snp,left_bound_bp,right_bound_snp,right_bound_bp))
}


bounds=as.data.frame(bounds,stringsAsFactors = F)
names(bounds)=c('phenotype','environment','CHR','SNP','BP','left_bound_snp','left_bound_bp','right_bound_snp','right_bound_bp')

bounds$alt_right_bound_bp=all_founder_blocks[match(bounds$right_bound_snp,all_founder_blocks$focal_snp),]$end


#fwrite(bounds,'QTL/all_adjusted_QTL_support_intervals.txt',row.names=F,quote=F,sep='\t')
fwrite(bounds,'QTL/MITE_only/all_MITE_only_QTL_support_intervals.txt',row.names=F,quote=F,sep='\t')



# For factor eqtl
#factoreqtl=fread('eqtl/results/all_residual_factor_fdr_peaks_FIXED.txt',data.table=F)

factoreqtl=fread('eqtl/results/all_factor_fdr_peaks_FIXED.txt',data.table=F)

bounds=c()
for(i in 1:nrow(factoreqtl)){
	row=factoreqtl[i,]
	chr=row$CHR
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
	snp=row$X_ID
	time=row$time
	factor=row$factor
	max_pos=pmap[pmap$marker==snp,]$pos
	results=fread(sprintf('eqtl/trans/results/%s_c%.0f_%s_trans_results_FIXED.txt',time,chr,factor),data.table=F)  	# Set window size to a 100MB region around the highest SNP

	#results=fread(sprintf('eqtl/trans/results/%s_residuals_c%.0f_%s_trans_results_FIXED.txt',time,chr,factor),data.table=F)  	# Set window size to a 100MB region around the highest SNP
 	window_range=c(max_pos-5e7, max_pos + 5e7)
  	max_logp=-log10(results[results$X_ID==snp,'p_value_ML'])
  	window=results[-log10(results$p_value_ML)>=(max_logp-2),c('X_ID','p_value_ML')]
  	window=merge(window,pmap,by.x='X_ID',by.y='marker')
  	# Grab only SNPs within that window
  	window=window[window$pos>=window_range[1] & window$pos<=window_range[2],]
  	rownames(window)=seq(1,dim(window)[1])
  	left_bound_snp=window[which.min(window$pos),]$X_ID
  	left_bound_bp=min(window$pos)
  	right_bound_snp=window[which.max(window$pos),]$X_ID
  	right_bound_bp=as.numeric(max(window$pos))
  	bounds=rbind(bounds,c(time,factor,chr,snp,left_bound_snp,left_bound_bp,right_bound_snp,right_bound_bp))
}


bounds=as.data.frame(bounds,stringsAsFactors = F)
names(bounds)=c('time','factor','CHR','SNP','left_bound_snp','left_bound_bp','right_bound_snp','right_bound_bp')

bounds$alt_right_bound_bp=all_founder_blocks[match(bounds$right_bound_snp,all_founder_blocks$focal_snp),]$end

fwrite(bounds,'eqtl/results/all_factor_fdr_SIs_FIXED.txt',row.names=F,quote=F,sep='\t')

#fwrite(bounds,'eqtl/results/all_residual_factor_fdr_SIs_FIXED.txt',row.names=F,quote=F,sep='\t')

########### trans-eQTL #################

#trans=fread('eqtl/results/all_trans_fdr_peaks_FIXED.txt',data.table=F)
trans=fread('eqtl/results/all_trans_fdr_hits_FIXED.txt',data.table=F)
trans$log10p=-log10(trans$p_value_ML)

all_founder_blocks=c()
for(chr in 1:10){#
  founder_blocks=fread(sprintf('eqtl/data/founder_recomb_blocks_c%s.txt',chr),data.table=F)
  all_founder_blocks=rbind(all_founder_blocks,founder_blocks)
}


bounds=c()
times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=1:10
for(chr in chroms){
	for(time in times){
		subtrans0=trans[trans$time==time & trans$CHR==chr,]
		tgenes=unique(subtrans0$Trait)
		#pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
		for(gene in tgenes){
			subtrans=subtrans0[subtrans0$Trait==gene,]
  			while(nrow(subtrans)!=0){
  				if(nrow(subtrans)==1){
  					startrow=subtrans
  					snp=startrow$X_ID
  					max_pos=startrow$BP
  					bounds=rbind(bounds,c(gene,time,chr,snp,max_pos,max_pos,max_pos,snp,snp))
  					subtrans=subtrans[subtrans$X_ID!=snp,]
  				}else{
  					startrow=subtrans[which.max(subtrans$log10p),]
					snp=startrow$X_ID
					max_pos=startrow$BP
					window_range=c(max_pos-2.5e7, max_pos + 2.5e7) # Within 50Mb of highest SNP
					max_logp=startrow$log10p
					window=subtrans[subtrans$BP>=window_range[1] & subtrans$BP<=window_range[2],]
					#window=subtrans[subtrans$log10p>=(max_logp-2),]
  					rownames(window)=seq(1,dim(window)[1])
  					left_bound_snp=window[which.min(window$BP),]$X_ID
  					left_bound_bp=min(window$BP)
  					right_bound_snp=window[which.max(window$BP),]$X_ID
  					right_bound_bp=as.numeric(max(window$BP))
  					bounds=rbind(bounds,c(gene,time,chr,snp,max_pos,left_bound_bp,right_bound_bp,left_bound_snp,right_bound_snp))
  					drop_snps=unique(c(snp,window$X_ID))
  					subtrans=subtrans[!(subtrans$X_ID %in% drop_snps),]
  				}
  			}
		}
	}
}



bounds=as.data.frame(bounds,stringsAsFactors = F)
names(bounds)=c('gene','time','CHR','SNP','BP','left_bound_bp','right_bound_bp','left_bound_snp','right_bound_snp')

bounds$alt_right_bound_bp=all_founder_blocks[match(bounds$right_bound_snp,all_founder_blocks$focal_snp),]$end
bounds$left_bound_bp=as.numeric(left_bound_bp)
fwrite(bounds,'eqtl/results/all_trans_fdr_SIs_FIXED.txt',row.names=F,quote=F,sep='\t')

bounds=fread('eqtl/results/all_trans_fdr_SIs_FIXED.txt',data.table=F)
for(chr in chroms){
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',chr),data.table=F)
	for(time in times){
		subtrans0=bounds[bounds$time==time & bounds$CHR==chr,]
		subtrans0=subtrans0[(subtrans0$SNP==subtrans0$left_bound_snp & subtrans0$SNP==subtrans0$right_bound_snp),]
		#tgenes=unique(subtrans0$gene)
		if(nrow(subtrans0)!=0){
			results=fread(sprintf('eqtl/trans/results/trans_eQTL_%s_c%.0f_weights_results_filtered_FIXED.txt',time,chr),data.table=F)
			results$log10p=-log10(results$p_value_ML)
			#results=merge(results,pmap,by.x='X_ID',by.y='marker')
			for(i in 1:nrow(subtrans0)){
				row=subtrans0[i,]
				snp=row$SNP
				gene=row$gene
				res=results[results$Trait==gene,]
				max_logp=res[res$X_ID==snp,]$log10p
				max_pos=row$BP
				window_range=c(max_pos-1e7, max_pos + 1e7) # Within 20Mb of highest SNP
				res=merge(res,pmap,by.x='X_ID',by.y='marker')
				window=res[res$pos>=window_range[1] & res$pos<=window_range[2],]
				window=window[window$log10p>=(max_logp-2),]
				rownames(window)=seq(1,dim(window)[1])
  				left_bound_snp=window[which.min(window$pos),]$X_ID
  				left_bound_bp=min(window$pos)
  				right_bound_snp=window[which.max(window$pos),]$X_ID
  				right_bound_bp=as.numeric(max(window$pos))
  			
  				bounds[bounds$gene==gene & bounds$SNP==snp,]$left_bound_bp=left_bound_bp
  				bounds[bounds$gene==gene & bounds$SNP==snp,]$left_bound_snp=left_bound_snp
  				bounds[bounds$gene==gene & bounds$SNP==snp,]$right_bound_bp=right_bound_bp
		  		bounds[bounds$gene==gene & bounds$SNP==snp,]$right_bound_snp=right_bound_snp
		  		bounds[bounds$gene==gene & bounds$SNP==snp,]$alt_right_bound_bp=all_founder_blocks[all_founder_blocks$focal_snp==snp,]$end
			}
		}	
	}
}

fwrite(bounds,'eqtl/results/all_trans_fdr_SIs_FIXED.txt',row.names=F,quote=F,sep='\t')

