#!/usr/bin/env Rscript

library('data.table')


sig=fread('eqtl/results/all_trans_fdr_hits_FIXED.txt',data.table=F)
sig$log10p=-log10(sig$p_value_ML)

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=unique(sig$CHR)
cutoff=0.5

ld=fread('eqtl/data/all_founder_ld.txt',data.table=F)


peaks=c()
for(chr in chroms){
	#ld=readRDS(sprintf('eqtl/data/founder_ld_c%.0f_c%.0f.rds',chr,chr))
	#ld=rbindlist(ld)
	#ld=as.data.frame(ld)
	sub0=sig[sig$CHR==chr,]
	for(time in times){
		sub1=sub0[sub0$time==time,]
		tgenes=unique(sub1$Trait)
		for(gene in tgenes){
			subsig=sub1[sub1$Trait==gene,]
			#subsig=subsig[order(subsig$value,reverse=T),]
			if(nrow(subsig)>1){
				snps=unique(subsig$X_ID)
				subld=ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]
				startrow=subsig[which.max(subsig$log10p),]
				startsnp=startrow$X_ID
				final=c()
				clumps=c()
				for(i in 1:nrow(subsig)){
					row=subsig[i,]
					snp2=row$X_ID
					tests=c(startsnp,snp2)
					r2=subld[subld$SNP_A==startsnp & subld$SNP_B==snp2,]$r2
					#pair=subld[(subld$SNP_A==startsnp & subld$SNP_B==snp2) | (subld$SNP_A==snp2 & subld$SNP_B==startsnp),]
					#print(r2)
					if(r2<cutoff){
						clumps=rbind(clumps,row)
					}
				}
				final=rbind(final,startrow)
				while(!is.null(clumps)){
					drop=c()
					clumprow=clumps[which.max(clumps$log10p),]
					final=rbind(final,clumprow)
					clumpsnp=clumprow$X_ID
					#tmp=clumps
					for(i in 1:nrow(clumps)){
						row=clumps[i,]
						snp2=row$X_ID
						r2=subld[subld$SNP_A==clumpsnp & subld$SNP_B==snp2,]$r2
						if(r2>=cutoff){
							drop=c(drop,i)
						}
					}
					clumps=clumps[-drop,]
					if(nrow(clumps)==0){
						clumps=c()
					}
				}
			}else{
				final=subsig
			}
			peaks=rbind(peaks,final)
		}
	}
}

fwrite(peaks,'eqtl/results/all_trans_fdr_peaks_FIXED.txt',row.names=F,quote=F,sep='\t')

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]
cutoff=0.5
thresh=0.10

hitfile='eqtl/results/all_trans_lfdr_hits_FIXED.txt'
peakfile='eqtl/results/all_trans_lfdr_peaks_FIXED.txt'
sig=fread(hitfile,data.table=F)
peaks=fread(peakfile,data.table=F)
chroms=unique(sig$CHR)
for(chr in chroms){
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	pmap=pmap[order(pmap$pos),]
	dropped=readRDS(sprintf('../genotypes/probabilities/geno_probs/dropped/bg%s_dropped_markers_genoprobs.rds',chr))
	subpeaks=peaks[peaks$CHR==chr,]
	subsig=sig[sig$CHR==chr,]
	if(nrow(subsig)>1){
		snps=unique(subsig$SNP)
		subld=ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]
		for(i in 1:nrow(subpeaks)){
			row=subpeaks[i,]
			snp2=row$SNP
			inld=subld[subld$SNP_A==snp2,]
			inld=inld[inld$r2>=0.5,]
			loc=which.min(inld$BP_B)
			left=inld[loc,]$BP_B
			left_SNP=inld[loc,]$SNP_B
			loc2=which.max(inld$BP_B)
			right=inld[loc2,]$BP_B
			right_SNP=inld[loc2,]$SNP_B
			peaks[peaks$SNP==snp2,]$leftmost=left
			peaks[peaks$SNP==snp2,]$leftmost_SNP=left_SNP
			peaks[peaks$SNP==snp2,]$rightmost=right
			peaks[peaks$SNP==snp2,]$rightmost_SNP=right_SNP
			r=which(unlist(lapply(dropped,function(x) x$marker==right_SNP)))
			linked=dropped[[r]]$linked
			last=which.max(pmap[match(linked,pmap$marker),]$pos)
			new_rightmost=pmap[pmap$marker==linked[last],]$pos
			peaks[peaks$SNP==snp2,]$alt_rightmost=new_rightmost
		}
	}			
}
fwrite(peak,peakfile,row.names=F,quote=F,sep='\t')




sig=fread('eqtl/results/all_factor_trans_eqtl_fdr_hits_FIXED.txt',data.table=F)

#sig=fread('eqtl/results/all_residuals_factor_trans_eqtl_fdr_hits_FIXED.txt',data.table=F)

ld=fread('eqtl/data/all_founder_ld.txt',data.table=F)

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=unique(sig$CHR)
cutoff=0.5

peaks=c()
for(chr in chroms){
	#ld=readRDS(sprintf('eqtl/data/founder_ld_c%.0f_c%.0f.rds',chr,chr))
	#ld=rbindlist(ld)
	#ld=as.data.frame(ld)
	times=unique(sig$time)
	sub0=sig[sig$CHR==chr,]
	for(time in times){
		sub1=sub0[sub0$time==time,]
		tgenes=unique(sub1$Trait)
		for(gene in tgenes){
			subsig=sub1[sub1$Trait==gene,]
			#subsig=subsig[order(subsig$value,reverse=T),]
			if(nrow(subsig)>1){
				snps=unique(subsig$X_ID)
				subld=ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]
				startrow=subsig[which.max(subsig$value),]
				startsnp=startrow$X_ID
				final=c()
				clumps=c()
				for(i in 1:nrow(subsig)){
					row=subsig[i,]
					snp2=row$X_ID
					tests=c(startsnp,snp2)
					r2=subld[subld$SNP_A==startsnp & subld$SNP_B==snp2,]$r2
					#pair=subld[(subld$SNP_A==startsnp & subld$SNP_B==snp2) | (subld$SNP_A==snp2 & subld$SNP_B==startsnp),]
					#print(r2)
					if(r2<cutoff){
						clumps=rbind(clumps,row)
					}
				}
				final=rbind(final,startrow)
				while(!is.null(clumps)){
					drop=c()
					clumprow=clumps[which.max(clumps$value),]
					final=rbind(final,clumprow)
					clumpsnp=clumprow$X_ID
					#tmp=clumps
					for(i in 1:nrow(clumps)){
						row=clumps[i,]
						snp2=row$X_ID
						r2=subld[subld$SNP_A==clumpsnp & subld$SNP_B==snp2,]$r2
						if(r2>=cutoff){
							drop=c(drop,i)
						}
					}
					clumps=clumps[-drop,]
					if(nrow(clumps)==0){
						clumps=c()
					}
				}
			}else{
				final=subsig
			}
			peaks=rbind(peaks,final)
		}
	}
}
# Gene counts
fwrite(peaks,'eqtl/results/all_factor_fdr_peaks_FIXED.txt',row.names=F,quote=F,sep='\t')

#fwrite(peaks,'eqtl/results/all_residual_factor_fdr_peaks_FIXED.txt',row.names=F,quote=F,sep='\t')


#hitfile='eqtl/results/all_residuals_factor_trans_eqtl_fdr_hits_FIXED.txt'
#peakfile='eqtl/results/all_residual_factor_fdr_peaks_FIXED.txt'

hitfile='eqtl/results/all_factor_trans_eqtl_fdr_hits_FIXED.txt'
peakfile='eqtl/results/all_factor_fdr_peaks_FIXED.txt'
sig=fread(hitfile,data.table=F)
peaks=fread(peakfile,data.table=F)
peaks$leftmost=0
peaks$rightmost=0
peaks$alt_rightmost=0
peaks$leftmost_SNP=""
peaks$rightmost_SNP=""
chroms=unique(sig$CHR)
for(chr in chroms){
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	pmap=pmap[order(pmap$pos),]
	dropped=readRDS(sprintf('../genotypes/probabilities/geno_probs/dropped/bg%s_dropped_markers_genoprobs.rds',chr))
	subpeaks=peaks[peaks$CHR==chr,]
	subsig=sig[sig$CHR==chr,]
	if(nrow(subsig)>1){
		snps=unique(subsig$X_ID)
		subld=ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]
		for(i in 1:nrow(subpeaks)){
			row=subpeaks[i,]
			snp2=row$X_ID
			inld=subld[subld$SNP_A==snp2,]
			inld=inld[inld$r2>=0.5,]
			loc=which.min(inld$BP_B)
			left=inld[loc,]$BP_B
			left_SNP=inld[loc,]$SNP_B	
			loc2=which.max(inld$BP_B)
			right=inld[loc2,]$BP_B
			right_SNP=inld[loc2,]$SNP_B
			peaks[peaks$X_ID==snp2,]$leftmost=left
			peaks[peaks$X_ID==snp2,]$leftmost_SNP=left_SNP
			peaks[peaks$X_ID==snp2,]$rightmost=right
			peaks[peaks$X_ID==snp2,]$rightmost_SNP=right_SNP			
			r=which(unlist(lapply(dropped,function(x) x$marker==right_SNP)))
			linked=dropped[[r]]$linked
			last=which.max(pmap[match(linked,pmap$marker),]$pos)
			new_rightmost=pmap[pmap$marker==linked[last],]$pos
			peaks[peaks$X_ID==snp2,]$alt_rightmost=new_rightmost
		}
	}else{
		snp2=subsig$X_ID
		left=subsig$BP
		left_SNP=snp2
		right=subsig$BP
		right_SNP=snp2
		peaks[peaks$X_ID==snp2,]$leftmost=left
		peaks[peaks$X_ID==snp2,]$leftmost_SNP=left_SNP
		peaks[peaks$X_ID==snp2,]$rightmost=right
		peaks[peaks$X_ID==snp2,]$rightmost_SNP=right_SNP
		r=which(unlist(lapply(dropped,function(x) x$marker==right_SNP)))
		linked=dropped[[r]]$linked
		last=which.max(pmap[match(linked,pmap$marker),]$pos)
		new_rightmost=pmap[pmap$marker==linked[last],]$pos
		peaks[peaks$X_ID==snp2,]$alt_rightmost=new_rightmost
	}
}


# Residuals
fwrite(peaks,'eqtl/results/all_factor_fdr_peaks_FIXED.txt',row.names=F,quote=F,sep='\t')

#fwrite(peaks,'eqtl/results/all_residual_factor_fdr_peaks_FIXED.txt',row.names=F,quote=F,sep='\t')

##### For all QTL updated
phenotypes=fread('phenotypes/phenotypes_all.csv',data.table=F)
ld=fread('eqtl/data/all_founder_ld.txt',data.table=F)

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]
cutoff=0.25
thresh=0.10
for(p in phenos){
	for(e in envs){
		hitfile=sprintf('QTL/MITE_only/%s_%s_QTL_scan_%.2f0000_hits.txt',p,e,thresh)
		if(file.exists(hitfile)){
			sig=fread(hitfile,data.table=F)
			chroms=unique(sig$CHR)
			peaks=c()
			for(chr in chroms){
				pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
				pmap=pmap[order(pmap$pos),]
				subsig=sig[sig$CHR==chr,]
				subsig=subsig[order(subsig$BP),]
				if(nrow(subsig)>1){
					snps=unique(subsig$SNP)
					subld=ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]
					startrow=subsig[which.max(subsig$value),]
					startsnp=startrow$SNP
					#leftmost=startrow$BP
					#rightmost=startrow$BP
					final=c()
					clumps=c()
					for(i in 1:nrow(subsig)){
						row=subsig[i,]
						snp2=row$SNP
						tests=c(startsnp,snp2)
						r2=subld[subld$SNP_A==startsnp & subld$SNP_B==snp2,]$r2
						if(r2<cutoff){
							clumps=rbind(clumps,row)
						}#else{
						#	leftmost=min(leftmost,row$BP)
						#	rightmost=max(rightmost,row$BP)
						#}
					}
					final=rbind(final,startrow)
					#final=cbind(final)#,leftmost,rightmost)
					while(!is.null(clumps)){
						drop=c()
						#leftmost=min(clumps$BP)
						#rightmost=max(clumps$BP)
						clumprow=clumps[which.max(clumps$value),]
						#clumprow=cbind(clumprow)#,leftmost,rightmost)
						final=rbind(final,clumprow)
						clumpsnp=clumprow$SNP
						#leftmost=clumprow$BP
						#rightmost=clumprow$BP
						for(i in 1:nrow(clumps)){
							row=clumps[i,]
							snp2=row$SNP
							r2=subld[subld$SNP_A==clumpsnp & subld$SNP_B==snp2,]$r2
							if(r2>=cutoff){
								drop=c(drop,i)
							}	
						}
						clumps=clumps[-drop,]
						if(nrow(clumps)==0){
							clumps=c()
						}
					}
				}else{
					final=subsig#,subsig$BP,subsig$BP)
					names(final)=c('CHR','BP','SNP','value')#,'leftmost','rightmost')
				}
				peaks=rbind(peaks,final)
			}
			fwrite(peaks,sprintf('QTL/MITE_only/%s_%s_QTL_scan_%.2f_peaks.txt',p,e,thresh),row.names=F,quote=F,sep='\t')
		}
	}
}


# Get bounds of peak
envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]
cutoff=0.25
thresh=0.10
for(p in phenos){
	for(e in envs){
		hitfile=sprintf('QTL/MITE_only/%s_%s_QTL_scan_%.2f0000_hits.txt',p,e,thresh)
		peakfile=sprintf('QTL/MITE_only/%s_%s_QTL_scan_%.2f_peaks.txt',p,e,thresh)
		if(file.exists(hitfile)){
			sig=fread(hitfile,data.table=F)
			peaks=fread(peakfile,data.table=F)
			peaks$leftmost=0
			peaks$rightmost=0
			peaks$alt_rightmost=0
			peaks$leftmost_SNP=""
			peaks$rightmost_SNP=""
			chroms=unique(sig$CHR)
			for(chr in chroms){
				pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
				pmap=pmap[order(pmap$pos),]
				dropped=readRDS(sprintf('../genotypes/probabilities/geno_probs/dropped/bg%s_dropped_markers_genoprobs.rds',chr))
				subpeaks=peaks[peaks$CHR==chr,]
				subsig=sig[sig$CHR==chr,]
				if(nrow(subsig)>1){
					snps=unique(subsig$SNP)
					subld=ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]
					for(i in 1:nrow(subpeaks)){
						row=subpeaks[i,]
						snp2=row$SNP
						inld=subld[subld$SNP_A==snp2,]
						inld=inld[inld$r2>=0.5,]
						loc=which.min(inld$BP_B)
						left=inld[loc,]$BP_B
						left_SNP=inld[loc,]$SNP_B
						
						loc2=which.max(inld$BP_B)
						right=inld[loc2,]$BP_B
						right_SNP=inld[loc2,]$SNP_B
						peaks[peaks$SNP==snp2,]$leftmost=left
						peaks[peaks$SNP==snp2,]$leftmost_SNP=left_SNP
						peaks[peaks$SNP==snp2,]$rightmost=right
						peaks[peaks$SNP==snp2,]$rightmost_SNP=right_SNP
						
						r=which(unlist(lapply(dropped,function(x) x$marker==right_SNP)))
						linked=dropped[[r]]$linked
						last=which.max(pmap[match(linked,pmap$marker),]$pos)
						new_rightmost=pmap[pmap$marker==linked[last],]$pos
						peaks[peaks$SNP==snp2,]$alt_rightmost=new_rightmost
					}
				}else{
					snp2=subsig$SNP
					left=subsig$BP
					left_SNP=snp2
					right=subsig$BP
					right_SNP=snp2
					peaks[peaks$SNP==snp2,]$leftmost=left
					peaks[peaks$SNP==snp2,]$leftmost_SNP=left_SNP
					peaks[peaks$SNP==snp2,]$rightmost=right
					peaks[peaks$SNP==snp2,]$rightmost_SNP=right_SNP
					r=which(unlist(lapply(dropped,function(x) x$marker==right_SNP)))
					linked=dropped[[r]]$linked
					last=which.max(pmap[match(linked,pmap$marker),]$pos)
					new_rightmost=pmap[pmap$marker==linked[last],]$pos
					peaks[peaks$SNP==snp2,]$alt_rightmost=new_rightmost
				}
			}
			fwrite(peaks,sprintf('QTL/MITE_only/%s_%s_QTL_scan_%.2f_peaks.txt',p,e,thresh),row.names=F,quote=F,sep='\t')
		}
	}
}

envs=unique(phenotypes$Loc.Year.Treat)
phenos=names(phenotypes)[3:9]
allq=c()
for(p in phenos){
	for(e in envs){
		peakfile=sprintf('QTL/MITE_only/%s_%s_QTL_scan_%.2f_peaks.txt',p,e,thresh)
		if(file.exists(peakfile)){
			#sig=fread(hitfile,data.table=F)
			peaks=fread(peakfile,data.table=F)
			peaks$phenotype=p
			peaks$environment=e
			allq=rbind(allq,peaks)
		}
	}
}
fwrite(allq,'QTL/MITE_only/all_MITE_only_QTL_peaks.txt',row.names=F,quote=F,sep='\t')


m1=intersect(pmap$marker,subld$SNP_A)
subld$SNP_A_f=factor(subld$SNP_A,levels=c(m1))
subld$SNP_B_f=factor(subld$SNP_B,levels=c(m1))

subpeaks$SNP_A_f=factor(subpeaks$SNP,levels=c(m1))
subpeaks$SNP_B_f=factor(subpeaks$SNP,levels=c(m1))

h1=ggplot(subld,aes(x=SNP_A_f,y=SNP_B_f,fill=r2)) + geom_tile() +
	theme_classic() + scale_fill_gradient(low="white", high="blue") +
	theme(axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank()) +
    geom_hline(yintercept=subpeaks$SNP_B_f,color="red") +
    geom_vline(xintercept=subpeaks$SNP_A_f,color="red") +
    xlab(sprintf("Chromosome %s",chr)) + ylab(sprintf("Chromosome %s",chr)) + 
    ggtitle("Founder LD around qDTA8")
    
pdf('QTL/qDTA8_ALL_founder_LD.pdf')
print(h1)
dev.off()
