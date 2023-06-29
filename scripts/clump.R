#!/usr/bin/env Rscript

library('data.table')


sig=fread('eqtl/results/all_trans_fdr_hits.txt',data.table=F)

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

fwrite(peaks,'eqtl/results/all_trans_fdr_peaks.txt',row.names=F,quote=F,sep='\t')

#Zm00001d048515 WD_0712


sig=fread('eqtl/results/all_factor_trans_eqtl_fdr_hits',data.table=F)

times=c("WD_0712","WD_0718","WD_0720","WD_0727")
chroms=unique(sig$CHR)
cutoff=0.5

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

fwrite(peaks,'eqtl/results/all_factor_fdr_peaks.txt',row.names=F,quote=F,sep='\t')


sig=fread('QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.100000_hits.txt',data.table=F)
chroms=unique(sig$CHR)
cutoff=0.5

ld=fread('eqtl/data/all_founder_ld.txt',data.table=F)

peaks=c()
for(chr in chroms){
	#ld=readRDS(sprintf('eqtl/data/founder_ld_c%.0f_c%.0f.rds',chr,chr))
	#ld=rbindlist(ld)
	#ld=as.data.frame(ld)
	#subld=ld[ld$CHR_A==chr & ld$CHR_B==chr,]
	subsig=sig[sig$CHR==chr,]
	if(nrow(subsig)>1){
		snps=unique(subsig$SNP)
		subld=ld[ld$SNP_A %in% snps & ld$SNP_B %in% snps,]
		startrow=subsig[which.max(subsig$value),]
		startsnp=startrow$SNP
		final=c()
		clumps=c()
		for(i in 1:nrow(subsig)){
			row=subsig[i,]
			snp2=row$SNP
			tests=c(startsnp,snp2)
			r2=subld[subld$SNP_A==startsnp & subld$SNP_B==snp2,]$r2
			if(r2<cutoff){
				clumps=rbind(clumps,row)
			}
		}
		final=rbind(final,startrow)
		while(!is.null(clumps)){
			drop=c()
			clumprow=clumps[which.max(clumps$value),]
			final=rbind(final,clumprow)
			clumpsnp=clumprow$SNP
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
		final=subsig
	}
	peaks=rbind(peaks,final)
}

fwrite(peaks,'QTL/male_flowering_d6_EXP_STPAUL_2017_WD_QTL_scan_0.10_peaks.txt',row.names=F,quote=F,sep='\t')



