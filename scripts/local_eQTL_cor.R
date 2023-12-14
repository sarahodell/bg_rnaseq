#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

local=fread('QTT/QTL_cis_eQTL_interval_overlap.txt',data.table=F)

candidates=fread('QTT/sig_candidate_genes.txt',data.table=F)


founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")


cand_df=c()
for(i in 1:nrow(candidates)){
	row1=candidates[i,]
	pei=row1$pei
	chr=row1$CHR
	time1=row1$time
	sublocal=local[local$pheno_env_ID==pei,]
	sublocal=sublocal[sublocal$time==time1,]
	sublocal$gene_snp=paste0(sublocal$Trait,'_',sublocal$X_ID)
	effect_sizes=fread(sprintf('eqtl/cis/results/eQTL_%s_c%s_weights_results_FIXED.txt',time1,chr),data.table=F)
	effect_sizes$gene_snp=paste0(effect_sizes$Trait,'_',effect_sizes$X_ID)
	effect_sizes=effect_sizes[effect_sizes$gene_snp %in% unique(sublocal$gene_snp),]
	n=nrow(effect_sizes)
	for(j in 1:(n-1)){
		rowj=effect_sizes[j,]
		genej=rowj$Trait
		snpj=rowj$X_ID
		betaj=unlist(rowj[,c(6,10:24)])
		wn=which(!is.na(betaj))[1]
		betaj[-wn]=betaj[-wn]+betaj[wn]
		for(k in (j+1):n){
			rowk=effect_sizes[k,]
			genek=rowk$Trait
			snpk=rowk$X_ID
			betak=unlist(rowk[,c(6,10:24)])
			wn2=which(!is.na(betak))[1]
			betak[-wn2]=betak[-wn2]+betak[wn2]
			
			ecor=cor(betaj,betak,use="complete.obs")
			
			line1=data.frame(pei=pei,time=time1,genej=genej,snpj=snpj,genek=genek,snpk=snpk,r=ecor)
			cand_df=rbind(cand_df,line1)
		}
	}
}


fwrite(cand_df,'eqtl/results/cand_local_eQTL_cors.txt',row.names=F,quote=F,sep='\t')


for(i in 1:nrow(candidates)){
	row1=candidates[i,]
	pei=row1$pei
	chr=row1$CHR
	subcand=cand_df[cand_df$pei==pei,]
	pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
	subcand$posj=pmap[match(subcand$snpj,pmap$marker),]$pos
	subcand$posk=pmap[match(subcand$snpk,pmap$marker),]$pos
	
 	genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)
	
	submap=genetable[genetable$Gene_ID %in% unique(c(subcand$genek,subcand$genej)),]
	submap=submap[order(submap$START),]
	rownames(submap)=seq(1,nrow(submap))
	
	subcand$genej_f=factor(subcand$genej,levels=c(submap$Gene_ID))
	subcand$genek_f=factor(subcand$genek,levels=c(submap$Gene_ID))
	
	heatp=ggplot(data=subcand,aes(x=genej,y=genek,fill=r)) + geom_tile() + xlab("local-eQTL") +
	ylab("local-eQTL") + scale_fill_gradientn(colors=c("red","white","blue"),limits=c(-1, 1), breaks=seq(-1,1,by=0.25)) +
	theme(axis.text.x=element_text(angle=-45,size=4),axis.text.y=element_text(size=4))
	
	pdf(sprintf('QTT/%s_local_eQLT_r_heatmap.pdf',pei))
	print(heatp)
	dev.off()	
}

