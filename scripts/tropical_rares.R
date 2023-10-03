#!/usr/bin/env Rscript

library('data.table')


all_rares=c()

for(chr in 1:10){
	bgrare=fread(sprintf('datasets/overlap/biogemma_%s_tropical_rare_alleles.txt',chr),data.table=F)
	hprare=fread(sprintf('datasets/hapmap/hapmap321_%s_tropical_rare_alleles.txt',chr),data.table=F)
	overlap=bgrare[bgrare$V2 %in% hprare$V2,]
	all_rares=rbind(all_rares,overlap)
}
all_rares=as.data.frame(all_rares,stringsAsFactors=F)
fwrite(all_rares,'datasets/overlap/biogemma_trop_rare_hapmap_common.txt',row.names=F,quote=F,sep='\t',col.names=F)

all_rares %>% group_by(V1) %>% count()
# A tibble: 10 × 2
# Groups:   V1 [10]
#      V1     n
#   <int> <int>
# 1     1    58
# 2     2    49
# 3     3    54
# 4     4    35
# 5     5    32
# 6     6    38
# 7     7    38
# 8     8    37
# 9     9    24
#10    10    36


# tropical dh allele probs
bimbam=fread(sprintf('datasets/chr%s_biogemma_rare_allele_probs.txt',chr),data.table=F)
bnames=names(bimbam)[1:347]

fullbim=c()

for(chr in 1:10){
	bim=fread(sprintf('datasets/overlap/chr%s_tropical_alleleprobs.txt',chr),data.table=F)
	fullbim=rbind(fullbim,bim)
}
fullbim=as.data.frame(fullbim,stringsAsFactors=F)
names(fullbim)=bnames

all_rares=fread('datasets/overlap/biogemma_trop_rare_hapmap_common.txt',data.table=F,header=F)
trop_m=paste0('S',all_rares$V1,'_',all_rares$V2)

trop_m=as.data.frame(trop_m)
fwrite(trop_m,'datasets/overlap/tropical_rares_markers.txt',row.names=F,quote=F,col.names=F,sep='\t')
founder=fread('datasets/overlap/biogemma_founder_tropicalrare_hmpcommon_genos.txt',data.table=F)
names(founder)=c("SNP","CHR","POS","REF","ALT","A632_usa","A654_inra","B73_inra","B96","C103_inra","CO255_inra","D105_inra","DK63","EP1_inra","F492", "FV252_inra", "FV2_inra","MBS847","ND245","OH43_inra","VA85","W117_inra")

founders=c("MBS847","B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
founder=founder[,c("SNP","CHR","POS","REF","ALT",founders)]

founder[,founders]=apply(founder[,founders],MARGIN=2,function(x) ifelse(x=="0/0",0,ifelse(x=="1/1",1,ifelse(x=="2/2",2,NA))))
# remove multiallelic sites?
#founder$ALT2=NA
#mult=which(grepl(',',founder$ALT))
#alts=founder$ALT[mult]
#alt1=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[1]])
#alt2=sapply(seq(1,length(alts)),function(x) strsplit(alts[x],',')[[1]][[2]])
#founder[mult,]$ALT=alt1
#founder[mult,]$ALT2=alt2

founder=founder[,c("SNP","CHR","POS","REF","ALT","ALT2",founders)]

founders2=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra", "A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
founder$f_freq=rowSums(founder[,founders2],na.rm=T)/16


# what is the frequency of these sites in the founders?
min_founder=readRDS(sprintf('datasets/hapmap/minor_allele_founders_c%s.rds',chr))
bimbam=fread(sprintf('../genotypes/probabilities/allele_probs/bg%s_wgs_alleleprobs.txt.gz',chr),data.table=F)
bimbam=bimbam[bimbam$marker %in% names(min_founder),]
bimbam$pos=sapply(seq(1,nrow(bimbam)),function(x) as.numeric(strsplit(bimbam$marker[x],'_')[[1]][[2]]))

# in the DH lines?