#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('tidyverse')
time="WD_0712"

lambda=fread(sprintf('MegaLMM/MegaLMM_%s_all_Lambda_means_FIXED.txt',time),data.table=F)
genetable=fread('eqtl/data/Zea_mays.B73_RefGen_v4.46_gene_list.txt',data.table=F)

#sub=lambda[,1:6]
factors=names(lambda)[-1]
sub=merge(x=lambda,y=genetable,by.x="V1",by.y="Gene_ID")


sub.tmp <- sub %>% group_by(CHROM) %>%
summarise(chr_len=max(END)) %>%
mutate(tot=cumsum(chr_len)-chr_len) %>%
#select(-chr_len) %>%
left_join(sub, ., by=c("CHROM"="CHROM")) %>%
arrange(CHROM, END) %>%
mutate( BPcum=as.numeric(END+tot)) %>%
gather(key, value,-START,-END,-V1,-CHROM,-BPcum,-tot,-all_of(factors))




# Only plot values that are greater than 0.2 (explain more than 4% of variance)
plot_list=vector("list",length(factors))
count=1
for(f in factors){
  p1=ggplot(aes_string(y=f),data=sub.tmp[which(abs(sub.tmp[,f])>0.2),])+ geom_point(aes(x=BPcum)) + ggtitle(sprintf("%s",f))
  plot_list[[count]]=p1
  count=count+1
}


pdf('images/lambda_by_pos.pdf')
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()
