#!/usr/bin/env Rscript

library('data.table')


times=c("WD_0712","WD_0718","WD_0720","WD_0727")

#wd12=fread('eqtl/normalized/WD_0712_voom_normalized_gene_counts_formatted.txt',data.table=F)
wd18=fread('eqtl/normalized/WD_0718_voom_normalized_gene_counts_formatted.txt',data.table=F)
wd20=fread('eqtl/normalized/WD_0720_voom_normalized_gene_counts_formatted.txt',data.table=F)
wd27=fread('eqtl/normalized/WD_0727_voom_normalized_gene_counts_formatted.txt',data.table=F)

#ind_inter=intersect(wd12$V1,wd18$V1)
#ind_inter=intersect(ind_inter,wd20$V1)
#ind_inter=intersect(ind_inter,wd27$V1)
#Only 24 if I do all 4 timepoints


ind_inter=intersect(wd18$V1,wd20$V1)
ind_inter=intersect(ind_inter,wd27$V1)

gene_inter=intersect(names(wd18)[-1],names(wd20)[-1])
gene_inter=intersect(gene_inter,names(wd27)[-1])
# 17233 genes
# 108 individuals for 18,20,and 27
rownames(wd18)=wd18$V1
wd18=wd18[,-1]
wd18=wd18[ind_inter,gene_inter]
rownames(wd20)=wd20$V1
wd20=wd20[,-1]
wd20=wd20[ind_inter,gene_inter]
rownames(wd27)=wd27$V1
wd27=wd27[,-1]
wd27=wd27[ind_inter,gene_inter]

print("Beginning to make dataset")

data=c()
for(i in ind_inter){
	print(i)
    line=data.frame(ind=i,gene=gene_inter,WD_0718=unlist(wd18[i,]),WD_0720=unlist(wd20[i,]),WD_0727=unlist(wd27[i,]),stringsAsFactors=F)
	data=rbind(data,line)
}
data=as.data.frame(data)
names(data)=c('ind','gene','WD_0718','WD_0720','WD_0727')
print("Calculating slopes")
xs=c(18,20,27)
slopes <- apply(data, 1, function(row) {
  lm(y ~ x, list(x = xs, y = as.numeric(row[3:5])))$coefficients['x']
})
data$slopes=slopes
print("Writing out dataset")
fwrite(data,'eqtl/expression_slopes.txt',row.names=F,quote=F,sep='\t')