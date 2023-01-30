#!/usr/bin/env Rscript

library('data.table')
library('lme4')
library('lme4qtl')
library('emmeans')
library('multcomp')
library('ggplot2')

time="WD_0718"
y="Zm00001d042196"

colorcodes=fread('founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

# Read in expression data
data=fread('test_Zm00001d042196_data.txt',data.table=F)
fkeep=names(data)[3:ncol(data)]
chr=3
snp="AX-90833615"
pos=155200685
  # Read in cis-eQTL effect sizes

  # Read in kinship matrix
K = fread(sprintf('K_matrix_chr%s.txt',chr),data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])

  # lmer model
f=as.formula(paste('y',paste(c(0,fkeep,"(1|ID)"),collapse=" + "),sep=" ~ "))
print(f)
m3=lmer(f,data=data,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",
     check.nobs.vs.nRE="ignore"))

# combine emmeans objects for each founder
for(f in fkeep){
  atlist=list()
  for(f2 in fkeep){
    if(f==f2){
      atlist[[f2]]=1
    }else{
      atlist[[f2]]=0
    }
  }
  refgrid=ref_grid(m3,at=atlist)
  em=emmeans(refgrid,specs=f,lmer.df="kenward-roger")

  if(f==fkeep[1]){
    allem=em
  }else{
    allem=rbind(allem,em)
  }
}
# run tukey cld
cld=cld(allem,Letters=letters)
cld=as.data.frame(cld,stringsAsFactors=F)

#sort by mean
cld=cld[order(cld$emmean),]
rownames(cld)=1:nrow(cld)
variable_f=sapply(seq(1,length(fkeep)),function(x) fkeep[which(cld[x,fkeep]==1)])
cld$variable_f=factor(variable_f,levels=variable_f)

fgroups=cld$.group



# plot founder effect sizes
p1=ggplot(cld,aes(x=variable_f,y=emmean,color=variable_f)) + geom_point() +
geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),data=cld) +
geom_text(aes(label=.group),vjust=-5,color="black",size=10)+
scale_color_manual(values=colorcodes[levels(cld$variable_f),]$hex_color,labels=levels(cld$variable_f))+
theme(axis.text.x=element_text(size=10)) +
xlab("Founder") + ylab("Expression (log2CPM)") +
labs(title=sprintf("%s cis-eQTL Founder Effect Sizes (lmer)",y))

png('lmer_test.png')
print(p1)
dev.off()
