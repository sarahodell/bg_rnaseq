library(data.table)
library(psych)
library(readxl)
library(parallel)
library(foreach)
library(doParallel)

env = system('env',intern=T)
#env = env[grep('SLURM_CPUS_PER_TASK',env)]
env=env[grep('SLURM_CPUS_ON_NODE',env)]
#ncores = as.numeric(sub('SLURM_CPUS_PER_TASK=','',env))
ncores = as.numeric(sub('SLURM_CPUS_ON_NODE=','',env))

# load fastq and bam file info
file_info = read.csv('updated_file_genotype_info_2.csv')
field = as.data.frame(readxl::read_excel('/home/sodell/projects/biogemma/expression/metadata/Limagrain_transcriptome_list_GS17.xlsx'))


# run samtools mpileup on each file
# skips files that already exist
samtools = '/share/apps/22.04/spack/spack-v0.19.1/opt/spack/linux-ubuntu22.04-x86_64_v2/gcc-9.5.0/samtools-1.14-v3mtfijiokswub5qb56xmzigfzwjcbdv/bin/samtools'
genome='/group/jrigrp/Share/assemblies/Zea_mays.B73_RefGen_v4.dna.toplevel.fa'
base=sprintf('%s mpileup -f %s --positions private_vcf.vcf --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends',samtools,genome)

out_dir = '/group/runciegrp2/Projects/Runcie/biogemma/samtools_output'
try(dir.create(out_dir,showWarnings=F))
# for(i in 1:nrow(bam_meta_expt)) {
parallel::mclapply(1:nrow(file_info),function(i) { #1:nrow(file_info)
	print(i)
	if(is.na(file_info$bam_file[i])) return(c())
	outfile = sprintf('%s/%s.mpileup',out_dir,file_info$sample_name[i])
	if(file.exists(outfile)) {
		if(file.info(outfile)$size > 0) return(c())
	}
	system(sprintf('%s -o %s %s',base,outfile,file_info$bam_file[i]))
},mc.cores=ncores)

# loads genotype calls for each line

geno_mats = list()
axio = fread('/home/sodell/projects/biogemma/genotypes/600K/Axiom600K_AGPv4.bed',data.table=F)

geno_info = fread('/home/sodell/projects/biogemma/expression/metadata/Lines_Names_BALANCE.txt',data.table=F)
geno_from_eh = tapply(geno_info$CODE_HYBRIDE,geno_info$ID,function(x) x[1])
names(geno_from_eh) = gsub('-','.',names(geno_from_eh),fixed=T)
gt_priv = as.matrix(read.csv('gt_priv.csv',row.names=1))

for(chr in 1:10) {
	genos = readRDS(sprintf('/home/sodell/projects/biogemma/genotypes/probabilities/geno_probs/bg%d_filtered_genotype_probs.rds',chr))		
	for(i in 1:length(genos)) {
		colnames(genos[[i]]) = axio$V2[match(colnames(genos[[i]]),axio$V4)]
	}
	geno_IDs = lapply(1:nrow(genos[[1]]),function(j) sapply(genos[colnames(gt_priv)],function(x) x[j,]))
	names(geno_IDs) = geno_from_eh[rownames(genos[[1]])]
	geno_mats[[chr]] = geno_IDs
}


# calculates correlation between RNAseq calls and each matrix of calls from the genotypes

thresh=10
registerDoParallel(ncores)
cor_results_allchr = foreach(i=1:nrow(file_info),.combine = c,.errorhandling='remove') %dopar% { 
	mp_file = sprintf('%s/%s.mpileup',out_dir,file_info$sample_name[i])
	experiment = file_info$experiment[i]
	genotype = file_info$genotype[i]
	X_Y = file_info$X_Y[i]
	if(!file.exists(mp_file)) return(c())
	mp = fread(mp_file,data.table=F,sep='\t')
	mp$dp0 = sapply(mp$V5,function(x) sum(strsplit(x,'')[[1]] %in% c('.',',')))
	mp$dp1 = sapply(mp$V5,function(x) sum(strsplit(x,'')[[1]] %in% c('a','c','g','t','A','C','G','T')))
	count = mp
	gt_1 = c()
	gt_geno = list()
	for(chr in 1:10) {
		count_chr = subset(count,V1==chr)
		count_chr = subset(count_chr,dp1 >= thresh & dp1/(dp1+dp0)>.1)
		count_chr$ID = paste(count_chr[,1],count_chr[,2],sep='_')
		count_chr = subset(count_chr,ID %in% rownames(gt_priv))
		pos = count_chr$V2
		if(length(pos) < 10) next
		gt_1 = rbind(gt_1,gt_priv[count_chr$ID,])

		geno_IDs = geno_mats[[chr]]
		f = approxfun(as.numeric(rownames(geno_IDs[[1]])),1:nrow(geno_IDs[[1]]),f=0.5,rule=2)
		match_pos = f(pos)
		for(geno in names(geno_IDs)) {
			gt_geno[[geno]] = rbind(gt_geno[[geno]],geno_IDs[[geno]][match_pos,])
		}
	}
	cors = sapply(gt_geno,function(x) cor(c(gt_1),c(x)))	
	cors = sort(cors,decreasing=T)
	list(
		data.frame(
			experiment=experiment,
			genotype=genotype,
			X_Y = X_Y,
			n = nrow(gt_1),
			geno_ID = names(cors),
			cor = cors))
}
saveRDS(cor_results_allchr,file = 'cor_results_all_chr_vs_genotypes.rds')
cor_results_allchr = readRDS('cor_results_all_chr_vs_genotypes.rds')


# pulls out the best match from each
matches_all_chr = do.call(rbind,lapply(cor_results_allchr,function(x) {
	rank_match = match(x$genotype[1],x$geno_ID)
	r = fisherz(x$cor)
	data.frame(x[1,],rank_match = rank_match, cor_match = x$cor[rank_match],
		Z_best=(r[1]-mean(r))/sd(r),Z2=(-diff(r)[1])/sd(r))
	}))

# classifier to decide which matches to genotypes are real
m = lm(Z_best~poly(log(n),4),subset(matches_all_chr,genotype == geno_ID))
resid = (predict(m,newdata = matches_all_chr) - matches_all_chr$Z_best)/sigma(m)
matches_all_chr$good_match = (resid < 4 & matches_all_chr$cor > 0.2) | matches_all_chr$genotype == matches_all_chr$geno_ID

# fixes genotype calls
# note: some genotypes end up multiple times per time point!

matches_all_chr_info = c()
for(expt in unique(matches_all_chr$experiment)) {
	matches_all_chr_expt = subset(matches_all_chr,experiment == expt)
	match = subset(matches_all_chr_expt,genotype == geno_ID)
	if(nrow(match)>0) {
		match$final_genotype = match$genotype
		match$final_X_Y = match$X_Y
		match$match_type = 'good'
	}

	no_match = subset(matches_all_chr_expt,genotype != geno_ID & !good_match)
	if(nrow(no_match) > 0) {
		no_match$final_genotype = NA
		no_match$final_X_Y = NA
		no_match$match_type = 'no_match'
	}

	bad_match = subset(matches_all_chr_expt,genotype != geno_ID & good_match)
	if(nrow(bad_match)>0) {
		bad_match$final_genotype = NA
		bad_match$final_X_Y = NA
		bad_match$match_type = NA

		for(i in 1:nrow(bad_match)) {
			geno_ID. = bad_match$geno_ID[i]
			if(geno_ID. %in% bad_match$genotype) {
				if(subset(bad_match,genotype==geno_ID.)$geno_ID == bad_match$genotype[i]) {
					bad_match$final_genotype[i] = geno_ID.
					bad_match$final_X_Y[i] = subset(bad_match,genotype==geno_ID.)$X_Y
					bad_match$match_type[i] = 'switch'
				} else {
					# print(c('fdas',i))
					bad_match$final_genotype[i] = geno_ID.					
					bad_match$match_type[i] = 'other'  # possible multi-way switch
				}	
			} else if(geno_ID. %in% match$genotype) {				
				bad_match$final_genotype[i] = geno_ID.					
				bad_match$match_type[i] = 'duplicate'	# duplicate sample because "correct" sample already exists
			} else if(geno_ID. %in% no_match$genotype) {	
				# print(c('bbbb',i))				
				bad_match$final_genotype[i] = geno_ID.					
				bad_match$match_type[i] = 'possible_duplicate'	# unclear if multi-way switch or duplicate because can't call other entry
			} else if(geno_ID. %in% subset(file_info,experiment == expt)$genotype) {
				# print(c('aaaa',i))	
				bad_match$final_genotype[i] = geno_ID.					
				bad_match$match_type[i] = 'other'    # maybe in other fq file?
			} else {
				# print(c('asdf',i))
				bad_match$final_genotype[i] = geno_ID.					
				bad_match$match_type[i] = 'other'
			}	
			if(!geno_ID. %in% field$GENOTYPE) {
				bad_match$final_genotype[i] = NA
				bad_match$match_type[i] = 'match not in field'
			}
		}
	}
	matches_all_chr_info = rbind(matches_all_chr_info,
		match,
		no_match,
		bad_match
		)
}


print(table(matches_all_chr_info$experiment,matches_all_chr_info$match_type,useNA='a'))

i = is.na(matches_all_chr_info$final_X_Y) & !is.na(matches_all_chr_info$final_genotype)
matches_all_chr_info$final_X_Y[i] = field$CLE[match(paste(ifelse(grepl('WD',matches_all_chr_info$experiment[i]),'WD','OPT'),matches_all_chr_info$final_genotype[i]),paste(field$TREATMENT,field$GENOTYPE))]


file_info$final_genotype = NA
file_info$match_type = NA
file_info$final_X_Y = NA
for(i in 1:nrow(matches_all_chr_info)) {
	expt = matches_all_chr_info$experiment[i]
	X_Y = matches_all_chr_info$X_Y[i]
	old_genotype = matches_all_chr_info$genotype[i]
	if(X_Y == 'NA_NA') next
	final_genotype = matches_all_chr_info$final_genotype[i]
	match_type = matches_all_chr_info$match_type[i]
	final_X_Y = matches_all_chr_info$final_X_Y[i]
	j = which(file_info$experiment == expt & file_info$X_Y == X_Y)
	# j = which(file_info$experiment == expt & file_info$genotype == old_genotype)
	if(length(j) == 0) {
		# stop('no hits')
		j = which(file_info$experiment == expt & file_info$genotype == old_genotype)
	}
	if(length(j) > 1) stop('multiple hits')
	if(file_info$genotype[j] != old_genotype) stop('wrong_genotype')
	file_info$final_genotype[j] = final_genotype
	file_info$match_type[j] = match_type
	file_info$final_X_Y[j] = final_X_Y
}

write.csv(file_info,file = 'updated_file_genotype_info_FIXED_genotypes_2.csv')
file_info = read.csv('updated_file_genotype_info_FIXED_genotypes_2.csv')

file_info$X = as.numeric(sapply(file_info$X_Y,function(x) strsplit(x,'_')[[1]][1]))
file_info$Y = as.numeric(sapply(file_info$X_Y,function(x) strsplit(x,'_')[[1]][2]))
file_info$final_X = as.numeric(sapply(file_info$final_X_Y,function(x) strsplit(x,'_')[[1]][1]))
file_info$final_Y = as.numeric(sapply(file_info$final_X_Y,function(x) strsplit(x,'_')[[1]][2]))


pdf('field_map.pdf')
expts = unique(file_info$experiment)
expts = sort(expts[grep('_07',expts)])
for(expt in expts) {
	try({
	a=subset(file_info,experiment==expt)
	print(expt)
	# a[order(a$plate_orig,a$loading_order_orig),]

	plot(a$X,a$Y,xlim = range(a$X,na.rm=T),ylim = range(a$Y,na.rm=T),col = ifelse(a$genotype==a$final_genotype,1,2),main=expt)
	with(subset(a,genotype != final_genotype),arrows(X,Y,final_X,final_Y,code=2))
	})
}
dev.off()

file_info = subset(file_info,experiment %in% expts)
table(file_info$experiment,file_info$genotype == file_info$final_genotype,useNA = 'a')
