#!/usr/bin/env Rscript

library('WGCNA')
library('data.table')

time="WD_0712"

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
datExpr=fread(sprintf('eqtl/normalized/%s_voom_normalized_gene_counts_formatted.txt',time),data.table=F)
rownames(datExpr)=datExpr$ID
datExpr=datExpr[,-1]

datTraits=fread(sprintf('eqtl/results/%s_vst_counts_x_phenotypes.txt',time),data.table=F)
phenotypes=names(datTraits)[grepl('_',names(datTraits))]
datTraits=datTraits[,c('ID',phenotypes)]
rownames(datTraits)=datTraits$ID
datTraits=datTraits[,-1]
#QC
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

png(sprintf('images/%s_wgcna_outliers.png',time))
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 15, col = "red");
dev.off()
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(sprintf('images/%s_wgcna_dendrogram.pdf',time))
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
dev.off()

save(datExpr, datTraits, file = sprintf("eqtl/results/WGCNA_%s_Input.RData",time))

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
png(sprintf('images/%s_wgcna_power.pdf',time))
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

lnames = load(file = "eqtl/results/WGCNA_WD_0712_Input.RData");
#The variable lnames contains the names of loaded variables.
lnames


net = blockwiseModules(datExpr, power = 5,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "WD_0712-1step",
verbose = 3)

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
pdf('images/wgnca_1stepnet_dendrogram.pdf')
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()


bwnet = blockwiseModules(datExpr, maxBlockSize = 20000,
power = 4, TOMType = "unsigned",
#minModuleSize = 30,
#reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE,
saveTOMs = TRUE,
saveTOMFileBase = "WD_0712-blockwise",
verbose = 3)

# Relabel blockwise modules
#bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
#bwModuleColors = labels2colors(bwLabels)


# Plot the dendrogram and the module colors underneath for block 1
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)
pdf('images/wgnca_bwnet_dendrogram.pdf')
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "BW_0712_bw_networkConstruction-auto.RData")

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

flowerpheno=phenotypes[grep('flowering',phenotypes)]
flowering=moduleTraitCor[,flowerpheno]
floweringpval=moduleTraitPvalue[,flowerpheno]
pdf('WGCNA/images/WD_0712_module_flowering_heatmap.pdf')

sizeGrWindow(20,12)
# Will display correlations and their p-values
textMatrix = paste(signif(flowering, 2), "\n(",
signif(floweringpval, 1), ")", sep = "");
dim(textMatrix) = dim(flowering)
par(mar = c(10, 8.5, 4, 4));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = flowering,
xLabels = c(flowerpheno),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Module-flowering trait relationships"))
dev.off()

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$"BLOIS_2014_OPT-female_flowering_d6");
names(weight) = "BLOIS_2014_female_flowering"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#Identifying high GS and MM
module = "thistle1"
column = match(module, modNames);
moduleGenes = moduleColors==module;
pdf('WGCNA/images/BLOIS_2014_thistle1.pdf')
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for BLOIS_2014_female_flowering_d6",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()


notflowerpheno=phenotypes[!phenotypes %in% flowerpheno]
notflowering=moduleTraitCor[,notflowerpheno]
notfloweringpval=moduleTraitPvalue[,notflowerpheno]
pdf('WGCNA/images/WD_0712_module_notflowering_heatmap.pdf')

sizeGrWindow(40,24)
# Will display correlations and their p-values
textMatrix = paste(signif(notflowering, 2), "\n(",
signif(notfloweringpval, 1), ")", sep = "");
dim(textMatrix) = dim(notflowering)
par(mar = c(10, 8.5, 4, 4));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = notflowering,
xLabels = c(notflowerpheno),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$"BLOIS_2014_OPT-grain_yield_15");
names(weight) = "BLOIS_2014_grain_yield"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#Identifying high GS and MM
module = "mediumpurple3"
column = match(module, modNames);
moduleGenes = moduleColors==module;
pdf('WGCNA/images/BLOIS_2014_gy_mediumpurple3.pdf')
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for BLOIS_2014_grain_yield_15",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# Read in annotation data and combine information
#geneInfo0 = data.frame(substanceBXH = probes,
#geneSymbol = annot$gene_symbol[probes2annot],
#LocusLinkID = annot$LocusLinkID[probes2annot],
#moduleColor = moduleColors,
#geneTraitSignificance,
#GSPvalue)
# Order modules by their significance for weight
#modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
#for (mod in 1:ncol(geneModuleMembership))
#{
#oldNames = names(geneInfo0)
#geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
#MMPvalue[, modOrder[mod]]);
#names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
#paste("p.MM.", modNames[modOrder[mod]], sep=""))
#}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
#geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
#geneInfo = geneInfo0[geneOrder, ]

#write.csv(geneInfo, file = "geneInfo.csv")
