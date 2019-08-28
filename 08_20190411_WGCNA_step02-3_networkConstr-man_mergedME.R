library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source("./config.R")
test_list = test_list

Sft = c(5,5,6,6,7,7,5,5)

inf_Expr = paste(       results_folder, 
                        "/05_WGCNA_step01_datainput/",
                        "datExpr_",
                         test_list,
                         ".csv",
                         sep="")

dir.create(     paste(results_folder,
                "/08_WGCNA_step02-3_mergedME",
                sep=""))

outf_ME = paste(
		results_folder,
                "/08_WGCNA_step02-3_mergedME/",
		test_list,
                sep="")

MEDissThres_list = rep(0.3, 8)


merged_ME <- function(sft, MEDissThres, inf_expr, outf_me)
{
	softPower = sft;

	datExpr   = read.csv(file=inf_expr, row.names=1)	

	adjacency = adjacency(datExpr, power = softPower);
	TOM = TOMsimilarity(adjacency);

	dissTOM = 1-TOM

	geneTree = hclust(as.dist(dissTOM), method = "average");


	minModuleSize = 30;	

	dynamicMods = cutreeDynamic(	dendro = geneTree, distM = dissTOM,

					deepSplit = 2, pamRespectsDendro = FALSE,

					minClusterSize = minModuleSize);

	dynamicColors = labels2colors(dynamicMods)
	
	MEList = moduleEigengenes(datExpr, colors = dynamicColors)
	MEs = MEList$eigengenes
	MEDiss = 1-cor(MEs);
	METree = hclust(as.dist(MEDiss), method = "average");

	merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	mergedColors = merge$colors;
	mergedMEs = merge$newMEs;

	png(filename = paste(outf_me, ".png", sep=""))
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
	dev.off()

	moduleColors = mergedColors
	colorOrder = c("grey", standardColors(50));
	moduleLabels = match(moduleColors, colorOrder)-1;
	MEs = mergedMEs;

	save(MEs, moduleLabels, moduleColors, geneTree, file = paste(outf_me, ".RData", sep=""))

}



for(i in 1:length(test_list))
{
	merged_ME(  sft = Sft[i],
		    MEDissThres = MEDissThres_list[i],
                    inf_expr = inf_Expr[i],
                    outf_me = outf_ME[i])
}
