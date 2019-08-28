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
                "/07_WGCNA_step02-2_clusteringME",
                sep=""))

outf_ME = paste(       results_folder,
	                "/07_WGCNA_step02-2_clusteringME/",
                         test_list,
			".png",
                         sep="")

clustering_ME <- function(sft, inf_expr, outf_me)
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

	png(filename = outf_me)
	plot(METree, main = "Clustering of module eigengenes",
 xlab = "", sub = "")
	abline(h=0.25, col = "red")
	abline(h=0.3, col = "blue")
	abline(h=0.4, col = "green")
	dev.off()
}



for(i in 1:length(test_list))
{
        clustering_ME(	sft = Sft[i],
			inf_expr = inf_Expr[i],
                        outf_me = outf_ME[i])
}
