library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source("./config.R")
test_list = test_list

inf_Expr = paste(       results_folder, 
                        "/05_WGCNA_step01_datainput/",
                        "datExpr_",
                         test_list,
                         ".csv",
                         sep="")

inf_Traits = paste(     results_folder,
                        "/05_WGCNA_step01_datainput/",
                        "datTraits_",
                         test_list,
                         ".csv",
                         sep="")

inf_ME_RData = paste(   results_folder,
                        "/08_WGCNA_step02-3_mergedME/",
                         test_list,
                         ".RData",
                         sep="")


dir.create(paste(results_folder,
                 "/09_WGCNA_step03_relating_to_traits",
                 sep=""))

Outf = paste(           results_folder,
			"/09_WGCNA_step03_relating_to_traits/",
                         test_list,
                         sep="")

library("biomaRt")
library("dplyr")
mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

TTG <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)
TTG <- dplyr::rename(TTG, ensembl_gene_id = ensembl_gene_id, geneSymbol = external_gene_name)


relating_to_traits <- function(inf_expr, inf_traits, inf_me_RData, ttg, outf)
{
	lnames = load(file= inf_me_RData)
	datExpr = read.csv(file=inf_expr, row.names=1)
	datTraits = read.csv(file=inf_traits, row.names=1)

	nGenes = ncol(datExpr);
	
	nSamples = nrow(datExpr);
	
	MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
	MEs = orderMEs(MEs0)

	moduleTraitCor = cor(MEs, datTraits, use = "p");
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


	textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");

	dim(textMatrix) = dim(moduleTraitCor)


	outf1 = paste(outf, "_module-trait-relationships.png", sep="")
	png(filename = outf1)
	par(mar = c(6, 8.5, 3, 3));
	labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
	dev.off()

	
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
modNames = substring(names(MEs), 3)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


trait_1 = as.data.frame(datTraits$pT);
names(trait_1) = "pT"
geneTraitSignificance_1 = as.data.frame(cor(datExpr, trait_1, use = "p"));
GSPvalue_1 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_1), nSamples));
names(geneTraitSignificance_1) = paste("GS.", names(trait_1), sep="");
names(GSPvalue_1) = paste("p.GS.", names(trait_1), sep="");

trait_2 = as.data.frame(datTraits$PNI);
names(trait_2) = "pNI"
geneTraitSignificance_2 = as.data.frame(cor(datExpr, trait_2, use = "p"));
GSPvalue_2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_2), nSamples));
names(geneTraitSignificance_2) = paste("GS.", names(trait_2), sep="");
names(GSPvalue_2) = paste("p.GS.", names(trait_2), sep="");

trait_3 = as.data.frame(datTraits$pN);
names(trait_3) = "pN"
geneTraitSignificance_3 = as.data.frame(cor(datExpr, trait_3, use = "p"));
GSPvalue_3 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_3), nSamples));
names(geneTraitSignificance_3) = paste("GS.", names(trait_3), sep="");
names(GSPvalue_3) = paste("p.GS.", names(trait_3), sep="");

trait_4 = as.data.frame(datTraits$pM);
names(trait_4) = "pM"
geneTraitSignificance_4 = as.data.frame(cor(datExpr, trait_4, use = "p"));
GSPvalue_4 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_4), nSamples));
names(geneTraitSignificance_4) = paste("GS.", names(trait_4), sep="");
names(GSPvalue_4) = paste("p.GS.", names(trait_4), sep="");

trait_5 = as.data.frame(datTraits$FU.months);
names(trait_5) = "FU.months"
geneTraitSignificance_5 = as.data.frame(cor(datExpr, trait_5, use = "p"));
GSPvalue_5 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_5), nSamples));
names(geneTraitSignificance_5) = paste("GS.", names(trait_5), sep="");
names(GSPvalue_5) = paste("p.GS.", names(trait_5), sep="");

trait_6 = as.data.frame(datTraits$FU.status);
names(trait_6) = "FU.status"
geneTraitSignificance_6 = as.data.frame(cor(datExpr, trait_6, use = "p"));
GSPvalue_6 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_6), nSamples));
names(geneTraitSignificance_6) = paste("GS.", names(trait_6), sep="");
names(GSPvalue_6) = paste("p.GS.", names(trait_6), sep="");




# Create the starting data frame
probes = names(datExpr)
probes2annot = match(probes, ttg$ensembl_gene_id)

geneInfo0 = data.frame(gene_id = probes,
		      geneSymbol  = ttg$geneSymbol[probes2annot],
                      moduleColor = moduleColors,
                      geneTraitSignificance_1,
                      GSPvalue_1,
                      geneTraitSignificance_2,
                      GSPvalue_2,
                      geneTraitSignificance_3,
                      GSPvalue_3,
                      geneTraitSignificance_4,
                      GSPvalue_4,
                      geneTraitSignificance_5,
                      GSPvalue_5,
                      geneTraitSignificance_6,
                      GSPvalue_6)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, trait_4, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pM));
geneInfo = geneInfo0[geneOrder, ]

outf2 = paste(outf, "_geneInfo.csv", sep="")
write.csv(geneInfo, file = outf2, row.names=FALSE)

}


for(i in 1:length(test_list))
{
	relating_to_traits(	inf_expr = inf_Expr[i],
				inf_traits = inf_Traits[i],
				inf_me_RData = inf_ME_RData[i],
				ttg = TTG,
				outf = Outf[i])
}
