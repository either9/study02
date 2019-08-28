library("clusterProfiler")
library("org.Hs.eg.db")
library("biomaRt")

options(stringsAsFactors=FALSE)


source("./config.R")

old_test_var = "11_survival_analysis"
new_test_var = "12_enrichment_analysis"

inf_old_GeneInfo = paste(
                        results_folder, "/",
                        old_test_var, "/",
                        test_list,
                        "_geneInfo.csv", sep="")

dir.create(paste(
                        results_folder, "/",
                        new_test_var, sep=""))

outf_Enrichment_analysis = paste(
                results_folder, "/",
                new_test_var, "/",
                test_list,
                sep="")

#pAdjustMethod = "none"
#pvalueCutoff  = 0.05
#qvalueCutoff  = 1


mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
all_gene = biomaRt::getBM(attributes = "ensembl_gene_id", mart=mart)
all_gene = dplyr::rename(all_gene, gene_id = ensembl_gene_id)
all_gene_entrez = bitr(all_gene[,1], fromType="ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")


enrichment_analysis <- function(inf, outf_d, all_gene_entrez)
{
	pAdjustMethod = "none"
	pvalueCutoff  = 0.05
	qvalueCutoff  = 1

	dir.create(outf_d)
	
	geneInfo = read.csv(inf)
	moduleColors = unique(geneInfo$moduleColor)
	
	df <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]	
	
	for(n in 1:length(moduleColors))
	{
		modcol = moduleColors[n]
		table = subset(geneInfo, geneInfo$moduleColor == modcol)
		sub_gene = as.character(table$gene_id)
		
		sub_gene_entrez = bitr(sub_gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
		
		ego_bp <- enrichGO(	gene = sub_gene_entrez[, 2],
							universe = all_gene_entrez[, 2],
							OrgDb = "org.Hs.eg.db",
							ont = "BP",
							pAdjustMethod = pAdjustMethod,
							pvalueCutoff = pvalueCutoff,
							qvalueCutoff = qvalueCutoff,
							readable = TRUE)
							
                ego_mf <- enrichGO(        gene = sub_gene_entrez[, 2],
                                                        universe = all_gene_entrez[, 2],
                                                        OrgDb = "org.Hs.eg.db",
                                                        ont = "MF",
                                                        pAdjustMethod = pAdjustMethod,
                                                        pvalueCutoff = pvalueCutoff,
                                                        qvalueCutoff = qvalueCutoff,
                                                        readable = TRUE)

                ego_cc <- enrichGO(        gene = sub_gene_entrez[, 2],
                                                        universe = all_gene_entrez[, 2],
                                                        OrgDb = "org.Hs.eg.db",
                                                        ont = "CC",
                                                        pAdjustMethod = pAdjustMethod,
                                                        pvalueCutoff = pvalueCutoff,
                                                        qvalueCutoff = qvalueCutoff,
                                                        readable = TRUE)

		ekeg <- enrichKEGG(	gene = sub_gene_entrez[, 2],
							universe = all_gene_entrez[, 2],
							organism = "hsa",
							keyType = "kegg",
							pAdjustMethod = pAdjustMethod,
							pvalueCutoff = pvalueCutoff,
							qvalueCutoff = qvalueCutoff)
		
		res_go_bp <- as.data.frame(ego_bp)
		res_go_mf <- as.data.frame(ego_mf)
		res_go_cc <- as.data.frame(ego_cc)
		res_keg <- as.data.frame(ekeg)
		
		if(length(res_go_bp$ID)>0)
		{
			file.name = paste(outf_d, "/GOenrichment_BP_", modcol, ".csv", sep="")
			write.csv(res_go_bp, file = file.name, row.names=FALSE)
		}
		if(length(res_go_mf$ID)>0)
                {
                        file.name = paste(outf_d, "/GOenrichment_MF_", modcol, ".csv", sep="")
                        write.csv(res_go_mf, file = file.name, row.names=FALSE)
                }
		if(length(res_go_cc$ID)>0)
                {
                        file.name = paste(outf_d, "/GOenrichment_CC_", modcol, ".csv", sep="")
                        write.csv(res_go_cc, file = file.name, row.names=FALSE)
                }

		
		if(length(res_keg$ID)>0)
		{
			file.name = paste(outf_d, "/KEGGenrichment_", modcol, ".csv", sep="")
			write.csv(res_keg, file = file.name, row.names=FALSE)
		}
		
		df <- as.data.frame(rbind(df, data.frame(modcol, length(res_go_bp$ID), length(res_go_mf$ID), length(res_go_cc$ID), length(res_keg$ID))))
		
	}
	
	names(df) = c("moduleColor", "num_GOenriched_term_BP", "num_GOenriched_term_MF", "num_GOenriched_term_CC", "num_KEGGenriched_term")
	file.name = paste(outf_d, "/summary_", "adj-", pAdjustMethod, "_p-", pvalueCutoff, ".csv", sep="")
	write.csv(df, file = file.name, row.names=FALSE)
}



for(i in 1:length(test_list))
{
	enrichment_analysis(
		inf = inf_old_GeneInfo[i],
		outf_d = outf_Enrichment_analysis[i],
		all_gene_entrez = all_gene_entrez
	)
}
