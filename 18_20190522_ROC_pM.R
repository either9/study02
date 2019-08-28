library("ROCR")
library("OptimalCutpoints")
options(stringsAsFactors=FALSE)

source("./config.R")

old_test_var = "17_CHAT"
new_test_var = "18_ROC_pM"

inf_old_GeneInfo = paste(
                        results_folder, "/",
                        old_test_var, "/",
                        test_list,
                        "_geneInfo.csv", sep="")

dir.create(paste(
                        results_folder, "/",
                        new_test_var, sep=""))

outf_new_GeneInfo = paste(
                        results_folder, "/",
                        new_test_var, "/",
                        test_list,
                        "_geneInfo.csv", sep="")

outf_d = paste(
                        results_folder, "/",
                        new_test_var, "/",
                        test_list,
	                sep="")


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


add_ROC_pM <- function(inf_geneInfo, inf_datExpr, inf_datTraits, outf_geneInfo, outf_d)
{
	status = "pM"

	dir.create(outf_d)

	geneInfo  = read.csv(inf_geneInfo)
	datExpr   = read.csv(inf_datExpr)
	datTraits = read.csv(inf_datTraits)

	
	rownames(datExpr) = datExpr[,1]
	datExpr = datExpr[,-1]
	datExpr = datExpr[!is.na(datTraits$pM),]
	
	rownames(datTraits) = datTraits[,1]
	datTraits = datTraits[,-1]
	datTraits = datTraits[!is.na(datTraits$pM),]

	moduleColors = unique(geneInfo$moduleColor)

	for(n in 1:length(moduleColors))
	{
		dir.create(paste(outf_d, "/", moduleColors[n], sep=""))
	}

	
	new_line_AUC  = rep(NA, nrow(geneInfo))
	new_line_cutoff  = rep(NA, nrow(geneInfo))
	new_line_median  = rep(NA, nrow(geneInfo))

	for(n in 1:length(geneInfo$gene_id))
	{
	
		gene = geneInfo$gene_id[n]
		gene_n = geneInfo$geneSymbol[n]
		modcol = geneInfo$moduleColor[n]

		num = match(gene, names(datExpr))
		data = as.data.frame(cbind(datExpr[,num], datTraits$pM))

		
		cutoff <- optimal.cutpoints(X = "V1", status = "V2", methods = "Youden", data = data, tag.healthy = 0)$Youden$Global$optimal.cutoff$cutoff

		if(length(cutoff)>1)
		{
			num2 = which(abs(cutoff-median)==min(abs(cutoff-median)))
			cutoff = cutoff[num2]

			if(length(cutoff)>1)
			cutoff = cutoff[2]
		}
		
		median = median(data[,1])

		
		pred <- prediction(data[,1], data[,2])
		perf <- performance(pred, "tpr", "fpr")
	
		
		auc <- as.numeric(performance(pred, "auc")@y.values)

			
		file_path = paste(outf_d, "/", modcol, "/", status, "_ROC_", gene, "_", gene_n, ".png", sep="")
		png(file_path)
		plot(perf, main=paste("ROC curve (" , status, ")\n", gene, " (", gene_n, ")\n", "AUC= ", auc, sep=""))
		dev.off()

		new_line_AUC[n]     = auc
		new_line_cutoff[n]  = cutoff
		new_line_median[n]  = median
	}

	new_geneInfo = as.data.frame(cbind(geneInfo[,1:22], new_line_AUC))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_cutoff))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_median))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, geneInfo[,23:ncol(geneInfo)]))

        names(new_geneInfo)[23] = paste("ROC_", status, "_AUC", sep="")
        names(new_geneInfo)[24] = paste("ROC_", status, "_Youden", sep="")
        names(new_geneInfo)[25] = paste("ROC_", status, "_median", sep="")
	
	write.csv(new_geneInfo, file=outf_geneInfo, row.names=FALSE)
	
}



for(i in 1:length(test_list))
{
	add_ROC_pM(
		inf_geneInfo  = inf_old_GeneInfo[i],
		inf_datExpr   = inf_Expr[i],
		inf_datTraits = inf_Traits[i],
		outf_geneInfo = outf_new_GeneInfo[i],
		outf_d        = outf_d[i]
	)
}
