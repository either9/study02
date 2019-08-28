library("survival")
options(stringsAsFactors=FALSE)

source("./config.R")

old_test_var = "10_PPI"
new_test_var = "11_survival_analysis"

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

outf_Survival_Curve = paste(
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

P_value_cutoff = 0.05


survival_analysis <- function(inf_expr, inf_traits, inf_old_geneInfo, outf_new_geneInfo, outf_survival_curve, p_value_cutoff)
{
	dir.create(outf_survival_curve)
	dir.create(paste(outf_survival_curve, "/signif", sep=""))
	
	datExpr   = read.csv(inf_expr)
	datTraits = read.csv(inf_traits)
	geneInfo  = read.csv(inf_old_geneInfo)
	
	# remove sample no.14
	# because it lack FU.months property
	datTraits = datTraits[-14, ]
	datExpr   = datExpr[-14, ]
	
	# remove another property in datTraits
	datTraits = datTraits[, -c(2:7)]
	
	new_line_pval   = rep(NA, nrow(geneInfo))
	new_line_median = rep(NA, nrow(geneInfo))
	
	for(n in 1:nrow(geneInfo))
	{
		gene_id = geneInfo$gene_id[n]
		gene_n  = geneInfo$geneSymbol[n]
		modcol  = geneInfo$moduleColor[n]
		
		num = match(geneInfo$gene_id[n], names(datExpr))
		data = cbind(datTraits, datExpr[num])
		names(data) = c("ID", "time", "status", "gene_counts")

		new_line_median[n] = median(datExpr[,num], na.rm = TRUE)
		
		for(k in 1:length(data$gene_counts))
		{
			if(data$gene_counts[k] >= new_line_median[n])	data$gene_counts[k] = 1
			else											data$gene_counts[k] = 0
		}
	
		data$time=as.numeric(data$time)
	
		if(length(unique(data$gene_counts))!=1)
		{
			Surv(data$time, data$status)
	
			data.sf = survfit(Surv(time, status)~gene_counts, data=data)
			dif <- survdiff(Surv(time, status==1)~gene_counts, data=data)
			new_line_pval[n] = pchisq(dif$chisq, 1, lower.tail=FALSE)
			
			png(paste(outf_survival_curve, "/", "Kaplan-Meier_plot_", modcol, "_", gene_id, "_", gene_n, ".png", sep=""))
			plot(data.sf, lty=1:2, mark.t=F, ylim=c(0,1.0), main=paste("KAplan-Meier_plot(median)\n", gene_id, "(", gene_n, ")\n", "p = " , new_line_pval[n]))
			dev.off()
		
				if(new_line_pval[n] <= p_value_cutoff)
				{
					png(paste(outf_survival_curve, "/signif/", "sig_Kaplan-Meier_plot_", modcol, "_", gene_id, "_", gene_n, ".png", sep=""))
					plot(data.sf, lty=1:2, mark.t=F, ylim=c(0,1.0), main=paste("KAplan-Meier_plot(median)\n", gene_id, "(", gene_n, ")\n", "p = " , new_line_pval[n]))
					dev.off()
				}
		}
	}
	
	new_geneInfo = as.data.frame(cbind(geneInfo[,1:5], new_line_median))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_pval))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, geneInfo[,6:ncol(geneInfo)]))

	names(new_geneInfo)[6] = "surv_curv_cutoff"
	names(new_geneInfo)[7] = "surv_curv_pval"
	
	write.csv(new_geneInfo, file=outf_new_geneInfo, row.names=FALSE)
}


for(i in 1:length(test_list))
{
	survival_analysis(
		inf_expr = inf_Expr[i],
		inf_traits = inf_Traits[i],
		inf_old_geneInfo = inf_old_GeneInfo[i],
		outf_new_geneInfo = outf_new_GeneInfo[i],
		outf_survival_curve = outf_Survival_Curve[i],
		p_value_cutoff = P_value_cutoff)
}
