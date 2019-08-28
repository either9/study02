library("data.table")
library("biomaRt")
options(stringsAsFactors=FALSE)

source("./config.R")

old_test_var = "11_survival_analysis"
new_test_var = "14_CMap"

inf_old_GeneInfo = paste(
                        results_folder, "/",
                        old_test_var, "/",
                        test_list,
                        "_geneInfo.csv", sep="")

dir.create(paste(
                        results_folder, "/",
                        new_test_var, sep=""))

outf_CMap = paste(
                results_folder, "/",
                new_test_var, "/",
                test_list,
                sep="")
			
outf_CMap_gmt_name = paste("CMap_", test_list, sep="")

Landmark = fread(ref_CMap_Landmark_gene)
Bing = fread(ref_CMap_BING_gene)

mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
e2e = biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart=mart)
e2e = dplyr::rename(e2e, gene_id = ensembl_gene_id)



generate_gmt_file <- function(inf, outf_d, outf_name, landmark, bing, e2e)
{
	dir.create(outf_d)
	
	geneInfo = read.csv(inf)
	moduleColors = unique(geneInfo$moduleColor)

	match2 = match(as.character(geneInfo$gene_id), e2e$gene_id)
	geneInfo = cbind(geneInfo, entrez_id = e2e$entrezgene[match2])
		
	match3 = match(geneInfo$entrez_id, landmark$'Entrez ID')
	geneInfo = cbind(geneInfo, cmap_landmark = landmark$'Entrez ID'[match3])
		
	match4 = match(geneInfo$entrez_id, bing$'Entrez ID')
	geneInfo = cbind(geneInfo, cmap_bing = bing$'Entrez ID'[match4])

	new_line_all = rep(NA, nrow(geneInfo))
	new_line_all[!is.na(geneInfo$cmap_landmark)] = geneInfo$cmap_landmark[!is.na(geneInfo$cmap_landmark)]
	new_line_all[!is.na(geneInfo$cmap_bing)] = geneInfo$cmap_bing[!is.na(geneInfo$cmap_bing)]
	geneInfo = cbind(geneInfo, cmap_all = new_line_all)
	
	summary = data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
	
	
	for(n in 1:length(moduleColors))
	{
		modcol = moduleColors[n]
		table = subset(geneInfo, geneInfo$moduleColor == modcol)
	
		match = match(paste("MM.", modcol, sep=""), names(table))
		order = order(table[,match], decreasing = TRUE)
		table = table[order, ]

		if(sum(!is.na(table$cmap_landmark)) >=150)
		{
	
			select_table1 = table[!is.na(table$cmap_landmark), ]
			select_table1 = select_table1[1:150, ]
	
		}else{

			select_table_land = table[!is.na(table$cmap_landmark), ]
			num = 150 - sum(!is.na(table$cmap_landmark))
			
			if(sum(!is.na(table$cmap_bing)) >= num){
				
				select_table_bing = table[!is.na(table$cmap_bing), ]
				select_table_bing = select_table_bing[1:num, ]
				
			}else{
			
				select_table_bing = table[!is.na(table$cmap_bing), ]
			}
			
			select_table1 = as.data.frame(rbind(select_table_land, select_table_bing))
		}	
		
		
		if(sum(!is.na(table$cmap_all)) >= 150)
		{
		
			select_table2 = table[!is.na(table$cmap_all), ]
			select_table2 = select_table2[1:150, ]
		
		}else{
		
			select_table2 = table[!is.na(table$cmap_all), ]
		}
		
		
		
		if(ncol(select_table1)>0)
		{
		
			filename = paste(outf_d, "/", outf_name, "_", modcol, "_pattern1.gmt", sep="")
			line = c(	paste(outf_name, "_", modcol, "_pattern1", sep=""),
						paste(outf_name, "_", modcol, "_pattern1", sep=""),
						select_table1$entrez_id)
			line = as.data.frame(t(line))
			write.table(line, file=filename, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
		
		}
		
		if(ncol(select_table2)>0)
		{
		
			filename = paste(outf_d, "/", outf_name, "_", modcol, "_pattern2.gmt", sep="")
			line = c(	paste(outf_name, "_", modcol, "_pattern1", sep=""),
						paste(outf_name, "_", modcol, "_pattern1", sep=""),
						select_table2$entrez_id)
			line = as.data.frame(t(line))
			write.table(line, file=filename, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
		
		}
	
	new_line_summary = data.frame(	moduleColor			= modcol,
									gene_num			= length(table$gene_id),
									entrez_num			= sum(!is.na(table$entrez_id)),
									CMap_Landmark_num	= sum(!is.na(table$cmap_landmark)),
									CMap_BING_num		= sum(!is.na(table$cmap_bing)),
									CMap_all			= sum(!is.na(table$cmap_all)),
									pattern1_gene_num	= length(select_table1$gene_id),
									pattern1_gene		= paste(select_table1$geneSymbol, collapse="|"),
									pattern2_gene_num	= length(select_table2$gene_id),
									pattern2_gene		= paste(select_table2$geneSymbol, collapse="|"))
	
	summary = as.data.frame(rbind(summary, new_line_summary))
	
	}
	
	write.csv(summary, file = paste(outf_d, "/summary.csv", sep=""), row.names=FALSE)

}




for(i in 1:length(test_list))
{
	generate_gmt_file(
		inf = inf_old_GeneInfo[i],
		outf_d = outf_CMap[i],
		outf_name = outf_CMap_gmt_name[i],
		landmark = Landmark,
		bing = Bing,
		e2e = e2e)
}
