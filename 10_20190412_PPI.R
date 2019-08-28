library("STRINGdb")
library("igraph")

source("./config.R")

old_test_var = "09_WGCNA_step03_relating_to_traits"
new_test_var = "10_PPI"


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


outf_PPI = paste(
		results_folder, "/",
		new_test_var, "/",
		test_list,
		sep="")


String_db <- STRINGdb$new(version="10", species = 9606, score_threshold = 0, input_directory="")


add_PPI <- function(string_db, inf_old_geneInfo, outf_new_geneInfo, outf_ppi)
{
	geneInfo0 = read.csv(inf_old_geneInfo)

	mod_list = unique(geneInfo0$moduleColor)
	new_line1 = rep(NA, nrow(geneInfo0))
	new_line2 = rep(NA, nrow(geneInfo0))
	
	for(n in 1:length(mod_list))
	{
		table = subset(geneInfo0, geneInfo0$moduleColor == mod_list[n])

		map = string_db$map(table, "gene_id", removeUnmappedRows = TRUE)

		if(length(map$STRING_id > 400))
		{
			match = match(paste("MM.", mod_list[n], sep=""), names(map))
			map2  = map[order(map[, match], decreasing = TRUE), ]
			hits  = map2$STRING_id[1:400]
		}else{
			hits = map$STRING_id
		}
		
		png(filename = paste(outf_ppi, "_", mod_list[n], ".png", sep=""))
		string_db$plot_network(hits)
		dev.off()
		
		ig <- string_db$get_subnetwork(hits)
		bw <- igraph::betweenness(ig, directed = TRUE)
		
		match1 = match(names(bw), map$STRING_id)
		bw = as.data.frame(cbind(bw, map$gene_id[match1]))
		names(bw)[2]="gene"
		
		match2 = match(bw$gene, geneInfo0$gene_id)
		new_line1[match2] = rownames(bw)
		new_line2[match2] = bw$bw
	}
	
	geneInfo = as.data.frame(cbind(geneInfo0[,1:3], new_line1))
	geneInfo = as.data.frame(cbind(geneInfo, new_line2))
	geneInfo = as.data.frame(cbind(geneInfo, geneInfo0[,4:ncol(geneInfo0)]))

	names(geneInfo)[4] = "STRING_id"
	names(geneInfo)[5] = "betweenness"
	
	write.csv(geneInfo, file=outf_new_geneInfo, row.names=FALSE)
}


for(i in 1:length(test_list))
{
	add_PPI(	string_db = String_db, 
			inf_old_geneInfo = inf_old_GeneInfo[i], 
			outf_new_geneInfo = outf_new_GeneInfo[i], 
			outf_ppi = outf_PPI[i])
}

