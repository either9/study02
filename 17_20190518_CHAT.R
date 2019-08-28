library("curl")
library("data.table")
options(stringsAsFactors=FALSE)


source("./config.R")

old_test_var = "16_DBIdb"
new_test_var = "17_CHAT"

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


add_CHAT <- function(inf, outf)
{
	geneInfo  = read.csv(inf)
	
	new_line_1  = rep(NA, nrow(geneInfo))
	new_line_2  = rep(NA, nrow(geneInfo))
	new_line_3  = rep(NA, nrow(geneInfo))
	new_line_4  = rep(NA, nrow(geneInfo))
	new_line_5  = rep(NA, nrow(geneInfo))
	new_line_6  = rep(NA, nrow(geneInfo))
	new_line_7  = rep(NA, nrow(geneInfo))
	new_line_8  = rep(NA, nrow(geneInfo))
	new_line_9  = rep(NA, nrow(geneInfo))
	new_line_10  = rep(NA, nrow(geneInfo))
	new_line_all  = rep(NA, nrow(geneInfo))
	
	
	for(n in 1:nrow(geneInfo))
	{
		gene_n  = geneInfo$geneSymbol[n]
		
		if(!is.na(gene_n))
		{
			order = paste("http://chat.lionproject.net/chartdata?q=", gene_n, "&measure=count&hallmarks=top", sep="")
			request = curl_fetch_memory(order)
			data  = rawToChar(request$content)
			data  = fread(data)
			
			new_line_1[n] = data$count[1]
			new_line_2[n] = data$count[2]
			new_line_3[n] = data$count[3]
			new_line_4[n] = data$count[4]
			new_line_5[n] = data$count[5]
			new_line_6[n] = data$count[6]
			new_line_7[n] = data$count[7]
			new_line_8[n] = data$count[8]
			new_line_9[n] = data$count[9]
			new_line_10[n] = data$count[10]
			new_line_all[n] = data$count[1] + data$count[2] + data$count[3] + data$count[4] + data$count[5]
							 + data$count[6] + data$count[7] + data$count[8] + data$count[9] + data$count[10]

		}
		
	}
		
	new_geneInfo = as.data.frame(cbind(geneInfo[,1:11], new_line_1))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_2))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_3))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_4))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_5))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_6))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_7))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_8))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_9))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_10))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_all))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, geneInfo[,12:ncol(geneInfo)]))

	names(new_geneInfo)[12]  = "CHAT_invasion_and_metastasis"
	names(new_geneInfo)[13]  = "CHAT_immune_destruction"
	names(new_geneInfo)[14]  = "CHAT_cellular_enegetics"
	names(new_geneInfo)[15]  = "CHAT_replicative_immortality"
	names(new_geneInfo)[16]  = "CHAT_evading_growth_suppressors"
	names(new_geneInfo)[17]  = "CHAT_genome_instability_and_mutation"
	names(new_geneInfo)[18]  = "CHAT_inducing_angiogenesis"
	names(new_geneInfo)[19]  = "CHAT_resisting_cell_death"
	names(new_geneInfo)[20]  = "CHAT_sustaining_proliferative_signaling"
	names(new_geneInfo)[21]  = "CHAT_tumor_promoting_inflammation"
	names(new_geneInfo)[22]  = "CHAT_all_counts"

	write.csv(new_geneInfo, file=outf, row.names=FALSE)
}


for(i in 7:length(test_list))
{
        add_CHAT(
                inf =  inf_old_GeneInfo[i],
                outf = outf_new_GeneInfo[i])
}
