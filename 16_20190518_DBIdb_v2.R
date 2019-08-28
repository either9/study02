library("jsonlite")
library("curl")
options(stringsAsFactors=FALSE)


source("./config.R")

old_test_var = "11_survival_analysis"
new_test_var = "16_DBIdb"

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


add_DBIdb <- function(inf, outf)
{
	geneInfo  = read.csv(inf)
	
	new_line_druglist  = rep(NA, nrow(geneInfo))
	new_line_scorelist = rep("", nrow(geneInfo))
	new_line_bestdrug  = rep("", nrow(geneInfo))
	new_line_bestscore = rep(NA, nrow(geneInfo))
	
	for(n in 1:nrow(geneInfo))
	
	{
		gene_n  = geneInfo$geneSymbol[n]
		
		if(!is.na(gene_n))
		{
			order = paste("http://dgidb.org/api/v2/interactions.json?genes=", gene_n, sep="")
		
			request = curl_fetch_memory(order)
			
			data = rawToChar(request$content)
			data = fromJSON(data)
			data = as.data.frame(data$matchedTerm$interactions)
			
			if(ncol(data)!=0)
			{
					new_line_druglist[n]   = paste(data$drugName, collapse = "|")
					new_line_scorelist[n]  = paste(data$score, collapse= "|")
					
					bestdrugs = data$drugName[(data$score == data$score[which.max(data$score)])]
					new_line_bestdrug[n]   = paste(bestdrugs, collapse = "|")
					new_line_bestscore[n]  = data$score[which.max(data$score)]
			}
		}
		
	}
		
	new_geneInfo = as.data.frame(cbind(geneInfo[,1:7], new_line_druglist))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_scorelist))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_bestdrug))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, new_line_bestscore))
	new_geneInfo = as.data.frame(cbind(new_geneInfo, geneInfo[,8:ncol(geneInfo)]))

	names(new_geneInfo)[8]  = "DGIdb_druglist"
	names(new_geneInfo)[9]  = "DGIdb_score"
	names(new_geneInfo)[10] = "DGIdb_bestdrugs"
	names(new_geneInfo)[11] = "DGIdb_bestscore"

	
	write.csv(new_geneInfo, file=outf, row.names=FALSE)
}


for(i in 1:length(test_list))
{
	add_DBIdb(
		inf =  inf_old_GeneInfo[i],
		outf = outf_new_GeneInfo[i])
}
