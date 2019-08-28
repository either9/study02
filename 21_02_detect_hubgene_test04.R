 library("STRINGdb")
library("igraph")
library("data.table")
options(stringsasFactors = FALSE)


source("./config.R")

old_test_var = "18_ROC_pM"


inf_miR_target_list = paste(    results_folder,
                                results_miRNA_predict_target, sep="")


inf_geneInfo = paste(
                        results_folder, "/",
                        old_test_var, "/",
                        test_list,
                        "_geneInfo.csv", sep="")

inf_mod = pM_related_module

dir.create(paste(results_folder, results_detect_hubgene, sep=""))

outf_name = paste(
                results_folder,
                results_detect_hubgene,
                "/test04_withGS.csv", sep="")

#################################
#
# load stringdb
# load miR data
#
#################################

string_db <- STRINGdb$new(version="10", species = 9606, score_threshold = 0, input_directory="")
miR_table_all = fread(inf_miR_target_list)

#################################
#
# generate cutoff points
#
#################################

cutoff_MM = seq(0.6, 0.9, by = 0.05)
cutoff_GS = as.numeric(c(0.2, 0.3))
cutoff_betweenness = as.numeric(c(3,5,10))

#cutoff_miR_TargetScan = seq(-0.5, -0.9, by = -0.05)
#cutoff_miR_miRDB = seq(80, 95, by = 5)
#cutoff_miR_miRTarbase = as.numeric(c(0,1))

cutoffs_gene = expand.grid(
	MM = cutoff_MM,
	GS = cutoff_GS,
	betweenness = cutoff_betweenness)

#cutoffs_miRNA = expand.grid(	
#	miR_TargetScan = cutoff_miR_TargetScan,
#	miR_miRDB = cutoff_miR_miRDB,
#	miR_miRTarbase = cutoff_miR_miRTarbase)

	#################################
	#
	# generate output data
	#
	#################################
	outf_data = data.frame(matrix(rep(NA, 9), nrow=1))[numeric(0),]

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------



for(i in 1:8){

	geneInfo0 = read.csv(inf_geneInfo[i])
	geneInfo0 = subset(geneInfo0, geneInfo0$moduleColor == inf_mod[i])

	print(paste("test ", i, " checking start", sep=""))

for(n in 1:nrow(cutoffs_gene)){

	print(paste("test ", i, ", num", n, " start", sep=""))
	
	###################
	# adapte MM cutoff
	###################
	
	geneInfo = subset(geneInfo0, geneInfo0[,(names(geneInfo0)==paste("MM.", inf_mod[i], sep=""))] >= cutoffs_gene$MM[n])
	if(nrow(geneInfo)==0) next
	
	###################
	# adapte GS.pM cutoff
	###################
	
	geneInfo = subset(geneInfo0, geneInfo0$GS.pM >= cutoffs_gene$GS[n])
	if(nrow(geneInfo)==0) next
	
	
	####################
	# perform PPI analysis
	####################

	# remove "STRING_id" colum in geneInfo
	geneInfo = geneInfo[, !(names(geneInfo)=="STRING_id")]
	
	map = string_db$map(geneInfo, "gene_id", removeUnmappedRows = TRUE)
	if(nrow(map)==0) next

	if(length(map$STRING_id > 400))
	{
		match = match(paste("MM.", inf_mod[i], sep=""), names(map))
		map2  = map[order(map[, match], decreasing = TRUE), ]
		hits  = map2$STRING_id[1:400]
	}else{
		hits = map$STRING_id
	}

	ig <- string_db$get_subnetwork(hits)
	bw <- igraph::betweenness(ig, directed = TRUE)
		
	match1 = match(names(bw), map$STRING_id)
	bw = data.frame(bw, gene_name = map$geneSymbol[match1])
	
	
	###################
	# adapte betweenness cutoff
	###################
	
	#cut_num = round(nrow(map) * cutoffs_gene$betweenness[n])
	cut_num = cutoffs_gene$betweenness[n]

	if(nrow(bw) <= cut_num){
		gene_list = bw$gene_name
	}else{
		gene_list = bw$gene_name[1:cut_num]
	}


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#
#	###################
#	# making new miRNA-mRNA network
#	###################
#
#	miR_table = subset(miR_table_all,
#				miR_table_all$TargetScan_score				<= as.numeric(cutoffs$miR_Targetscan[n])|
#				miR_table_all$miRDB_score					>= as.numeric(cutoffs$miR_miRDB[n])|
#				miR_table_all$miRTarBase_strong_evidence	>= as.numeric(cutoffs$miR_miRTarbase[n]))
#
#	miR_table = miR_table[is.element(miR_table$target_gene, gene_list), ]
#
#
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


	###################
	# checking by survival curve and log-lank test
	###################

	geneInfo_sub = subset(geneInfo, is.element(geneInfo$geneSymbol, gene_list))
	
	t = sum((geneInfo_sub$surv_curv_pval <= 0.05))/nrow(geneInfo_sub)
	good = geneInfo_sub$geneSymbol[geneInfo_sub$surv_curv_pval <= 0.05]
	
			new_line = line(1:9)
			new_line[1] = i
			new_line[2] = n
			new_line[3] = cutoffs_gene$MM[n]
			new_line[4] = cutoffs_gene$GS[n]
			new_line[5] = cutoffs_gene$betweenness[n]
			new_line[6] = t
			new_line[7] = length(gene_list)
			new_line[8] = length(good)
			new_line[9] = paste(good, collapse="|")
		
			names(new_line) = c(
				"test_num",
				"condition_num",
				"cutoff_MM",
				"cutoff_GS",
				"cutoff_bw",
				"ok_ratio",
				"selected_gene_num",
				"good_gene_num",
				"good_gene_list")
						
			outf_data = rbind(outf_data, new_line)	
			outf_data$good_gene_list = as.character(outf_data$good_gene_list)
}
}

outf_data = outf_data[order(outf_data$ok_ratio, decreasing = TRUE), ]

write.csv(outf_data, file = outf_name, row.names=FALSE)

