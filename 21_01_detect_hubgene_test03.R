library("STRINGdb")
library("igraph")
library("data.table")
options(stringsasFactors = FALSE)

source("./config.R")

old_test_var = "18_ROC_pM"


inf_miR_target_list = paste(	results_folder,
		                results_miRNA_predict_target, sep="")


inf_geneInfo = paste(
                        results_folder, "/",
                        old_test_var, "/",
                        test_list,
                        "_geneInfo.csv", sep="")

inf_mod = pM_related_module

dir.create(paste(results_folder, results_detect_hubgene, sep=""))

outf_name_1 = paste(
		results_folder,
		results_detect_hubgene,
		"/test03_1_withGS.csv", sep="")

outf_name_2 = paste(
                results_folder,
                results_detect_hubgene,
                "/test03_2_withGS.csv", sep="")


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
cutoff_betweenness = as.numeric(c(0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.75, 1))

cutoff_miR_TargetScan = c(0, seq(-0.5, -0.9, by = -0.05))
cutoff_miR_miRDB = c(49, seq(80, 95, by = 5))
cutoff_miR_miRTarbase = as.numeric(c(0,1))

cutoffs_gene = expand.grid(
	MM = cutoff_MM,
	GS = cutoff_GS,
	betweenness = cutoff_betweenness)

cutoffs_miRNA = expand.grid(	
	miR_TargetScan = cutoff_miR_TargetScan,
	miR_miRDB = cutoff_miR_miRDB,
	miR_miRTarbase = cutoff_miR_miRTarbase)

	#################################
	#
	# generate output data
	#
	#################################
	outf_data_1 = data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0),]
	outf_data_2 = data.frame(matrix(rep(NA, 14), nrow=1))[numeric(0),]


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# check all miRNA-mRNA networks in modules
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


for(i in 1:8)
{

	geneInfo0 = read.csv(inf_geneInfo[i])
	geneInfo0 = subset(geneInfo0, geneInfo0$moduleColor == inf_mod[i])

	group_ok = subset(geneInfo0, geneInfo0$surv_curv_pval <= 0.05)

	for(n in 1:nrow(cutoffs_miRNA))
	{
	
		miR_table = subset(miR_table_all,
				miR_table_all$TargetScan_score <= cutoffs_miRNA$miR_TargetScan[n]|
				miR_table_all$miRDB_score >= cutoffs_miRNA$miR_miRDB[n]|
				miR_table_all$miRTarBase_strong_evidence >= cutoffs_miRNA$miR_miRTarbase[n])
		
		miR_table_b = miR_table[is.element(miR_table$target_gene, geneInfo0$geneSymbol), ]
		if(nrow(miR_table_b)==0) next
		
		print(paste("test ", i, ", cutoffs_miRNA ", n, " start", sep=""))
		
		miR_list = unique(miR_table_b$miRNA)
		
		
		for(j in 1:length(miR_list))
		{
			miR_table_sub = subset(miR_table_b, miR_table_b$miRNA == miR_list[j])
			
			good = intersect(miR_table_sub$target_gene, group_ok$geneSymbol)
			t =length(good)/nrow(miR_table_sub)
			
			new_line = list(1:10)
			new_line[1] = i
			new_line[2] = n
			new_line[3] = cutoffs_miRNA$miR_TargetScan[n]
			new_line[4] = cutoffs_miRNA$miR_miRDB[n]
			new_line[5] = cutoffs_miRNA$miR_miRTarbase[n]
			new_line[6] = miR_list[j]
			new_line[7] = t
			new_line[8] = nrow(miR_table_sub)
			new_line[9] = length(good)
			new_line[10] = paste(good, collapse="|")
			
			names(new_line) = c(
				"test_num",
				"condition_num",
				"cutoff_TargetScan",
				"cutoff_miRDB",
				"cutoff_miRTarbase",
				"miRNA",
				"score",
				"target_gene_num",
				"good_gene_num",
				"good_gene_list")
				
			outf_data_1 = rbind(outf_data_1, new_line)

			outf_data_1$miRNA = as.character(outf_data_1$miRNA)
			outf_data_1$good_gene_list = as.character(outf_data_1$good_gene_list)			

		}
	}
}

names(outf_data_1) = c(
				"test_num",
				"condition_num",
				"cutoff_TargetScan",
				"cutoff_miRDB",
				"cutoff_miRTarbase",
				"miRNA",
				"score",
				"target_gene_num",
				"good_gene_num",
				"good_gene_list")

outf_data_1 = outf_data_1[order(outf_data_1$score, decreasing = TRUE), ]

write.csv(outf_data_1, file = outf_name_1, row.names=FALSE)








#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# check all miRNA-mRNA networks in MM-PPI filterd genes
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


for(i in 1:8){

	geneInfo0 = read.csv(inf_geneInfo[i])
	geneInfo0 = subset(geneInfo0, geneInfo0$moduleColor == inf_mod[i])

for(n in 1:nrow(cutoffs_gene)){

	
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
	
	gene_cutoff_nums = round(nrow(map) * cutoffs_gene$betweenness[n])
	
	if(nrow(bw) <= gene_cutoff_nums){
		gene_list = bw$gene_name
	}else{
		gene_list = bw$gene_name[1:gene_cutoff_nums]
	}

	geneInfo = subset(geneInfo, is.element(geneInfo$geneSymbol, gene_list))
	group_ok = subset(geneInfo0, geneInfo$surv_curv_pval <= 0.05)
	if(nrow(group_ok)==0) next
	
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

	for(j in 1:nrow(cutoffs_miRNA))
	{
	
		###################
		# making new miRNA-mRNA network
		###################
	
	if(cutoffs_miRNA$miR_TargetScan[j] == 0)
	{
		if(cutoffs_miRNA$miR_miRDB[j] == 49)
		{
			if(cutoffs_miRNA$miR_miRTarbase[j] == 0)
			{
				
			}else{
			
				miR_table = subset(miR_table_all,
					miR_table_all$miRTarBase_strong_evidence >= cutoffs_miRNA$miR_miRTarbase[j])
			}			
		
		}else{
		
			if(cutoffs_miRNA$miR_miRTarbase[j] == 0)
			{
				miR_table = subset(miR_table_all,
					miR_table_all$miRDB_score >= cutoffs_miRNA$miR_miRDB[j])
			}else{
			
				miR_table = subset(miR_table_all,
					miR_table_all$miRDB_score >= cutoffs_miRNA$miR_miRDB[j]&
					miR_table_all$miRTarBase_strong_evidence >= cutoffs_miRNA$miR_miRTarbase[j])
			}		
		}
	}else{
	
		if(cutoffs_miRNA$miR_miRDB[j] == 49)
		{
			if(cutoffs_miRNA$miR_miRTarbase[j] == 0)
			{
				miR_table = subset(miR_table_all,
					miR_table_all$TargetScan_score <= cutoffs_miRNA$miR_TargetScan[j])
			}else{
			
				miR_table = subset(miR_table_all,
					miR_table_all$TargetScan_score <= cutoffs_miRNA$miR_TargetScan[j]&
					miR_table_all$miRTarBase_strong_evidence >= cutoffs_miRNA$miR_miRTarbase[j])
			}			
		
		}else{
		
			if(cutoffs_miRNA$miR_miRTarbase[j] == 0)
			{
				miR_table = subset(miR_table_all,
					miR_table_all$TargetScan_score <= cutoffs_miRNA$miR_TargetScan[j]&
					miR_table_all$miRDB_score >= cutoffs_miRNA$miR_miRDB[j])
			}else{
			
				miR_table = subset(miR_table_all,
					miR_table_all$TargetScan_score <= cutoffs_miRNA$miR_TargetScan[j]&
					miR_table_all$miRDB_score >= cutoffs_miRNA$miR_miRDB[j]&
					miR_table_all$miRTarBase_strong_evidence >= cutoffs_miRNA$miR_miRTarbase[j])
			}		
		}	
	}



		miR_table = miR_table[is.element(miR_table$target_gene, gene_list), ]
		if(nrow(miR_table)==0) next
		
		miR_list = unique(miR_table$miRNA)
		
		
		for(k in 1:length(miR_list))
		{
			miR_table_sub = subset(miR_table, miR_table$miRNA == miR_list[k])
			
			good = intersect(miR_table_sub$target_gene, group_ok$geneSymbol)
			t =length(good)/nrow(miR_table_sub)
			
			new_line = list(1:14)
			new_line[1] = i
			new_line[2] = n
			new_line[3] = j
			new_line[4] = cutoffs_gene$MM[n]
			new_line[5] = cutoffs_gene$GS[n]
			new_line[6] = cutoffs_gene$betweenness[n]
			new_line[7] = cutoffs_miRNA$miR_TargetScan[j]
			new_line[8] = cutoffs_miRNA$miR_miRDB[j]
			new_line[9] = cutoffs_miRNA$miR_miRTarbase[j]
			new_line[10] = miR_list[k]
			new_line[11] = t
			new_line[12] = nrow(miR_table_sub)
			new_line[13] = length(good)
			new_line[14] = paste(good, collapse="|")
			
			names(new_line) = c(
				"test_num",
				"condition_gene_num",
				"condition_miRNA_num",
				"cutoff_MM",
				"cutoff_GS",
				"cutoff_PPI_ratio",
				"cutoff_TargetScan",
				"cutoff_miRDB",
				"cutoff_miRTarbase",
				"miRNA",
				"score",
				"target_gene_num",
				"good_gene_num",
				"good_gene_list")
				
			outf_data_2 = rbind(outf_data_2, new_line)
			outf_data_2$miRNA = as.character(outf_data_2$miRNA)
			outf_data_2$good_gene_list = as.character(outf_data_2$good_gene_list)
		}
	}
}
}

names(outf_data_2) = c(
				"test_num",
				"condition_gene_num",
				"condition_miRNA_num",
				"cutoff_MM",
				"cutoff_GS",
				"cutoff_PPI_ratio",
				"cutoff_TargetScan",
				"cutoff_miRDB",
				"cutoff_miRTarbase",
				"miRNA",
				"score",
				"target_gene_num",
				"good_gene_num",
				"good_gene_list")

outf_data_2 = outf_data_2[order(outf_data_2$score, decreasing = TRUE), ]

write.csv(outf_data_2, file = outf_name_2, row.names=FALSE)

