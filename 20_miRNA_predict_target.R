library("stringr")
library("data.table")
library("biomaRt")
library("readxl")
options(stringsAsFactors = FALSE)

source("./config.R")

# targetscan_70 results import
db_ts <- fread(ref_Targetscan)
db_ts <- db_ts[(db_ts$'Gene Tax ID' == 9606), ]
db_ts$miRNA = str_replace(db_ts$miRNA, "[.].", "")

# miRDB results import and colname add
db_RDB <- fread(ref_miRDB)
db_RDB <- db_RDB[str_detect(db_RDB$V1, "hsa"), ]
mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
refseq2g = biomaRt::getBM(attributes = c("refseq_mrna", "external_gene_name"), mart=mart)
db_RDB <- cbind(db_RDB, refseq2g$external_gene_name[match(db_RDB$V2, refseq2g$refseq_mrna)])
names(db_RDB) = c("miRNA", "refseq_ID", "score", "gene_name")
db_RDB <- db_RDB[!is.na(db_RDB$gene_name), ]

# miRTarBase results import
db_tar <- read_excel(ref_miRTarbase)
db_tar <- subset(db_tar, db_tar$`Species (miRNA)`=="Homo sapiens")


# import miRNA_list
miRNA_list = read.csv(
		paste(
                results_folder,
                results_DEmiRNA, "/",
                results_DEmiRNA_down, sep=""))
miRNA_list = as.character(miRNA_list$X)


outf = paste(	results_folder,
		results_miRNA_predict_target, sep="")

# --------------------------------------------------------------------


all_table = data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]

for(n in 1:length(miRNA_list))
{
	
	table_ts  = data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
	table_tar = data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
	table_RDB  = data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]

	mir = miRNA_list[n]
	print(paste("miRNA: ", mir, " is searched", sep=""))

	# Targetscan
	db_sub = subset(db_ts, db_ts$miRNA==mir)

	if(length(db_sub$miRNA)!=0)
	{
		unique_gene = unique(db_sub$'Gene Symbol')

		for(i in 1:length(unique_gene))
		{
			gene = unique_gene[i]

			db_subsub = subset(db_sub, db_sub$'Gene Symbol'==gene)
			score = min(db_subsub$`weighted context++ score`)

			new_line = data.frame(mir, gene, score)
			table_ts = rbind(table_ts, new_line)	
		}
	print("TargetScan databes, complete")
	}



	# miRTarbase
	db_sub = subset(db_tar, db_tar$miRNA == mir)

	if(length(db_sub$miRNA)!=0)
	{
		unique_gene = unique(db_sub$`Target Gene`)

		for(i in 1:length(unique_gene))
		{
			gene = unique_gene[i]

			db_subsub = subset(db_sub, db_sub$`Target Gene`==gene)

			exp = db_subsub$Experiments

			check = 0		

			score_1 = 0
			score_2 = 0
			score_3 = 0
			score_4 = 0
			score_5 = 0
			score_6 = 0	

			if(sum(grepl("Reporter assay", exp))>0)	{score_1 = 1; check=1}
			if(sum(grepl("Western blot", exp))>0)	{score_2 = 1; check=1}
			if(sum(grepl("q-PCR", exp))>0)		{score_3 = 1; check=1}
			if(sum(grepl("NGS", exp))>0) 		{score_4 = 1; check=1}
			if(sum(grepl("pSILAC", exp))>0)		{score_5 = 1; check=1}
			if(check == 0)				{score_6 = 1}

			score_7 = length(exp)

			score_8 = 0
			score_9 = 0

			if((score_1 == 1) | (score_2 == 1) | (score_3 == 1)){score_8 = 1}
			else{score_9 = 1}

			new_line = data.frame(mir, gene, score_1, score_2, score_3, score_4, score_5, score_6, score_7, score_8, score_9)
			table_tar = rbind(table_tar, new_line)	
		}
	print("miRTarBase databes, complete")
	}



	db_sub = subset(db_RDB, db_RDB$miRNA==mir)
	if(length(db_sub$miRNA)!=0)
	{
		unique_gene = unique(db_sub$gene_name)

		for(i in 1:length(unique_gene))
		{
			gene = unique_gene[i]

			db_subsub = subset(db_sub, db_sub$gene_name==gene)
			score = max(db_subsub$score)

			new_line = data.frame(mir, gene, score)
			table_RDB = rbind(table_RDB, new_line)	
		}

		table_RDB = table_RDB[!is.na(table_RDB$gene),]
	
	print("miRDB databes, complete")
	}



	merged_table = data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]

	if(length(table_ts[,1])!=0)
	{
		merged_table = data.frame(matrix(rep(NA, 10), nrow=length(table_ts$mir), ncol=10))
		merged_table = cbind(table_ts, merged_table)
		
		new_table = data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
		
		if(length(table_RDB[,1])!=0)
		{
			for(p in 1:length(table_RDB$mir))
			{
				mir_p  = table_RDB$mir[p]
				gene_p = table_RDB$gene[p]
				score  = table_RDB$score[p]

				check = 0

				for(m in 1:length(merged_table$mir))
				{
					if(mir_p == merged_table$mir[m] & gene_p == merged_table$gene[m])
					{
						merged_table[m, 4] = score
						check = 1
					}
					if(check == 1) break
				}

		
				if(check == 0)
				{
					new_line = data.frame(mir_p, gene_p, NA, score, NA, NA, NA, NA, NA, NA, NA, NA, NA)
					new_table = rbind(new_table, new_line)
				}	
			}		

		names(new_table) = names(merged_table)
		merged_table = rbind(merged_table, new_table)
		new_table = data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
		}
		
		
		if(length(table_tar[,1])!=0)
		{
			for(p in 1:length(table_tar$mir))
			{
				mir_p  = table_tar$mir[p]
				gene_p = table_tar$gene[p]
				score_1  = table_tar$score_1[p]
				score_2  = table_tar$score_2[p]
				score_3  = table_tar$score_3[p]
				score_4  = table_tar$score_4[p]
				score_5  = table_tar$score_5[p]
				score_6  = table_tar$score_6[p]
				score_7  = table_tar$score_7[p]
				socre_8  = table_tar$score_8[p]
				score_9  = table_tar$score_9[p]				

				check = 0
	
				for(m in 1:length(merged_table$mir))
				{
					if(mir_p == merged_table$mir[m] & gene_p == merged_table$gene[m])
					{
						merged_table[m, 5] = score_1
						merged_table[m, 6] = score_2
						merged_table[m, 7] = score_3
						merged_table[m, 8] = score_4
						merged_table[m, 9] = score_5
						merged_table[m, 10] = score_6
						merged_table[m, 11] = score_7
						merged_table[m, 12] = score_8
						merged_table[m, 13] = score_9
					
						check = 1
					}
					if(check == 1) break
				}
		
			
				if(check == 0)
				{
					new_line = data.frame(mir_p, gene_p, NA, NA, score_1, score_2, score_3, score_4, score_5, score_6, score_7, score_8, score_9)
					new_table = rbind(new_table, new_line)
				}	
			}
		names(new_table) = names(merged_table)
		merged_table = rbind(merged_table, new_table)
		}
		
	}else{
	
		if(length(table_RDB[,1])!=0){
		
			merged_table = data.frame(matrix(rep(NA, 9), nrow=length(table_RDB$mir), ncol=9))
			merged_table = cbind(table_RDB[,c(1:2)], merged_table[,1] , table_RDB[,3], merged_table)

			if(length(table_tar[,1])!=0)
			{
			new_table = data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
			
				for(p in 1:length(table_tar$mir))
				{
					mir_p  = table_tar$mir[p]
					gene_p = table_tar$gene[p]
					score_1  = table_tar$score_1[p]
					score_2  = table_tar$score_2[p]
					score_3  = table_tar$score_3[p]
					score_4  = table_tar$score_4[p]
					score_5  = table_tar$score_5[p]
					score_6  = table_tar$score_6[p]
					score_7  = table_tar$score_7[p]
					score_8  = table_tar$score_8[p]
					socre_9  = table_tar$score_9[p]

					check = 0
	
					for(m in 1:length(merged_table$mir))
					{
						if(mir_p == merged_table$mir[m] & gene_p == merged_table$gene[m])
						{
							merged_table[m, 5] = score_1
							merged_table[m, 6] = score_2
							merged_table[m, 7] = score_3
							merged_table[m, 8] = score_4
							merged_table[m, 9] = score_5
							merged_table[m, 10] = score_6
							merged_table[m, 11] = score_7
							merged_table[m, 12] = score_8
							merged_table[m, 13] = score_9
						
							check = 1
						}
						if(check == 1) break
					}
		
			
					if(check == 0)
					{
						new_line = data.frame(mir_p, gene_p, NA, NA, score_1, score_2, score_3, score_4, score_5, score_6, score_7, score_8, score_9)
						new_table = rbind(new_table, new_line)
					}	
				}
			names(new_table) = names(merged_table)
			merged_table = rbind(merged_table, new_table)
			}
		}else{
		
			if(length(table_tar[,1])!=0)
			{
				merged_table = data.frame(matrix(rep(NA, 4), nrow=length(table_tar$mir), ncol=4))
				merged_table = cbind(table_tar[,c(1:2)], merged_table, table_tar[,c(3:11)])			
			}
		}
	}
	
	if(length(merged_table[,1])!=0)
	{
		names(merged_table) = names(all_table)
		all_table = rbind(all_table, merged_table)
	}
}

names(all_table) = c(
	 "miRNA",
	 "target_gene",
	 "TargetScan_score",
	 "miRDB_score",
	 "miRTarBase_Reporter.assay",
	 "miRTarBase_Western.blot",
	 "miRTarBase_q-PCR",
	 "miRTarBase_NGS",
	 "miRTarBase_pSILAC",
	 "miRTarBase_only-others",
	 "miRTarBase_sum_of_papers",
	 "miRTarBase_strong_evidence",
	 "miRTarBase_weak_evidence")

write.csv(all_table, file = outf, row.names= FALSE)

