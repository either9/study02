source("./config.R")



merge_featureCounts <- function(seq_list, in_dir, out_dir)
{
	merged_table = 0

	for(num in 1:length(seq_list))
	{
		new_table = read.table(paste(in_dir, "/counts_", seq_list[num], ".txt", sep=""), header = TRUE)	

		if(num==1)
		{
			names(new_table)[7] = seq_list[num]
			merged_table = new_table
		}else{
			sample_n = seq_list[num]
			add = as.data.frame(new_table[,7])
			names(add) = sample_n
			merged_table = cbind(merged_table, add)
		}		
	}
	
	write.csv(merged_table, file= out_dir, row.names=FALSE)
}



merge_featureCounts(	seq_list = seq_list, 
			in_dir = paste(results_folder, 
				       results_featureCounts_STAR_100,
				       sep = ""),
			out_dir = paste(results_folder,
				        results_merged_featureCounts_STAR_100,
					sep = ""))

merge_featureCounts(    seq_list = seq_list,
                        in_dir = paste(results_folder,
                                       results_featureCounts_STAR_50,
                                       sep = ""),
                        out_dir = paste(results_folder,
                                        results_merged_featureCounts_STAR_50,
                                        sep = ""))

