source("./config.R")

STAR_table = 0

for(n in 1:length(miRNA_seq_list))
{
        s_table = read.table(paste(
		results_folder,
		results_miRNA,
		"/counts_",
		miRNA_seq_list[n],
		".txt", sep=""), header = TRUE)        
        
        if(n==1)
        {
                names(s_table)[7] = miRNA_seq_list[n]
                STAR_table = s_table
        }else{
                sample_n = miRNA_seq_list[n]
                STAR_table = cbind(STAR_table, sample_n = s_table[,7])
                names(STAR_table)[n+6] = sample_n
        }               
}

write.csv(STAR_table, file= paste(results_folder, results_miRNA, "/", results_miRNA_merged, sep=""), row.names=FALSE)
