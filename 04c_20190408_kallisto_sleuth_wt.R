#options(mc.cores = 0)
options(stringsAsFactors = FALSE)

source("./config.R")

inf = paste(results_folder, results_kallisto, sep="")
con_f = c(rep("acc", 42), rep("normal", 5))

outf_so  = paste(results_folder, "/04_DEG/sleuth_test08.so", sep="")
outf_deg = paste(results_folder, "/04_DEG/test08.csv", sep="")
outf_counts = paste(results_folder, "/04_DEG/sleuth_matrix_test08.csv", sep="")

#----------------------------------------
#
# test07 kallisto, sleuth likelihood ratio test
#
#----------------------------------------


library("sleuth")
library("biomaRt")

mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
  "ensembl_gene_id", "external_gene_name", "description",
  "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name) 
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene')) 


sample_id <- dir(inf)
kal_dirs <- file.path(inf, sample_id)
s2c <- data.frame(sample = seq_list, condition = con_f)
s2c <- dplyr::mutate(s2c, path= kal_dirs)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = ttg, aggregation_column='ens_gene', gene_mode=TRUE, num_cores = 4)

so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
#so <- sleuth_lrt(so, 'reduced', 'full')
so <- sleuth_wt(so, 'conditionnormal', 'full')

#res <- sleuth_results(so, 'reduced:full', test_type='lrt')
res <- sleuth_results(so, 'conditionnormal', test_type='wt')
sig_res <- subset(res, qval <= 0.05)

table = sleuth_to_matrix(so, "obs_norm", "scaled_reads_per_base")
table = as.data.frame(table)

sleuth_save(so, outf_so)
write.csv(sig_res, file=outf_deg, row.names=FALSE)
write.csv(table, file=outf_counts)
