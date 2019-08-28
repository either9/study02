library(DESeq2)

source("./config.R")

dir.create(paste(results_folder, results_DEmiRNA, sep=""))

inf = paste(results_folder, results_miRNA, "/", results_miRNA_merged, sep="")
outf_up   = paste(results_folder, results_DEmiRNA, "/", results_DEmiRNA_up, sep="")
outf_down = paste(results_folder, results_DEmiRNA, "/", results_DEmiRNA_down, sep="")

Group = miRNA_group

#--------------------------------------

table <- read.csv(inf, header=TRUE)
rownames(table) = table$Geneid
table = table[,-c(1:6)]

group <- data.frame(con = factor(Group))

dds <- DESeqDataSetFromMatrix(countData = as.matrix(table), colData = group, design = ~ con)
dds <- DESeq(dds)
res <- results(dds)

sig <- subset(res, padj <= 0.05)
up   <- subset(sig, log2FoldChange > 0)
down <- subset(sig, log2FoldChange < 0)

write.csv(up, file = outf_up)
write.csv(down, file = outf_down)
