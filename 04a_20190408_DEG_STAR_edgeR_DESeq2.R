source("./config.R")

inf1 = paste(results_folder, results_merged_featureCounts_STAR_100, sep="")
inf2 = paste(results_folder, results_merged_featureCounts_STAR_50, sep="")

dir.create(paste(results_folder, "/04_DEG", sep=""))

outf1 = paste(results_folder, "/04_DEG/test01.csv", sep="")
outf2 = paste(results_folder, "/04_DEG/test02.csv", sep="")
outf3 = paste(results_folder, "/04_DEG/test03.csv", sep="")
outf4 = paste(results_folder, "/04_DEG/test04.csv", sep="")
outf5 = paste(results_folder, "/04_DEG/test05.csv", sep="")
outf6 = paste(results_folder, "/04_DEG/test06.csv", sep="")

library(edgeR)
library(DESeq2)
library(magrittr)

#--------------------------------------------------------
#
# test1 : STAR_results01_100 edgeR likelihood ratio test
# test2 : STAR_results02_50  edgeR likelihood ratio test
# test3 : STAR_results01_100 DESeq2 likelihood ratio test
# test4 : STAR_results02_50  DESeq2 likelihood ratio test
# test5 : STAR_results01_100 DESeq2 Wald test
# test6 : STAR_results02_50  DESeq2 Wald test
#
#--------------------------------------------------------

# import STAR_results01 to table1
# import STAR_results02 to table2
table1 = read.csv(inf1)
rownames(table1) = table1$Geneid
table1 = table1[,-c(1:6)]
table1 = as.matrix(table1)

table2 = read.csv(inf2)
rownames(table2) = table2$Geneid
table2 = table2[,-c(1:6)]
table2 = as.matrix(table2)


# define group design
group_e = factor(c(rep("ACC", 42), rep("normal", 5)))
group_d = data.frame(con = factor(c(rep("ACC", 42), rep("normal", 5))))
design = model.matrix(~ group_e)

# making edgeR object for table1(STAR_results01)
d1 <- DGEList(counts = table1, group = group_e)

	#ok_c = (d1$counts > 5) %>% rowSums %>% {. >0}
	#ok_cpm = d1 %>% cpm %>% {. >1} %>% rowSums %>% {. > 0}
	#d1 = d1[ok_c & ok_cpm, , keep.lib.sizes=FALSE]


# making edgeR object for table2(STAR_results02)
d2 <- DGEList(counts = table2, group = group_e)

	#ok_c = (d2$counts > 5) %>% rowSums %>% {. >0}
	#ok_cpm = d2 %>% cpm %>% {. >1} %>% rowSums %>% {. > 0}
	#d2 = d2[ok_c & ok_cpm, , keep.lib.sizes=FALSE]


# making DESeq2 object for table1, table2
dds1 <- DESeqDataSetFromMatrix(countData = table1, colData = group_d, design = ~ con)
dds2 <- DESeqDataSetFromMatrix(countData = table2, colData = group_d, design = ~ con)


for(i in 1:2)
{
	if(i==1)
	{
		d = d1
		dds_lrt = dds1
		dds_wald = dds1
	}
	if(i==2)
	{
		d = d2
		dds_lrt = dds2
		dds_wald = dds2
	}

	#edgeR
	d <- calcNormFactors(d)
	d <- estimateGLMCommonDisp(d, design)
	d <- estimateGLMTrendedDisp(d, design)
	d <- estimateGLMTagwiseDisp(d, design)

	fit <- glmFit(d, design)
	lrt <- glmLRT(fit, coef = 2)
	results <- as.data.frame(topTags(lrt, n = 100000000))
	sig_e <- subset(results, results$FDR <= 0.05 & (results$logFC >= 1 | results$logFC <= -1))
	

	#DESeq2 likelihood ratio test
	dds_lrt <- DESeq(dds_lrt, test = "LRT", reduced = ~1)
	res_lrt <- results(dds_lrt)
	sig_lrt <- subset(res_lrt, padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1))

	#DESeq2 wald test
        dds_wald <- DESeq(dds_wald, test = "Wald")
        res_wald <- results(dds_wald)
        sig_wald <- subset(res_wald, padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1))


	if(i==1)
	{
		write.csv(sig_e, file=outf1)
                write.csv(sig_lrt, file=outf3)
                write.csv(sig_wald, file=outf5)
	}

	if(i==2)
	{
		write.csv(sig_e, file=outf2)
                write.csv(sig_lrt, file=outf4)
                write.csv(sig_wald, file=outf6)
	}
}
