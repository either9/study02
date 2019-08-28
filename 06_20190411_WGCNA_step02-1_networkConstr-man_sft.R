library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source("./config.R")
test_list = test_list

inf_Expr = paste(	results_folder, 
			"/05_WGCNA_step01_datainput/",
			"datExpr_",
			 test_list,
			 ".csv",
			 sep="")

dir.create(	paste(results_folder,
		"/06_WGCNA_step02-1_sft",
		sep=""))

outf_Sft = paste(       results_folder,
                        "/06_WGCNA_step02-1_sft/",
                         test_list,
                         sep="")


calculate_sft <- function(inf_expr, outf_sft)
{
	datExpr   = read.csv(file=inf_expr, row.names=1)

	powers = c(c(1:10), seq(from = 12, to=20, by=2))
	sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)



        png(filename=paste(outf_sft, "_sft.png", sep=""))

	par(mfrow = c(1,2))
	
	cex1 = 0.9

	plot(	sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
		main = paste("Scale independence"));

	text(	sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		labels=powers,cex=cex1,col="red");

	abline(h=0.90,col="red")

	dev.off()

	png(filename=paste(outf_sft, "_con.png", sep=""))
	plot(	sft$fitIndices[,1], sft$fitIndices[,5],
	     	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	     	main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

	dev.off()
}



for(i in 1:length(test_list))
{
	calculate_sft(	inf_expr = inf_Expr[i],
			outf_sft = outf_Sft[i])
}
