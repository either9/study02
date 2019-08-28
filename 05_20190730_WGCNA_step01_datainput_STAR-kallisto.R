library(WGCNA)
options(stringsAsFactors = FALSE)

source("./config.R")

test_list = test_list

inf_Counts = c(	paste(results_folder, results_merged_featureCounts_STAR_100, sep=""),
		paste(results_folder, results_merged_featureCounts_STAR_50, sep=""),
		paste(results_folder, results_merged_featureCounts_STAR_100, sep=""),
                paste(results_folder, results_merged_featureCounts_STAR_50, sep=""),
		paste(results_folder, results_merged_featureCounts_STAR_100, sep=""),
                paste(results_folder, results_merged_featureCounts_STAR_50, sep=""),
		paste(results_folder, "/04_DEG/sleuth_matrix_test07.csv", sep=""),
		paste(results_folder, "/04_DEG/sleuth_matrix_test08.csv", sep="")
		)

which_STAR_or_kallisto = as.character(c(
				"STAR",
				"STAR",
				"STAR",
				"STAR",
				"STAR",
				"STAR",
				"kallisto",
				"kallisto"
			))

inf_DEG = paste(results_folder, "/04_DEG/", test_list, ".csv", sep="")

inf_Trait = ref_clinical_data

dir.create( paste(results_folder, "/05_WGCNA_step01_datainput", sep=""))
outf_Expr   = paste(results_folder, "/05_WGCNA_step01_datainput/", "datExpr_", test_list, ".csv", sep="")
outf_Traits = paste(results_folder, "/05_WGCNA_step01_datainput/", "datTraits_", test_list, ".csv", sep="")



datainput <- function(STAR_or_kallisto, inf_counts, inf_deg, inf_traits, outf_expr, outf_traits)
{

	if(STAR_or_kallisto == "STAR")
	{
	# for STAR, featureCounts, edgeR, DESeq2 pipeline
	# (featureCounts data)
		datExpr0 = read.csv(inf_counts)
		rownames(datExpr0)=datExpr0[,1]
		datExpr0 = datExpr0[, -c(1:6)]
		datExpr0 = as.data.frame(t(datExpr0))
	}
	if(STAR_or_kallisto == "kallisto")
	{
	# for kallsito-sleuth pipeline
	# (sleuth_to_matrix data)
		datExpr0 = read.csv(inf_counts)
		rownames(datExpr0)=datExpr0[,1]
		datExpr0 = datExpr0[, -c(1)]
		datExpr0 = as.data.frame(t(datExpr0))
	}

	#exclude normal sample data
	datExpr0 = datExpr0[-c(43:47), ]


	deg = read.csv(inf_deg)
	bool = is.element(names(datExpr0), deg[,1])
	datExpr0 = datExpr0[, bool]


#WGCNA goodSample, goodGene selection
	gsg = goodSamplesGenes(datExpr0, verbose = 3)
	if (!gsg$allOK)
	{
		if (sum(!gsg$goodGenes)>0)
			printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
		if (sum(!gsg$goodSamples)>0)
			printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));

		datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
	}

	datTraits0 = read.csv(inf_traits)
	rownames(datTraits0) = datTraits0[,1]

	#exclude not-using data in datTraits
	datTraits0 = datTraits0[,-c(1, 2, 5, 12:17)]
	datTraits0$FU.status = as.numeric(datTraits0$FU.status == "DOD")

	 if (!gsg$allOK)
        {
                datTraits0 = datTraits0[gsg$goodSamples, ]
        }


	write.csv(datExpr0, outf_expr)
	write.csv(datTraits0, outf_traits)
}



for(i in 1:length(test_list))
{
	datainput(	STAR_or_kallisto = which_STAR_or_kallisto[i], 
			inf_counts  = inf_Counts[i],
			inf_deg     = inf_DEG[i],
			inf_traits  = inf_Trait,
			outf_expr   = outf_Expr[i],
			outf_traits = outf_Traits[i]
		)
}
