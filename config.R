seq="../sequence/SRP067524"

seq_list = c("SRR3081601", "SRR3081602", "SRR3081603", "SRR3081604",
"SRR3081605", "SRR3081606", "SRR3081607", "SRR3081608", "SRR3081609",
"SRR3081610", "SRR3081611", "SRR3081612", "SRR3081613", "SRR3081614",
"SRR3081615", "SRR3081616", "SRR3081617", "SRR3081618", "SRR3081619",
"SRR3081620", "SRR3081621", "SRR3081622", "SRR3081623", "SRR3081624",
"SRR3081625", "SRR3081626", "SRR3081627", "SRR3081628", "SRR3081629",
"SRR3081630", "SRR3081631", "SRR3081632", "SRR3081633", "SRR3081634",
"SRR3081635", "SRR3081636", "SRR3081637", "SRR3081638", "SRR3081639",
"SRR3081640", "SRR3081641", "SRR3081642", "SRR3081643", "SRR3081644",
"SRR3081645", "SRR3081646", "SRR3081647")

miRNA_seq="../sequence/SRP154282"

miRNA_seq_list=c("SRR7535281", "SRR7535282", "SRR7535283",
"SRR7535284", "SRR7535285", "SRR7535286")

miRNA_group = c("Tumor", "normal", "normal", "Tumor", "normal", "Tumor")

pM_related_module = c(
	"brown", "brown", "turquoise", "blue", "blue",
	"blue", "blue", "blue")


ref_genome_fa="../reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ref_cdna_fa="../reference/Homo_sapiens.GRCh38.cdna.all.fa"
ref_gtf="../reference/Homo_sapiens.GRCh38.97.gtf"

ref_kallisto_index="../reference/kallisto-index-Ens97"
ref_STAR_index_100="../reference/star-index_100"
ref_STAR_index_50="../reference/star-index_50"

ref_clinical_data = "../data/clinical_data_modified.csv"

ref_CMap_Landmark_gene = "../data/CMap_Landmark_gene.txt"
ref_CMap_BING_gene = "../data/CMap_BING_gene.txt"

ref_Targetscan = "../data/Predicted_Targets_Context_Scores.default_predictions.txt"
ref_miRDB      = "../data/miRDB_v6.0_prediction_result.txt"
ref_miRTarbase = "../data/hsa_MTI.xlsx"

test_list = c(	"test01",
		"test02",
		"test03",
		"test04",
		"test05",
		"test06",
		"test07",
		"test08")

results_folder="./results"

results_kallisto="/01_kallsito_results"
results_pizzly="/01_pizzly_results"
results_STAR_100="/01_STAR_results01_100"
results_STAR_50="/01_STAR_results02_50"

results_featureCounts_STAR_100="/02_featureCounts_STAR_100"
results_featureCounts_STAR_50="/02_featureCounts_STAR_50"

results_merged_featureCounts_STAR_100 = "/03_merged_featureCounts_STAR_100.csv"
results_merged_featureCounts_STAR_50 = "/03_merged_featureCounts_STAR_50.csv"

results_miRNA = "/19_02_miRNA_mapping"
results_miRNA_merged = "miRNA_merged_counts.csv"

results_DEmiRNA      = "/19_04_DEmiRNA"
results_DEmiRNA_up   = "DEmiRNA_up.csv"
results_DEmiRNA_down = "DEmiRNA_down.csv"

results_miRNA_predict_target = "/20_miRNA_predict_target.csv"

results_detect_hubgene = "/21_detect_hubgene"
