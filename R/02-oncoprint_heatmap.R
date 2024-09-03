#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(time_point)) %>%
	   dplyr::left_join(clinical, by = "patient_id") %>%
	   dplyr::rename(Tumor_Sample_Barcode = sample_id) %>%
	   dplyr::mutate(time_point = factor(time_point, levels = c("pre-treatment", "on-treatment", "post-treatment", "follow-up recurrence"), ordered = TRUE))

readr::write_tsv(manifest, path = "../res/clinical.txt", append = FALSE, col_names = TRUE)

maf = maftools::read.maf(maf = url_mutation_smry_ft,
			 clinicalData = "../res/clinical.txt",
			 verbose = TRUE)

maf@data = maf@data %>%
	   dplyr::mutate(is_hotspot = ifelse(is.na(is_hotspot), FALSE, is_hotspot)) %>%
	   dplyr::mutate(is_loh = ifelse(is.na(is_loh), FALSE, is_loh)) %>%
	   dplyr::mutate(Variant_Classification = ifelse(is_hotspot, "Hotspot", as.character(Variant_Classification))) %>%
	   dplyr::mutate(Hugo_Symbol = paste0(Hugo_Symbol, "\n", gsub("p.", "", HGVSp_Short, fixed = TRUE)))


################################################################################################################################
## BC01
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC01") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC01-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC01__.pdf", width = 4/3 * nrow(manifest_ft), height = 2.25/2 * 2)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC01-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nS241F", "IRS2\nL806V"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5)
dev.off()

################################################################################################################################
## BC02
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC02") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC02-|BZ340691", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode, perl = TRUE)
pdf(file = "../res/__BC02__.pdf", width = 4/3 * nrow(manifest_ft), height = 2.25/2 * 2)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC02-|BZ340691", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]], perl = TRUE),
	 genes = c("CDKN1B\nP137Rfs*8", "NOTCH4\nL16del"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5)
dev.off()

################################################################################################################################
## BC03
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC03") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC03-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC03__.pdf", width = 4/3 * nrow(manifest_ft) * 1.4 / 1.1, height = 2.25/2 * 3 / 1.1)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC03-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nY126*", "NOTCH3\nR2031C", "RAD51\nL91F"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5)
dev.off()

################################################################################################################################
## BC04
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC04") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC04-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC04__.pdf", width = 4/3 * nrow(manifest_ft) * 1.4 / 1.32, height = 2.25/2 * 10 / 1.45)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC04-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nR280G", "PIK3CA\nH1047R", "NF1\nL899Cfs*3", "PIK3CA\nV125E", "KRAS\nK147N", "TP53\nL194F", "PTPRS\nV1300D", "TP53\nH179R", "RICTOR\nF1072L", "DOCK8\nX1385_splice"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC06
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC06") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC06-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC06__.pdf", width = 4/3 * nrow(manifest_ft) * 1.55 / 1.25, height = 2.25/2 * 6 / 1.35)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC06-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nN239Kfs*25", "RUNX1\nR346Pfs*254", "ERBB3\nD297Y", "SF3B1\nK213*", "IL7R\nK120N", "TET2\nE1973K"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC07
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC07") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC07-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "LE984537-T", replacement = "-T3", x = maf_ft@data$Tumor_Sample_Barcode)
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "UO030693-T", replacement = "-T4", x = maf_ft@data$Tumor_Sample_Barcode)
manifest_ft$Tumor_Sample_Barcode = gsub(pattern = "BC07-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
manifest_ft$Tumor_Sample_Barcode = gsub(pattern = "LE984537-T", replacement = "-T3", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
manifest_ft$Tumor_Sample_Barcode = gsub(pattern = "UO030693-T", replacement = "-T4", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
pdf(file = "../res/__BC07__.pdf", width = 4/3 * nrow(manifest_ft) * .835, height = 2.25/2 * 2)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = (manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nR282W", "NTRK3\nN714K"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC08
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC08") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC08-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC08__.pdf", width = 4/3 * nrow(manifest_ft), height = 2.25/2 * 2)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC08-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("ARID1A\nN935Kfs*71", "GATA3\nP409Afs*99"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC09
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC09") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC09-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC09__.pdf", width = 4/3 * nrow(manifest_ft) * 1.15 / 1.15, height = 2.25/2 * 3 / 1.19)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC09-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("PIK3CA\nH1047R", "TP53\nP278R", "PRDM1\nV811L"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC10
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC10") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC10-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "WA450553-T", replacement = "-T3", x = maf_ft@data$Tumor_Sample_Barcode)
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "OI612646-T", replacement = "-T4", x = maf_ft@data$Tumor_Sample_Barcode)
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "YW750609-T", replacement = "-T5", x = maf_ft@data$Tumor_Sample_Barcode)
manifest_ft$Tumor_Sample_Barcode = gsub(pattern = "BC10-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
manifest_ft$Tumor_Sample_Barcode = gsub(pattern = "WA450553-T", replacement = "-T3", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
manifest_ft$Tumor_Sample_Barcode = gsub(pattern = "OI612646-T", replacement = "-T4", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
manifest_ft$Tumor_Sample_Barcode = gsub(pattern = "YW750609-T", replacement = "-T5", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
pdf(file = "../res/__BC10__.pdf", width = 4/3 * nrow(manifest_ft) * 1.17 / 1.65, height = 2.25/2 * 9 / 1.75)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = (manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nR248Q", "ERBB3\nG854V", "ARID1B\nS52L", "ERCC5\nP909L", "TP53\nR306*", "KMT2D\nQ3750*", "NOTCH4\nP938A", "MYCN\nP196T", "RB1\nD507A"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC12
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC12") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC12-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC12__.pdf", width = 4/3 * nrow(manifest_ft) * 1.55 / 1.35, height = 2.25/2 * 5 / 1.6)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC12-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nX307_splice", "EPHB1\nT21M", "GNAS\nP640T", "MED12\nM959I", "PTCH1\nL575F"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC13
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC13") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC13-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC13__.pdf", width = 4/3 * nrow(manifest_ft) * 1.45 / 1.20, height = 2.25/2 * 4 / 1.30)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC13-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nR273C", "NBN\nK236*", "CASP8\nD514N", "TP63\nA273S"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC15
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC15") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC15-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC15__.pdf", width = 4/3 * nrow(manifest_ft) * 1.25 / 1.01, height = 2.25/2 * 2 / 1)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC15-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("FAT1\nN81Tfs*16", "TET1\nR1753H"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC16
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC16") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC16-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC16__.pdf", width = 4/3 * nrow(manifest_ft) * 1.45 / 1.20, height = 2.25/2 * 4 / 1.30)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC16-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nR342Efs*3", "GNAQ\nY101*", "GNAS\nD448A", "MED12\nQ2100_Q2101ins*"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC17
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC17") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC17-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC17__.pdf", width = 4/3 * nrow(manifest_ft) * 1.45 / 1.20, height = 2.25/2 * 4 / 1.30)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC17-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("ERBB2\nV777L", "NBN\nR28K", "TP53\nR282G", "MITF\nN454K"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC18
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC18") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC18-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC18__.pdf", width = 4/3 * nrow(manifest_ft) * 1.45 / 1.20, height = 2.25/2 * 4 / 1.30)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC18-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("TP53\nP250R", "GATA3\nT419Hfs*89", "GATA3\nR353K", "GRIN2A\nL845M"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC20
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC20") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC20-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC20__.pdf", width = 4/3 * nrow(manifest_ft) * 1.4 / 1.35, height = 2.25/2 * 6 / 1.35)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC20-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("ARID1A\nQ1512Rfs*15", "GATA3\nL355*", "CTCF\nL366S", "KDM6A\nL1288F", "PIK3R1\nE443_A444del", "INSR\nR410W"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

################################################################################################################################
## BC21
################################################################################################################################
manifest_ft = manifest %>%
	      dplyr::filter(patient_id == "BC21") %>%
	      dplyr::filter(is_tumor_yes_no == "Yes") %>%
	      dplyr::arrange(time_point)
maf_ft = subsetMaf(maf = maf,
	   	   tsb = manifest_ft %>% .[["Tumor_Sample_Barcode"]])
maf_ft@data$Tumor_Sample_Barcode = gsub(pattern = "BC21-", replacement = "", x = maf_ft@data$Tumor_Sample_Barcode)
pdf(file = "../res/__BC21__.pdf", width = 4/3 * nrow(manifest_ft) * 1.4 / 1.1, height = 2.25/2 * 3 / 1.1)
oncoplot(maf = maf_ft,
	 minMut = 0,
 	 drawRowBar = FALSE,
 	 drawColBar = FALSE,
	 showTumorSampleBarcodes = TRUE,
	 barcode_mar = 4,
	 barcodeSrt = 0,
	 sampleOrder = gsub(pattern = "BC21-", replacement = "", x = manifest_ft %>% .[["Tumor_Sample_Barcode"]]),
	 genes = c("GATA3\nR331Lfs*24", "MDM4\nS13C", "RPS6KA4\nF519L"),
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 colors = maftools::vcColors(),
 	 showTitle = FALSE,
 	 showPct = FALSE,
 	 legend_height = 0,
 	 sepwd_genes = 1.5,
 	 sepwd_samples = 1.5,
 	 additionalFeature = c("is_loh", TRUE),
 	 additionalFeaturePch = 3,
 	 additionalFeatureCol = "white",
 	 additionalFeatureCex = 1.5,
	 removeNonMutated = FALSE)
dev.off()

pdf(file = "../res/__Legend__.pdf", width = 2, height = 2.5)
mafLegend()
dev.off()
