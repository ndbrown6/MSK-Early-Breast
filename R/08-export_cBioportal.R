#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

# clinical
manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(is_tumor_yes_no == "Yes") %>%
	   dplyr::select(patient_id, sample_id) %>%
	   dplyr::mutate(sample_id = case_when(
		   sample_id == "BZ340691-T" ~ "BC02-M",
		   sample_id == "LE984537-T" ~ "BC07-M1",
		   sample_id == "UO030693-T" ~ "BC07-M2",
		   sample_id == "WA450553-T" ~ "BC10-M1",
		   sample_id == "OI612646-T" ~ "BC10-M2",
		   sample_id == "YW750609-T" ~ "BC10-M3",
		   TRUE ~ sample_id
	   ))

clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(ER_status_pre_CT = factor(ER_status_pre_CT, levels = c("Positive", "Negative"), ordered = TRUE)) %>%
	   dplyr::mutate(PR_status_pre_CT = factor(PR_status_pre_CT, levels = c("Positive", "Negative"), ordered = TRUE)) %>%
	   dplyr::mutate(HER2_status_pre_CT = factor(HER2_status_pre_CT, levels = c("Positive", "Negative"), ordered = TRUE)) %>%
	   dplyr::mutate(BC_subtype = case_when(
		HER2_status_pre_CT == "Positive" ~ "HER2+",
		HER2_status_pre_CT == "Negative" & ER_status_pre_CT == "Negative" & PR_status_pre_CT == "Negative" ~ "TN",
		TRUE ~ "ER+"
	   )) %>%
	   dplyr::mutate(BC_subtype = factor(BC_subtype, levels = c("HER2+", "TN", "ER+"), ordered = TRUE)) %>%
	   dplyr::mutate(tumor_size_pre_CT = gsub(pattern = " cm", replacement = "", x = tumor_size_pre_CT, fixed = TRUE)) %>%
	   readr::type_convert() %>%
	   dplyr::select(patient_id,
			 age_at_diagnosis,
			 tumor_size = tumor_size_pre_CT,
			 tumor_stage = tumor_stage_pre_CT,
			 histological_grade,
			 ER_status = ER_status_pre_CT,
			 PR_status = PR_status_pre_CT,
			 HER2_status = HER2_status_pre_CT,
			 pCR_status = pCR_status_yes_no,
			 RCB) %>%
	   dplyr::mutate(`IMPACT Version` = "341",
			 `Matched Normal` = TRUE)

manifest %>%
dplyr::left_join(clinical, by = "patient_id") %>%
dplyr::select(-patient_id) %>%
readr::write_tsv(path = "../res/cBioportal_clinical.txt", append = FALSE, col_names = TRUE)

# maf
manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(time_point)) %>%
	   dplyr::rename(Tumor_Sample_Barcode = sample_id) %>%
	   dplyr::mutate(time_point = case_when(
		time_point == "pre-treatment" ~ "Baseline\nTumor",
		time_point == "on-treatment" ~ "On-treat\n-ment",
		time_point == "post-treatment" ~ "Residual\nDisease",
		time_point == "follow-up recurrence" ~ "Local/ Distant\nRelapse",
	   )) %>%
	   dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nTumor", "On-treat\n-ment", "Residual\nDisease", "Local/ Distant\nRelapse"), ordered = TRUE))

maf_ft = readr::read_tsv(file = url_mutation_smry_ft, comment = "#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	 readr::type_convert() %>%
	 dplyr::left_join(manifest, by = "Tumor_Sample_Barcode") %>%
	 dplyr::filter(is_tumor_yes_no == "Yes") %>%
	 dplyr::filter(!is.na(HGVSp_Short),
		       !is.na(Hugo_Symbol),
		       Variant_Classification != "Silent") %>%
	 dplyr::mutate(af = 100 * t_alt_count / t_depth) %>%
	 ################################################################################################################################
	 ## BC02
	 ################################################################################################################################
	 dplyr::filter(!(patient_id == "BC02" & Hugo_Symbol == "NOTCH4")) %>%
	 ################################################################################################################################
	 ## BC06
	 ################################################################################################################################
	 dplyr::filter(!(patient_id == "BC06" & Hugo_Symbol == "NOTCH4")) %>%
	 ################################################################################################################################
	 ## BC16
	 ################################################################################################################################
	 dplyr::filter(!(patient_id == "BC16" & grepl("APC|ASXL|GRIN|MDC1|MITF|RNF4|SMO|TRAF", Hugo_Symbol, perl=TRUE))) %>%
	 dplyr::filter(!(patient_id == "BC16" & Hugo_Symbol == "MED1")) %>%
	 ################################################################################################################################
	 ## BC17
	 ################################################################################################################################
	 dplyr::filter(!(patient_id == "BC17" & grepl("AMER|AR|EPHA|FAT1|GRIN|MDC|POLE|RFWD|TSC1", Hugo_Symbol, perl=TRUE))) %>%
	 ################################################################################################################################
	 ## BC21
	 ################################################################################################################################
	 dplyr::filter(!(patient_id == "BC21" & Hugo_Symbol == "ARAF")) %>%
	 dplyr::filter(!(patient_id == "BC21" & Variant_Classification == "In_Frame_Del")) %>%
	 dplyr::select(`Tumor Sample` = Tumor_Sample_Barcode,
		       Chromosome,
		       Position = Start_Position,
		       `Reference Allele` = Reference_Allele,
		       `Alternate Allele` = Tumor_Seq_Allele2,
		       Gene = Hugo_Symbol,
		       `Mutation Type` = Variant_Classification,
		       `AA change` = HGVSp_Short,
		       t_ref_count,
		       t_alt_count,
		       n_ref_count,
		       n_alt_count) %>%
	 dplyr::mutate(`Tumor Sample` = case_when(
		   `Tumor Sample` == "BZ340691-T" ~ "BC02-M",
		   `Tumor Sample` == "LE984537-T" ~ "BC07-M1",
		   `Tumor Sample` == "UO030693-T" ~ "BC07-M2",
		   `Tumor Sample` == "WA450553-T" ~ "BC10-M1",
		   `Tumor Sample` == "OI612646-T" ~ "BC10-M2",
		   `Tumor Sample` == "YW750609-T" ~ "BC10-M3",
		   TRUE ~ `Tumor Sample`
	   ))

gene_list = readr::read_tsv(file = url_gene_list, col_names = FALSE, col_types = cols(.default = col_character())) %>%
	    readr::type_convert()

maf_ft %>%
dplyr::filter(Gene %in% gene_list$X1) %>%
readr::write_tsv(path = "../res/cBioportal_variants.txt", append = FALSE, col_names = TRUE)

# copy number
manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(time_point)) %>%
	   dplyr::rename(Tumor_Sample_Barcode = sample_id) %>%
	   dplyr::mutate(time_point = case_when(
		time_point == "pre-treatment" ~ "Baseline\nTumor",
		time_point == "on-treatment" ~ "On-treat\n-ment",
		time_point == "post-treatment" ~ "Residual\nDisease",
		time_point == "follow-up recurrence" ~ "Local/ Distant\nRelapse",
	   )) %>%
	   dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nTumor", "On-treat\n-ment", "Residual\nDisease", "Local/ Distant\nRelapse"), ordered = TRUE))

bygene = readr::read_tsv(file = url_by_gene, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	 readr::type_convert() %>%
	 dplyr::select(-chrom, -start, -end, -band) %>%
	 dplyr::select(Hugo_Symbol = hgnc,
		       contains("_LRR_threshold")) %>%
	 reshape2::melt(id.vars = "Hugo_Symbol", variable.name = "sample_id", value.name = "call") %>%
	 dplyr::mutate(sample_id = unlist(lapply(as.character(sample_id), function(x) {strsplit(x, split = "_", fixed = TRUE)[[1]][1]}))) %>%
	 dplyr::mutate(sample_id = case_when(
		   sample_id == "BZ340691-T" ~ "BC02-M",
		   sample_id == "LE984537-T" ~ "BC07-M1",
		   sample_id == "UO030693-T" ~ "BC07-M2",
		   sample_id == "WA450553-T" ~ "BC10-M1",
		   sample_id == "OI612646-T" ~ "BC10-M2",
		   sample_id == "YW750609-T" ~ "BC10-M3",
		   TRUE ~ sample_id
	 ))

gene_list = readr::read_tsv(file = url_gene_list, col_names = FALSE, col_types = cols(.default = col_character())) %>%
	    readr::type_convert()

bygene %>%
dplyr::filter(Hugo_Symbol %in% gene_list$X1) %>%
reshape2::dcast(Hugo_Symbol ~ sample_id, value.var = "call", fun.aggregate = function(x) { max(x) }, fill = 0) %>%
readr::write_tsv(path = "../res/cBioportal_copynumber.txt", append = FALSE, col_names = TRUE)



