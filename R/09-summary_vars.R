#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(stage_at_diagnosis = gsub(pattern = "A|B|C", replacement = "", x = stage_at_diagnosis, perl = TRUE)) %>%
	   dplyr::mutate(stage_at_diagnosis = factor(stage_at_diagnosis, levels = c("II", "III"), ordered = TRUE)) %>%
	   dplyr::mutate(tumor_size_pre_CT = gsub(pattern = " cm", replacement = "", x = tumor_size_pre_CT, fixed = TRUE)) %>%
	   dplyr::mutate(tumor_size = case_when(
		tumor_size_pre_CT < 3 ~ "<3 cm",
		tumor_size_pre_CT >= 3 ~ ">=3 cm",
	   )) %>%
	  dplyr::mutate(tumor_size = factor(tumor_size, levels = c("<3 cm", ">=3 cm"), ordered = TRUE)) %>%
	  dplyr::mutate(node_stage_pre_CT = factor(node_stage_pre_CT, levels = c("N1", "N2", "N3"), ordered = TRUE)) %>%
	  dplyr::mutate(ER_status_pre_CT = factor(ER_status_pre_CT, levels = c("Positive", "Negative"), ordered = TRUE)) %>%
	  dplyr::mutate(PR_status_pre_CT = factor(PR_status_pre_CT, levels = c("Positive", "Negative"), ordered = TRUE)) %>%
	  dplyr::mutate(HER2_status_pre_CT = factor(HER2_status_pre_CT, levels = c("Positive", "Negative"), ordered = TRUE)) %>%
	  dplyr::mutate(pCR_status_yes_no = factor(pCR_status_yes_no, levels = c("Yes", "No"), ordered = TRUE)) %>%
	  dplyr::mutate(any_recurrence_yes_no = factor(pCR_status_yes_no, levels = c("Yes", "No"), ordered = TRUE)) %>%
	  dplyr::mutate(BC_subtype = case_when(
		HER2_status_pre_CT == "Positive" ~ "HER2+",
		HER2_status_pre_CT == "Negative" & ER_status_pre_CT == "Negative" & PR_status_pre_CT == "Negative" ~ "TN",
		TRUE ~ "ER+"
	  )) %>%
	  dplyr::mutate(BC_subtype = factor(BC_subtype, levels = c("HER2+", "TN", "ER+"), ordered = TRUE)) %>%
	  dplyr::mutate(histological_grade = factor(histological_grade, levels = c("2", "3"), ordered = TRUE)) %>%
	  dplyr::mutate(age_at_diagnosis = case_when(
		age_at_diagnosis < 50 ~ "<50 yrs",
		age_at_diagnosis >= 50 ~ ">=50 yrs"
	  )) %>%
	  dplyr::mutate(age_at_diagnosis = factor(age_at_diagnosis, levels = c("<50 yrs", ">=50 yrs"), ordered = TRUE)) %>%
	  readr::type_convert()

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(time_point)) %>%
	   dplyr::rename(Tumor_Sample_Barcode = sample_id) %>%
	   dplyr::mutate(time_point = case_when(
		time_point == "pre-treatment" ~ "Baseline",
		time_point == "on-treatment" ~ "On-treatment",
		time_point == "post-treatment" ~ "Post-treatment",
		time_point == "follow-up recurrence" ~ "Local/Distant Relapse",
	   )) %>%
	   dplyr::mutate(time_point = factor(time_point, levels = c("Baseline", "On-treatment", "Post-treatment", "Local/Distant Relapse"), ordered = TRUE))

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
	 dplyr::filter(!(patient_id == "BC21" & Variant_Classification == "In_Frame_Del"))


any_vars_cfdna = maf_ft %>%
		 dplyr::filter(is_plasma_yes_no == "Yes") %>%
		 dplyr::group_by(patient_id) %>%
		 dplyr::summarize(n_variants_in_plasma = n(),
				  mean_af = mean(af, na.rm=TRUE),
				  max_af = max(af, na.rm=TRUE)) %>%
		 dplyr::left_join(maf_ft %>%
				  dplyr::filter(is_plasma_yes_no == "No") %>%
				  dplyr::filter(time_point == "Baseline") %>%
				  dplyr::group_by(patient_id) %>%
				  dplyr::summarize(n = n()), by = "patient_id") %>%
		 dplyr::mutate(`%_variants_in_plasma` = 100*n_variants_in_plasma/n) %>%
		 dplyr::mutate(time_point = "Baseline") %>%
		 dplyr::select(-n)

any_vars_ddpcr = readr::read_tsv(file = url_ddpcr, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		 dplyr::filter(`is_tumor?` == "T") %>%
		 readr::type_convert() %>%
		 dplyr::mutate(af = 100 * droplet_count_mutant / (droplet_count_mutant + droplet_count_wt)) %>%
		 dplyr::mutate(af = ifelse(is.na(af), 0, af)) %>%
		 dplyr::filter(!(patient_id == "BC18" & Hugo_Symbol == "GATA3")) %>%
		 dplyr::mutate(time_point = case_when(
			grepl("_1", sample_type) ~ "Baseline",
			grepl("_2", sample_type) ~ "On-treatment",
			grepl("_3", sample_type) ~ "Post-treatment"
		 )) %>%
		 dplyr::group_by(patient_id, time_point) %>%
		 dplyr::summarize(`ctDNA_fraction_ddPCR_%` = max(af, na.rm = TRUE)) %>%
		 dplyr::ungroup() %>%
		 dplyr::mutate(ctDNA_fraction_ddPCR_yes_no = ifelse(`ctDNA_fraction_ddPCR_%` == 0, 0, 1))

manifest %>%
dplyr::filter(time_point != "Local/Distant Relapse") %>%
dplyr::left_join(any_vars_cfdna, by = c("patient_id", "time_point")) %>%
dplyr::left_join(any_vars_ddpcr, by = c("patient_id", "time_point")) %>%
dplyr::left_join(clinical, by = "patient_id") %>%
dplyr::select(-Tumor_Sample_Barcode, -is_tumor_yes_no, -is_plasma_yes_no, -cmo_or_dmp_id) %>%
readr::write_tsv(path = "../res/summary_ddpcr_impact.txt", append = FALSE, col_names = TRUE)
