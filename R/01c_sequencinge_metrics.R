#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::rename(SAMPLE_NAME = sample_id)
	   
gc_metrics = readr::read_tsv(file = url_gc_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert()

idx_metrics = readr::read_tsv(file = url_idx_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert()

insert_metrics = readr::read_tsv(file = url_insert_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		 readr::type_convert()

hs_metrics = readr::read_tsv(file = url_hs_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert()


manifest %>%
dplyr::left_join(hs_metrics, by = "SAMPLE_NAME") %>%
dplyr::filter(is_tumor_yes_no == "Yes") %>%
dplyr::filter(is_plasma_yes_no == "No") %>%
dplyr::summarize(min_target_converage = min(MEAN_TARGET_COVERAGE, na.rm = TRUE),
		 median_target_converage = median(MEAN_TARGET_COVERAGE, na.rm = TRUE),
		 max_target_converage = max(MEAN_TARGET_COVERAGE, na.rm = TRUE)) %>%
pander::pander(caption = "Tumor sample sequencing metrics")


manifest %>%
dplyr::left_join(hs_metrics, by = "SAMPLE_NAME") %>%
dplyr::filter(is_tumor_yes_no == "No") %>%
dplyr::summarize(min_target_converage = min(MEAN_TARGET_COVERAGE, na.rm = TRUE),
		 median_target_converage = median(MEAN_TARGET_COVERAGE, na.rm = TRUE),
		 max_target_converage = max(MEAN_TARGET_COVERAGE, na.rm = TRUE)) %>%
pander::pander(caption = "Normal sample sequencing metrics")


manifest %>%
dplyr::left_join(hs_metrics, by = "SAMPLE_NAME") %>%
dplyr::filter(is_tumor_yes_no == "Yes") %>%
dplyr::filter(is_plasma_yes_no == "Yes") %>%
dplyr::summarize(min_target_converage = min(MEAN_TARGET_COVERAGE, na.rm = TRUE),
		 median_target_converage = median(MEAN_TARGET_COVERAGE, na.rm = TRUE),
		 max_target_converage = max(MEAN_TARGET_COVERAGE, na.rm = TRUE)) %>%
pander::pander(caption = "Plasma sample sequencing metrics")
