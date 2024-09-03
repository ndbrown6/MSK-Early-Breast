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


smry_ = manifest %>%
	dplyr::left_join(hs_metrics, by = "SAMPLE_NAME") %>%
	dplyr::filter(is_tumor_yes_no == "Yes") %>%
	dplyr::filter(is_plasma_yes_no == "No")

smry_ %>%
.[["MEAN_TARGET_COVERAGE"]] %>%
median(., na.rm = TRUE)

smry_ %>%
.[["MEAN_TARGET_COVERAGE"]] %>%
range(., na.rm = TRUE)

smry_ = manifest %>%
	dplyr::left_join(hs_metrics, by = "SAMPLE_NAME") %>%
	dplyr::filter(is_tumor_yes_no == "No") 

smry_ %>%
.[["MEAN_TARGET_COVERAGE"]] %>%
median(., na.rm = TRUE)

smry_ %>%
.[["MEAN_TARGET_COVERAGE"]] %>%
range(., na.rm = TRUE)

smry_ = manifest %>%
	dplyr::left_join(hs_metrics, by = "SAMPLE_NAME") %>%
	dplyr::filter(is_tumor_yes_no == "Yes") %>%
	dplyr::filter(is_plasma_yes_no == "Yes")

smry_ %>%
.[["MEAN_TARGET_COVERAGE"]] %>%
median(., na.rm = TRUE)

smry_ %>%
.[["MEAN_TARGET_COVERAGE"]] %>%
range(., na.rm = TRUE)