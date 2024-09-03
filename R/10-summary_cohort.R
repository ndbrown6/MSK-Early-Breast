#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

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
	   readr::type_convert()
	  

# Age at disgnosis
pander::pander(clinical %>%
	       dplyr::summarize(Q1 = round(quantile(age_at_diagnosis, probs = .25)),
				Q2 = round(median(age_at_diagnosis)),
				Q3 = round(quantile(age_at_diagnosis, probs = .75))))

# Molecular subtype
pander::pander(clinical %>%
	       dplyr::group_by(BC_subtype) %>%
	       dplyr::summarize(n = n()) %>%
	       dplyr::ungroup() %>%
	       dplyr::mutate(`%` = 100*n/20))

# Tumor size
pander::pander(clinical %>%
	       dplyr::summarize(Q1 = signif(quantile(tumor_size_pre_CT, probs = .25), 2),
				Q2 = signif(median(tumor_size_pre_CT), 2),
				Q3 = signif(quantile(tumor_size_pre_CT, probs = .75), 2)))

# Lymph node status
pander::pander(clinical %>%
	       dplyr::group_by(node_stage_pre_CT) %>%
	       dplyr::summarize(n = n()) %>%
	       dplyr::ungroup() %>%
	       dplyr::mutate(`%` = 100*n/20))


# PCR status
pander::pander(clinical %>%
	       dplyr::group_by(pCR_status_yes_no) %>%
	       dplyr::summarize(n = n()) %>%
	       dplyr::ungroup() %>%
	       dplyr::mutate(`%` = 100*n/20))

# Recurrence status
pander::pander(clinical %>%
	       dplyr::group_by(any_recurrence_yes_no) %>%
	       dplyr::summarize(n = n()) %>%
	       dplyr::ungroup() %>%
	       dplyr::mutate(`%` = 100*n/20))