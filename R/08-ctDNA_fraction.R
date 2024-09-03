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
	   dplyr::rename(tumor_id = sample_id) %>%
	   dplyr::mutate(time_point = case_when(
		time_point == "pre-treatment" ~ "Baseline\nTumor",
		time_point == "on-treatment" ~ "On-treat\n-ment",
		time_point == "post-treatment" ~ "Residual\nDisease",
		time_point == "follow-up recurrence" ~ "Local/ Distant\nRelapse",
	   )) %>%
	   dplyr::left_join(readr::read_tsv(file = url_purity_ploidy, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     		    readr::type_convert(), by = "tumor_id") %>%
	   dplyr::mutate(purity = ifelse(is.na(purity), 0, purity)) %>%
	   dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nTumor", "On-treat\n-ment", "Residual\nDisease", "Local/ Distant\nRelapse"), ordered = TRUE)) %>%
	   dplyr::arrange(time_point)

plot_ = manifest %>%
	dplyr::filter(is_plasma_yes_no == "Yes") %>%
	dplyr::full_join(clinical, by = "patient_id") %>%
	dplyr::filter(!is.na(purity)) %>%
	reshape2::melt(id.vars = c("patient_id", "purity"),
		       measure.vars = c("tumor_size", "stage_at_diagnosis", "ER_status_pre_CT", "PR_status_pre_CT",
					"HER2_status_pre_CT", "histological_grade", "age_at_diagnosis", "pCR_status_yes_no")) %>%
	dplyr::mutate(variable = case_when(
		variable == "tumor_size" ~ "Tumor size",
		variable == "stage_at_diagnosis" ~ "Stage at diagnosis",
		variable == "ER_status_pre_CT" ~ "ER status",
		variable == "PR_status_pre_CT" ~ "PR status",
		variable == "HER2_status_pre_CT" ~ "HER2 status",
		variable == "histological_grade" ~ "Histological grade",
		variable == "age_at_diagnosis" ~ "Age at diagnosis",
		variable == "pCR_status_yes_no" ~ "pCR status"
	)) %>%
	dplyr::mutate(variable = factor(variable, levels = c("Age at diagnosis", "Stage at diagnosis", "Tumor size", "Histological grade",
							     "ER status", "PR status", "HER2 status", "pCR status"), ordered = TRUE)) %>%
	ggplot(aes(x = value, y = 100*purity)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 100)) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("") +
	ylab(expression("ctDNA fraction (%)")) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Age at diagnosis"),
		    comparisons = list(c("<50 yrs", ">=50 yrs")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "ER status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "HER2 status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "PR status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Histological grade"),
		    comparisons = list(c("2", "3")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "pCR status"),
		    comparisons = list(c("No", "Yes")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Stage at diagnosis"),
		    comparisons = list(c("II", "III")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Tumor size"),
		    comparisons = list(c("<3 cm", ">=3 cm")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 87.5,
		    tip_length = 0.01) +
	facet_wrap(~variable, scales = "free_x", nrow = 2, ncol = 4)

pdf(file = "../res/ctDNA_Fraction.pdf", width = 9, height = 7)
print(plot_)
dev.off()