#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

preanalytical_conditions = readr::read_tsv(file = url_preanalytical_conditions, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			   readr::type_convert()

smry_ = preanalytical_conditions %>%
	dplyr::group_by(patient_id, sample_id) %>%
	dplyr::summarize(concentration_pg_muL = mean(concentration_pg_muL),
			 cfDNA_total_yield_ng = mean(cfDNA_total_yield_ng)) %>%
	dplyr::ungroup() %>%
	dplyr::left_join(clinical, by = "patient_id") %>%
        dplyr::mutate(time_point = case_when(
		grepl("_1", sample_id) ~ "Baseline",
		grepl("_2", sample_id) ~ "On-treat\n-ment",
		grepl("_3", sample_id) ~ "Post-treat\n-ment"
	)) %>%
	dplyr::mutate(stage_at_diagnosis = gsub(pattern = "A|B|C", replacement = "", x = stage_at_diagnosis, perl = TRUE)) %>%
	dplyr::mutate(time_point = factor(time_point, levels = c("Baseline", "On-treat\n-ment", "Post-treat\n-ment"), ordered = TRUE)) %>%
	dplyr::mutate(stage_at_diagnosis = factor(stage_at_diagnosis, levels = c("II", "III"), ordered = TRUE)) %>%
	dplyr::mutate(tumor_size_pre_CT = gsub(pattern = " cm", replacement = "", x = tumor_size_pre_CT, fixed = TRUE)) %>%
	readr::type_convert() %>%
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

ddPCR = readr::read_tsv(file = url_ddpcr, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	dplyr::filter(`is_tumor?` == "T") %>%
	readr::type_convert() %>%
	dplyr::mutate(af = droplet_count_mutant / (droplet_count_mutant + droplet_count_wt)) %>%
	dplyr::mutate(af = ifelse(is.na(af), 0, af)) %>%
	dplyr::filter(!(patient_id == "BC18" & Hugo_Symbol == "GATA3")) %>%
	dplyr::mutate(time_point = case_when(
			grepl("_1", sample_type) ~ "Baseline",
			grepl("_2", sample_type) ~ "On-treat\n-ment",
			grepl("_3", sample_type) ~ "Post-treat\n-ment"
	)) %>%
	dplyr::group_by(patient_id, time_point) %>%
	dplyr::summarize(max_af = max(af),
			 mean_af = mean(af)) %>%
	dplyr::ungroup()

smry_ = smry_ %>%
	dplyr::left_join(ddPCR, by = c("patient_id", "time_point")) %>%
	dplyr::mutate(Max_ctDNA_concentration_pg_muL = concentration_pg_muL * max_af) %>%
	dplyr::mutate(Mean_ctDNA_concentration_pg_muL = concentration_pg_muL * mean_af)

################################################################################################################################
## Baseline max AF cfDNA concentration
################################################################################################################################
plot_ = smry_ %>%
	dplyr::mutate(molecular_subtype = case_when(
				 ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
				 ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
				 HER2_status_pre_CT=="Positive" ~ "HER2+"
	)) %>%
	dplyr::filter(time_point == "Baseline") %>%
	reshape2::melt(id.vars = c("sample_id", "Max_ctDNA_concentration_pg_muL"),
		       measure.vars = c("tumor_size", "stage_at_diagnosis", "ER_status_pre_CT", "PR_status_pre_CT",
					"HER2_status_pre_CT", "histological_grade", "age_at_diagnosis", "pCR_status_yes_no", "molecular_subtype")) %>%
	dplyr::mutate(variable = case_when(
		variable == "tumor_size" ~ "Tumor size",
		variable == "stage_at_diagnosis" ~ "Stage at diagnosis",
		variable == "molecular_subtype" ~ "Molecular subtype",
		variable == "ER_status_pre_CT" ~ "ER status",
		variable == "PR_status_pre_CT" ~ "PR status",
		variable == "HER2_status_pre_CT" ~ "HER2 status",
		variable == "histological_grade" ~ "Histological grade",
		variable == "age_at_diagnosis" ~ "Age at diagnosis",
		variable == "pCR_status_yes_no" ~ "pCR status"
	)) %>%
	dplyr::filter(!(variable %in% c("ER status", "PR status", "HER2 status"))) %>%
	dplyr::mutate(variable = factor(variable, levels = c("Age at diagnosis", "Stage at diagnosis", "Tumor size", "Histological grade",
							     "Molecular subtype", "pCR status"), ordered = TRUE)) %>%
	ggplot(aes(x = value, y = Max_ctDNA_concentration_pg_muL)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 25),
			   labels = scientific_10) +
	xlab("") +
	ylab(expression("ctDNA (pg"%.%mu~"L"^-1~")")) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Age at diagnosis"),
		    comparisons = list(c("<50 yrs", ">=50 yrs")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 20,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Molecular subtype"),
		    comparisons = list(c("HR+", "HER2+"),
				       c("HER2+", "TN"),
				       c("HR+", "TN")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = c(20, 22.5, 25)-2,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Histological grade"),
		    comparisons = list(c("2", "3")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 20,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "pCR status"),
		    comparisons = list(c("No", "Yes")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 20,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Stage at diagnosis"),
		    comparisons = list(c("II", "III")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 20,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Tumor size"),
		    comparisons = list(c("<3 cm", ">=3 cm")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 20,
		    tip_length = 0.01) +
	facet_wrap(~variable, scales = "free_x", nrow = 2, ncol = 3) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      strip.background = element_blank())

pdf(file = "../res/Baseline_cfDNA_Concentration_by_AF.pdf", width = 9, height = 7)
print(plot_)
dev.off()

################################################################################################################################
## On-treatment max AF cfDNA concentration
################################################################################################################################
plot_ = smry_ %>%
	dplyr::mutate(molecular_subtype = case_when(
				 ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
				 ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
				 HER2_status_pre_CT=="Positive" ~ "HER2+"
	)) %>%
	dplyr::filter(time_point == "On-treat\n-ment") %>%
	reshape2::melt(id.vars = c("sample_id", "Max_ctDNA_concentration_pg_muL"),
		       measure.vars = c("tumor_size", "stage_at_diagnosis", "ER_status_pre_CT", "PR_status_pre_CT",
					"HER2_status_pre_CT", "histological_grade", "age_at_diagnosis", "pCR_status_yes_no", "molecular_subtype")) %>%
	dplyr::mutate(variable = case_when(
		variable == "tumor_size" ~ "Tumor size",
		variable == "stage_at_diagnosis" ~ "Stage at diagnosis",
		variable == "molecular_subtype" ~ "Molecular subtype",
		variable == "ER_status_pre_CT" ~ "ER status",
		variable == "PR_status_pre_CT" ~ "PR status",
		variable == "HER2_status_pre_CT" ~ "HER2 status",
		variable == "histological_grade" ~ "Histological grade",
		variable == "age_at_diagnosis" ~ "Age at diagnosis",
		variable == "pCR_status_yes_no" ~ "pCR status"
	)) %>%
	dplyr::filter(!(variable %in% c("ER status", "PR status", "HER2 status"))) %>%
	dplyr::mutate(variable = factor(variable, levels = c("Age at diagnosis", "Stage at diagnosis", "Tumor size", "Histological grade",
							     "Molecular subtype", "pCR status"), ordered = TRUE)) %>%
	ggplot(aes(x = value, y = Max_ctDNA_concentration_pg_muL)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 5),
			   labels = scientific_10) +
	xlab("") +
	ylab(expression("ctDNA (pg"%.%mu~"L"^-1~")")) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Age at diagnosis"),
		    comparisons = list(c("<50 yrs", ">=50 yrs")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Molecular subtype"),
		    comparisons = list(c("HR+", "HER2+"),
				       c("HER2+", "TN"),
				       c("HR+", "TN")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = c(3, 3.5, 4),
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Histological grade"),
		    comparisons = list(c("2", "3")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "pCR status"),
		    comparisons = list(c("No", "Yes")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Stage at diagnosis"),
		    comparisons = list(c("II", "III")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Tumor size"),
		    comparisons = list(c("<3 cm", ">=3 cm")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	facet_wrap(~variable, scales = "free_x", nrow = 2, ncol = 3) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      strip.background = element_blank())

pdf(file = "../res/On-treatment_cfDNA_Concentration_by_AF.pdf", width = 9, height = 7)
print(plot_)
dev.off()

################################################################################################################################
## Post-treatment max AF cfDNA concentration
################################################################################################################################
plot_ = smry_ %>%
	dplyr::mutate(molecular_subtype = case_when(
				 ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
				 ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
				 HER2_status_pre_CT=="Positive" ~ "HER2+"
	)) %>%
	dplyr::filter(time_point == "Post-treat\n-ment") %>%
	reshape2::melt(id.vars = c("sample_id", "Max_ctDNA_concentration_pg_muL"),
		       measure.vars = c("tumor_size", "stage_at_diagnosis", "ER_status_pre_CT", "PR_status_pre_CT",
					"HER2_status_pre_CT", "histological_grade", "age_at_diagnosis", "pCR_status_yes_no", "molecular_subtype")) %>%
	dplyr::mutate(variable = case_when(
		variable == "tumor_size" ~ "Tumor size",
		variable == "stage_at_diagnosis" ~ "Stage at diagnosis",
		variable == "molecular_subtype" ~ "Molecular subtype",
		variable == "ER_status_pre_CT" ~ "ER status",
		variable == "PR_status_pre_CT" ~ "PR status",
		variable == "HER2_status_pre_CT" ~ "HER2 status",
		variable == "histological_grade" ~ "Histological grade",
		variable == "age_at_diagnosis" ~ "Age at diagnosis",
		variable == "pCR_status_yes_no" ~ "pCR status"
	)) %>%
	dplyr::filter(!(variable %in% c("ER status", "PR status", "HER2 status"))) %>%
	dplyr::mutate(variable = factor(variable, levels = c("Age at diagnosis", "Stage at diagnosis", "Tumor size", "Histological grade",
							     "Molecular subtype", "pCR status"), ordered = TRUE)) %>%
	ggplot(aes(x = value, y = Max_ctDNA_concentration_pg_muL)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 5),
			   labels = scientific_10) +
	xlab("") +
	ylab(expression("ctDNA (pg"%.%mu~"L"^-1~")")) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Age at diagnosis"),
		    comparisons = list(c("<50 yrs", ">=50 yrs")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Molecular subtype"),
		    comparisons = list(c("HR+", "HER2+"),
				       c("HER2+", "TN"),
				       c("HR+", "TN")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = c(3, 3.5, 4),
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Histological grade"),
		    comparisons = list(c("2", "3")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "pCR status"),
		    comparisons = list(c("No", "Yes")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Stage at diagnosis"),
		    comparisons = list(c("II", "III")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Tumor size"),
		    comparisons = list(c("<3 cm", ">=3 cm")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 4,
		    tip_length = 0.01) +
	facet_wrap(~variable, scales = "free_x", nrow = 2, ncol = 3) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      strip.background = element_blank())

pdf(file = "../res/Post-treatment_cfDNA_Concentration_by_AF.pdf", width = 9, height = 7)
print(plot_)
dev.off()