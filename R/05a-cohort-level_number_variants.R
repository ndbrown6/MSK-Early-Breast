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
	 dplyr::filter(!(patient_id == "BC21" & Variant_Classification == "In_Frame_Del"))



################################################################################################################################
## Any vars in plasma
################################################################################################################################
any_vars_cfdna = maf_ft %>%
		 dplyr::filter(is_plasma_yes_no == "Yes") %>%
		 dplyr::group_by(patient_id) %>%
		 dplyr::summarize(n = n())

plot_ = any_vars_cfdna %>%
	dplyr::full_join(clinical, by = "patient_id") %>%
	dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
	reshape2::melt(id.vars = c("patient_id", "n"),
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
	ggplot(aes(x = value, y = n)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 15)) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("") +
	ylab(expression("Number of variants in cfDNA")) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Age at diagnosis"),
		    comparisons = list(c("<50 yrs", ">=50 yrs")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "ER status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "HER2 status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "PR status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Histological grade"),
		    comparisons = list(c("2", "3")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "pCR status"),
		    comparisons = list(c("No", "Yes")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Stage at diagnosis"),
		    comparisons = list(c("II", "III")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Tumor size"),
		    comparisons = list(c("<3 cm", ">=3 cm")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 12.5,
		    tip_length = 0.01) +
	facet_wrap(~variable, scales = "free_x", nrow = 2, ncol = 4)

pdf(file = "../res/cfDNA_Number_Variants.pdf", width = 9, height = 7)
print(plot_)
dev.off()

plot_ = any_vars_cfdna %>%
	dplyr::full_join(clinical %>%
			 dplyr::mutate(molecular_subtype = case_when(
				 ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
				 ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
				 HER2_status_pre_CT=="Positive" ~ "HER2+"
			 )),
			 by = "patient_id") %>%
	dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
	reshape2::melt(id.vars = c("patient_id", "n"),
		       measure.vars = "molecular_subtype") %>%
	dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)) %>%
	ggplot(aes(x = value, y = n)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 15)) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("") +
	ylab(expression("Number of variants in cfDNA")) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    comparisons = list(c("HR+", "HER2+"),
				       c("HR+", "TN"),
				       c("HER2+", "TN")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = c(12, 13, 14),
		    tip_length = 0.01)

pdf(file = "../res/cfDNA_Number_Variants_by_Molecular_Subtype.pdf", width = 4, height = 5)
print(plot_)
dev.off()

################################################################################################################################
## % vars in plasma
################################################################################################################################
any_vars_cfdna = maf_ft %>%
		 dplyr::filter(is_plasma_yes_no == "Yes") %>%
		 dplyr::group_by(patient_id) %>%
		 dplyr::summarize(n = n())

any_vars_tumor = maf_ft %>%
		 dplyr::filter(is_plasma_yes_no == "No") %>%
		 dplyr::filter(time_point == "Baseline\nTumor") %>%
		 dplyr::group_by(patient_id) %>%
		 dplyr::summarize(n_0 = n())

any_vars_cfdna = any_vars_cfdna %>%
		 dplyr::left_join(any_vars_tumor, by = "patient_id") %>%
		 dplyr::mutate(n_0 = ifelse(is.na(n_0), 0, n_0)) %>%
		 dplyr::mutate(n_1 = n/(n_0+n))

plot_ = any_vars_cfdna %>%
	dplyr::full_join(clinical, by = "patient_id") %>%
	dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19", "BC10"))) %>%
	reshape2::melt(id.vars = c("patient_id", "n_1"),
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
	ggplot(aes(x = value, y = 100*n_1)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 200)) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("") +
	ylab(expression("Fraction of variants in cfDNA (%)")) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Age at diagnosis"),
		    comparisons = list(c("<50 yrs", ">=50 yrs")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "ER status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "HER2 status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "PR status"),
		    comparisons = list(c("Negative", "Positive")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Histological grade"),
		    comparisons = list(c("2", "3")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "pCR status"),
		    comparisons = list(c("No", "Yes")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Stage at diagnosis"),
		    comparisons = list(c("II", "III")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	geom_signif(stat = "signif",
		    data = . %>%
		    	   dplyr::filter(variable == "Tumor size"),
		    comparisons = list(c("<3 cm", ">=3 cm")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	facet_wrap(~variable, scales = "free_x", nrow = 2, ncol = 4)

pdf(file = "../res/cfDNA_Fraction_Variants.pdf", width = 9, height = 7)
print(plot_)
dev.off()

plot_ = any_vars_cfdna %>%
	dplyr::full_join(clinical %>%
			 dplyr::mutate(molecular_subtype = case_when(
				 ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
				 ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
				 HER2_status_pre_CT=="Positive" ~ "HER2+"
			 )),
			 by = "patient_id") %>%
	dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
	reshape2::melt(id.vars = c("patient_id", "n_1"),
		       measure.vars = "molecular_subtype") %>%
	dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)) %>%
	ggplot(aes(x = value, y = 100*n_1)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = value, y = 100*n_1, color = value, fill = value),
		    data = any_vars_cfdna %>%
		    	   dplyr::full_join(clinical %>%
					    dplyr::mutate(molecular_subtype = case_when(
						    ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
						    ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
						    HER2_status_pre_CT=="Positive" ~ "HER2+"
					    )),
					    by = "patient_id") %>%
		    	   dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
		    	   dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
		    	   reshape2::melt(id.vars = c("patient_id", "n_1"),
					  measure.vars = "molecular_subtype") %>%
			   dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)),
		    stat = "identity", width = .1, height = 0, shape = 21, alpha = .75, size = 4.5, inherit.aes = FALSE) +
	scale_y_continuous(limits = c(0, 200)) +
	scale_fill_brewer(type = "qual", palette = 6) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("") +
	ylab(expression("Fraction of variants in cfDNA (%)")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 18),
 	      axis.title.y = element_text(margin = margin(r = 20), size = 18),
	      axis.text.x = element_text(size = 16),
	      axis.text.y = element_text(size = 16),
	      strip.background = element_blank()) +
	geom_signif(stat = "signif",
		    comparisons = list(c("HR+", "HER2+"),
				       c("HR+", "TN"),
				       c("HER2+", "TN")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = c(170, 185, 200),
		    tip_length = 0.01,
		    size = .5,
		    textsize = 5) +
	guides(color = FALSE, fill = FALSE)

pdf(file = "../res/cfDNA_Fraction_Variants_by_Molecular_Subtype.pdf", width = 5.5, height = 5.5)
print(plot_)
dev.off()

smry_ = any_vars_cfdna %>%
	dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19", "BC10"))) %>%
	dplyr::summarize(Q1 = 100*signif(quantile(n_1, probs = .25), 2),
			 Q2 = 100*signif(median(n_1), 2),
			 Q3 = 100*signif(quantile(n_1, probs = .75), 2))

pander::pander(smry_)

SAVE = any_vars_cfdna

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
	 dplyr::filter(!(patient_id == "BC21" & Variant_Classification == "In_Frame_Del"))



################################################################################################################################
## Any vars in plasma
################################################################################################################################
any_vars_cfdna = maf_ft %>%
		 dplyr::filter(is_plasma_yes_no == "Yes") %>%
		 dplyr::group_by(patient_id) %>%
		 dplyr::summarize(mean_af = mean(af, na.rm=TRUE),
				  max_af = max(af, na.rm=TRUE))


plot_ = SAVE %>%
	dplyr::full_join(any_vars_cfdna) %>%
	dplyr::full_join(clinical %>%
			 dplyr::mutate(molecular_subtype = case_when(
				 ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
				 ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
				 HER2_status_pre_CT=="Positive" ~ "HER2+"
			 )),
			 by = "patient_id") %>%
	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
	dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
	dplyr::mutate(max_af = ifelse(is.na(max_af), 0, max_af)) %>%
	dplyr::mutate(mean_af = ifelse(is.na(mean_af), 0, mean_af)) %>%
	reshape2::melt(id.vars = c("patient_id", "max_af", "mean_af", "n", "n_0", "n_1"),
		       measure.vars = "molecular_subtype") %>%
	dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)) %>%
	ggplot(aes(x = max_af, y = n_1*100, color = value, fill = value)) +
	geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .85, size = 1.5) +
	geom_smooth(data = SAVE %>%
		    	   dplyr::full_join(any_vars_cfdna) %>%
		    	   dplyr::full_join(clinical %>%
					    dplyr::mutate(molecular_subtype = case_when(
						    ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
						    ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
						    HER2_status_pre_CT=="Positive" ~ "HER2+"
					    )),
					    by = "patient_id") %>%
		    	   dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
		    	   dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
		    	   dplyr::mutate(max_af = ifelse(is.na(max_af), 0, max_af)) %>%
		    	   dplyr::mutate(mean_af = ifelse(is.na(mean_af), 0, mean_af)) %>%
		    	   reshape2::melt(id.vars = c("patient_id", "max_af", "mean_af", "n", "n_0", "n_1"),
					  measure.vars = "molecular_subtype") %>%
		    	   dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)),
		    mapping = aes(x = max_af, y = n_1*100),
		    stat = "smooth", formula = y ~ x, method = "lm", color = "goldenrod3", alpha = .55, se = FALSE, fullrange = TRUE, size = 1.5, inherit.aes = FALSE) +
	geom_point(stat = "identity", shape = 21, alpha = .75, size = 4.5) +
	stat_cor(data = SAVE %>%
		    	dplyr::full_join(any_vars_cfdna) %>%
		    	dplyr::full_join(clinical %>%
					 dplyr::mutate(molecular_subtype = case_when(
						 		ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
						 		ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
						 		HER2_status_pre_CT=="Positive" ~ "HER2+"
					 )),
					 by = "patient_id") %>%
		 	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
		 	dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
		 	dplyr::mutate(max_af = ifelse(is.na(max_af), 0, max_af)) %>%
		 	dplyr::mutate(mean_af = ifelse(is.na(mean_af), 0, mean_af)) %>%
		 	reshape2::melt(id.vars = c("patient_id", "max_af", "mean_af", "n", "n_0", "n_1"),
				       measure.vars = "molecular_subtype") %>%
		 	dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)),
		 mapping = aes(x = max_af, y = n_1*100),
		 method = "spearman", size = 6, inherit.aes = FALSE) +
	scale_x_sqrt(limits = c(0, 100),
		     breaks = c(.1, 1, 5, 10, 20, 30, 50, 100),
		     labels = c(.1, 1, 5, 10, 20, 30, 50, 100)) +
	scale_y_sqrt(limits = c(0, 200),
		     breaks = c(.1, 1, 5, 10, 20, 30, 50, 100, 200),
		     labels = c(.1, 1, 5, 10, 20, 30, 50, 100, 200)) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_fill_brewer(type = "qual", palette = 6) +
	xlab(expression("Max AF in cfDNA (%)")) +
	ylab(expression("Fraction of variants in cfDNA (%)")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 18),
 	      axis.title.y = element_text(margin = margin(r = 20), size = 18),
	      axis.text.x = element_text(size = 16),
	      axis.text.y = element_text(size = 16),
	      strip.background = element_blank()) +
	guides(color = guide_legend(title = "Molecular\nsubtype"),
	       fill = guide_legend(title = "Molecular\nsubtype"))

pdf(file = "../res/cfDNA_Max_AF_Fraction_Variants_Baseline_by_Molecular_Subtype.pdf", width = 6.5, height = 5.5)
print(plot_)
dev.off()

plot_ = SAVE %>%
	dplyr::full_join(any_vars_cfdna) %>%
	dplyr::full_join(clinical %>%
			 dplyr::mutate(molecular_subtype = case_when(
				 ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
				 ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
				 HER2_status_pre_CT=="Positive" ~ "HER2+"
			 )),
			 by = "patient_id") %>%
	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
	dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
	dplyr::mutate(max_af = ifelse(is.na(max_af), 0, max_af)) %>%
	dplyr::mutate(mean_af = ifelse(is.na(mean_af), 0, mean_af)) %>%
	reshape2::melt(id.vars = c("patient_id", "max_af", "mean_af", "n", "n_0", "n_1"),
		       measure.vars = "molecular_subtype") %>%
	dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)) %>%
	ggplot(aes(x = mean_af, y = n_1*100, color = value, fill = value)) +
	geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .85, size = 1.5) +
	geom_smooth(data = SAVE %>%
		    	   dplyr::full_join(any_vars_cfdna) %>%
		    	   dplyr::full_join(clinical %>%
					    dplyr::mutate(molecular_subtype = case_when(
						    ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
						    ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
						    HER2_status_pre_CT=="Positive" ~ "HER2+"
					    )),
					    by = "patient_id") %>%
		    	   dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
		    	   dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
		    	   dplyr::mutate(max_af = ifelse(is.na(max_af), 0, max_af)) %>%
		    	   dplyr::mutate(mean_af = ifelse(is.na(mean_af), 0, mean_af)) %>%
		    	   reshape2::melt(id.vars = c("patient_id", "max_af", "mean_af", "n", "n_0", "n_1"),
					  measure.vars = "molecular_subtype") %>%
		    	   dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)),
		    mapping = aes(x = mean_af, y = n_1*100),
		    stat = "smooth", formula = y ~ x, method = "lm", color = "goldenrod3", alpha = .55, se = FALSE, fullrange = TRUE, size = 1.5, inherit.aes = FALSE) +
	geom_point(stat = "identity", shape = 21, alpha = .75, size = 4.5) +
	stat_cor(data = SAVE %>%
		    	dplyr::full_join(any_vars_cfdna) %>%
		    	dplyr::full_join(clinical %>%
					 dplyr::mutate(molecular_subtype = case_when(
						 		ER_status_pre_CT=="Positive" & HER2_status_pre_CT=="Negative" ~ "HR+",
						 		ER_status_pre_CT=="Negative" & PR_status_pre_CT=="Negative" & HER2_status_pre_CT=="Negative" ~ "TN",
						 		HER2_status_pre_CT=="Positive" ~ "HER2+"
					 )),
					 by = "patient_id") %>%
		 	dplyr::filter(!(patient_id %in% c("BC05", "BC11", "BC14", "BC19"))) %>%
		 	dplyr::mutate(n_1 = ifelse(is.na(n_1), 0, n_1)) %>%
		 	dplyr::mutate(max_af = ifelse(is.na(max_af), 0, max_af)) %>%
		 	dplyr::mutate(mean_af = ifelse(is.na(mean_af), 0, mean_af)) %>%
		 	reshape2::melt(id.vars = c("patient_id", "max_af", "mean_af", "n", "n_0", "n_1"),
				       measure.vars = "molecular_subtype") %>%
		 	dplyr::mutate(value = factor(value, levels = c("HR+", "HER2+", "TN"), ordered = TRUE)),
		 mapping = aes(x = mean_af, y = n_1*100),
		 method = "spearman", size = 6, inherit.aes = FALSE) +
	scale_x_sqrt(limits = c(0, 100),
		     breaks = c(.1, 1, 5, 10, 20, 30, 50, 100),
		     labels = c(.1, 1, 5, 10, 20, 30, 50, 100)) +
	scale_y_sqrt(limits = c(0, 200),
		     breaks = c(.1, 1, 5, 10, 20, 30, 50, 100, 200),
		     labels = c(.1, 1, 5, 10, 20, 30, 50, 100, 200)) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_fill_brewer(type = "qual", palette = 6) +
	xlab(expression("Mean AF in cfDNA (%)")) +
	ylab(expression("Fraction of variants in cfDNA (%)")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 18),
 	      axis.title.y = element_text(margin = margin(r = 20), size = 18),
	      axis.text.x = element_text(size = 16),
	      axis.text.y = element_text(size = 16),
	      strip.background = element_blank()) +
	guides(color = guide_legend(title = "Molecular\nsubtype"),
	       fill = guide_legend(title = "Molecular\nsubtype"))

pdf(file = "../res/cfDNA_Mean_AF_Fraction_Variants_Baseline_by_Molecular_Subtype.pdf", width = 6.5, height = 5.5)
print(plot_)
dev.off()

