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


ddPCR = readr::read_tsv(file = url_ddpcr, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	dplyr::filter(`is_tumor?` == "T") %>%
	readr::type_convert() %>%
	dplyr::mutate(af = 100 * droplet_count_mutant / (droplet_count_mutant + droplet_count_wt)) %>%
	dplyr::mutate(af = ifelse(is.na(af), 0, af)) %>%
	dplyr::filter(!(patient_id == "BC18" & Hugo_Symbol == "GATA3")) %>%
	dplyr::filter(grepl("_1", sample_type)) %>%
	dplyr::rename(HGVSp_Short = HGVSp)


################################################################################################################################
## AF vars in plasma
################################################################################################################################
cfdna_vars = maf_ft %>%
	     dplyr::filter(is_plasma_yes_no == "Yes") %>%
	     dplyr::right_join(ddPCR, by = c("patient_id", "Hugo_Symbol", "HGVSp_Short")) %>%
	     dplyr::mutate(af.x = ifelse(is.na(af.x), 0, af.x)) %>%
	     dplyr::mutate(af.y = ifelse(is.na(af.y), 0, af.y))

plot_ = cfdna_vars %>%
	ggplot(aes(x = af.x+.1, y = af.y+.1)) +
	geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .85, size = 1.5) +
	geom_smooth(stat = "smooth", formula = y ~ x +0, method = "lm", color = "goldenrod3", alpha = .55, se = FALSE, fullrange = TRUE, size = 1.5) +
	geom_point(stat = "identity", shape = 21, size = 3.5, color = "black", fill = "white", alpha = .75) +
	stat_cor(method = "spearman", size = 6) +
	scale_x_sqrt(limits = c(.1, 33),
		     breaks = c(0, 1, 2, 4, 8, 16, 32)+.1,
		     labels = c("ND", 1, 2, 4, 8, 16, 32)) +
	scale_y_sqrt(limits = c(.1, 33),
		     breaks = c(0, 1, 2, 4, 8, 16, 32)+.1,
		     labels = c("ND", 1, 2, 4, 8, 16, 32)) +
	xlab("AF by MSK-IMPACT (%)") +
	ylab("AF by ddPCR (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 18),
 	      axis.title.y = element_text(margin = margin(r = 20), size = 18),
	      axis.text.x = element_text(size = 16),
	      axis.text.y = element_text(size = 16),
	      strip.background = element_blank())

pdf(file = "../res/Figure_2A.pdf", width = 5.5, height = 5.5)
print(plot_)
dev.off()
