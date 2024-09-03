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
	   dplyr::mutate(purity = ifelse(is.na(purity), 1, purity)) %>%
	   dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nTumor", "On-treat\n-ment", "Residual\nDisease", "Local/ Distant\nRelapse"), ordered = TRUE)) %>%
	   dplyr::arrange(time_point)
	   

log2_ratio = readr::read_tsv(file = url_log2_ratio, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert()

total_copy = readr::read_tsv(file = url_total_copy, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert()

outlier = foreach (i=1:length(blacklist_samples)) %dopar% {
	outlier = log2_ratio %>%
		  dplyr::filter(tumor_id == blacklist_samples[i]) %>%
		  dplyr::mutate(p = mean(Log2_Ratio, na.rm=TRUE),
				q = sd(Log2_Ratio, na.rm = TRUE)) %>%
		  dplyr::filter(Log2_Ratio > p+1.5*q | Log2_Ratio < p-1.5*q) %>%
		  dplyr::select(Chromosome, Position)
	return(invisible(outlier))
}
outlier = do.call(rbind, outlier) %>%
	  dplyr::mutate(UUID = paste0(Chromosome, ":", Position)) %>%
	  dplyr::mutate(Is_Outlier = TRUE) %>%
	  dplyr::filter(!duplicated(UUID)) %>%
	  dplyr::select(-UUID)

Log2_Ratio = foreach (i=1:(length(unique(manifest$patient_id)))) %dopar% {
	smoothed_log2 = log2_ratio %>%
			dplyr::filter(tumor_id %in% (manifest %>%
						     dplyr::filter(patient_id == (unique(manifest$patient_id))[i]) %>%
						     .[["tumor_id"]])) %>%
			dplyr::filter(Chromosome != "Y") %>%
			dplyr::arrange(tumor_id, Chromosome, Position) %>%
			dplyr::left_join(outlier, by = c("Chromosome", "Position")) %>%
			dplyr::mutate(Is_Outlier = ifelse(is.na(Is_Outlier), FALSE, Is_Outlier)) %>%
			dplyr::filter(!Is_Outlier) %>%
			reshape2::dcast(formula = Chromosome + Position ~ tumor_id, value.var = "Log2_Ratio") %>%
			copynumber::winsorize(method = "mad", , tau = 2.0, k = 50, verbose = FALSE) %>%
			dplyr::rename(Chromosome = chrom, Position = pos) %>%
			reshape2::melt(id.vars = c("Chromosome", "Position"), variable.name = "sample_name", value.name = "Log2_Ratio")

	segmented_log2 = smoothed_log2 %>%
			 reshape2::dcast(formula = Chromosome + Position ~ sample_name, fill = 0, value.var = "Log2_Ratio") %>%
			 copynumber::multipcf(gamma = 5, normalize = FALSE, fast = TRUE, verbose = FALSE) %>%
			 dplyr::select(-arm, -n.probes) %>%
			 dplyr::rename(Chromosome = chrom, Start_Position = start.pos, End_Position = end.pos) %>%
			 dplyr::mutate(Chromosome = case_when(
				 		Chromosome == 23 ~ "X",
				 		TRUE ~ as.character(Chromosome)
			 )) %>%
			 dplyr::mutate(Chromosome = factor(Chromosome, levels = c(1:22, "X"), ordered = TRUE)) %>%
			 reshape2::melt(id.vars = c("Chromosome", "Start_Position", "End_Position"), variable.name = "sample_name", value.name = "Log2_Ratio") %>%
			 dplyr::left_join(manifest %>% dplyr::rename(sample_name = tumor_id), by = "sample_name") %>%
			 dplyr::mutate(absolute_copies = ((2^Log2_Ratio)*((purity*ploidy)+(2-2*purity)) - 2 + 2*purity)/purity) %>%
		 	 dplyr::mutate(absolute_copies = round(absolute_copies))

	segmented_log2 = segmented_log2 %>%
		 	 dplyr::left_join(segmented_log2 %>%
					  dplyr::group_by(sample_name, absolute_copies) %>%
					  dplyr::summarize(Log2_Ratio_C = mean(Log2_Ratio)) %>%
					  dplyr::ungroup(),
					  by = c("sample_name", "absolute_copies")) %>%
			 dplyr::mutate(Log2_Ratio_C = ifelse(Log2_Ratio_C>4, 3.95, Log2_Ratio_C))

	plot_ = smoothed_log2 %>%
		dplyr::mutate(Chromosome_Color = Chromosome %% 2) %>%
		dplyr::mutate(Chromosome_Color = as.character(Chromosome_Color)) %>%
		dplyr::mutate(Chromosome = case_when(
			Chromosome == 23 ~ "X",
			TRUE ~ as.character(Chromosome)
		)) %>%
		dplyr::mutate(Chromosome = factor(Chromosome, levels = c(1:22, "X"), ordered = TRUE)) %>%
		ggplot(aes(x = Position, y = Log2_Ratio, color = Chromosome_Color)) +
		geom_point(stat = "identity", fill = NA, shape = 16, size = .85, alpha = 1) +
		scale_color_manual(values = c("0" = "grey", "1" = "lightblue")) +
		geom_segment(data = segmented_log2,
			     mapping = aes(x = Start_Position, y = Log2_Ratio_C, xend = End_Position, yend = Log2_Ratio_C),
			     color = "#e41a1c", inherit.aes = FALSE, size = 1, lineend = "round") +
		geom_line(data = segmented_log2 %>%
				 tidyr::pivot_longer(c(Start_Position, End_Position)),
			     mapping = aes(x = value, y = Log2_Ratio_C),
			     color = "#e41a1c", inherit.aes = FALSE, size = .1, lineend = "round") +
		geom_hline(yintercept = 0, color = "grey50", linetype = 2, size = .5) +
		xlab("\n\n") +
		ylab(expression(Log[2]~"Ratio")) +
		scale_x_continuous(labels = scientific_10) +
		scale_y_continuous(expand = c(0, 0),
				   breaks = c(-2, -1, 0, 1, 2, 3, 4),
				   labels = c("-2", "", "0", "", "2", "", "4"),
				   limits = c(-2, 4)) +
		theme_minimal() +
		theme(axis.text.x = element_blank(),
		      axis.ticks.x = element_blank(),
		      axis.text.y = element_text(size = 9),
		      axis.title.y = element_text(size = 11, margin = margin(r = 30)),
		      strip.text.x = element_text(size = 8),
		      strip.text.y = element_text(size = 9),
		      panel.spacing.x = unit(0.05, 'lines'),
		      panel.spacing.y = unit(1, 'lines'),
		      panel.grid.minor = element_blank()) +
		facet_grid(sample_name~Chromosome, scales = "free_x", space = "free_x", switch = "x") +
		guides(color = FALSE)

	h = manifest %>%
	    dplyr::filter(patient_id == (unique(manifest$patient_id))[i]) %>%
	    nrow() * 2


	pdf(file = paste0("../res/Log2_", unique(manifest$patient_id)[i], ".pdf"), width = 9, height = h)
	print(plot_)
	dev.off()
	
}
