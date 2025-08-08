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
		grepl("_2", sample_id) ~ "On-treatment",
		grepl("_3", sample_id) ~ "Post-treatment"
	)) %>%
	readr::type_convert()


################################################################################################################################
## cfDNA concentration
################################################################################################################################
plot_ = smry_ %>%
	dplyr::mutate(rcb_classification = case_when(
		RCB == "0" | RCB == "I" ~ "0/I",
		RCB == "II" | RCB == "III" ~ "II/III"
	)) %>%
	reshape2::melt(id.vars = c("sample_id", "rcb_classification", "time_point"),
		       measure.vars = "concentration_pg_muL",
		       value.name = "concentration_pg_muL") %>%
	dplyr::select(-variable) %>%
	ggplot(aes(x = rcb_classification, y = concentration_pg_muL)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "grey20", fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", shape = 21, alpha = .75, size = 3.5) +
	scale_y_continuous(limits = c(0, 3500),
			   labels = scientific_10) +
	xlab("Residual Cancer Burden") +
	ylab(expression("cfDNA (pg"%.%mu~"L"^-1~")")) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_signif(stat = "signif",
		    comparisons = list(c("0/I", "II/III")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 3000,
		    tip_length = 0.02,
		    vjust = -.2) +
	facet_wrap(~time_point, scales = "free_x", nrow = 1, ncol = 3) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      strip.background = element_blank())

pdf(file = "../res/Supplementary_Figure_S7.pdf", width = 9, height = 7/1.75)
print(plot_)
dev.off()

