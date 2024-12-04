#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

ddPCR = readr::read_tsv(file = url_ddpcr, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	dplyr::filter(`is_tumor?` == "T") %>%
	readr::type_convert() %>%
	dplyr::mutate(af = 100 * droplet_count_mutant / (droplet_count_mutant + droplet_count_wt)) %>%
	dplyr::mutate(af = ifelse(is.na(af), 0, af)) %>%
	dplyr::mutate(af = ifelse(af==0, 0.02, af)) %>%
	dplyr::mutate(HGVSp_Short = gsub("p.", "", HGVSp, fixed = TRUE)) %>%
	dplyr::mutate(Hugo_Symbol = substr(Hugo_Symbol, start = 1, stop = 4)) %>%
        dplyr::mutate(HGVSp_Short = substr(HGVSp_Short, start = 1, stop = 3)) %>%
	dplyr::mutate(HGVSp_Short = paste0(Hugo_Symbol, " ", HGVSp_Short)) %>%
	dplyr::filter(!(patient_id == "BC18" & Hugo_Symbol == "GATA"))

patient_ids = unique(ddPCR$patient_id)

for (i in 1:length(patient_ids)) {
	ddPCR_ft = ddPCR %>%
		   dplyr::filter(patient_id == patient_ids[i])
	
	colourCount = length(unique(ddPCR_ft$HGVSp_Short))
	getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	
	plot_ = ddPCR_ft %>%
		dplyr::mutate(time_point = case_when(
			grepl("_1", sample_type) ~ "Baseline\nPlasma",
			grepl("_2", sample_type) ~ "On-Treat\n-ment",
			grepl("_3", sample_type) ~ "Post-Treat\n-ment"
		)) %>%
		dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nPlasma", "On-Treat\n-ment", "Post-Treat\n-ment"), ordered = TRUE)) %>%
		ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
		geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
		geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
		xlab("") +
		ylab("\nAllele Fraction (%)\n\n") +
		scale_color_manual(values = getPalette(colourCount)) +
		scale_x_discrete(expand = c(0.1, 0),
				 breaks = c("Baseline\nPlasma", "On-Treat\n-ment", "Post-Treat\n-ment"),
				 labels = c("Baseline\nPlasma", "On-Treat\n-ment", "Post-Treat\n-ment"),
				 drop = FALSE) +
		scale_y_log10(limits = c(0.02, 100),
			      breaks = c(0.02, 0.1, 1, 10, 100),
			      labels = c("ND   ", ".1   ", "1   ", "10   ", "100   ")) +
		coord_cartesian(clip = "off") +
		theme(legend.background = element_rect(fill = NA, size = 0.5),
		      legend.key = element_rect(fill = NA),
		      legend.text = element_text(size = 7)) +
		guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
		theme_cowplot()

	pdf(file = paste0("../res/Figure_3-", patient_ids[i], ".pdf"), width = 5, height = 3)
	print(plot_)
	dev.off()
}
