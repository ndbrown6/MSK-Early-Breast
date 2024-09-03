#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(time_point)) %>%
	   dplyr::left_join(clinical, by = "patient_id") %>%
	   dplyr::rename(Tumor_Sample_Barcode = sample_id) %>%
	   dplyr::mutate(time_point = case_when(
		time_point == "pre-treatment" ~ "Baseline\nTumor",
		time_point == "on-treatment" ~ "On-treat\n-ment",
		time_point == "post-treatment" ~ "Residual\nDisease",
		time_point == "follow-up recurrence" ~ "Local/ Distant\nRelapse",
	   )) %>%
	   dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nTumor", "On-treat\n-ment", "Residual\nDisease", "Local/ Distant\nRelapse"), ordered = TRUE))

maf = readr::read_tsv(file = url_mutation_smry, comment = "#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
      readr::type_convert() %>%
      dplyr::left_join(manifest, by = "Tumor_Sample_Barcode") %>%
      dplyr::filter(is_tumor_yes_no == "Yes") %>%
      dplyr::filter(is_plasma_yes_no == "No") %>%
      dplyr::filter(!is.na(HGVSp_Short),
		    !is.na(Hugo_Symbol),
		    Variant_Classification != "Silent") %>%
      dplyr::mutate(HGVSp_Short = gsub(pattern = "p.", replacement = "", x = HGVSp_Short, fixed = TRUE)) %>%
      dplyr::mutate(af = 100 * t_alt_count / t_depth) %>%
      dplyr::mutate(af = ifelse(t_depth==0, 0.1, af)) %>%
      dplyr::mutate(Hugo_Symbol = substr(Hugo_Symbol, start = 1, stop = 4)) %>%
      dplyr::mutate(HGVSp_Short = substr(HGVSp_Short, start = 1, stop = 3)) %>%
      dplyr::mutate(HGVSp_Short = paste0(Hugo_Symbol, " ", HGVSp_Short))

################################################################################################################################
## BC01
################################################################################################################################
i = 1
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i])

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	 
plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC02
################################################################################################################################
i = 2
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::mutate(time_point = case_when(
		 Tumor_Sample_Barcode == "BC02-T" ~ "Baseline\nTumor",
		 Tumor_Sample_Barcode == "BZ340691-T" ~ "Brain\nMetastasis"
	 )) %>%
	 dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nTumor", "Brain\nMetastasis"), ordered = TRUE)) %>%
	 dplyr::filter(Hugo_Symbol != "NOTC")

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC03
################################################################################################################################
i = 3
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::bind_rows(maf %>%
		   	  dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
		   	  dplyr::mutate(time_point = "Residual\nDisease",
				 	af = 0.1))

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC04
################################################################################################################################
i = 4
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i])

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC06
################################################################################################################################
i = 6
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::bind_rows(maf %>%
		   	  dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
		   	  dplyr::mutate(time_point = "Residual\nDisease",
				 	af = 0.1)) %>%
	 dplyr::filter(!grepl("NOTC", Hugo_Symbol, perl = TRUE))

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(
	      legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC07
################################################################################################################################
i = 7
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::filter(Tumor_Sample_Barcode != "LE984537-T") %>%
	 dplyr::filter(Tumor_Sample_Barcode != "UO030693-T") %>%
	 dplyr::mutate(time_point = case_when(
		 Tumor_Sample_Barcode == "BC07-T" ~ "Baseline\nTumor",
		 Tumor_Sample_Barcode == "BC07-T2" ~ "Residual\nDisease"
	 )) %>%
	 dplyr::mutate(time_point = factor(time_point, levels = c("Baseline\nTumor", "Residual\nDisease"), ordered = TRUE))

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC08
################################################################################################################################
i = 8
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i])

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC09
################################################################################################################################
i = 9
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i])

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC12
################################################################################################################################
i = 12
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::bind_rows(maf %>%
		   	  dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
		   	  dplyr::mutate(time_point = "Residual\nDisease",
				 	af = 0.1))

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC13
################################################################################################################################
i = 13
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::bind_rows(maf %>%
		   	  dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
		   	  dplyr::mutate(time_point = "Residual\nDisease",
				 	af = 0.1))

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC16
################################################################################################################################
i = 16
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::bind_rows(maf %>%
		   	  dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
		   	  dplyr::mutate(time_point = "Residual\nDisease",
				 	af = 0.1)) %>%
	 dplyr::filter(!grepl("APC|ASXL|GRIN|MDC1|MITF|RNF4|SMO|TRAF|MED1", Hugo_Symbol, perl=TRUE))

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC17
################################################################################################################################
i = 17
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
	 dplyr::bind_rows(maf %>%
		   	  dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i]) %>%
		   	  dplyr::mutate(time_point = "Residual\nDisease",
				 	af = 0.1)) %>%
	 dplyr::filter(!grepl("AMER|AR|EPHA|FAT1|GRIN|MDC|POLE|RFWD|TSC1", Hugo_Symbol, perl=TRUE))

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()

################################################################################################################################
## BC20
################################################################################################################################
i = 19
maf_ft = maf %>%
	 dplyr::filter(patient_id == (clinical %>% .[["patient_id"]])[i])

colourCount = length(unique(maf_ft$HGVSp_Short))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot_ = maf_ft %>%
	ggplot(aes(x = time_point, y = af, group = HGVSp_Short, color = HGVSp_Short)) +
	geom_line(stat = "identity", alpha = .65, size = .95, show.legend = FALSE) +
	geom_point(stat = "identity", shape = 21, size = 2, fill = "white") +
	xlab("") +
	ylab("\nAllele Fraction (%)\n\n") +
	scale_color_manual(values = getPalette(colourCount)) +
	scale_x_discrete(expand = c(0.1, 0)) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 20, 40, 60, 80, 100),
			   labels = c("ND", 20, 40, 60, 80, 100)) +
	theme(legend.background = element_rect(fill = NA, size = 0.5),
	      legend.key = element_rect(fill = NA),
	      legend.text = element_text(size = 7)) +
	guides(color = guide_legend(title = "", override.aes = list(size = 2.5))) +
	theme_cowplot()
pdf(file = paste0("../res/AF__", (clinical %>% .[["patient_id"]])[i], "__.pdf"), width = 5, height = 3)
print(plot_)
dev.off()
