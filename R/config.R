#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
suppressPackageStartupMessages(library("gdata"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("superheat"))
suppressPackageStartupMessages(library("ggsignif"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggforce"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("preseqR"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("fuzzyjoin"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggpubr"))

registerDoMC(8)

url_manifest <- "../data/manifest.txt"
url_clinical <- "../data/clinical.txt"
url_mutation_smry_ft <- "../data/mutation_summary_ft.maf"
url_mutation_smry <- "../data/mutation_summary.maf"
url_ddpcr <- "../data/ddPCR.txt"
url_ddpcr_manifest <- "../data/ddPCR_manifest.txt"
url_preanalytical_conditions <- "../data/preanalytical_conditions.txt"
url_log2_ratio <- "../data/Log2_Ratio.txt"
url_total_copy <- "../data/Total_Copy.txt"
url_purity_ploidy <- "../data/Purity_Ploidy.txt"

url_aln_metrics <- "../data/aln_metrics.txt"
url_gc_metrics <- "../data/gc_metrics.txt"
url_hs_metrics <- "../data/hs_metrics.txt"
url_idx_metrics <- "../data/idx_metrics.txt"
url_insert_metrics <- "../data/insert_metrics.txt"
url_oxog_metrics <- "../data/oxog_metrics.txt"

blacklist_samples <- c("BC01-P1", "BC02-P1", "BC03-P1", "BC04-P1", "BC05-P1", "BC06-P1", "BC09-P1")

'vc_colors' <- function() {
	init_cols = maftools:::get_vcColors(websafe = FALSE)
	init_cols = c(init_cols, "Hotspot" = as.character(init_cols["Nonsense_Mutation"]))
	init_cols["Nonsense_Mutation"] = as.character(init_cols["Frame_Shift_Del"])
	init_cols["Frame_Shift_Ins"] = as.character(init_cols["Frame_Shift_Del"])
	init_cols["In_Frame_Ins"] = as.character(init_cols["In_Frame_Del"])
	return(invisible(init_cols))
}

'vc_legend' <- function() {
	init_cols = dplyr::tibble(Variant_Classification = names(vc_colors()),
				  Variant_Color = as.character(vc_colors())) %>%
		    dplyr::filter(!(Variant_Classification %in% c("IGR", "Silent", "RNA", "Intron", "ITD", "Multi_Hit", "Amp", "Del", "Complex_Event", "pathway"))) %>%
		    dplyr::mutate(Variant_Classification = case_when(
			    Variant_Classification == "Nonsense_Mutation" ~ "Nonsense SNV",
			    Variant_Classification == "Frame_Shift_Del" ~ "Frameshifting indel",
			    Variant_Classification == "Frame_Shift_Ins" ~ "Frameshifting indel",
			    Variant_Classification == "In_Frame_Del" ~ "In-frame indel",
			    Variant_Classification == "In_Frame_Ins" ~ "In-frame indel",
			    Variant_Classification == "Splice_Site" ~ "Splice site",
			    Variant_Classification == "Missense_Mutation" ~ "Missense SNV",
			    Variant_Classification == "Translation_Start_Site" ~ "Translation start site",
			    Variant_Classification == "Nonstop_Mutation" ~ "Nonstop SNV",
			    TRUE ~ Variant_Classification
		    )) %>%
		   dplyr::filter(!duplicated(Variant_Classification)) %>%
		   dplyr::bind_rows(dplyr::tibble(Variant_Classification = "Loss of Heterozygosity",
					 	  Variant_Color = "#000000FF")) %>%
		   dplyr::mutate(Variant_Classification = factor(Variant_Classification, levels = c("Hotspot", "Frameshifting indel", "Nonstop SNV", "Nonsense SNV", "Missense SNV", "In-frame indel", "Splice site", "Translation start site", "Loss of Heterozygosity"))) %>%
		   dplyr::mutate(Variant_Shape = case_when(
			   Variant_Classification == "Loss of Heterozygosity" ~ 1,
			   TRUE ~ 0
		   )) %>%
		   dplyr::mutate(Variant_Shape = factor(Variant_Shape))
		  
		

	manual_colors = init_cols$Variant_Color
	names(manual_colors) = init_cols$Variant_Classification
	plot_ = ggplot(data = init_cols, aes(x = Variant_Classification, y = Variant_Color, fill = Variant_Classification)) +
		geom_point(stat = "identity", shape = 22, size = 5, color = "white") +
		scale_fill_manual(values = manual_colors) +
		guides(fill = guide_legend(title = "Somatic mutation type")) +
		theme_classic()
	legend_ = cowplot::get_legend(plot_)
	grid.newpage()
	grid.draw(legend_)
}
	
'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}
