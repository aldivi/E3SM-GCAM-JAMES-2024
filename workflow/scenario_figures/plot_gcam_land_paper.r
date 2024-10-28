# plot_gcam_land_paper.r

# panel plot of global crop/forest gcam outputs for 4 simulations

# also plot the gcam-elm crop forest comparison figure

# read in the data from the diagnostic plots and replot them

library(ggplot2)

indir = "./"
outdir = "./scenario_outputs/"

# create the output directory
dir.create(outdir)

#### gcam outputs

# the control data are in each file
gcam_fdbk_file = paste0(indir, "compare_gcam_cropforest_land_full_fdbk.csv")
gcam_agfdbk_file = paste0(indir, "compare_gcam_cropforest_land_ag_fdbk.csv")
gcam_cfdbk_file = paste0(indir, "compare_gcam_cropforest_land_c_fdbk.csv")

gcam_fdbk_in = read.csv(gcam_fdbk_file, stringsAsFactors = FALSE)
gcam_agfdbk_in = read.csv(gcam_agfdbk_file, stringsAsFactors = FALSE)
gcam_cfdbk_in = read.csv(gcam_cfdbk_file, stringsAsFactors = FALSE)

plot_data = rbind(gcam_fdbk_in, gcam_agfdbk_in[gcam_agfdbk_in$scenario == "AG_FDBK",], gcam_cfdbk_in[gcam_cfdbk_in$scenario == "C_FDBK",])
plot_data$X = NULL
plot_data = plot_data[plot_data$region == "Global",]
panlabs = c("a) Bioenergy crop", "b) Harvested forest", "c) Non-bioenergy crop", "d) Non-harvested forest", "e) All crop", "f) All forest")
scenlabs = c("C_FDBK", "AG_FDBK", "FULL_FDBK", "CONTROL")
plot_data$scenario = factor(plot_data$scenario, levels = scenlabs)

p1 <- ggplot(data=plot_data, aes(x=year, y=gcam_cover_area_thouskm2, color=scenario)) +
	geom_line(linewidth=1) +
	#scale_x_continuous(breaks=seq(from=2015, to=2090, by=15), expand=c(0.01,0.01)) +
	xlim(2015, 2100) +
	labs(y=paste("Thousand",expression(km^2)), x="Year") +
	theme_bw() +   
	facet_wrap(~(factor(gcam_cover, levels=c("gcam_biomass", "gcam_man_forest", "gcam_nb_crop",
	                                         "gcam_unman_forest", "gcam_crop", "gcam_forest"), labels = panlabs)), scales="free_y", ncol=2) +
	scale_color_manual(values=c(C_FDBK = "purple", AG_FDBK = "red", FULL_FDBK = "black", CONTROL = "orange"))
	
print(p1)

ggsave(paste0(outdir, "compare_cropforest_paper.pdf"), plot=p1, device="pdf")
write.csv(plot_data, file=paste0(outdir, "compare_cropforest_paper.csv"), row.names=FALSE)




#### gcam elm comparison
# out values in thousand km^2
# plot changes only for forest and crop
# normalize to elm area

# elm

elm_file = paste0(indir, "compare_surfdata_iESM_dyn_fdbk_nofdbk.csv")
elm_in = read.csv(elm_file, stringsAsFactors = FALSE)
	
forest_fdbk_df = elm_in[, c("year", "forest_with_fdbks", "units")]
names(forest_fdbk_df) = c("Year", "value", "units")
forest_fdbk_df$scenario = "FULL_FDBK"
forest_fdbk_df$type = "b) Forest"
forest_fdbk_df$model = "ELM"
forest_nofdbk_df = elm_in[, c("year", "forest", "units")]
names(forest_nofdbk_df) = c("Year", "value", "units")
forest_nofdbk_df$scenario = "CONTROL"
forest_nofdbk_df$type = "b) Forest"
forest_nofdbk_df$model = "ELM"

crop_fdbk_df = elm_in[, c("year", "crop_with_fdbks", "units")]
names(crop_fdbk_df) = c("Year", "value", "units")
crop_fdbk_df$scenario = "FULL_FDBK"
crop_fdbk_df$type = "a) Crop"
crop_fdbk_df$model = "ELM"
crop_nofdbk_df = elm_in[, c("year", "crop", "units")]
names(crop_nofdbk_df) = c("Year", "value", "units")
crop_nofdbk_df$scenario = "CONTROL"
crop_nofdbk_df$type = "a) Crop"
crop_nofdbk_df$model = "ELM"

elm_data = rbind(forest_fdbk_df, forest_nofdbk_df, crop_fdbk_df, crop_nofdbk_df)

# convert elm to thousand km^2
elm_data$value = elm_data$value / 1000
elm_data$units = "area_thous_km2"

# gcam

gcam_df = gcam_fdbk_in[gcam_fdbk_in$region == "Global" & (gcam_fdbk_in$gcam_cover == "gcam_forest" | gcam_fdbk_in$gcam_cover == "gcam_crop"), c("year", "gcam_cover", "gcam_cover_area_thouskm2", "scenario")]
names(gcam_df) = c("Year", "type", "value", "scenario")
gcam_df$units = "area_thous_km2"
gcam_df$model = "GCAM"
#gcam_df = gcam_df[,names(plot_data)]
gcam_df$type[gcam_df$type == "gcam_crop"] = "a) Crop"
gcam_df$type[gcam_df$type == "gcam_forest"] = "b) Forest"

# shift gcam crop area to align with elm crop area
gcam_df$gcam_glbarea_shifted = NA
gcam_df$gcam_glbarea_shifted[gcam_df$type == "a) Crop" & gcam_df$scenario == "CONTROL"] =
	gcam_df$value[gcam_df$type == "a) Crop" & gcam_df$scenario == "CONTROL"] - gcam_df$value[gcam_df$type == "a) Crop" & gcam_df$scenario == "CONTROL" & gcam_df$Year == 2015] +
	crop_nofdbk_df$value[crop_nofdbk_df$Year == 2015] / 1000
gcam_df$gcam_glbarea_shifted[gcam_df$type == "a) Crop" & gcam_df$scenario == "FULL_FDBK"] =
	gcam_df$value[gcam_df$type == "a) Crop" & gcam_df$scenario == "FULL_FDBK"] - gcam_df$value[gcam_df$type == "a) Crop" & gcam_df$scenario == "FULL_FDBK" & gcam_df$Year == 2015] +
	crop_fdbk_df$value[crop_fdbk_df$Year == 2015] / 1000
gcam_df$gcam_glbarea_shifted[gcam_df$type == "b) Forest" & gcam_df$scenario == "CONTROL"] =
	gcam_df$value[gcam_df$type == "b) Forest" & gcam_df$scenario == "CONTROL"] - gcam_df$value[gcam_df$type == "b) Forest" & gcam_df$scenario == "CONTROL" & gcam_df$Year == 2015] +
	forest_nofdbk_df$value[forest_nofdbk_df$Year == 2015] / 1000
gcam_df$gcam_glbarea_shifted[gcam_df$type == "b) Forest" & gcam_df$scenario == "FULL_FDBK"] =
	gcam_df$value[gcam_df$type == "b) Forest" & gcam_df$scenario == "FULL_FDBK"] - gcam_df$value[gcam_df$type == "b) Forest" & gcam_df$scenario == "FULL_FDBK" & gcam_df$Year == 2015] +
	forest_fdbk_df$value[forest_fdbk_df$Year == 2015] / 1000

gcam_out = gcam_df
gcam_out$value = gcam_out$gcam_glbarea_shifted
gcam_out$gcam_glbarea_shifted = NULL

# this will show gcam in steps
if(FALSE){
for (s in unique(gcam_out$scenario)) {
	for (t in unique(gcam_out$type)) {
		for (y in plot_data$Year) {
			if (!(y %in% gcam_out$Year)) {
				# only add record if year does not exist
				temp_df = gcam_out[gcam_out$scenario == s & gcam_out$type == t & gcam_out$Year < y,]
				temp_df = temp_df[which(temp_df$Year == max(temp_df$Year)),]
				temp_df$Year = y
				gcam_out = rbind(gcam_out, temp_df)
			} # end if add record
		} # end y loop over all years
	} # end t loop over type
} # end s loop over scenario
}

plot2_data = rbind(elm_data, gcam_out)

p2 <- ggplot(data=plot2_data, aes(x=Year, y=value, color=scenario, linetype = model)) +
	geom_line(linewidth=1) +
	#scale_x_continuous(breaks=seq(from=2015, to=2090, by=15), expand=c(0.01,0.01)) +
	xlim(2015, 2100) +
	labs(y=paste("Thousand",expression(km^2))) +
	theme_bw() +   
	facet_wrap(~(factor(type, levels=c("a) Crop", "b) Forest"))), scales="free_y", ncol=1) +
	scale_color_manual(values=c(FULL_FDBK = "black", CONTROL = "orange"))
	
print(p2)

ggsave(paste0(outdir, "compare_elm_gcam_land_paper.pdf"), plot=p2, device="pdf")
write.csv(plot_data, file=paste0(outdir, "compare_elm_gcam_land_paper.csv"), row.names=FALSE)





