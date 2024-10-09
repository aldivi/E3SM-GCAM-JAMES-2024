require(devtools)
library(rgcam)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

indir = "./"
outdir = "./"

# gcam output data file with detailed land allocation and ag yield scalar data (with feedbacks)
gcam_output_name_fdbk = paste0(indir,"20240730_SSP245_ZATM.dat")

# gcam output data file with detailed land allocation and ag yield scalar data (without feedbacks)
gcam_output_name_nofdbk = paste0(indir,"20240730_SSP245_ZATM_without_feedbacks.dat")

start_year=2015
end_year=2100

oname = paste0(outdir, "crop_prices_fdbk_nofdbk")

if(substr(outdir,nchar(outdir), nchar(outdir)) != "/") { outdir = paste0(outdir, "/") }
dir.create(outdir, recursive=TRUE)

scheme_basic <- theme_bw() +
  theme(legend.text = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 18, face = "bold")) +
  theme(plot.title = element_text(size = 15, face = "bold", vjust = 1)) +
  theme(plot.subtitle = element_text(size = 9, face = "bold", vjust = 1))+
  theme(strip.text = element_text(size = 7))+
  theme(strip.text.x = element_text(size = 10, face = "bold"))+
  theme(strip.text.y = element_text(size = 15, face = "bold"))+
  theme(legend.position = "right")+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_text(size = 12,color = "black",face="bold"))+
  theme(axis.text.x= element_text(hjust=1,angle=90))+
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
        

prj0 <- loadProject(gcam_output_name_nofdbk)
gcam_agcom_in0 <-getQuery(prj0, query="ag commodity prices")
agcom0 = as.data.frame(gcam_agcom_in0)
agcom0$scenario = "CONTROL"

prj <- loadProject(gcam_output_name_fdbk)
gcam_agcom_in <-getQuery(prj, query="ag commodity prices")
agcom = as.data.frame(gcam_agcom_in)
agcom$scenario = "FULL_FDBK"

crop_prices = rbind(agcom0, agcom)

write.csv(crop_prices, paste0(oname, ".csv"))

crop_prices_baseline<- crop_prices %>% filter(scenario=="CONTROL") %>% rename(value_baseline=value) %>% dplyr::select(-scenario)

# these need to be in order
plot_crops = c("BioenergyCrop","Corn","OilPalm","Rice","Soybean","SugarCrop","Wheat")

crop_prices %>% 
  filter(scenario != "CONTROL") %>% 
  filter(year>2015) %>% 
  left_join(crop_prices_baseline) %>% 
  filter(sector %in% c(plot_crops, "biomass", "Biomass")) %>% 
  mutate(diff=((value- value_baseline)/value_baseline)*100)->crop_marginal_cost

crop_marginal_cost$sector[crop_marginal_cost$sector=="biomass"] = "Biomass"
crop_marginal_cost$sector[crop_marginal_cost$sector=="Biomass"] = "BioenergyCrop"

g <- ggplot(data=crop_marginal_cost %>% filter(year %in% c(2030,2050,2075,2100)) , aes(x=sector,y=diff))+
  geom_boxplot(aes(color=scenario), color="black")+
  scale_shape_manual(values= c(0))+
  scale_alpha_manual(values = c(0.1, 1)) +
  facet_wrap(~year, labeller=labeller(year = c('2030'="a) 2030", '2050'="b) 2050", '2075'="c) 2075", '2100'="d) 2100"))) +
  xlab("Commodity")+
  ylab("% difference between FULL_FDBK and CONTROL")

g+scheme_basic

ggsave( paste0(oname,'_relative_change_subbox.png'),width = 10.5, height = 10)

# write these boxplot stats to a file - format it for easier access
# colour is the scenario, PANEL is the year, group denotes the scenario and commodity
plot_df = data.frame(ggplot_build(g)[[1]])

# panel ids are the indices of panel_names
panel_names = c(2030, 2050, 2075, 2100)
# map colour to scenarios
scenario_names = c("FULL_FDBK")
# map groups to scenario and year
group_vals = plot_df$group
group_scen_names = rep(c("FULL_FDBK"), 7)
group_commodity = c(rep(plot_crops[1],1), rep(plot_crops[2],1), rep(plot_crops[3],1), rep(plot_crops[4],1), rep(plot_crops[5],1), rep(plot_crops[6],1), rep(plot_crops[7],1))

plot_df = plot_df[,c("colour", "PANEL", "group", "ymin_final", "ymin", "lower", "middle", "upper", "ymax", "ymax_final")]
   
# ymin/max_final is the true min/max; ymin/max is the 1.5*IQR value; lower/upper is the 25%/75$% qurtile; middle is the median
names(plot_df) <- c("Scenario", "Year", "commodity", "min", "low_1.5*IQR", "Q25", "median", "Q75", "hi_1.5*IQR", "max")
plot_df$Year <- as.character(plot_df$Year)
for (i in 1:length(group_vals)){
  	plot_df$Scenario[which(plot_df$commodity == group_vals[i])] = group_scen_names[i]
   	plot_df$commodity[which(plot_df$commodity == group_vals[i])] = group_commodity[i]
}
for (i in 1:length(panel_names)){
   	plot_df$Year[which(plot_df$Year == i)] = panel_names[i]
}
plot_df = plot_df[order(plot_df$Year, plot_df$commodity),]

# write the region box plot summary
write.csv(plot_df, paste0(oname,'_relative_change_subbox.csv'))


### make line plots of all crops in each region just to visualize

crop_prices %>% 
  filter(scenario != "CONTROL") %>% 
  filter(year>2015) %>% 
  left_join(crop_prices_baseline) %>% 
  mutate(diff=((value- value_baseline)/value_baseline)*100)->all_crop_marginal_cost

all_crop_marginal_cost$sector[all_crop_marginal_cost $sector=="biomass"] = "Biomass"
all_crop_marginal_cost$sector[all_crop_marginal_cost $sector=="Biomass"] = "BioenergyCrop"

write.csv(all_crop_marginal_cost, paste0(oname, "_relative_change.csv"))

pdf(file = paste0(oname, "_relative_change.pdf"), width = 10.5, height = 10)

for (r in unique(all_crop_marginal_cost$region)) {

	g <- ggplot(data= all_crop_marginal_cost %>% filter(year >= 2015, region == r) , aes(x=year,y=diff))+
	  geom_line(aes(color=sector))+
	  xlab("Year")+
	  ylab("% difference between FULL_FDBK and CONTROL")+
	  ggtitle(paste(r, "change in crop prices"))
  
	print(g+scheme_basic)

} # end for loop over region

dev.off()


#### now do the production-allocation-price correlations


plot_crops = c("BioenergyCrop", "Corn", "FiberCrop", "FodderGrass", "FodderHerb", "Fruits", "Legumes", "MiscCrop", "NutsSeeds", "OilCrop", "OilPalm",
               "OtherGrain", "Rice", "RootTuber", "Soybean", "SugarCrop", "Vegetables", "Wheat")

gcam_crop_names = c("biomassGrass", "biomassTree", "CornC4", "FiberCrop", "FodderGrass", "FodderHerb", "FodderHerbC4", "Fruits", "FruitsTree",
					"Legumes", "MiscCrop", "MiscCropTree", "NutsSeeds", "NutsSeedsTree", "OilCrop", "OilCropTree", "OilPalmTree", "OtherGrain", "OtherGrainC4",
					"Rice", "RootTuber", "Soybean", "SugarCrop", "SugarCropC4", "Vegetables", "Wheat")

plot_map_names = c("BioenergyCrop", "BioenergyCrop", "Corn", "FiberCrop", "FodderGrass", "FodderHerb", "FodderHerb", "Fruits", "Fruits",
                   "Legumes", "MiscCrop", "MiscCrop", "NutsSeeds", "NutsSeeds", "OilCrop", "OilCrop", "OilPalm", "OtherGrain", "OtherGrain",
                   "Rice", "RootTuber", "Soybean", "SugarCrop", "SugarCrop", "Vegetables", "Wheat")

ag_production <- getQuery(prj0, "ag production by crop type") %>%
  mutate(scenario="CONTROL") %>% bind_rows(getQuery(prj, "ag production by crop type") %>%mutate(scenario="FULL_FDBK")) %>%
  filter(year >= 2015 & sector %in% c(plot_crops, "biomass"))
ag_production$sector[ag_production$sector=="biomass"] = "BioenergyCrop"
ag_production$output = NULL
names(ag_production)[which(names(ag_production) == "value")] = "prod_Mt_EJ"
ag_production$Units = NULL

ag_prices <- getQuery(prj0, "ag commodity prices") %>%
  mutate(scenario="CONTROL") %>% bind_rows(getQuery(prj, "ag commodity prices") %>%mutate(scenario="FULL_FDBK")) %>%
  filter(year >= 2015 & sector %in% c(plot_crops, "biomass"))
ag_prices$sector[ag_prices$sector=="biomass"] = "BioenergyCrop"
names(ag_prices)[which(names(ag_prices) == "value")] = "price_1975_usd_perkg_perGJ"
ag_prices$Units = NULL

ag_allocation <- getQuery(prj0, "detailed land allocation") %>%
  mutate(scenario="CONTROL") %>% bind_rows(getQuery(prj, "detailed land allocation") %>%mutate(scenario="FULL_FDBK")) %>%
  filter(year >= 2015)

# separate the land type and basin names and determine the glu code
ag_allocation$LT_GCAM = sapply(strsplit(ag_allocation$landleaf,"_"),"[[",1)
ag_allocation$gcam_glu_abbr = sapply(strsplit(ag_allocation$landleaf,"_"),"[[",2)
ag_allocation$water = NA
ag_allocation$water[ag_allocation$LT_GCAM %in% gcam_crop_names] = 
			sapply(strsplit(ag_allocation$landleaf[ag_allocation$LT_GCAM %in% gcam_crop_names],"_"),"[[",3)
ag_allocation$fert = NA
ag_allocation$fert[ag_allocation$LT_GCAM %in% gcam_crop_names] =
			sapply(strsplit(ag_allocation$landleaf[ag_allocation$LT_GCAM %in% gcam_crop_names],"_"),"[[",4)
ag_allocation$landleaf = NULL
filter(ag_allocation, LT_GCAM %in% gcam_crop_names)
# aggregate to region
ag_allocation = aggregate(value ~ Units + scenario + region + year + LT_GCAM, data = ag_allocation, FUN=sum, na.rm=TRUE)
# aggregate to production/price types
ag_allocation$sector = NA
for (i in 1:length(gcam_crop_names)) {
  ag_allocation$sector[ag_allocation$LT_GCAM == gcam_crop_names[i]] = plot_map_names[i]
}
ag_allocation = aggregate(value ~ Units + scenario + region + sector + year, data = ag_allocation, FUN=sum, na.rm=TRUE)
names(ag_allocation)[which(names(ag_allocation) == "value")] = "alloc_thous_km2"
ag_allocation$Units = NULL

# merge the three data streams into one df
plot_df = merge(ag_allocation, ag_production, by=c("scenario", "region", "sector", "year"), all.x = TRUE)
plot_df = merge(plot_df, ag_prices, by=c("scenario", "region", "sector", "year"), all.x = TRUE)

plot_df[is.na(plot_df)] <- 0

# add global region
# only price needs weighting; Mt to kg and GJ to EJ are both *1e9
plot_df$tot_price = plot_df$price_1975_usd_perkg_perGJ * plot_df$prod_Mt_EJ * 1e9
global = aggregate(cbind(alloc_thous_km2, prod_Mt_EJ, tot_price) ~ scenario + sector + year, data = plot_df, FUN = sum, na.rm = TRUE)
global$price_1975_usd_perkg_perGJ = global$tot_price / 1e9
global$region = "Global"
global = global[,names(plot_df)]
all_df = rbind(plot_df, global)



# first do the correlation calcs by scenario and sector and region
# include data from all years in each group
# include the globe as a region?
plot_df = all_df

for(s in c(unique(plot_df$scenario))){
  
  # write separate file for each for each scenario
  # include each region

  pdf(paste0(oname, "_correlations_", s, ".pdf"), width=8.5, height=4)
  
  for(r in c(unique(plot_df$region))){
  for(c in c(unique(plot_df$sector))){

    # production to allocation
    prod_alloc_error = FALSE
    tr <- tryCatch( {lmout=lm(prod_Mt_EJ ~ alloc_thous_km2, data = filter(plot_df, sector == c & scenario == s & region == r))},
                    error = function(e) {prod_alloc_error <<- TRUE})
    if(!prod_alloc_error){
      slm = summary(lmout)
      r2 = slm$r.squared
      r2adj = slm$adj.r.squared
      pval = slm$coefficients[8]
      min_x = min(plot_df$alloc_thous_km2, na.rm=TRUE)
      max_x = max(plot_df$alloc_thous_km2, na.rm=TRUE)
      coef=slm$coefficients
      lint = coef[1]
      lslp = coef[2]
      # form: y = a + b*x
      plot_df$prod_alloc_a[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(lint,2)
      plot_df$prod_alloc_b[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(lslp,3)
      plot_df$prod_alloc_r2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(r2,2)
      plot_df$prod_alloc_slp_pval[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(pval,5)
      plot_df$prod_alloc_y[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = 
        predict(lmout, list(alloc_thous_km2 = plot_df$alloc_thous_km2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r]))
    }

    # price to production
    price_prod_error = FALSE
    tr <- tryCatch( {lmout=lm(price_1975_usd_perkg_perGJ ~ prod_Mt_EJ, data = filter(plot_df, sector == c & scenario == s & region == r))},
                    error = function(e) {price_prod_error <<- TRUE})
    if(!price_prod_error){
      slm = summary(lmout)
      r2 = slm$r.squared
      r2adj = slm$adj.r.squared
      pval = slm$coefficients[8]
      min_x = min(plot_df$alloc_thous_km2, na.rm=TRUE)
      max_x = max(plot_df$alloc_thous_km2, na.rm=TRUE)
      coef=slm$coefficients
      lint = coef[1]
      lslp = coef[2]
      # form: y = a + b*x
      plot_df$price_prod_a[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(lint,2)
      plot_df$price_prod_b[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(lslp,5)
      plot_df$price_prod_r2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(r2,2)
      plot_df$price_prod_slp_pval[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(pval,5)
      plot_df$price_prod_y[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = 
        predict(lmout, list(prod_Mt_EJ = plot_df$prod_Mt_EJ[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r]))
    }
  
    # price to allocation
    price_alloc_error = FALSE
    tr <- tryCatch( {lmout=lm(price_1975_usd_perkg_perGJ ~ alloc_thous_km2, data = filter(plot_df, sector == c & scenario == s & region == r))},
                    error = function(e) {price_alloc_error <<- TRUE})
    if(!price_alloc_error){
      slm = summary(lmout)
      r2 = slm$r.squared
      r2adj = slm$adj.r.squared
      pval = slm$coefficients[8]
      min_x = min(plot_df$alloc_thous_km2, na.rm=TRUE)
      max_x = max(plot_df$alloc_thous_km2, na.rm=TRUE)
      coef=slm$coefficients
      lint = coef[1]
      lslp = coef[2]
      # form: y = a + b*x
      plot_df$price_alloc_a[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(lint,2)
      plot_df$price_alloc_b[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(lslp,5)
      plot_df$price_alloc_r2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(r2,2)
      plot_df$price_alloc_slp_pval[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = round(pval,5)
      plot_df$price_alloc_y[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r] = 
        predict(lmout, list(alloc_thous_km2 = plot_df$alloc_thous_km2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r]))
  }
  
    ## now make plots
    
    # production to allocation
  
    pstitle=paste("r2 =", plot_df$prod_alloc_r2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1], "pval =", 
                plot_df$prod_alloc_slp_pval[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1])
    modtext = paste0("\n", "y = ", plot_df$prod_alloc_a[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1], " + ", 
                   plot_df$prod_alloc_b[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1]," * x")
    pstitle = paste(pstitle, modtext)
    
    if(!prod_alloc_error){
      g <- ggplot(data=plot_df %>% filter(sector==c & scenario == s & plot_df$region == r) ) +
        geom_point(data=plot_df %>% filter(sector==c & scenario == s & plot_df$region == r), aes(x=alloc_thous_km2, y=prod_Mt_EJ, color=year)) +
        ylab("Mt (EJ for BioeneryCrop)") +
        xlab("thousand km^2") +
		    labs(title = paste(r, c, "production vs allocation"), subtitle = pstitle, color = "year") +
	  	  geom_line(aes(x=alloc_thous_km2, y = prod_alloc_y))

      g+scheme_basic
      print(g)
    }
    
    # price to production
  
    pstitle=paste("r2 =", plot_df$price_prod_r2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1], "pval =", 
                plot_df$price_prod_slp_pval[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1])
    modtext = paste0("\n", "y = ", plot_df$price_prod_a[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1], " + ", 
                   plot_df$price_prod_b[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1]," * x")
    pstitle = paste(pstitle, modtext)
    
    if(!price_prod_error){
      g <- ggplot(data=plot_df %>% filter(sector==c & scenario == s & plot_df$region == r) ) +
        geom_point(data=plot_df %>% filter(sector==c & scenario == s & plot_df$region == r), aes(x=prod_Mt_EJ, y=price_1975_usd_perkg_perGJ, color=year)) +
        ylab("1975$ per kg (per GJ for BioeneryCrop)") +
        xlab("Mt (EJ for BioeneryCrop)") +
		    labs(title = paste(r, c, "price vs production"), subtitle = pstitle, color = "year") +
		    geom_line(aes(x=prod_Mt_EJ, y = price_prod_y))

      g+scheme_basic
      print(g)
    }
    
    # price to allocation
  
    pstitle=paste("r2 =", plot_df$price_alloc_r2[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1], "pval =", 
                plot_df$price_alloc_slp_pval[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1])
    modtext = paste0("\n", "y = ", plot_df$price_alloc_a[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1], " + ", 
                   plot_df$price_alloc_b[plot_df$sector == c & plot_df$scenario == s & plot_df$region == r][1]," * x")
    pstitle = paste(pstitle, modtext)
    
    if(!price_alloc_error){
      g <- ggplot(data=plot_df %>% filter(sector==c & scenario == s & plot_df$region == r) ) +
        geom_point(data=plot_df %>% filter(sector==c & scenario == s & plot_df$region == r), aes(x=alloc_thous_km2, y=price_1975_usd_perkg_perGJ, color=year)) +
        ylab("1975$ per kg (per GJ for BioeneryCrop)") +
        xlab("thousand km^2") +
		    labs(title = paste(r, c, "price vs allocation"), subtitle = pstitle, color = "year") +
		    geom_line(aes(x=alloc_thous_km2, y = price_alloc_y))

      g+scheme_basic
      print(g)
    }
  
  } # end c loop over sectors

    
  } # end r loop over regions
  
  dev.off()

  
} # end s loop over scenarios

write.csv(plot_df, paste0(oname, "_correlations.csv"))