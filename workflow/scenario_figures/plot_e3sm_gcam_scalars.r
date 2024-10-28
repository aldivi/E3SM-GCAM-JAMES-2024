# plot_e3sm_gcam_scalars.r

# make timeseries plots of above and below ground scalers
# make maps of above and below ground scalers

# do this at different levels of regional and plant aggregation (area weighted)
#    globe or region or regionXbasin
#    gcam land types or aggregated land types

# the input files are the scaler_yyyy.csv files output by the iac
# the above scalar is the fourth column and the below scalar is the fifth column
# these are the values that are sent to gcam


library(raster)
library(rasterVis)
library(sf)     # this replaces rgdal
library(RColorBrewer)
library(sp)
#library(plyr)
library(ggplot2)
library(viridis)
library(dplyr)
# these are to get the gcam data directly from the output data file
require(devtools)
library(rgcam)

# with feedbacks
indir = "./"
outdir = "./scenario_outputs/"
# gcam output data file with detailed land allocation and ag yield scalar data (with feedbacks)
gcam_output_name = paste0(indir,"20240730_SSP245_ZATM.dat")

# without feedbacks
#indir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_20240718_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3/no_feedbacks/"
#outdir = "~/projects/e3sm/iesmv2/test_diagnostics/eva_20240718_SSP245_ZATM_BGC_ne30pg2_f09_oEC60to30v3/no_feedbacks/"
# gcam output data file with detailed land allocation and ag yield scalar data (with feedbacks)
#gcam_output_name = paste0(indir,"20240718_SSP245_ZATM_without_feedbacks.dat")

map_year = 2050

start_year = 2025
end_year = 2100



gcam_region_info_file = "./iac_region_glu_codes.csv"

gcam_boundaries_dir = "./gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/"
gcam_reg_glu_shape_fname = "reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp"
gcam_region_shape_fname = "region_boundaries_moirai_combined_3p1_0p5arcmin.shp"

# land type groups for aggregation and plotting

# enter some detailed ones, then combine them for the generic types
# these match the land allocation plot script
gcam_forest_names = c("Forest", "ProtectedUnmanagedForest", "UnmanagedForest")
gcam_pasture_names = c("Pasture", "UnmanagedPasture", "ProtectedUnmanagedPasture")
gcam_grassland_names = c("Grassland", "ProtectedGrassland")
gcam_crop_names = c("CornC4", "FiberCrop", "FodderGrass", "FodderHerb", "FodderHerbC4", "Fruits", "FruitsTree", "Legumes", "MiscCrop", "MiscCropC4",
	"MiscCropTree", "NutsSeeds", "NutsSeedsTree", "OilCrop", "OilCropTree", "OilPalmTree", "OtherGrain", "OtherGrainC4", "Rice", "RootTuber", "Soybean",
	"SugarCrop", "SugarCropC4", "Vegetables", "Wheat", "biomassGrass", "biomassTree")
gcam_other_arable_name = "OtherArableLand"		# this is considered a managed cropland type
gcam_shrubland_names = c("Shrubland", "ProtectedShrubland")
# these are constant
gcam_other_names = c("RockIceDesert", "Tundra")
gcam_urban_name = "UrbanLand"

crops = c(gcam_crop_names, gcam_other_arable_name)
static = c(gcam_other_names, gcam_urban_name)
forest = c(gcam_forest_names)
pasture = c(gcam_pasture_names)
shrubland = c(gcam_shrubland_names)
grassland = c(gcam_grassland_names)


plot_e3sm_gcam_scalars <- function(indir = "./",
                                   outdir = "./scenario_outputs/",
                                   gcam_region_info_file = "./iac_region_glu_codes.csv",
                                   gcam_boundaries_dir = "./gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/",
								   gcam_reg_glu_shape_fname = "reg_glu_boundaries_moirai_combined_3p1_0p5arcmin.shp",
								   gcam_region_shape_fname = "region_boundaries_moirai_combined_3p1_0p5arcmin.shp",
								   gcam_output_name = paste0(indir,"20240730_SSP245_ZATM.dat"),
								   map_year = 2050,
                                   start_year = 2025,
                                   end_year = 2100) {
                                   	
	cat("started plot_e3sm_gcam_scalars.r", date())                                	
        
    # hardcode these output file names          
    regglu_map_outname = "regXglu_scalar_maps"
	region_map_outname = "region_scalar_maps"
        
    # create the output directory
    dir.create(outdir)
                                   	
    # read the region and basin info
    region_info = read.csv(gcam_region_info_file, stringsAsFactors=FALSE)
    
    # load the output data file
    prj <- loadProject(gcam_output_name)
    
    # process the area data
    
    # get the GCAM areas with default protected area
    gcam_area_in <-getQuery(prj, query="detailed land allocation")
    gcam_area = as.data.frame(gcam_area_in)
    gcam_area$scenario = NULL
    gcam_area$Units = NULL
    
    # split the basin, land type, and techs
    # then aggregate the techs for each type
    tl = strsplit(gcam_area$landleaf, "_")
    gcam_area$basin_abr = sapply( tl, "[[", 2 )
    gcam_area$land_type = sapply( tl, "[[", 1 )
    gcam_area$water_tech = NA
    gcam_area$fert_tech = NA
    wt_inds = which(sapply( tl, length ) > 2)
    ft_inds = which(sapply( tl, length ) > 3)
    gcam_area$water_tech[wt_inds] = sapply( tl[wt_inds], "[[", 3 )
    gcam_area$fert_tech[ft_inds] = sapply( tl[ft_inds], "[[", 4 )
    gcam_area$landleaf = NULL
    
    gcam_area = aggregate(value ~ year + region + basin_abr + land_type, data = gcam_area, FUN = sum)
    colnames(gcam_area)[ncol(gcam_area)] = "area_thous_sqkm"
    
    # not any more because the spaces have been added
    # need to remove whitespace from gcam_area region
    # gcam_area$region = gsub("\\s", "", gcam_area$region)
    
    # read in the GCAM ag vegetation scalars and split technology to match the EHC files
    gcam_ag_scalars_in <-getQuery(prj, query="yield-scaler")
    gcam_ag_scalars = as.data.frame(gcam_ag_scalars_in)
    gcam_ag_scalars$Units = NULL
    gcam_ag_scalars$scenario = NULL
    gcam_ag_scalars$sector = NULL
    gcam_ag_scalars$subsector = NULL
    colnames(gcam_ag_scalars) = c("Region", "Year", "technology", "AgYield_scalar")
    tl = strsplit(gcam_ag_scalars$technology, "_")
    land_type = sapply( tl, "[[", 1 )
    basin_abr = sapply( tl, "[[", 2 )
    gcam_ag_scalars$LandType_Basin = paste0(land_type, "_", basin_abr)
    gcam_ag_scalars$technology = NULL
    gcam_ag_scalars = unique(gcam_ag_scalars)
    
	# the scalars are available at 5-year intervals
	num_years = (end_year - start_year)/5
	year_array = start_year + 5*c(0:num_years)
    
    # read in the EHC scalars and append to one data frame
    all_df = NULL
    for(y in year_array) {
    	cat("processing year",y,"\n")
    		
        # read EHC scalar files
    	scalar_df = read.csv(paste0(indir,"scalars_",y,".csv"), stringsAsFactors=FALSE)
    		
    	# include gcam otuput for comparison with the EHC file
		temp_df = merge(scalar_df, gcam_ag_scalars[gcam_ag_scalars$Year == y,], all.x=TRUE)
    		
    	#temp_df = scalar_df
    	#temp_df = scalar_df[c(1:10,1000:1010,8000:8010),]
    		
    	# split the basin and land type
    	tl = strsplit(temp_df$LandType_Basin,"_")
    	temp_df$basin_abr = sapply( tl, "[[", 2 )
    	temp_df$land_type = sapply( tl, "[[", 1 )
    		
    	# need to do the area merging by year
    	
    	# if there is area but no scalar, need to set scalar to 1
        # if scalar but no area is fine - just drop it
        
        temp_df = merge(temp_df, gcam_area[gcam_area$year == y,],
        	by.x = c("Year", "Region", "basin_abr", "land_type"), by.y = c("year", "region", "basin_abr", "land_type"), all.x=TRUE)
       
       #### delete, i think 
        # need different line when area file has multiple years
        #main_df = merge(gcam_area[gcam_area$year == y,], all_df,
        #                by.x = c("year", "region", "basin", "land_type"), by.y = c("year", "gcam_reg_name", "basin_abr", "land_type"), all.x=TRUE)
        # need gcam 6 allocation file!
        #temp_df = merge(gcam_area[,c("region", "basin", "land_type", "area_thous_sqkm")], temp_df,
        #                by = c("region", "basin", "land_type"), all.x=TRUE)
        # for now just keep the land types that match
        #temp_df = merge(gcam_area[,c("region", "basin", "land_type", "area_thous_sqkm")], temp_df,
        #               by = c("region", "basin", "land_type"))
       #### 
        
        # drop NA area records
        temp_df = temp_df[!is.na(temp_df$area_thous_sqkm),]
        # drop zero area records
        temp_df = temp_df[temp_df$area_thous_sqkm > 0,]

		# there shouldn't be any missing years or scalars because the files have been trimmed to valid records
		temp_df$Year[is.na(temp_df$Year)] = y
        # set missing scalars to 1
        temp_df$Vegetation[is.na(temp_df$Vegetaiton)] = 1
        temp_df$Soil[is.na(temp_df$Soil)] = 1

		# check gcam output against EHC - so far these are identical
		temp_df$ag_diff = temp_df$AgYield_scalar - temp_df$Vegetation
		diff_inds = which(temp_df$ag_diff != 0 & !is.na(temp_df$AgYield_scalar))
		#if(length(diff_inds) > 0) {
		#	cat("Ag yield scalars don't match EHC veg scalars in year", y, "num diff values is", length(diff_inds), "\n")
		#	stop("Error!")
		#}
		temp_df$ag_diff = NULL

        # determine the raster id for this regionXbasin
    	temp_df = merge(temp_df, region_info[,c("gcam_reg_name", "basin_abr", "gcam_reg_code", "gcam_basin_code")],
    				by.x = c("Region", "basin_abr"), by.y = c("gcam_reg_name","basin_abr"), all.x = TRUE)			
        na_rinds = which(is.na(temp_df$gcam_reg_code))
        na_binds = which(is.na(temp_df$gcam_basin_code))
        if (length(na_rinds)>0 | length(na_binds) > 0) {
        		cat("regions or basins in mapping file not matching those in area/scalar files\n")
        		stop("Error!")	
        }
        temp_df$raster_id = temp_df$gcam_reg_code * 1000 + temp_df$gcam_basin_code
        
        # get some useful values for aggregating
        temp_df$above_weighted = temp_df$Vegetation * temp_df$area_thous_sqkm
        temp_df$below_weighted = temp_df$Soil * temp_df$area_thous_sqkm
        
    	# add this year records to all_df
    	all_df = rbind(all_df,temp_df)
    } # end for y loop over years
    
    # make these column names consistent with the following data frames
    colnames(all_df) = c("region", "basin_abr", "year", "land_type", "LandType_Basin", "above_scalar", "below_scalar", "AgYield_scalar", "area_thous_sqkm",
    	"gcam_reg_code", "gcam_basin_code", "raster_id", "above_weighted", "below_weighted")
    
	# aggregate each land type to region
	
	# land type area per region
    region_df = aggregate(cbind(area_thous_sqkm, above_weighted, below_weighted) ~ year + region + land_type + gcam_reg_code, data = all_df, FUN = sum)
    colnames(region_df) = c("year", "region", "land_type", "gcam_reg_code", "reg_lt_area_thous_sqkm", "above_wsum", "below_wsum")
    # if there isn't any area can drop it
    region_df = region_df[region_df$reg_lt_area_thous_sqkm > 0,]
    region_df$above_scalar = region_df$above_wsum / region_df$reg_lt_area_thous_sqkm
    region_df$below_scalar = region_df$below_wsum / region_df$reg_lt_area_thous_sqkm
	
	# aggregate each land type to globe
	globe_df = aggregate(cbind(area_thous_sqkm, above_weighted, below_weighted) ~ year + land_type, data = all_df, FUN = sum)
	colnames(globe_df) = c("year", "land_type", "glb_lt_area_thous_sqkm", "above_wsum", "below_wsum")
	# restore the region column to label it as globe
	globe_df$region = "Globe"
	globe_df = globe_df[,c("year", "region", "land_type", "glb_lt_area_thous_sqkm", "above_wsum", "below_wsum")]
    # if there isn't any area can drop it
    globe_df = globe_df[globe_df$glb_lt_area_thous_sqkm > 0,]
    globe_df$above_scalar = globe_df$above_wsum / globe_df$glb_lt_area_thous_sqkm
    globe_df$below_scalar = globe_df$below_wsum / globe_df$glb_lt_area_thous_sqkm
	
	# aggregate to assigned types: crop, forest, grassland, shrubland, pasture
	# use only crops here because nothing is being produced on other arable land
	# get the max and min as well as the average
	
	# aggregate to assigned types within regions
	temp_df = region_df
	temp_df$agglt = NA
	temp_df$agglt[temp_df$land_type %in% gcam_crop_names] = "Crop"
	temp_df$agglt[temp_df$land_type %in% gcam_forest_names] = "Forest"
	temp_df$agglt[temp_df$land_type %in% gcam_grassland_names] = "Grassland"
	temp_df$agglt[temp_df$land_type %in% gcam_shrubland_names] = "Shrubland"
	temp_df$agglt[temp_df$land_type %in% gcam_pasture_names] = "Pasture"
	# weighted average
	region_agglt_df = aggregate(cbind(reg_lt_area_thous_sqkm, above_wsum, below_wsum) ~ year + region + agglt + gcam_reg_code, data = temp_df, FUN = sum)
    colnames(region_agglt_df) = c("year", "region", "agglt", "gcam_reg_code", "reg_agglt_area_thous_sqkm", "above_wsum", "below_wsum")
    # if there isn't any area can drop it
    region_agglt_df = region_agglt_df[region_agglt_df$reg_agglt_area_thous_sqkm > 0,]
    region_agglt_df$above_scalar = region_agglt_df$above_wsum / region_agglt_df$reg_agglt_area_thous_sqkm
    region_agglt_df$below_scalar = region_agglt_df$below_wsum / region_agglt_df$reg_agglt_area_thous_sqkm
    # max and min
    temp_df = do.call(data.frame, aggregate(cbind(above_scalar, below_scalar) ~ year + region + agglt + gcam_reg_code, data = temp_df,
    										FUN = function(x) c(mn=min(x), mx=max(x)) ) )
    colnames(temp_df) = c("year", "region", "agglt", "gcam_reg_code", "above_scalar_min", "above_scalar_max", "below_scalar_min", "below_scalar_max")
    region_agglt_df = merge(region_agglt_df, temp_df, by = c("year", "region", "agglt", "gcam_reg_code"), all.x=TRUE)

	
	# aggregate to assigned types for the globe
	temp_df = globe_df
	temp_df$agglt = NA
	temp_df$agglt[temp_df$land_type %in% gcam_crop_names] = "Crop"
	temp_df$agglt[temp_df$land_type %in% gcam_forest_names] = "Forest"
	temp_df$agglt[temp_df$land_type %in% gcam_grassland_names] = "Grassland"
	temp_df$agglt[temp_df$land_type %in% gcam_shrubland_names] = "Shrubland"
	temp_df$agglt[temp_df$land_type %in% gcam_pasture_names] = "Pasture"
	globe_agglt_df = aggregate(cbind(glb_lt_area_thous_sqkm, above_wsum, below_wsum) ~ year + region + agglt, data = temp_df, FUN = sum)
	colnames(globe_agglt_df) = c("year", "region", "agglt", "glb_agglt_area_thous_sqkm", "above_wsum", "below_wsum")
    # if there isn't any area can drop it
    globe_agglt_df = globe_agglt_df[globe_agglt_df$glb_agglt_area_thous_sqkm > 0,]
    globe_agglt_df$above_scalar = globe_agglt_df$above_wsum / globe_agglt_df$glb_agglt_area_thous_sqkm
    globe_agglt_df$below_scalar = globe_agglt_df$below_wsum / globe_agglt_df$glb_agglt_area_thous_sqkm
    # max and min
    temp_df = do.call(data.frame, aggregate(cbind(above_scalar, below_scalar) ~ year + region + agglt, data = temp_df,
    										FUN = function(x) c(mn=min(x), mx=max(x)) ) )
    colnames(temp_df) = c("year", "region", "agglt", "above_scalar_min", "above_scalar_max", "below_scalar_min", "below_scalar_max")
    globe_agglt_df = merge(globe_agglt_df, temp_df, by = c("year", "region", "agglt"), all.x=TRUE)
	
	# aggregate to assigned types within subregion
    temp_df = all_df
	temp_df$agglt = NA
	temp_df$agglt[temp_df$land_type %in% gcam_crop_names] = "Crop"
	temp_df$agglt[temp_df$land_type %in% gcam_forest_names] = "Forest"
	temp_df$agglt[temp_df$land_type %in% gcam_grassland_names] = "Grassland"
	temp_df$agglt[temp_df$land_type %in% gcam_shrubland_names] = "Shrubland"
	temp_df$agglt[temp_df$land_type %in% gcam_pasture_names] = "Pasture"
	# weighted average
	all_agglt_df = aggregate(cbind(area_thous_sqkm, above_weighted, below_weighted) ~ year + region + basin_abr + agglt + gcam_reg_code + gcam_basin_code, data = temp_df, FUN = sum)
    colnames(all_agglt_df) = c("year", "region", "basin_abr", "agglt", "gcam_reg_code", "gcam_basin_code", "agglt_area_thous_sqkm", "above_wsum", "below_wsum")
    # if there isn't any area can drop it
    all_agglt_df = all_agglt_df[all_agglt_df$agglt_area_thous_sqkm > 0,]
    all_agglt_df$above_scalar = all_agglt_df$above_wsum / all_agglt_df$agglt_area_thous_sqkm
    all_agglt_df$below_scalar = all_agglt_df$below_wsum / all_agglt_df$agglt_area_thous_sqkm
    # max and min
    temp_df = do.call(data.frame, aggregate(cbind(above_scalar, below_scalar) ~ year + region + basin_abr + agglt + gcam_reg_code + gcam_basin_code, data = temp_df,
    										FUN = function(x) c(mn=min(x), mx=max(x)) ) )
    colnames(temp_df) = c("year", "region", "basin_abr", "agglt", "gcam_reg_code", "gcam_basin_code", "above_scalar_min", "above_scalar_max", "below_scalar_min", "below_scalar_max")
    all_agglt_df = merge(all_agglt_df, temp_df, by = c("year", "region", "basin_abr", "agglt", "gcam_reg_code", "gcam_basin_code"), all.x=TRUE)
    
    
    
	# plot scalar time series of multiple land types on one or more plots
    
    # need some info to store plots and write them later into separate files

	regions = unique(region_df$region)
	num_reg = length(regions)
	lunits = unique(all_df[,c("region", "basin_abr")])
	num_lunits = nrow(lunits)
	basin_counts = aggregate(basin_abr ~ region, data = lunits, FUN=length)
	colnames(basin_counts) = c("reigon", "num_basins")
	
	FIGURE_DIMS <- list(dpi=300, width=2560/300, height=1440/300)
	theme_set(theme_bw())
	
	#########################################################################
	# global plots
	
	outname = paste0(outdir,"globe_scalar_plots.pdf")
	
	num_globe_plots = 6
	globe_plots = vector(num_globe_plots,mode="list")
	
	# crops above
	title = "Global crop vegetation scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	plot_df = globe_df[globe_df$land_type%in%crops,]
	p <- ( ggplot(plot_df, aes(year, above_scalar, color=land_type))
			+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= land_type), size = 1.5)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	globe_plots[[1]] <- recordPlot()
	
	# crops below
	title = "Global crop soil scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	plot_df = globe_df[globe_df$land_type%in%crops,]
	p <- ( ggplot(plot_df, aes(year, below_scalar, color=land_type))
			+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= land_type), size = 1.5)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	globe_plots[[2]] <- recordPlot()
	
	# non-crops and non-static above
	title = "Global non-crop vegetation scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	plot_df = globe_df[globe_df$land_type%in%forest | globe_df$land_type%in%pasture | globe_df$land_type%in%shrubland | globe_df$land_type%in%grassland,]
	p <- ( ggplot(plot_df, aes(year, above_scalar, color=land_type))
			+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= land_type), size = 1.5)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	globe_plots[[3]] <- recordPlot()
	
	# non-crops and non-static below
	title = "Global non-crop soil scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	plot_df = globe_df[globe_df$land_type%in%forest | globe_df$land_type%in%pasture | globe_df$land_type%in%shrubland | globe_df$land_type%in%grassland,]
	p <- ( ggplot(plot_df, aes(year, below_scalar, color=land_type))
			+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= land_type), size = 1.5)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	globe_plots[[4]] <- recordPlot()


	# global average scalar values for aggregated types
	
	# vegetation
	title = "Global average vegetation scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	plot_df = globe_agglt_df
	p <- ( ggplot(plot_df, aes(year, above_scalar, color=agglt))
			+ scale_shape_manual(values=1:length(unique(plot_df$agglt)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= agglt), size = 1.5)
			+ geom_ribbon(aes(ymin = above_scalar_min, ymax = above_scalar_max, fill=agglt, linetype=NA), alpha=0.1, show.legend=FALSE)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			+labs(color="land type", shape="land type")
			+scale_x_continuous(breaks=unique(plot_df$year), labels=as.character(unique(plot_df$year)))
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	globe_plots[[5]] <- recordPlot()

	# soil
	title = "Global average soil scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	plot_df = globe_agglt_df
	p <- ( ggplot(plot_df, aes(year, below_scalar, color=agglt))
			+ scale_shape_manual(values=1:length(unique(plot_df$agglt)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= agglt), size = 1.5)
			+ geom_ribbon(aes(ymin = below_scalar_min, ymax = below_scalar_max, fill=agglt, linetype=NA), alpha=0.1, show.legend=FALSE)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			+labs(color="land type", shape="land type")
			+scale_x_continuous(breaks=unique(plot_df$year), labels=as.character(unique(plot_df$year)))
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	globe_plots[[6]] <- recordPlot()
	


	
	# write the plots to pdf
	pdf(outname, width=6.5,height=3.25)
	
	for (i in 1:num_globe_plots) {
		replayPlot(globe_plots[[i]])
	}
	
	dev.off()
		
    ####### global timeseries scalar plots for the paper
    
    glb_above = globe_agglt_df[,c("year", "region", "agglt", "above_scalar", "above_scalar_min", "above_scalar_max")]
    names(glb_above) = c("year", "region", "agglt", "scalar", "scalar_min", "scalar_max")
    glb_above$type = "a) Vegetation"
    glb_below = globe_agglt_df[,c("year", "region", "agglt", "below_scalar", "below_scalar_min", "below_scalar_max")]
    names(glb_below) = c("year", "region", "agglt", "scalar", "scalar_min", "scalar_max")
    glb_below$type = "b) Soil"
    plot_df = rbind(glb_above, glb_below)
    
    title = "Global average scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	p <- ( ggplot(plot_df, aes(year, scalar, color=agglt))
			+ scale_shape_manual(values=1:length(unique(plot_df$agglt)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= agglt), size = 1.5)
			+ geom_ribbon(aes(ymin = scalar_min, ymax = scalar_max, fill=agglt, linetype=NA), alpha=0.1, show.legend=FALSE)
			+ facet_wrap(~type, ncol=1)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			#+ ggtitle(title)
			+labs(color="land type", shape="land type")
			+scale_x_continuous(breaks=unique(plot_df$year), labels=as.character(unique(plot_df$year)))
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	ggsave(paste0(outdir, "globe_scalar_plots_paper.pdf"), plot=p, device="pdf")

    #######	
	
	#########################################################################
	# regional plots	
	
	reg_outname = paste0(outdir,"region_scalar_plots.pdf")
	subreg_outname = paste0(outdir,"subregion_scalar_plots.pdf")
	
	num_region_plots = num_reg * 6
	region_plots = vector(num_region_plots, mode="list")
	
	num_subregion_plots = num_lunits * 6
	subregion_plots = vector(num_subregion_plots, mode="list")
	
	# loop over the regions to make separate plots
	rp_ind = 0
	sp_ind = 0
	for (r in regions) {
		
		# crops above
		title = paste("Region", r, "crop vegetation scalars")
		xlabel = "Year"
		ylabel = "Scalar value (unitless)"
		plot_df = region_df[region_df$land_type%in%crops & region_df$region == r,]
		p<- ( ggplot(plot_df, aes(year, above_scalar, color=land_type))
				+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
				+ geom_line(size = 0.3)
				+ geom_point(aes(shape= land_type), size = 1.5)
				+ ylab( ylabel )
				+ xlab( xlabel )
				+ theme(legend.key.size = unit(0.4,"cm"))
				+ ggtitle(title)
				)
		p$save_args <- FIGURE_DIMS
		print(p)
		rp_ind = rp_ind + 1
		region_plots[[rp_ind]] <- recordPlot()
	
		# crops below
		title = paste("Region", r, "crop soil scalars")
		xlabel = "Year"
		ylabel = "Scalar value (unitless)"
		plot_df = region_df[region_df$land_type%in%crops & region_df$region == r,]
		p <- ( ggplot(plot_df, aes(year, below_scalar, color=land_type))
				+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
				+ geom_line(size = 0.3)
				+ geom_point(aes(shape= land_type), size = 1.5)
				+ ylab( ylabel )
				+ xlab( xlabel )
				+ theme(legend.key.size = unit(0.4,"cm"))
				+ ggtitle(title)
				)
		p$save_args <- FIGURE_DIMS
		print(p)
		rp_ind = rp_ind + 1
		region_plots[[rp_ind]] <- recordPlot()
	
		# non-crops and non-static above
		title = paste("Region", r, "non-crop vegetation scalars") 
		xlabel = "Year"
		ylabel = "Scalar value (unitless)"
		plot_df = region_df[(region_df$land_type%in%forest | region_df$land_type%in%pasture |
		                    region_df$land_type%in%shrubland | region_df$land_type%in%grassland) & region_df$region == r,]
		p <- ( ggplot(plot_df, aes(year, above_scalar, color=land_type))
				+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
				+ geom_line(size = 0.3)
				+ geom_point(aes(shape= land_type), size = 1.5)
				+ ylab( ylabel )
				+ xlab( xlabel )
				+ theme(legend.key.size = unit(0.4,"cm"))
				+ ggtitle(title)
				)
		p$save_args <- FIGURE_DIMS
		print(p)
		rp_ind = rp_ind + 1
		region_plots[[rp_ind]] <- recordPlot()
	
		# non-crops and non-static below
		title = paste("Region", r, "non-crop soil scalars")
		xlabel = "Year"
		ylabel = "Scalar value (unitless)"
		plot_df = region_df[(region_df$land_type%in%forest | region_df$land_type%in%pasture |
		                    region_df$land_type%in%shrubland | region_df$land_type%in%grassland) & region_df$region == r,]
		p <- ( ggplot(plot_df, aes(year, below_scalar, color=land_type))
				+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
				+ geom_line(size = 0.3)
				+ geom_point(aes(shape= land_type), size = 1.5)
				+ ylab( ylabel )
				+ xlab( xlabel )
				+ theme(legend.key.size = unit(0.4,"cm"))
				+ ggtitle(title)
				)
		p$save_args <- FIGURE_DIMS
		print(p)
		rp_ind = rp_ind + 1
		region_plots[[rp_ind]] <- recordPlot()
		
		
		# region average scalar values for aggregated types
	
		# vegetation
		title = paste("Region", r, "average vegetation scalars")
		xlabel = "Year"
		ylabel = "Scalar value (unitless)"
		plot_df = region_agglt_df[region_agglt_df$region == r,]
		p <- ( ggplot(plot_df, aes(year, above_scalar, color=agglt))
			+ scale_shape_manual(values=1:length(unique(plot_df$agglt)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= agglt), size = 1.5)
			+ geom_ribbon(aes(ymin = above_scalar_min, ymax = above_scalar_max, fill=agglt, linetype=NA), alpha=0.1, show.legend=FALSE)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			+labs(color="land type", shape="land type")
			+scale_x_continuous(breaks=unique(plot_df$year), labels=as.character(unique(plot_df$year)))
			)
		p$save_args <- FIGURE_DIMS
		print(p)
		rp_ind = rp_ind + 1
		region_plots[[rp_ind]] <- recordPlot()

		# soil
		title = paste("Region", r, "average soil scalars")
		xlabel = "Year"
		ylabel = "Scalar value (unitless)"
		plot_df = region_agglt_df[region_agglt_df$region == r,]
		p <- ( ggplot(plot_df, aes(year, below_scalar, color=agglt))
			+ scale_shape_manual(values=1:length(unique(plot_df$agglt)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= agglt), size = 1.5)
			+ geom_ribbon(aes(ymin = below_scalar_min, ymax = below_scalar_max, fill=agglt, linetype=NA), alpha=0.1, show.legend=FALSE)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			+ ggtitle(title)
			+labs(color="land type", shape="land type")
			+scale_x_continuous(breaks=unique(plot_df$year), labels=as.character(unique(plot_df$year)))
			)
		p$save_args <- FIGURE_DIMS
		print(p)
		rp_ind = rp_ind + 1
		region_plots[[rp_ind]] <- recordPlot()
		
		
		
		
		
		##################################################
		# loop over the subregions in each region to make more plots

		
		subregions = lunits$basin_abr[lunits$region == r]
		for (s in subregions) {
			
			# crops above
			title = paste("Region", r, "basin", s, "crop vegetation scalars")
			xlabel = "Year"
			ylabel = "Scalar value (unitless)"
			plot_df = all_df[all_df$land_type%in%crops & all_df$region == r & all_df$basin_abr == s,]
			p<- ( ggplot(plot_df, aes(year, above_scalar, color=land_type))
					+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
					+ geom_line(size = 0.3)
					+ geom_point(aes(shape= land_type), size = 1.5)
					+ ylab( ylabel )
					+ xlab( xlabel )
					+ theme(legend.key.size = unit(0.4,"cm"))
					+ ggtitle(title)
					)
			p$save_args <- FIGURE_DIMS
			print(p)
			sp_ind = sp_ind + 1
			subregion_plots[[sp_ind]] <- recordPlot()
	
			# crops below
			title = paste("Region", r, "basin", s, "crop soil scalars")
			xlabel = "Year"
			ylabel = "Scalar value (unitless)"
			plot_df = all_df[all_df$land_type%in%crops & all_df$region == r & all_df$basin_abr == s,]
			p <- ( ggplot(plot_df, aes(year, below_scalar, color=land_type))
					+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
					+ geom_line(size = 0.3)
					+ geom_point(aes(shape= land_type), size = 1.5)
					+ ylab( ylabel )
					+ xlab( xlabel )
					+ theme(legend.key.size = unit(0.4,"cm"))
					+ ggtitle(title)
					)
			p$save_args <- FIGURE_DIMS
			print(p)
			sp_ind = sp_ind + 1
			subregion_plots[[sp_ind]] <- recordPlot()
	
			# non-crops and non-static above
			title = paste("Region", r, "basin", s, "non-crop vegetation scalars")
			xlabel = "Year"
			ylabel = "Scalar value (unitless)"
			plot_df = all_df[(all_df$land_type%in%forest | all_df$land_type%in%pasture |
			                    all_df$land_type%in%shrubland | all_df$land_type%in%grassland) & all_df$region == r & all_df$basin_abr == s,]
			p <- ( ggplot(plot_df, aes(year, above_scalar, color=land_type))
					+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
					+ geom_line(size = 0.3)
					+ geom_point(aes(shape= land_type), size = 1.5)
					+ ylab( ylabel )
					+ xlab( xlabel )
					+ theme(legend.key.size = unit(0.4,"cm"))
					+ ggtitle(title)
					)
			p$save_args <- FIGURE_DIMS
			print(p)
			sp_ind = sp_ind + 1
			subregion_plots[[sp_ind]] <- recordPlot()
	
			# non-crops and non-static below
			title = paste("Region", r, "basin", s, "non-crop soil scalars")
			xlabel = "Year"
			ylabel = "Scalar value (unitless)"
			plot_df = all_df[(all_df$land_type%in%forest | all_df$land_type%in%pasture |
			                    all_df$land_type%in%shrubland | all_df$land_type%in%grassland) & all_df$region == r & all_df$basin_abr == s,]
			p <- ( ggplot(plot_df, aes(year, below_scalar, color=land_type))
					+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
					+ geom_line(size = 0.3)
					+ geom_point(aes(shape= land_type), size = 1.5)
					+ ylab( ylabel )
					+ xlab( xlabel )
					+ theme(legend.key.size = unit(0.4,"cm"))
					+ ggtitle(title)
					)
			p$save_args <- FIGURE_DIMS
			print(p)
			sp_ind = sp_ind + 1
			subregion_plots[[sp_ind]] <- recordPlot()
			
			
			# subregion average scalar values for aggregated types
	
			# vegetation
			title = paste("Region", r, "basin", s, "average vegetation scalars")
			xlabel = "Year"
			ylabel = "Scalar value (unitless)"
			plot_df = all_agglt_df[all_agglt_df$region == r & all_agglt_df$basin_abr == s,]
			p <- ( ggplot(plot_df, aes(year, above_scalar, color=agglt))
				+ scale_shape_manual(values=1:length(unique(plot_df$agglt)))
				+ geom_line(linewidth = 0.3)
				+ geom_point(aes(shape= agglt), size = 1.5)
				+ geom_ribbon(aes(ymin = above_scalar_min, ymax = above_scalar_max, fill=agglt, linetype=NA), alpha=0.1, show.legend=FALSE)
				+ ylab( ylabel )
				+ xlab( xlabel )
				+ theme(legend.key.size = unit(0.4,"cm"))
				+ ggtitle(title)
				+labs(color="land type", shape="land type")
				+scale_x_continuous(breaks=unique(plot_df$year), labels=as.character(unique(plot_df$year)))
			)
			p$save_args <- FIGURE_DIMS
			print(p)
			sp_ind = sp_ind + 1
			subregion_plots[[sp_ind]] <- recordPlot()

			# soil
			title = paste("Region", r, "basin", s, "average soil scalars")
			xlabel = "Year"
			ylabel = "Scalar value (unitless)"
			plot_df = all_agglt_df[all_agglt_df$region == r & all_agglt_df$basin_abr == s,]
			p <- ( ggplot(plot_df, aes(year, below_scalar, color=agglt))
				+ scale_shape_manual(values=1:length(unique(plot_df$agglt)))
				+ geom_line(linewidth = 0.3)
				+ geom_point(aes(shape= agglt), size = 1.5)
				+ geom_ribbon(aes(ymin = below_scalar_min, ymax = below_scalar_max, fill=agglt, linetype=NA), alpha=0.1, show.legend=FALSE)
				+ ylab( ylabel )
				+ xlab( xlabel )
				+ theme(legend.key.size = unit(0.4,"cm"))
				+ ggtitle(title)
				+labs(color="land type", shape="land type")
				+scale_x_continuous(breaks=unique(plot_df$year), labels=as.character(unique(plot_df$year)))
			)
			p$save_args <- FIGURE_DIMS
			print(p)
			sp_ind = sp_ind + 1
			subregion_plots[[sp_ind]] <- recordPlot()
			
			
		} # end for loop over subregions
		
	} # end for loop over regions
    
	
	# write the regional plots to pdf
    pdf(reg_outname, width=6.5,height=3.25)
	
	for (i in 1:num_region_plots) {
		replayPlot(region_plots[[i]])
	}
	
	dev.off()
	
	# write the subregional plots to pdf
	pdf(subreg_outname, width=6.5,height=3.25)
	
	for (i in 1:num_subregion_plots) {
		replayPlot(subregion_plots[[i]])
	}
	
	dev.off()
    
    ####### reginoal example timeseries crop scalar plots for the paper
    
    usa_above = region_agglt_df[region_agglt_df$region == r, c("year", "region", "agglt", "above_scalar", "above_scalar_min", "above_scalar_max")]
    names(glb_above) = c("year", "region", "agglt", "scalar", "scalar_min", "scalar_max")
    glb_above$type = "a) Vegetation"
    glb_below = globe_agglt_df[,c("year", "region", "agglt", "below_scalar", "below_scalar_min", "below_scalar_max")]
    names(glb_below) = c("year", "region", "agglt", "scalar", "scalar_min", "scalar_max")
    glb_below$type = "b) Soil"
    
    plot_df = region_df[region_df$land_type%in%crops & (region_df$region == "USA" | region_df$region == "Japan"),]
    plot_df$region[plot_df$region=="USA"] = "a) USA"
    plot_df$region[plot_df$region=="Japan"] = "b) Japan"
    x_breaks = unique(plot_df$year)[( (1:(length(unique(plot_df$year))/2))*2 ) - 1]
    
    title = "Average regional crop scalars"
	xlabel = "Year"
	ylabel = "Scalar value (unitless)"
	p <- ( ggplot(plot_df, aes(year, above_scalar, color=land_type))
			+ scale_shape_manual(values=1:length(unique(plot_df$land_type)))
			+ geom_line(linewidth = 0.3)
			+ geom_point(aes(shape= land_type), size = 1.5)
			+ facet_wrap(~region, ncol=1)
			+ ylab( ylabel )
			+ xlab( xlabel )
			+ theme(legend.key.size = unit(0.4,"cm"))
			#+ ggtitle(title)
			+labs(color="crop", shape="crop")
			+scale_x_continuous(breaks=x_breaks, labels=as.character(x_breaks))
			)
	p$save_args <- FIGURE_DIMS
	print(p)
	ggsave(paste0(outdir, "region_scalar_plots_paper.pdf"), plot=p, device="pdf")

    #######	
    
    ############## plot scalar maps
    
    ##### scheme for plotting maps
	scheme_basic <- theme_bw() +
  		theme(legend.text = element_text(size = 15)) +
  		theme(legend.title = element_text(size = 15)) +
  		theme(axis.text = element_text(size = 18)) +
  		theme(axis.title = element_text(size = 18, face = "bold")) +
  		theme(plot.title = element_text(size = 15, face = "bold", vjust = 1)) +
  		theme(plot.subtitle = element_text(size = 9, face = "bold", vjust = 1))+ 
  		theme(strip.text = element_text(size = 7))+
  		theme(strip.text.x = element_text(size = 18, face = "bold"))+
  		theme(strip.text.y = element_text(size = 15, face = "bold"))+
  		theme(legend.position = "right")+
  		theme(legend.text = element_text(size = 12))+
  		theme(legend.title = element_text(size = 12,color = "black",face="bold"))+
  		theme(axis.text.x= element_text(hjust=1,angle=90))+
  		theme(legend.background = element_blank(),
    	    legend.box.background = element_rect(colour = "black"))
    
    #### plot the aggregated land type maps
    
    ## region X glu
    
    pdf(paste0(outdir, regglu_map_outname, "_", map_year, ".pdf"), width=6.5,height=3.25)
    
    shapefile <- sf::st_read(paste0(gcam_boundaries_dir, gcam_reg_glu_shape_fname))
	plot_map_df = left_join( shapefile, all_agglt_df, join_by(reg_id == gcam_reg_code, glu_id == gcam_basin_code))
	plot_map_df = plot_map_df[plot_map_df$year == map_year,]
    
    # same scale for veg across types and same scale for soil across types
    amax = max(plot_map_df$above_scalar, na.rm=TRUE)
    amin = min(plot_map_df$above_scalar, na.rm=TRUE)
    abreaks = seq(from=round(amin,2), to=round(amax,2), by=round((round(amax,2) - round(amin,2))/5, 2))
    
    bmax = max(plot_map_df$below_scalar, na.rm=TRUE)
    bmin = min(plot_map_df$below_scalar, na.rm=TRUE)
    bbreaks = seq(from=round(bmin,2), to=round(bmax,2), by=round((round(bmax,2) - round(bmin,2))/5, 2))
    
    ## loop over agglt
    alts = unique(plot_map_df$agglt[!is.na(plot_map_df$agglt)])
    for(l in alts){
    
	    plot_df = plot_map_df[plot_map_df$agglt == l,]
    
		g <- ggplot()+
	     	geom_sf(data= plot_df, aes(fill=above_scalar))+
	     	scale_fill_gradient2(low="red", high="blue", mid="white", midpoint=1, limits=c(amin,amax), breaks=abreaks, labels=abreaks, na.value="white")+
	     	#scale_fill_viridis(limits=c(amin,amax), breaks=abreaks, labels=abreaks, na.value="grey")+
	     	ggtitle(paste(map_year, l, "vegetation scalars"))+
	     	scheme_basic

		print(g)

	
		g <- ggplot()+
	     	geom_sf(data= plot_df, aes(fill=below_scalar))+
	     	scale_fill_gradient2(low="red", high="blue", mid="white", midpoint=1, limits=c(bmin,bmax), breaks=bbreaks, labels=bbreaks, na.value="white")+
	     	#scale_fill_viridis(limits=c(bmin,bmax), breaks=bbreaks, labels=bbreaks, na.value="grey")+
	     	ggtitle(paste(map_year, l, "soil scalars"))+
	     	scheme_basic

		print(g)

    } # end regXglu loop over agglt
    
    dev.off()
    
    ############ agglt veg scalar maps for paper
    
    plot_lts = c("Crop", "Forest", "Grassland", "Shrubland")
    plot_df = plot_map_df[plot_map_df$agglt %in% plot_lts,]
    plot_df$agglt[plot_df$agglt=="Crop"] = "a) Crop"
    plot_df$agglt[plot_df$agglt=="Forest"] = "b) Forest"
    plot_df$agglt[plot_df$agglt=="Grassland"] = "c) Grassland"
    plot_df$agglt[plot_df$agglt=="Shrubland"] = "d) Shrubland"
    g <- ggplot()+
	    geom_sf(data= plot_df, aes(fill=above_scalar))+
	   	scale_fill_gradient2(low="red", high="blue", mid="white", midpoint=1, limits=c(amin,amax), breaks=abreaks, labels=abreaks, na.value="white")+
	    #ggtitle(paste(map_year, l, "vegetation scalars"))+
	    facet_wrap(~agglt, ncol=1) +
	    labs(fill = "Value") +
	    scheme_basic +
	    theme(axis.text.x= element_blank()) +
	    theme(axis.ticks= element_blank()) +
	    theme(strip.text.x = element_text(size = 12, face="plain"))

	print(g)
    ggsave(paste0(outdir, regglu_map_outname, "_", map_year, "_paper.pdf"), plot=g, device="pdf")
    
    
    ############
    
    
    ## region
    
    pdf(paste0(outdir, region_map_outname, "_", map_year, ".pdf"), width=6.5,height=3.25)
    
    shapefile <- sf::st_read(paste0(gcam_boundaries_dir, gcam_region_shape_fname))
	plot_map_df = left_join( shapefile, region_agglt_df, join_by(reg_id == gcam_reg_code) )
	plot_map_df = plot_map_df[plot_map_df$year == map_year,]
    
    # same scale for veg across types and same scale for soil across types
    amax = max(plot_map_df$above_scalar, na.rm=TRUE)
    amin = min(plot_map_df$above_scalar, na.rm=TRUE)
    abreaks = seq(from=round(amin,2), to=round(amax,2), by=round((round(amax,2) - round(amin,2))/5, 2))
    
    bmax = max(plot_map_df$below_scalar, na.rm=TRUE)
    bmin = min(plot_map_df$below_scalar, na.rm=TRUE)
    bbreaks = seq(from=round(bmin,2), to=round(bmax,2), by=round((round(bmax,2) - round(bmin,2))/5, 2))
    
    ## loop over agglt
    alts = unique(plot_map_df$agglt)
    for(l in alts){
    
    plot_df = plot_map_df[plot_map_df$agglt == l,]
    
	g <- ggplot()+
     	geom_sf(data= plot_df, aes(fill=above_scalar), na.rm=TRUE)+
     	scale_fill_gradient2(low="red", high="blue", mid="white", midpoint=1, limits=c(amin,amax), breaks=abreaks, labels=abreaks, na.value="white")+
     	#scale_fill_viridis(limits=c(amin,amax), breaks=abreaks, labels=abreaks, na.value="grey")+
     	ggtitle(paste(map_year, l, "vegetation scalars"))+
     	scheme_basic

	print(g)
	
	g <- ggplot()+
     	geom_sf(data= plot_df, aes(fill=below_scalar))+
     	scale_fill_gradient2(low="red", high="blue", mid="white", midpoint=1, limits=c(bmin,bmax), breaks=bbreaks, labels=bbreaks, na.value="white")+
     	#scale_fill_viridis(limits=c(bmin,bmax), breaks=bbreaks, labels=bbreaks, na.value="grey")+
     	ggtitle(paste(map_year, l, "soil scalars"))+
     	scheme_basic

	print(g)

    } # end region loop over agglt
    
    dev.off()
    
    #### write the data to files
    write.csv(all_df, paste0(outdir, "regXglu_scalar_data.csv"), row.names=FALSE)
    write.csv(region_df, paste0(outdir, "region_scalar_data.csv"), row.names=FALSE)
    write.csv(globe_df, paste0(outdir, "globe_scalar_data.csv"), row.names=FALSE)
    write.csv(all_agglt_df, paste0(outdir, "regXglu_scalar_data_agglt.csv"), row.names=FALSE)
    write.csv(region_agglt_df, paste0(outdir, "region_scalar_data_agglt.csv"), row.names=FALSE)
    write.csv(globe_agglt_df, paste0(outdir, "globe_scalar_data_agglt.csv"), row.names=FALSE)
    
                                   	
	cat("finished plot_e3sm_gcam_scalars.r", date())
                                   	
}