library(raster)
library(dplyr)
library(sf)
library(lidR)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(beepr)
library(terra)

aoi<-shapefile(".../dati/shp/area_polygon/area_ledro.shp")
chm_als<-raster(vrt(x = list.files(path = ".../dati/ALS/ledro/chm/",  full.names = TRUE)))
dtm_als<-raster(vrt(x = list.files(path = ".../dati/ALS/ledro/dtm/",  full.names = TRUE)))
chm_meta<-raster(".../dati/Meta/ledro/CanopyHeight_AOI_ledro.tif")



# Function to count the number of pixels per unique raster value
count_pixels_per_value <- function(raster_layer) {
  values <- getValues(raster_layer)
  values <- values[!is.na(values)]
  unique_values <- unique(values)
  
  counts <- sapply(unique_values, function(value) {
    sum(values == value, na.rm = TRUE)
  })
  
  data.frame(Value = unique_values, Count = counts)
}



# Load 400 random points
random_points <- shapefile("...qgis/shp/random_points/random_points_ledro.shp")


# Convert random points to sf object
random_points_sf <- st_as_sf(random_points)

# Define the buffer size (100x100 meters)
buffer_size <- 50 # Since the buffer radius is half the side length

# Initialize lists to store results
results_list <- list()

#delete what inside the folder that will be filled
unlink(".../dati_results/ledro/shp_point_trees/*")
unlink(".../results/ledro/chm_trees/*")
unlink(".../results/ledro/rds/final_results_ledro_400_lm4.rds")


# Function to compute window size based on mean canopy height
variable_window_size <- function(mean_height) {
  return(mean_height * 0.1 + 3)
}

# Function to perform tree top detection using a variable window size based on mean CHM
perform_tree_detection <- function(chm, mean_chm) {
  # Determine the window size based on the mean CHM
  window_size <- variable_window_size(mean_chm)
  
  # Perform tree top detection using the calculated window size
  ttops <- locate_trees(chm, lmf(window_size))
  
  # Segment trees using the detected tree tops
  trees <- dalponte2016(chm = chm, treetops = ttops)()
  
  return(list(ttops = ttops, trees = trees))
}


# Loop through each random point
for (i in 1:nrow(random_points_sf)) {
  point <- random_points_sf[i, ]
  
  # Create a 100x100m buffer around the point
  buffer <- st_buffer(point, dist = buffer_size)
  
  # Crop the CHM and DTM rasters to the buffer and set value between 0 and 45
  chm_crop_als <- crop(chm_als, buffer)
  chm_crop_als[chm_crop_als < 0 | chm_crop_als > 45] <- NA
  
  dtm_crop <- crop(dtm_als, buffer)
  chm_crop_meta <- crop(chm_meta, buffer)
  
  # Calculate canopy cover for ALS
  canopy_pixels_als <- sum(getValues(chm_crop_als) > 2, na.rm = TRUE)
  total_pixels_als <- sum(!is.na(getValues(chm_crop_als)))
  canopy_cover_als <- (canopy_pixels_als / total_pixels_als) * 100
  
  # Calculate canopy cover for Meta
  canopy_pixels_meta <- sum(getValues(chm_crop_meta) > 2, na.rm = TRUE)
  total_pixels_meta <- sum(!is.na(getValues(chm_crop_meta)))
  canopy_cover_meta <- (canopy_pixels_meta / total_pixels_meta) * 100
  
  
  # Calculate mean, max and min CHM for ALS
  mean_chm_als <- mean(getValues(chm_crop_als), na.rm = TRUE)
  max_chm_als <- max(getValues(chm_crop_als), na.rm = TRUE)
  min_chm_als <- min(getValues(chm_crop_als), na.rm = TRUE)
  
  # Calculate mean, max and min CHM for Meta
  mean_chm_meta <- mean(getValues(chm_crop_meta), na.rm = TRUE)
  max_chm_meta <- max(getValues(chm_crop_meta), na.rm = TRUE)
  min_chm_meta <- min(getValues(chm_crop_meta), na.rm = TRUE)
  
  
  # Calculate tree tops using dalponte2016 function for chm_als
  tree_detection_als <- perform_tree_detection(chm_crop_als, mean_chm_als)
  ttops_als <- tree_detection_als$ttops
  trees_als <- tree_detection_als$trees
  
  # Calculate tree tops using dalponte2016 function for chm_meta
  tree_detection_meta <- perform_tree_detection(chm_crop_meta, mean_chm_meta)
  ttops_meta <- tree_detection_meta$ttops
  trees_meta <- tree_detection_meta$trees
  
  # Skip analysis if no trees are detected in CHM ALS
  if (nrow(ttops_als) == 0) {
    message("No trees detected for CHM ALS at point ", i, ", skipping.")
    next
  }
  # Skip analysis if no trees are detected in CHM Meta
  if (nrow(ttops_meta) == 0) {
    message("No trees detected for CHM Meta at point ", i, ", skipping.")
    next
  }
  
  
  
  # Calculate average slope, exposition, and altitude
  slope <- terrain(dtm_crop, opt = "slope", unit = "degrees")
  aspect <- terrain(dtm_crop, opt = "aspect", unit = "degrees")
  altitude <- getValues(dtm_crop)
  
  avg_slope <- mean(getValues(slope), na.rm = TRUE)
  avg_aspect <- mean(getValues(aspect), na.rm = TRUE)
  avg_altitude <- mean(altitude, na.rm = TRUE)
  
  # Count the number of pixels for each unique value in trees_als and trees_meta
  pixel_counts_als <- count_pixels_per_value(trees_als)
  pixel_counts_meta <- count_pixels_per_value(trees_meta)
  
  # Add crown area information to ttops_als and ttops_meta
  ttops_als$crown <- pixel_counts_als$Count[match(ttops_als$treeID, pixel_counts_als$Value)]
  ttops_meta$crown <- pixel_counts_meta$Count[match(ttops_meta$treeID, pixel_counts_meta$Value)]
  
  ttops_als$cc <- canopy_cover_als
  ttops_meta$cc <- canopy_cover_meta
  
  
  # Save ttops_als and ttops_meta as shapefiles
  setwd(".../dati_results/ledro/shp_point_trees/")
  st_write(ttops_als, paste0("ttops_als_", i, ".shp"), layer_options = "SHPT=POINTZ")
  st_write(ttops_meta, paste0("ttops_meta_", i, ".shp"), layer_options = "SHPT=POINTZ")
  
  setwd(".../dati_results/ledro/chm_trees/")
  # Save trees_als and trees_meta rasters
  writeRaster(trees_als, paste0("trees_als_", i, ".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(trees_meta, paste0("trees_meta_", i, ".tif"), format = "GTiff", overwrite = TRUE)
  
  # Store results in a list
  results_list[[i]] <- data.frame(
    Point_ID = i,
    Canopy_Cover_als = canopy_cover_als,
    Canopy_Cover_meta = canopy_cover_meta,
    Avg_Slope = avg_slope,
    Avg_Aspect = avg_aspect,
    Avg_Altitude = avg_altitude,
    Trees_als = nrow(ttops_als),
    Trees_meta = nrow(ttops_meta),
    Mean_chm_als = mean_chm_als,
    Max_chm_als = max_chm_als,
    Min_chm_als = min_chm_als,
    Mean_chm_meta = mean_chm_meta,
    Max_chm_meta = max_chm_meta,
    Min_chm_meta = min_chm_meta
  )
}


# Combine results into a single data frame
results_df <- bind_rows(results_list)

saveRDS(results_df, ".../dati_results/ledro/rds/final_results_ledro_400_lm4.rds")
