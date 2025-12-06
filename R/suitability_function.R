#' Function: Mapping Aquaculture Suitability Using Sea Surface Temperature and Bathymetry as Constraints  
#'
#' @param tmin : species' suitable minimum temperature in °Celsius
#' @param tmax : species' suitable maximum temperature in °Celsius
#' @param dmin : species' suitable minimum depth in meters (negative below sea level)
#' @param dmax : species' suitable depth in meters (negative below sea level)
#' @param species_name : character string of the species name, either common or scientific, for the map title
#'
#' @returns : a table of suitable area (km²) per Exclusive Economic Zone (EEZ) and a map
#' @export
#'
#' @examples
map_suitability <- function(tmin, tmax, 
                            dmin, dmax, 
                            species_name = "species") {
  # Load and process SST
  sst_files <- list.files(here::here("data"), 
                          pattern = "average_annual_sst", 
                          full.names = TRUE)
  
  # Stack SST rasters (2008–2012)
  
  sst_stack <- c(SST_2008, SST_2009, SST_2010, SST_2011, SST_2012)

  # Compute mean SST over time (per pixel)
  
  mean_sst <- terra::mean(sst_stack)
  
  # Convert temperature units from Kelvin to Celsius
  
  mean_sst_c <- mean_sst - 273.15 
  
  # Prepare bathymetry (depth) raster
  
  # Load the bathymetry raster
  bathy <- rast(here("data", "depth.tif"))
  
  # Make CRS match SST:
  
  bathy_proj <- project(bathy, crs(mean_sst_c))
  
  
  # Crop bathymetry to SST's extent
  
  bathy_crop <- crop(bathy_proj, mean_sst_c)
  
  
  # Resample bathymetry to SST resolution (so each cell lines up):
  
  depth_res <- resample(bathy_crop, mean_sst_c, method = "near")
  
  
  #...............................................................................
  #                                                                              .
  #  Find Suitable Locations                                                     .
  #                                                                              .
  #...............................................................................
  
  
  # Reclassify SST into suitable / unsuitable
  # The oyster thresholds for SST are between 11–30 Celsius
  
  sst_rcl <- matrix(c(-Inf, tmin,  0,
                      tmin, tmax,  1,
                      tmax, Inf, 0), ncol = 3, byrow = TRUE)
  
  sst_suit <- classify(mean_sst_c, rcl = sst_rcl)
  
  
  # Reclassify depth into suitable or unsuitable (1, 0)
  # The depth thresholds for oysters range is 0-70 meters b.s.l.
  
  depth_rcl <- matrix(c(-Inf, dmin, 0,
                        dmin, dmax, 1,
                        dmax, Inf, 0), ncol = 3, byrow = TRUE)
  
  depth_suit <- classify(depth_res, depth_rcl)
  
  # Combine SST and depth suitability
  
  # Multiply to get 1 only where both are suitable
  suitability <- sst_suit * depth_suit # 1 = suitable, 0 = unsuitable
  
  # Prepare EEZ vector and mask
  eez_sf <- st_read(here("data", "wc_regions_clean.shp"), quiet = TRUE)
  
  # convert to SpatVector format
  eez <- vect(eez_sf)
  
  # Re-project with mean sst's CRS
  eez_proj <- project(eez, mean_sst_c)
  
  # Mask suitability by EEZ
  suitability_eez <- mask(suitability, eez_proj)
  
  # Count suitable cells per EEZ
  
  # a. Compute area in squared kilometers of each raster cell
  cell_area_km <- terra::cellSize(suitability_eez, unit = "km")
  
  # b. Multiply suitability (0/1) by area to get "suitable area per cell"
  suitable_area_raster <- suitability_eez * cell_area_km
  
  # c. Extract total suitable area (km²) for each EEZ polygon
  eez_suitable_area <- terra::extract(suitable_area_raster,
                               eez_proj,
                               fun = sum,
                               na.rm = TRUE)
  
  # d. Add the EEZ region names to the output
  # Replace "region" with the correct column name from your EEZ shapefile
  eez_suitable_area$region <- eez_proj$rgn
  
  # Rename the "mean"column
  names(eez_suitable_area)[2] <- "suitable_km2"
  
  # Format the table for EEZ suitable area totals
  eez_suitable_area_table <- eez_suitable_area %>%
    select(region, suitable_km2) %>%
    kable(col.names = c("EEZ Region", "Suitable Area (km²)"),
          caption = paste("Total Suitable Area for ",
                          species_name, 
                          " by EEZ Region"),
          align = "lr", digits = 1) %>%
    kable_styling(
      full_width = FALSE,
      bootstrap_options = c("striped", "hover"))
  
  # Merge table back into EEZ SpatVector
  eez_map <- eez_proj
  eez_map$suitable_km2 <- eez_suitable_area$suitable_km2
  
  # Create the final map
  map_final <- tm_basemap("Esri.OceanBasemap") +
    
    tm_shape(eez_map) +
    
    tm_polygons(fill = "suitable_km2", 
                fill.scale = tm_scale_intervals(values = "brewer.blues", 
                                                style = "quantile"),
                fill.legend = tm_legend(bg.color = "white",
                                        bg.alpha = 0.2,
                                        title = "Suitable Area (km²)",
                                        frame = FALSE,
                                        position = tm_pos_out("right", "center")),
                col = "black",
                lwd = 0.2) +
    
    tm_text("rgn", col = "black", fontface = "bold", size = 0.7) + 
    
    tm_title(size = 0.8,
             text = paste("Suitable Area for ",
                          species_name, "\n",
                          "Aquaculture in West Coast EEZ")) +
    
    tm_layout(frame = TRUE, frame.lwd = 0.4, frame.r = 0) +
    
    tm_scalebar(breaks = c(0, 250, 500), 
                position = tm_pos_in("left", "bottom")) +
    
    tm_compass(type = "arrow", position = tm_pos_in("right", "top")) +
    
    tmap_options(component.autoscale = FALSE)
  
  return(list(map = map_final,
              table = eez_suitable_area_table))
}
