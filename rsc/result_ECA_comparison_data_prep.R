library("terra")
library("tidyterra")

dir_dst <- "L:/poppman/shared/dami/dat/top9"
dir_res_maps <- "L:/poppman/shared/dami/dat/res"
dir_top_9 <- "L:/ogundipe/shared/Manuel/top9"
dir_pat_all <- "L:/poppman/shared/dami/dat/pat"

f_pat_all <- file.path(dir_pat_all, "patches_1.3_ha_dj.tif")
f_pa_existing <- file.path(dir_top_9, "PA_.shp")

f_scenario_maps <- list.files(dir_top_9, pattern = ".tif$", full.names = TRUE)

protected_areas <- terra::vect(f_pa_existing)

protected_areas_buff <- terra::buffer(
  protected_areas, width = 25
)

raster_template <- terra::rast(f_pat_all)

patches <- terra::rast(f_pat_all) %>%
  terra::as.polygons(
    dissolve = TRUE,
    values = TRUE
    ) %>%
  dplyr::mutate(
    area = terra::expanse(., unit = "ha")
  )

for (f_scenario_map in c(f_scenario_maps, "baseline")) {
  if (f_scenario_map == "baseline") {
    scenario <- protected_areas %>%
      terra::rasterize(
        raster_template
      )
  } else {
    scenario <- terra::rast(f_scenario_map) %>%
      terra::project(raster_template, method = "near")
  }
  
  pixelarea <- do.call("*", as.list(terra::res(scenario))) / 1e4
  
  scenario_area <- terra::freq(scenario) %>%
    dplyr::filter(value == 1) %>%
    dplyr::pull(count) * pixelarea
    
  
  # Create simple grid-based raster of top 9 percent areas plus existing PAs
  terra::rasterize(protected_areas, raster_template) %>%
    c(scenario) %>%
    terra::app(fun = "max", na.rm = TRUE) %>%
    terra::writeRaster(
      filename = file.path(
        dir_dst,
        paste0(
          tools::file_path_sans_ext(basename(f_scenario_map)),
          "_gridbased.tif"
          )
        ),
      datatype = "INT2S",
      overwrite = TRUE
      )
  
  if (f_scenario_map == "baseline") {
    file.copy(
      file.path(
        dir_dst,
        paste0(
          tools::file_path_sans_ext(basename(f_scenario_map)),
          "_gridbased.tif"
        )
      ),
      file.path(
        dir_dst,
        paste0(
          tools::file_path_sans_ext(basename(f_scenario_map)),
          "_patchbased.tif"
        )
      )
    )
    next
  }
  # Extract patch-level overlap with top 9 percent suitability
  patch_suitabilities <- terra::extract(
    scenario, patches, fun = "mean", exact = TRUE
    )
  
  patch_suitabilities[which(is.na(patch_suitabilities[, 2])), 2] <- 0
  suited <- patches[order(patch_suitabilities[, 2], decreasing = TRUE), ]
  suited$cumArea <- cumsum(suited$area)
  
  patch_selection <- suited[which(suited$cumArea < scenario_area),] %>%
    terra::rasterize(raster_template)
  
  terra::rasterize(protected_areas, raster_template) %>%
    c(patch_selection) %>%
    terra::app(fun = "max", na.rm = TRUE) %>%
    terra::writeRaster(
      filename = file.path(
        dir_dst,
        paste0(
          tools::file_path_sans_ext(basename(f_scenario_map)),
          "_patchbased.tif"
        )
      ),
      datatype = "INT2S",
      overwrite = TRUE
    )
}

# Calculate area effect of buffering PAs
# This is obsolete, since we decided to use the unbiffered patches and just
# assume a new PA will not immediately touch an existing one or, at least,
# that the effect of touching vs 25 m distance is negligible
area_pa <- sum(terra::expanse(protected_areas, unit = "ha"), na.rm = TRUE)
area_pa_buff <- sum(
  terra::expanse(protected_areas_buff, unit = "ha"),
  na.rm = TRUE
  )
diff_ha <- area_pa_buff - area_pa
