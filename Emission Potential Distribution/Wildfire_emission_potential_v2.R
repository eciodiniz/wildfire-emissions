library(data.table)
library(sf)
library(terra)
library(ggplot2)
library(scales)

# --------------------------------------------------------
# 1. Configuration: input/output paths
# --------------------------------------------------------
base_path <- "D:/AXA-Chile/Emissions_tests/Output"           # Main output folder (maps are saved here)
rds_path  <- file.path(base_path, "Pollutant_emissions_RDS")  # Folder with the pollutant emissions RDS file
shp_base  <- "D:/Burn_probability_results/KDE map"            # Folder with the study-area boundary shapefiles


pollutant <- "PM10"   # Pollutant to map (must match a name in the RDS list)

out_png <- file.path(base_path, paste0("Map_", pollutant, "_Regional.jpeg")) # Output map file path

areas <- c("Chillan", "Concepcion", "Los_Angeles", "SB (Curanilahue)", # Study areas: site names as in the RDS 'Site' column
           "San_Antonio", "Santiago", "Talca", "Valparaiso")

shp_file_map <- c(                         # Study area -> boundary shapefile name (without .shp)
  "Chillan"          = "Chillan",
  "Concepcion"       = "concepcion",
  "Los_Angeles"      = "los_angeles",
  "SB (Curanilahue)" = "curanilahue",
  "San_Antonio"      = "san_antonio",
  "Santiago"         = "santiago",
  "Talca"            = "Talca",
  "Valparaiso"       = "valparaiso"
)

crs_wgs84 <- "EPSG:4326"   # Coordinate reference system (WGS84, lon/lat)
res_grid  <- 0.01          # Raster cell size in degrees (~1 km)

# --------------------------------------------------------
# 2. Load pollutant data
# --------------------------------------------------------
pollutant_list <- readRDS(file.path(rds_path, "All_pollutants_list_corrected.rds"))  # List of emission tables, one per pollutant

if (!pollutant %in% names(pollutant_list)) {  # Stop if the selected pollutant is not in the RDS
  stop("Pollutant not found in RDS: ", pollutant)
}

dt_all <- copy(pollutant_list[[pollutant]]) # Working copy of the selected pollutant table

# Keep only spatially usable records (valid coordinates and emission value)
dt_all <- dt_all[!is.na(lon) & !is.na(lat) & !is.na(Emi_kton)]

# --------------------------------------------------------
# 3. Boundary shapefiles of the 8 study areas
# --------------------------------------------------------
read_area_shp <- function(area) {                              # Read and standardize one study-area boundary
  shp_name <- shp_file_map[[area]]                             # Shapefile name for this area
  shp_file <- file.path(shp_base, paste0(shp_name, ".shp"))    # Full path to the shapefile
  if (!file.exists(shp_file)) stop("File not found: ", shp_file)

  g <- st_read(shp_file, quiet = TRUE)                         # Read shapefile
  g <- st_transform(g, crs_wgs84)                              # Reproject to WGS84
  g <- st_zm(g, drop = TRUE)                                   # Drop Z/M dimensions (2D geometry only)
  g$Area <- area                                               # Add study-area name
  st_geometry(g) <- "geometry"                                 # Standardize geometry column name
  g[, "Area"]                                                  # Keep only area name + geometry
}

shp_all <- do.call(rbind, lapply(areas, read_area_shp))        # Combine all study-area boundaries
shp_all <- st_make_valid(shp_all)                              # Fix invalid geometries


# --------------------------------------------------------
# 4. Rasterization
# --------------------------------------------------------
ext_total <- ext(min(dt_all$lon), max(dt_all$lon),             # Spatial extent covering all emission points
                 min(dt_all$lat), max(dt_all$lat))

r_template <- rast(ext_total, resolution = res_grid, crs = crs_wgs84)  # Raster template

pts <- vect(dt_all, geom = c("lon", "lat"), crs = crs_wgs84)   # Emission points as a spatial vector

r_emis <- rasterize(pts, r_template, field = "Emi_kton",       # Total emissions per cell (kton)
                    fun = "sum", background = NA)              # Cells without points stay NA

rm(pts, dt_all)                                                # Free memory
gc()

# --------------------------------------------------------
# 5. Final map
# --------------------------------------------------------
r_df <- as.data.frame(r_emis, xy = TRUE, na.rm = TRUE)   # Raster to data frame for ggplot
names(r_df)[3] <- "Emi_kton"

# Color-scale cap (visualization only; emission values are not changed)
cap_map <- c("PM10" = 0.04)   # Fixed cap for PM10; other pollutants use the 97th percentile
cap <- if (pollutant %in% names(cap_map)) cap_map[[pollutant]] else quantile(r_df$Emi_kton, 0.97)
r_df$Emi_kton <- pmin(r_df$Emi_kton, cap)

p <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = Emi_kton)) +
  scale_fill_viridis_c(
    name = paste0(pollutant, " (kton)"),
    option = "inferno",
    na.value = "transparent",
    trans = "sqrt",     # Square-root scale to enhance low values
    labels = label_number(accuracy = 0.001)
  ) +
  geom_sf(data = shp_all, fill = NA, color = "white", linewidth = 0.6) +   # Study-area boundaries
  scale_x_continuous(breaks = seq(-74, -70, by = 0.5)) +
  coord_sf(expand = FALSE) +
  labs(title = paste(pollutant, "Emissions - Central Chile"),
       subtitle = " ", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.height = unit(1.6, "cm"),
    legend.key.width = unit(0.9, "cm")
  )
p

# --------------------------------------------------------
# 6. Export
# --------------------------------------------------------
ggsave(out_png, plot = p, width = 12, height = 14, dpi = 300, bg = "white")   # Save map (300 dpi)
cat("Map saved to:", out_png, "\n")

